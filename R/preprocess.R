#' Prepare chunks for signal loading
#'
#' Function to process signal and read mappings,
#' ready for signal loading
#'
#' @param nat_mapping Path to NAT sorted bam file
#' @param pcr_mapping Path to PCR sorted bam file
#' @param hdf5 list of nat and pcr hdf5 objects
#' @param chunk_size Size to chunk reference into (default 1e5)
#' @param debug logical, print additional information
#' @param pvp logical, load only PCR and split it into 1:1 labeled PCR and NAT
#' @return metainfo list with nested bacthes and types
#' @export
prepare_metainfo <- function(
                       nat_mapping,
                       pcr_mapping,
                       nat_hdf5,
                       pcr_hdf5,
                       chunk_size = 1e5,
                       threads = 1,
                       debug = FALSE,
                       pvp = FALSE) {
  if (debug) logger::log_threshold(logger::TRACE)

  logger::log_debug("Loading read mapping")
  pcr <- load_read_mapping(pcr_mapping)
  pcr[, type := "pcr"]

  pcr_hdf5_obj <- hdf5r::H5File$new(pcr_hdf5, mode = "r")
  if (pvp) {
    read_mapping <- pcr[
        sample(seq_len(.N), .N / 2), type := "nat"
      ]
    nat_hdf5_obj <- pcr_hdf5_obj
  } else {
    nat <- load_read_mapping(nat_mapping)
    nat[, type := "nat"]
    nat_hdf5_obj <- hdf5r::H5File$new(nat_hdf5, mode = "r")
    read_mapping <- rbind(pcr, nat)
  }

  read_mapping[
      , chunk := pos %/% chunk_size
    ]
  read_mapping <- rbind(
      read_mapping,
      process_multichunk_reads(read_mapping, chunk_size)
  )
  read_mapping <- downsample(read_mapping,  chunk_size = chunk_size)


  metainfo_nat <- load_metainfo_all(nat_hdf5_obj)
  metainfo_pcr <- load_metainfo_all(pcr_hdf5_obj)

  metainfo <- rbind(
    metainfo_nat[
      read_mapping[type == "nat"],
      on = c("read_id" = "qname"),
      nomatch = NULL
    ],
    metainfo_pcr[
      read_mapping[type == "pcr"],
      on = c("read_id" = "qname"),
      nomatch = NULL
    ]
  )
  data.table::setnames(metainfo, c("rname"), c("reference"))
  metainfo[
    , chunk_ref := paste(chunk, reference, sep = "_")
  ]
  rm(read_mapping)
  gc()
  no_cov_chunk <- metainfo[
      length(unique(type)) < 2, chunk_ref, by = chunk_ref
    ]$chunk_ref %>%
    unique()
  if (length(no_cov_chunk) > 0) {
    logger::log_warn(
      paste0("NAT or PCR missing: ", paste(no_cov_chunk, collapse = ", "))
    )
    metainfo <- metainfo[!(chunk_ref %in% no_cov_chunk)]
  }

  metainfo_chunk <- split(
    metainfo,
    by = c("chunk_ref", "batch", "type"),
    flatten = FALSE
  )
  return(metainfo_chunk)
}

################################################################################

#' Load sorted read mappings
#'
#' Load read mapping required for correct assignment of signal mapping reference.
#'
#' @param path_mapping Preprocessed signal mappings
#' @return read mapping data.table
#' @import data.table
#' @export
load_read_mapping <- function(path_mapping) {
  if (!file.exists(glue::glue("{path_mapping}.bai"))) {
    Rsamtools::indexBam(path_mapping)
  }
  what <- c("qname", "rname", "pos", "strand", "qwidth")
  param <- Rsamtools::ScanBamParam(what = what)
  read_mapping <-
    Rsamtools::scanBam(
      path_mapping,
      param = param,
      index = glue::glue("{path_mapping}.bai")
    ) %>%
    unlist(recursive = FALSE)
  setDT(read_mapping)
}

################################################################################

#' Get read mappings that overextent their chunk
#'
#' Gets the read that overextent their assigned chunk and return dt with new 
#' chunk assignments
#'
#'
#' @param read_mapping Read mapping dt
#' @param chunk_size Size of chunks
#' @return read mapping data.table
#' @import data.table
#' @export
process_multichunk_reads <- function(read_mapping, chunk_size) {
  cnames <- colnames(read_mapping)[!(colnames(read_mapping) == "chunk")]
  max_chunk <- max(read_mapping$chunk)
  corrected_read_mappings <- read_mapping[
      pos + qwidth > (1 + chunk) * chunk_size,
    ][
      , over_extension := pos + qwidth - (1 + chunk) * chunk_size
    ][
      , .(chunk = 1:(1 + over_extension %/% chunk_size) + chunk),
      by = mget(cnames)
    ][
      chunk <= max_chunk
    ]
  unnest_dt(corrected_read_mappings, "chunk", cnames)
}

################################################################################

#' Downsample either nat or pcr to least abundant type
#'
#' Removes read of the most bundant mapping type, such there are almost equal mapped bases of each type to each chunk.
#'
#'
#' @param read_mapping Preprocessed signal mappings
#' @param chunk_size size of chunks
#' @param min_cov Minimum coverage of type
#' @return read mapping data.table
#' @import data.table
#' @export
downsample <- function(read_mapping, chunk_size, min_cov = 20) {
  min_chunk_cov <- min_cov * chunk_size
  read_mapping <- read_mapping[sample(seq_len(nrow(read_mapping)))]
  read_mapping[
      , read_chunk_cov := fcase(
          pos + qwidth > (chunk + 1) * chunk_size,  as.integer((chunk + 1) * chunk_size - pos),
          pos < chunk * chunk_size, as.integer(pos + qwidth - chunk * chunk_size),
          rep(TRUE, length(pos)), qwidth
        )
    ][
      , coverage := cumsum(as.numeric(read_chunk_cov)), by = .(type, chunk, strand)
    ][
      , max_coverage := max(coverage), by = .(type, chunk, strand)
    ][
      , max_allowed := min(max_coverage), by = .(chunk, strand)
    ][
      (coverage < max_allowed) | (coverage < min_chunk_cov),
    ][
      , `:=`(max_coverage = NULL, max_allowed = NULL)
    ]
}

################################################################################

#' Read metainfo from signal_mappings.hdf
#'
#' This function load all  metainformation in a singal mapping hdf5 file.
#'
#' @param hdf5 Open hdf5 object
#' @return list of data.tables with metainformation of all batches
#' @export
load_metainfo_all <- function(hdf5) {
  batches <- get_batches(hdf5)
  rbindlist(
    lapply(
      batches,
      function(batch) {
        signal_mapping <- load_metainfo(hdf5, batch)
        signal_mapping[, batch := batch]
        return(signal_mapping)
      }
    )
  )
}

################################################################################

#' See which batches are present in signal_mappings.hdf
#'
#'
#'
#' @param hdf5 Open signal_mappings.hdf5 hdf5 object
#' @return vector with batch names
#' @export
get_batches <- function(hdf5) {
  hdf5r::list.datasets(hdf5) %>%
    stringr::str_extract("Batch_[\\d]*") %>%
    na.omit() %>%
    c() %>%
    unique()
}

################################################################################

#' Read metainfo from signal_mappings.hdf
#'
#' This function load metainformation of each read in a batch in a singal mapping hdf5 file.
#'
#' @param hdf5 Open hdf5 object
#' @param batch Batch to read from
#' @return data.table with metainformation
#' @export
load_metainfo <- function(hdf5, batch) {
  metainfo <- data.table::data.table(
    Ref_to_signal_lengths = hdf5[[glue::glue("Batches/{batch}/Ref_to_signal_lengths")]][],
    Dacs_lengths = hdf5[[glue::glue("Batches/{batch}/Dacs_lengths")]][],
    digitisation = hdf5[[glue::glue("Batches/{batch}/digitisation")]][],
    offset = hdf5[[glue::glue("Batches/{batch}/offset")]][],
    range = hdf5[[glue::glue("Batches/{batch}/range")]][],
    scale_frompA = hdf5[[glue::glue("Batches/{batch}/scale_frompA")]][],
    shift_frompA = hdf5[[glue::glue("Batches/{batch}/shift_frompA")]][],
    read_id = hdf5[[glue::glue("Batches/{batch}/read_id")]][]
  )
  # Add indices for Ref_to_signal and Dacs
  metainfo[
      , dacs_end := cumsum(as.numeric(Dacs_lengths))
    ][
      , dacs_start := dacs_end - Dacs_lengths + 1
    ][
      , ref_to_signal_end := cumsum(as.numeric(Ref_to_signal_lengths))
    ][
      , ref_to_signal_start := ref_to_signal_end - Ref_to_signal_lengths + 1
    ]
}
