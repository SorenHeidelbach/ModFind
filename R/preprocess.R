#' Run preprocess 
#' 
#' Function to process signal and read mappings into preprocessed .tsv files, ready for current difference calculation
#' 
#' @param nat_mapping Path to NAT sorted bam file
#' @param pcr_mapping Path to PCR sorted bam file
#' @param nat_signal HDF5 object of NAT signal mappings
#' @param pcr_signal HDF5 object of PCR signal mappings
#' @param out Output path for preprocessed signal mappings
#' @param chunk_size Size to chunk reference into (default 1e4 bases)
#' @param threads Number of jobs to run parallel (default 1; >1 not supported on windows)
#' @return metainfo with signal added to each leaf
#' @export
add_signal_chunk <- function(metainfo_chunk, hdf5) {
  batches <- names(metainfo_chunk)
  successfully_read <- lapply(batches,
  function(batch) {
    logger::log_debug(glue::glue("\tProcessing {batch}"))
    types <- names(metainfo_chunk[[batch]])
    lapply(types,
    function(type) {
      logger::log_debug(glue::glue("\t\t{type}"))
      tryCatch(
        add_signal(
          metainfo_chunk[[batch]][[type]],
          hdf5,
          batch,
          ref_to_signal
        ),
        error = function(e) FALSE
      )
      return(TRUE)
    } # type
    )
  } # batch
  )
}


#' Run preprocess 
#' 
#' Function to process signal and read mappings into preprocessed .tsv files, ready for current difference calculation
#' 
#' @param nat_mapping Path to NAT sorted bam file
#' @param pcr_mapping Path to PCR sorted bam file
#' @param nat_signal HDF5 object of NAT signal mappings
#' @param pcr_signal HDF5 object of PCR signal mappings
#' @param out Output path for preprocessed signal mappings
#' @param chunk_size Size to chunk reference into (default 1e4 bases)
#' @param threads Number of jobs to run parallel (default 1; >1 not supported on windows)
#' @return metainfo list with nesated bacthes and types
#' @export
prepare_metainfo <- function(
                       nat_mapping,
                       pcr_mapping,
                       nat_signal,
                       pcr_signal,
                       chunk_size = 1e5,
                       threads = 1,
                       debug = FALSE,
                       pvp = FALSE) {
  
  if (debug) {logger::log_threshold(logger::TRACE)}

  hdf5 <- list(
    pcr = pcr_signal,
    nat = nat_signal
  )

  logger::log_debug("Loading read mapping")
  pcr <- load_mapping(pcr_mapping)
  pcr[, type := "pcr"]

  if (pvp) {
    read_mapping <- pcr[
        sample(1:.N, .N/2), type := "nat"
      ]
  } else {
    nat <- load_mapping(nat_mapping)
    nat[, type := "nat"]
    read_mapping <- rbind(pcr, nat)
  }

  read_mapping[
      , chunk := pos %/% chunk_size
    ]
  read_mapping <- rbind(
      read_mapping,
      get_overextending_reads(read_mapping, chunk_size)
  )
  read_mapping <- downsample(read_mapping,  chunk_size = chunk_size)


  metainfo_nat <- read_metainfo_all(hdf5[["nat"]])
  metainfo_pcr <- read_metainfo_all(hdf5[["pcr"]])

  metainfo <- rbind(
    metainfo_nat[
      read_mapping[type == "nat"],
      on = c('read_id'='qname')
    ],
    metainfo_pcr[
      read_mapping[type == "pcr"],
      on = c('read_id'='qname')
    ]
  )
  data.table::setnames(metainfo, c("rname"), c("reference"))
  metainfo[
    , chunk_ref := paste(chunk, reference, sep = "_")
  ]
  rm(read_mapping)
  gc()

  metainfo_chunk <- split(
    metainfo,
    by = c("chunk_ref", "batch", "type"),
    flatten = FALSE
  )
  return(metainfo_chunk)
}

#' Read metainfo from signal_mappings.hdf
#' 
#' This function load metainformation of each read in a batch in a singal mapping hdf5 file.
#' 
#' @param hdf5 Open hdf5 object
#' @param batch Batch to read from
#' @return data.table with metainformation
#' @export
read_metainfo = function(hdf5, batch) {
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
      , dacs_end := cumsum(Dacs_lengths)
    ][
      , dacs_start := dacs_end - Dacs_lengths + 1
    ][
      , ref_to_signal_end := cumsum(Ref_to_signal_lengths)
    ][
      , ref_to_signal_start := ref_to_signal_end - Ref_to_signal_lengths + 1
    ]
}

#' Read metainfo from signal_mappings.hdf
#' 
#' This function load all  metainformation in a singal mapping hdf5 file.
#' 
#' @param hdf5 Open hdf5 object
#' @return list of data.tables with metainformation of all batches
#' @export
read_metainfo_all <- function(hdf5) {
  batches <- get_batches(hdf5)
  rbindlist(
    lapply(
      batches,
      function(batch){
        signal_mapping = read_metainfo(hdf5, batch)
        signal_mapping[, batch := batch]
        return(signal_mapping)
      }
    )
  )
}

#' See which batches are present in signal_mappings.hdf
#' 
#' 
#' 
#' @param hdf5 Open signal_mappings.hdf5 hdf5 object
#' @return vector with batch names
#' @export
get_batches <- function(hdf5){
  hdf5r::list.datasets(hdf5) %>%
    stringr::str_extract("Batch_[\\d]*") %>%
    na.omit() %>%
    c() %>%
    unique()
}

#' Add signal mapping to metainfo
#' 
#' Loads the signal mappings associated with a batch
#' 
#' @param metainfo data.table of metainfo loaded with read_metainfo
#' @param hdf5 Open hdf5 object
#' @param hdf_batch Batch to load signal mappings from
#' @return Nothing, signal mappings will be added to the metainfo object in memory
#' @export
add_signal = function(metainfo, hdf5, batch, ref_to_sig) {
  # Input checks
  assert::assert({
    required_columns <- c(
      "ref_to_signal_start",
      "ref_to_signal_end",
      "dacs_start",
      "dacs_end",
      "read_id",
      "offset",
      "range",
      "digitisation",
      "shift_frompA",
      "scale_frompA")
    all(required_columns %in% names(metainfo))
    },
    msg = paste0(
      "Following required columns are missing: ",
      paste0(
        required_columns[!required_columns %in% names(metainfo)],
        collapse = ", "
      )
    )
  )

  assert::assert({
    duplicated_read_ids <- which(duplicated(metainfo$read_id))
    length(duplicated_read_ids) == 0
    },
    msg = paste0(
      "Following rows have duplicate read_ids: ",
      paste0(duplicated_read_ids, collapse = ", ")
    )
  )
  logger::log_trace(glue::glue("\t\t\tAdding ref to signal"))
  metainfo[
      , ref_to_signal := list(
              list(hdf5[[glue::glue('/Batches/{batch}/Ref_to_signal')]][ref_to_signal_start:ref_to_signal_end])
              ),
      by = read_id
    ]
  logger::log_trace(glue::glue("\t\t\tAdding dacs"))
  metainfo[
      , dacs := list(
              list(hdf5[[glue::glue('/Batches/{batch}/Dacs')]][dacs_start:dacs_end])
              ),
      by = read_id
    ]
  logger::log_trace(glue::glue("\t\t\tDacs to current"))
  metainfo[
      , current := list(
          list(((unlist(dacs) + offset) * range) / digitisation)
        ),
      by = read_id
    ][
      , `:=`(
        dacs = NULL,
        offset = NULL,
        range = NULL,
        digitisation = NULL
        )
    ]
  logger::log_trace(glue::glue("\t\t\tNormalising current"))
  metainfo[
      , current_norm := list(
          list((unlist(current) - shift_frompA) / scale_frompA)
        ),
      by = read_id
    ][
      , `:=`(
        current = NULL,
        shift_frompA = NULL,
        scale_frompA = NULL
        )
    ]
  logger::log_trace(glue::glue("\t\t\tAdding reference position to signal"))
  metainfo[
      , current_index := list(
          list(add_index_to_vector(unlist(ref_to_signal), unlist(current_norm)))
        ),
      by = read_id
    ][
      , `:=`(ref_to_signal = NULL, current_norm = NULL)
    ]

    return(NULL)
}


#' Load signal mappings of one batch
#' 
#' Loads the signal mappings associated with a batch
#' 
#' @param signal data.table of metainfo loaded with read_metainfo
#' @return Nothing, signal mappings will be added to the metainfo object in memory
#' @export
get_reference_context <- function(signal) {
  signal[
      , unlist(current_index, recursive = FALSE),
      by = .(read_id, type, chunk, reference, pos, strand)
    ][
      ,
      pos_ref := ind + pos - 1
    ][
      pos_ref <= (chunk + 1) * chunk_size
    ][
      pos_ref >= (chunk) * chunk_size
    ]
}

