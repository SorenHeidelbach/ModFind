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
#' @return preprocessed chunk saved to 'out'
#' @export
preprocess <- function(nat_mapping,
                       pcr_mapping,
                       nat_signal,
                       pcr_signal,
                       out,
                       chunk_size = 1e4,
                       threads = 1,
                       debug = FALSE,
                       contigs = "all",
                       pvp = FALSE) {
    # Forking not possible on windows systems
  if (.Platform$OS.type == "windows") {
    threads <- 1
  }
  if (debug) {
    logger::log_threshold(logger::TRACE)
  }

  logger::log_debug("Loading read mapping")
  hdf5 <- list(
    pcr = pcr_signal,
    nat = nat_signal
  )

  pcr <- load_mapping(pcr_mapping)
  pcr[, type := "pcr"]
  nat <- load_mapping(nat_mapping)
  nat[, type := "nat"]
  read_mapping <- rbind(pcr, nat)
  read_mapping[
      , chunk := pos %/% chunk_size
    ]
  read_mapping <- rbind(
      read_mapping,
      get_overextending_reads(read_mapping, chunk_size)
  )

  if(pvp){
    read_mapping <- read_mapping[
      type == "pcr"
    ][
      sample(1:.N, .N/2), type := "nat"
    ]
  }

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

  metainfo_list <- split(metainfo, by = c("chunk_ref", "batch", "type"), flatten = FALSE)

  chunks <- names(metainfo_list)

  if (contigs != "all"){
    chunks <- chunks[contigs %>%
      lapply(function(contig) grepl(contig, chunks)) %>%
      transpose() %>%
      lapply(any) %>%
      unlist()]
    logger::log_info(paste0("Only processing chunks: ", paste0(chunks, collapse = ", ")))
  }

  lapply(
    chunks,
    function(chunk) {
      logger::log_debug(glue::glue("Processing {chunk}"))
      batches <- names(metainfo_list[[chunk]])
      successfully_read <- lapply(
        batches,
        function(batch) {
          logger::log_debug(glue::glue("\tProcessing {batch}"))
          types <- names(metainfo_list[[chunk]][[batch]])
          lapply(
            types,
            function(type) {
              logger::log_debug(glue::glue("\t\t{type}"))
              dacs <- hdf5[[type]][[glue::glue('/Batches/{batch}/Dacs')]][]
              ref_to_signal <-  hdf5[[type]][[glue::glue('/Batches/{batch}/Ref_to_signal')]][]
              tryCatch(
                add_signal(
                  metainfo_list[[chunk]][[batch]][[type]],
                  dacs,
                  ref_to_signal
                ),
                error = function(e) FALSE
              )
              rm("dacs", "ref_to_signal")
              gc()
            } # type
          )
        } # batch
      )
      signal <- rbindlist(unlist(metainfo_list[[chunk]], recursive = FALSE))
      metainfo_list[chunk] <- NULL
      signal <- get_reference_context(signal)

      # Current difference
      logger::log_debug(glue::glue("Calculating statistics"))
      signal_stats <- calculate_statistics(signal)


      fwrite(
          signal_stats,
          file = glue::glue("{out}/{chunk}.tsv"),
          sep = "\t")
      rm(signal); gc()
    } # chunk
  )
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
add_signal = function(metainfo, dacs_vec, ref_to_sig) {
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
      paste0(
        duplicated_read_ids,
        collapse = ", "
      )
    )
  )



  metainfo[
    # Add ref_to_signal
      , ref_to_signal := list(
          list(ref_to_sig[ref_to_signal_start:ref_to_signal_end])
        ),
      by = read_id
    ][
    # Add Dacs
      , dacs := .(
          list(dacs_vec[dacs_start:dacs_end])
        ),
      by = read_id
    ][
    # Convert Dacs to current
      , current := .(
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
    ][
    # Normalise current
      , current_norm := .(
          list((unlist(current) - shift_frompA) / scale_frompA)
        ),
      by = read_id
    ][
      , `:=`(
        current = NULL,
        shift_frompA = NULL,
        scale_frompA = NULL
        )
    ][
    # Index current with respective ref position
      , current := .(list(add_index_to_vector(unlist(ref_to_signal), unlist(current_norm)))),
      by = read_id
    ][
      , `:=`(ref_to_signal = NULL, current_norm = NULL)
    ]
  return(TRUE)
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
      , unlist(current, recursive = FALSE),
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


read_ref_to_signal = function(ref_to_signal_start, ref_to_signal_end, hdf5, batch) {
  ref_to_sig <- hdf5[[glue::glue('/Batches/{batch}/Ref_to_signal')]][]
  return(
    mapply(
      function(start, end) {
        ref_to_sig[start:end]
      },
      ref_to_signal_start,
      ref_to_signal_end,
      SIMPLIFY = FALSE
    )
  )
}

read_dacs = function(dacs_start, dacs_end, hdf5, batch) {
  dacs <- hdf5[[glue::glue('/Batches/{batch}/Dacs')]][]
  length_dacs <- length(dacs)
  return(
    mapply(
      function(start, end) {
        if (length_dacs < max(end)) {
          print(batch)
        }
        dacs[start:end]
      },
      dacs_start,
      dacs_end,
      SIMPLIFY = FALSE
    )
  )
}