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
                       contigs = "all") {
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

  read_mapping <- preprocess_mappings(nat_mapping, pcr_mapping, chunk_size)

  logger::log_debug("Loading metainfo")
  metainfo <- rbind(
    read_metainfo_all(hdf5[["nat"]]),
    read_metainfo_all(hdf5[["pcr"]])
    )

  add_mapping(metainfo, read_mapping)
  rm(read_mapping)
  gc()

  logger::log_debug("Splitting metainfo")
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
      lapply(
        batches,
        function(batch) {
          logger::log_debug(glue::glue("    Processing {batch}"))
          types <- names(metainfo_list[[chunk]][[batch]])
          lapply(
            types,
            function(type) {
              add_signal(metainfo_list[[chunk]][[batch]][[type]], hdf5[[type]], batch)
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
  add_index(metainfo)
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
add_signal = function(metainfo, hdf5, hdf_batch) {
  logger::log_trace("Loading ref to signal")
  metainfo[
      , ref_to_signal := read_ref_to_signal(ref_to_signal_start, ref_to_signal_end, hdf5 = hdf5, batch = hdf_batch)
    ][
      , `:=`(ref_to_signal_start = NULL, ref_to_signal_end = NULL, Ref_to_signal_lengths = NULL)
    ]
  logger::log_trace("Loading dacs")
  metainfo[
      , dacs := read_dacs(dacs_start, dacs_end, hdf5 = hdf5, batch = hdf_batch)
    ][
      , `:=`(dacs_start = NULL, dacs_end = NULL, Dacs_lengths = NULL)
    ]
  logger::log_trace("Normalising dacs")
  calculate_current_from_dacs(metainfo)
  normalise_current(metainfo)
  logger::log_trace("Batching dacs")
  metainfo[
    , current := batch_current_by_reference(ref_to_signal, current_norm)
  ][
    , `:=`(ref_to_signal = NULL, current_norm = NULL)
  ]
  gc()
  return("Added signal")
}

#' Add read mappings to metainfo
#' 
#' Add necessary read mapping information to metainfo
#' 
#' @param metainfo data.table of metainfo loaded with read_metainfo
#' @param read_mapping read mappings
#' @return Nothing, signal mappings will be added to the metainfo object in memory
#' @export
add_mapping <- function(metainfo, read_mapping) {
  metainfo[
    read_mapping,
    on = c('read_id'='qname'),
    `:=`(chunk = i.chunk, reference = i.rname, type = i.type, pos = i.pos, strand = i.strand)
  ][
    , chunk_ref := paste0(chunk, "_", reference)
  ]
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
      pos_ref := pos_read + pos - 1
    ][
      pos_ref <= (chunk + 1) * chunk_size
    ][
      pos_ref >= (chunk) * chunk_size
    ]
}

add_index = function(signal_mapping) {
  signal_mapping[
      , dacs_end := cumsum(Dacs_lengths)
    ][
      , dacs_start := dacs_end - Dacs_lengths + 1
    ][
      , ref_to_signal_end := cumsum(Ref_to_signal_lengths)
    ][
      , ref_to_signal_start := ref_to_signal_end - Ref_to_signal_lengths + 1
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
  return(
    mapply(
      function(start, end) {
        dacs[start:end]
      },
      dacs_start,
      dacs_end,
      SIMPLIFY = FALSE
    )
  )
}

calculate_current_from_dacs = function(signal_mapping) {
  signal_mapping[, current := mapply(
    function(dacs, offset, range, digitisation) {
      ((dacs + offset) * range) / digitisation
    },
    dacs,
    offset,
    range,
    digitisation,
    SIMPLIFY = FALSE
  )]
  signal_mapping[, `:=`(
    dacs = NULL,
    offset = NULL,
    range = NULL,
    digitisation = NULL
  )]
}

normalise_current = function(signal_mapping) {
  signal_mapping[, current_norm := mapply(
    function(current, shift, scale) {
      (current - shift) / scale
    },
    current,
    shift_frompA,
    scale_frompA,
    SIMPLIFY = FALSE
  )]
  signal_mapping[, `:=`(
    current = NULL,
    shift_frompA = NULL,
    scale_frompA = NULL
  )]
}

batch_current_by_reference = function(ref_to_signal, current_norm) {
  current_norm_batch = mapply(
    function(start, current) {
      end = data.table::shift(start, type = "lead") %>%
        na.omit() %>%
        c() - 1
      start = start[1:(length(start) - 1)]
      list(
        read_pos = rep(1:length(start), end - start + 1),
        signal = current[min(start):max(end)]
        )
    },
    ref_to_signal,
    current_norm,
    SIMPLIFY = FALSE
  )
}
