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
                       threads = 1) {
  if (.Platform$OS.type == "windows") {
    threads <- 1
  }
  check_hdf5(nat_signal)
  check_hdf5(pcr_signal)

  hdf5 <- list(
    pcr = pcr_signal,
    nat = nat_signal
  )
  print("Loading read mapping")
  read_mapping <- preprocess_mappings(nat_mapping, pcr_mapping, chunk_size)

  print("Loading metainfo")
  metainfo <- rbind(
    read_metainfo_all(hdf5[["nat"]]),
    read_metainfo_all(hdf5[["pcr"]])
    )
  add_mapping(metainfo, read_mapping)

  print("Splitting metainfo")
  metainfo_list <- split(metainfo, by = c("chunk_ref", "batch", "type"), flatten = FALSE)

  print("Loading singal")
  chunks <- names(metainfo_list)
  parallel::mclapply(
    mc.cores = threads,
    chunks,
    function(chunk) {
      print(chunk)
      batches <- names(metainfo_list[[chunk]])
      lapply(
        batches,
        function(batch) {
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
      fwrite(
        signal,
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
  check_hdf5(hdf5)
  metainfo <- data.table::data.table(
    Ref_to_signal_lengths = rhdf5::h5read(hdf5, glue::glue(
      "Batches/{batch}/Ref_to_signal_lengths"
    )),
    Dacs_lengths = rhdf5::h5read(hdf5, glue::glue("Batches/{batch}/Dacs_lengths")),
    digitisation = rhdf5::h5read(hdf5, glue::glue("Batches/{batch}/digitisation")),
    offset = rhdf5::h5read(hdf5, glue::glue("Batches/{batch}/offset")),
    range = rhdf5::h5read(hdf5, glue::glue("Batches/{batch}/range")),
    scale_frompA = rhdf5::h5read(hdf5, glue::glue("Batches/{batch}/scale_frompA")),
    shift_frompA = rhdf5::h5read(hdf5, glue::glue("Batches/{batch}/shift_frompA")),
    read_id = rhdf5::h5read(hdf5, glue::glue("Batches/{batch}/read_id"))
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
  rhdf5::h5ls(hdf5) %>%
    dplyr::pull(group) %>%
    stringr::str_extract("Batch_[\\d]*") %>%
    na.omit() %>%
    c() %>%
    unique()
}

#' Add signal mapping to metainfo
#' 
#' Loads the signal mappings associated with a batch
#' 
#' 
#' @param metainfo data.table of metainfo loaded with read_metainfo
#' @param hdf5 Open hdf5 object
#' @param hdf_batch Batch to load signal mappings from
#' @return Nothing, signal mappings will be added to the metainfo object in memory
#' @export
add_signal = function(metainfo, hdf5, hdf_batch) {
  metainfo[
      , ref_to_signal := read_ref_to_signal(ref_to_signal_start, ref_to_signal_end, hdf5 = hdf5, batch = hdf_batch)
    ][
      , `:=`(ref_to_signal_start = NULL, ref_to_signal_end = NULL, Ref_to_signal_lengths = NULL)
    ]
  metainfo[
      , dacs := read_dacs(dacs_start, dacs_end, hdf5 = hdf5, batch = hdf_batch)
    ][
      , `:=`(dacs_start = NULL, dacs_end = NULL, Dacs_lengths = NULL)
    ]
  calculate_current_from_dacs(metainfo)
  normalise_current(metainfo)
  metainfo[
    , current_batched := batch_current_by_reference(ref_to_signal, current_norm)
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
#' @param mapping Read mappings dt
#' @return Nothing, signal mappings will be added to the metainfo object in memory
#' @export
get_reference_context <- function(signal) {
  signal[
      , list(unlist(current_batched, recursive = FALSE)),
      by = .(read_id, type, chunk, reference, pos, strand)
    ][
      , pos_read := 1:.N, by = .(read_id, type, chunk)
    ][
      ,
      pos_ref := pos_read + pos - 1
    ][
      pos_ref <= (chunk + 1) * chunk_size
    ][
      pos_ref >= (chunk) * chunk_size
    ][
      ,
      unlist(V1),
      by = .(read_id, pos_ref, type, chunk, strand)
    ][
      ,
      pos_current := 1:.N,
      by = .(read_id, pos_ref)
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
  mapply(
    function(start, end) {
      rhdf5::h5read(
        file = hdf5,
        name = glue::glue(
          '/Batches/{batch}/Ref_to_signal'
        ),
        index = list(start:end)
      )
    },
    ref_to_signal_start,
    ref_to_signal_end, 
    SIMPLIFY = FALSE
  )
}

read_dacs = function(dacs_start, dacs_end, hdf5, batch) {
  mapply(
    function(start, end) {
      rhdf5::h5read(
        file = hdf5,
        name = glue::glue('/Batches/{batch}/Dacs'),
        index = list(start:end)
      )
    },
    dacs_start,
    dacs_end,
    SIMPLIFY = FALSE
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
      mapply(
        function(i, j) {
          current[i:j]
        },
        start,
        end,
        SIMPLIFY = FALSE
      )
    },
    ref_to_signal,
    current_norm,
    SIMPLIFY = FALSE
  )
}

check_hdf5 <- function(hdf5_obj) {
  assert::assert(
    class(hdf5_obj) == "H5IdComponent",
    msg = "hdf5 arguments is not a proper HDF5 object, please load with rhdf5::H5Fopen(hdf5_path)"
  )
}
