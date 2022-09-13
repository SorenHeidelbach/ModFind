#' Unnests data.table
#'
#' This function unnest list column of data.table. 
#'
#' @param dt data.table containing listed column
#' @param list_col Column that contain list
#' @param keep_col Columns that are kept in the returned data.table
#' @return Unnested data.table
#' @import data.table
#' @export
unnest_dt = function(dt, list_col, keep_col) {
  stopifnot(data.table::is.data.table(dt))
  dt[, get(list_col), by = mget(keep_col)]
}

#' Read metainfo from signal_mappings.hdf
#' 
#' This function load metainformation of each read in a batch in a singal mapping hdf5 file.
#' 
#' @param hdf5 Open hdf5 object
#' @param batch Batch to read from
#' @return data.table with metainformation
#' @export
metainfo_load_batch = function(hdf5, batch) {
  data.table::data.table(
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
}

#' Read metainfo from signal_mappings.hdf
#' 
#' This function load all  metainformation in a singal mapping hdf5 file.
#' 
#' @param hdf5 Open hdf5 object
#' @return list of data.tables with metainformation of all batches
#' @export
metainfo_load_all <- function(hdf5) {
  batches <- get_hdf_batches(hdf5)
  lapply(
    batches,
    function(batch){
      signal_mapping = metainfo_load_batch(hdf5, batch)
      return(signal_mapping)
    }
  )
}

#' See which batches are present in signal_mappings.hdf
#' 
#' 
#' 
#' @param hdf5 Open hdf5 object
#' @return vector with batch names
#' @export
get_hdf_batches <- function(hdf5){
  rhdf5::h5ls(hdf5) %>%
    dplyr::pull(group) %>%
    stringr::str_extract("Batch_[\\d]*") %>%
    na.omit() %>%
    c() %>%
    unique()
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

#' Load signal mappings of one batch
#' 
#' Loads the signal mappings associated with a batch
#' 
#' 
#' @param metainfo data.table of metainfo loaded with metainfo_load_batch
#' @param hdf5 Open hdf5 object
#' @param hdf5 Open hdf5 object
#' @return Nothing, signal mappings will be added to the metainfo object in memory
#' @export
load_signal_mappings_batch = function(metainfo, hdf5, batch) {
  add_index(metainfo)
  metainfo[
      , ref_to_signal := read_ref_to_signal(ref_to_signal_start, ref_to_signal_end, hdf5 = hdf5, batch = unique(batch))
    ][
      , `:=`(ref_to_signal_start = NULL, ref_to_signal_end = NULL, Ref_to_signal_lengths = NULL)
    ]
  metainfo[
      , dacs := read_dacs(dacs_start, dacs_end, hdf5 = hdf5, batch = unique(batch))
    ][
      , `:=`(dacs_start = NULL, dacs_end = NULL, Dacs_lengths = NULL)
    ]
  calculate_current_from_dacs(metainfo)
  normalise_current(metainfo)
  metainfo[
    , current_batched := batch_current_by_reference(ref_to_signal, current_norm)
  ]
  gc()
  return("Loaded signal mappings")
}
