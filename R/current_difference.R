#' Processes all chunks in metainfo list
#'
#' Loads signal and calculates current difference
#'
#' @param chunks_list list of batches with nat and pcr metainfo
#' @param h5_list list of nat and pcr h5 objects
#' @param out output folder of preprocessed chunks
#' @param chunk_size size of chunks
#' @param plot_path path to output generated plot
#' @return dt with calculated statistics
#' @import data.table
#' @export
preprocess_all_chunks <- function(
  metainfo,
  nat_hdf5,
  pcr_hdf5,
  out,
  chunk_size = 1e5
) {
  future.apply::future_lapply(
    metainfo,
    function(x) {
      # For future, worker nodes to copy in memory to acces in nested functions
      x <- copy(x)
      chunk_reference <- unique(x[[1]][[1]]$chunk_ref)
      logger::log_debug(chunk_reference)
      chunk_stats <- process_chunk(
          x,
          nat_hdf5 = nat_hdf5,
          pcr_hdf5 = pcr_hdf5,
          chunk_size = chunk_size,
          plot_path = out
        )
      rm(x)
      chunk_stats[, chunk_ref := chunk_reference]
      chunk_out_path <- paste_path(out, "/chunks/", chunk_reference)
      dir.create(chunk_out_path, recursive = TRUE)
      fwrite(chunk_stats, paste_path(chunk_out_path, "chunk.tsv"))
      rm(chunk_stats)
    }
  )
}

################################################################################

#' Process chunk
#'
#' Loads signal and calculates current difference
#'
#' @param chunk_list list of batches with nat and pcr metainfo
#' @param h5_list list of nat and pcr h5 objects
#' @param chunk_size size of chunks
#' @param plot_path path to output generated plot
#' @return dt with calculated statistics
#' @import data.table
#' @export
process_chunk <- function(
  chunk_list,
  nat_hdf5,
  pcr_hdf5,
  chunk_size,
  plot_path
) {
  hdf5_list <- list(
    nat = hdf5r::H5File$new(nat_hdf5, mode = "r"),
    pcr = hdf5r::H5File$new(pcr_hdf5, mode = "r")
  )
  add_signal_chunk(chunk_list, hdf5_list = hdf5_list)
  logger::log_debug("Reformating chunk dt")
  chunk <- rbindlist(unlist(chunk_list, recursive = FALSE))
  chunk <- get_reference_context(chunk, chunk_size)

  ## signal viz
  logger::log_debug("Plotting signal subset")
  plot_chunk_current(chunk, paste0(plot_path, "/current_plots/", unique(chunk_list[[1]][[1]]$chunk_ref), "/"))

  ## difference
  logger::log_debug("Calculaitng statistics")
  chunk_stats <- calculate_statistics(chunk)
  return(chunk_stats)
}

#' Adds signal to a complete chunk
#' Iterator for adding signal to all bacthes and type (nat and pcr) of a
#' chunk of the chunk list outputted by prepare_metainfo
#' @param metainfo_chunk a list of batches, with NAt and PCR metainfo
#' @param hdf5 a list of nat and pcr hdf5 objects
#' @return logical of success status (signal is loaded in memory)
#' @export
add_signal_chunk <- function(metainfo_chunk, hdf5_list) {
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
          hdf5_list[[type]],
          batch
        ),
        error = function(e) FALSE
      )
    } # type
    )
  } # batch
  )
}

################################################################################

#' Add signal mapping to metainfo
#'
#' Loads the signal mappings associated with a batch
#'
#' @param metainfo data.table of metainfo loaded with load_metainfo
#' @param hdf5_obj Open hdf5 object
#' @param batch Batch to load signal mappings from
#' @return Nothing, signal mappings will be added to the metainfo object in memory
#' @export
add_signal <- function(metainfo, hdf5_obj, batch) {
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
  # +1 as R, compared to python, indexes from 1
  metainfo[
      , ref_to_signal := list(
              list(hdf5_obj[[glue::glue("/Batches/{batch}/Ref_to_signal")]][ref_to_signal_start:ref_to_signal_end] + 1)
              ),
      by = read_id
    ]
  logger::log_trace(glue::glue("\t\t\tAdding dacs"))
  metainfo[
      , dacs := list(
              list(hdf5_obj[[glue::glue("/Batches/{batch}/Dacs")]][dacs_start:dacs_end])
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
      , signal := list(
          list(add_index_to_vector(unlist(ref_to_signal), unlist(current_norm)))
        ),
      by = read_id
    ][
      , `:=`(current_norm = NULL)
    ]

}

################################################################################

#' Load signal mappings of one batch
#'
#' Loads the signal mappings associated with a batch
#'
#' @param signal data.table of metainfo loaded with load_metainfo
#' @param chunk_size integer, size of chunk
#' @return data.table
#' @export
get_reference_context <- function(signal_dt, chunk_size) {
  signal_unlisted <- signal_dt[
      , list(unlist(signal, recursive = FALSE)), by = .(read_id, type, chunk, reference, pos, strand)
    ]
  setnames(signal_unlisted, "V1", "signal")
  signal_unlisted  <- signal_unlisted[
      strand == "+", pos_read := seq_len(.N), by = .(read_id, chunk)
    ][
      strand == "-", pos_read := rev(seq_len(.N)), by = .(read_id, chunk)
    ][
      , pos_ref := pos + pos_read - 1
    ][
      pos_ref <= (chunk + 1) * chunk_size
    ][
      pos_ref >= (chunk) * chunk_size
    ]
  return(signal_unlisted)
}

################################################################################

#' Calculates statistics between NAT and PCR signal
#'
#' Calculates statistics between NAT and PCR signal
#'
#' @param signal Preprocessed signal mappings
#' @param window Number of position to in frame include at each side of event
#' @param value_col Statistic to include in event frame
#' @return current difference data.table
#' @import data.table
#' @export
calculate_statistics <- function(
  signal,
  frame_size = 8,
  value_col = "mean_diff"
) {
  chunk_diff <- dcast(
    signal,
    strand + pos_ref ~ type,
    value.var = "signal", fun.agg = function(x) list(unlist(x))
  )
  chunk_diff[
    , mean_diff := mean(unlist(nat)) - mean(unlist(pcr), na.rm = TRUE),
    by = .(pos_ref, strand)
  ][
    , median_diff := median(unlist(nat)) - median(unlist(pcr), na.rm = TRUE),
    by = .(pos_ref, strand)
  ][
    , paste0(value_col, -frame_size:frame_size) := get_moving_cunks(get(value_col), frame_size = frame_size),
    by = strand
  ]
  chunk_diff <- chunk_diff[
    (lengths(nat) > 1) & (lengths(pcr) > 1)
  ][
    , p_val := wilcox.test(unlist(nat), unlist(pcr), exact = FALSE)$p.value,
    by = .(pos_ref, strand)
  ]
  chunk_diff[
    , `:=`(nat = NULL, pcr = NULL)
  ]
  return(chunk_diff)
}

################################################################################

get_moving_cunks <- function(vector, frame_size = 8) {
  # Return list of 2*window+1 vectors.
  # The list is transposed such the first value of each list belong to the same
  # frame
  outer(
      -frame_size:frame_size,
      seq_along(vector),
      FUN = function(x, y) {
        z <- x + y + 1
        fifelse(z < 0, NA_real_, z)
      }) %>%
    apply(MARGIN = 2, FUN = function(x) vector[x]) %>%
    data.table::transpose()
  }
