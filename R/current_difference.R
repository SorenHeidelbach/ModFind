
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
          chunk_size = chunk_size
        )
      rm(x)
      gc()
      chunk_stats[, chunk_ref := chunk_reference]
      chunk_out_path <- paste_path(out, "/chunks/", chunk_reference)
      dir.create(chunk_out_path, recursive = TRUE)
      fwrite(chunk_stats, paste_path(chunk_out_path, "chunk.tsv"))
      rm(chunk_stats)
    }
  )
}
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
    , mean_diff := mean(unlist(nat)) - mean(unlist(pcr), na.rm=TRUE),
    by = .(pos_ref, strand)
  ][
    , median_diff := median(unlist(nat)) - median(unlist(pcr), na.rm=TRUE),
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

#' Extract features from event statistics
#'
#' Selects events based on p-value and extract feature from signal dt 
#' with calculated statistics, which is then ready for embedding and clustering
#'
#' @param signal Preprocessed signal mappings
#' @param p_value_threshold Threshold for event selection (based on two sided wilcox test)
#' @return data.table with feature columns and reference position
#' @export
get_event_features <- function(signal, p_value_threshold = 1e-6) {
  signal[
    , event := min(p_val) < p_value_threshold, by = .(pos_ref)
    ]
  
  events <- gather_plus_and_minus_strand(signal)
  return(events)
}

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
  plot_path = paste0("./current_plots/", unique(chunk_list[[1]][[1]]$chunk_ref), "/")
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
  plot_chunk_current(chunk, plot_path)

  ## difference
  logger::log_debug("Calculaitng statistics")
  chunk_stats <- calculate_statistics(chunk)
  return(chunk_stats)
}


get_moving_cunks = function(vector, frame_size = 8) {
  # Return list of 2*window+1 vectors.
  # The list is transposed such the first value of each list belong to the 
  # same frame
  outer(
      -frame_size:frame_size,
      seq_along(vector),
      FUN = function(x, y) {
        z = x + y + 1
        fifelse(z < 0, NA_real_, z)
      }) %>%
    apply(MARGIN = 2, FUN = function(x) vector[x]) %>%
    data.table::transpose()
  }


gather_plus_and_minus_strand <- function(chunk) {
    chunk_merged = merge(
    chunk[
        event == TRUE & strand == "-",
      ][
        , .SD, .SDcols = names(chunk) %like% "diff|pos_ref|chunk_ref"
      ],
    chunk[
        event == TRUE & strand == "+",
      ][
        , .SD, .SDcols = names(chunk) %like% "diff|pos_ref|chunk_ref"
      ],
    by = c("pos_ref", "chunk_ref"),
    all = TRUE,
    suffixes = c("_minus", "_plus")
  )
  return(chunk_merged)
}

add_reference <- function(chunk, sequence) {
  chunk[
      strand == "+", base := toupper(sequence[pos_ref])
    ][
      strand == "-", base := get_complement_sequence(sequence[pos_ref])
    ]
}
