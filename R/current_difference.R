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
calculate_statistics <- function(signal,
                               window = 8,
                               value_col = "mean_diff"
                               ) {
    chunk_diff <- signal %>%
      gather_nat_and_pcr() %>%
      calculate_mean_dif() %>%
      calculate_median_dif() %>%
      add_frame_statistics(frame_size = window, value_col = value_col) %>%
      remove_when_signal_missing() %>%
      calculate_wilcox_two_sided()
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
get_event_features <- function(signal, p_value_threshold = 1e-10) {
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
process_chunk <- function(chunk_list, h5_list, chunk_size, plot_path = paste0("./current_plots/", unique(metainfo[[2]][[1]][[1]]$chunk_ref), "/")) {
  add_signal_chunk(chunk_list, h5_list)
  chunk <- rbindlist(unlist(chunk_list, recursive = FALSE))
  chunk <- get_reference_context(chunk, chunk_size)

  ## signal viz
  plot_chunk_current(chunk, plot_path)

  ## difference
  chunk_stats <- calculate_statistics(chunk)
  return(chunk_stats)
}

gather_nat_and_pcr <- function(chunk){
  dcast(
    chunk,
    strand + pos_ref ~ type,
    value.var = "signal", fun.agg = function(x) list(unlist(x))
  )
}

remove_when_signal_missing <- function(chunk){
  chunk[
    (lengths(nat) > 0) & (lengths(pcr) > 0)
  ]
}

calculate_mean_dif <- function(chunk){
  chunk[
    ,
    mean_diff := mean(unlist(nat)) - mean(unlist(pcr)),
    by = .(pos_ref, strand)
  ]
}

calculate_median_dif <- function(chunk){
  chunk[
    ,
    median_diff := median(unlist(nat)) - median(unlist(pcr)),
    by = .(pos_ref, strand)
  ]
}

calculate_wilcox_two_sided <- function(chunk){
  chunk[
    ,
    p_val := wilcox.test(unlist(nat), unlist(pcr), exact = FALSE)$p.value,
    by = .(pos_ref, strand)
  ]
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

add_frame_statistics <- function(chunk, frame_size, value_col) {
  chunk[
      ,
      paste0(value_col, "_dif_", -frame_size:frame_size) := get_moving_cunks(get(value_col), frame_size = frame_size),
      by = strand
    ]
}

gather_plus_and_minus_strand <- function(chunk) {
    chunk_merged = merge(
    chunk[
        event == TRUE & strand == "-",
      ][
        , .SD, .SDcols = names(chunk) %like% "diff_|pos_ref"
      ],
    chunk[
        event == TRUE & strand == "+",
      ][
        , .SD, .SDcols = names(chunk) %like% "diff_|pos_ref"
      ],
    by = "pos_ref",
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

