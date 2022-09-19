#' Calculates statistics between NAT and PCR signal
#'
#' Complete processing from difference calculation to event selection
#'
#' @param signal Preprocessed signal mappings
#' @param window Number of position to in frame include at each side of event
#' @param value_col Statistic to include in event frame
#' @param fill_na Value to fill NA values with
#' @return current difference data.table
#' @import data.table
#' @export
calculate_statistics <- function(signal,
                               window = 8,
                               value_col = "mean_diff",
                               fill_na = 0
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

#' Calculates current difference from preprocess output
#'
#' Complete processing from difference calculation to event selection
#'
#' @param signal Preprocessed signal mappings
#' @param p_value_threshold Threshold for event selection (based on two sided wilcox test)
#' @return current difference data.table
#' @import data.table
#' @export
get_event_features <- function(signal, p_value_threshold = 1e-10) {
  events <- signal %>%
    identify_events(p_value_threshold) %>%
    gather_plus_and_minus_strand()
    return(events)
}

gather_nat_and_pcr <- function(chunk){
  dcast(
    chunk,
    strand + pos_ref ~ type,
    value.var = "signal", fun.agg = function(x) list(x)
  )
}

remove_when_signal_missing <- function(chunk){
  chunk[
    (lengths(nat) > 0) & (lengths(pcr) > 0)
  ]
}

identify_events <- function(chunk, p_value_threshold) {
  chunk[
    , event := min(p_val) < p_value_threshold, by = .(pos_ref)
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
        , .SD, .SDcols = names(chunk) %like% "dif_|pos_ref"
      ],
    chunk[
        event == TRUE & strand == "+",
      ][
        , .SD, .SDcols = names(chunk) %like% "dif_|pos_ref"
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

