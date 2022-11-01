
unnest_to_row_per_position = function(signal_mapping) {
  signal_mapping = signal_mapping[
      , list(unlist(current_batch, recursive = FALSE)), by = read_id
    ][
      , read_position := 1:.N, by = read_id
    ]
  signal_mapping %>% setnames("V1", "current")
}

unnest_to_row_per_current_value = function(signal_mapping) {
  signal_mapping[, list(unlist(V1, recursive = FALSE)), by = .(read_id, position)][
    , current_position := 1:.N, by = .(read_id, position)]
}

filter_signal_mapping = function(signal_mapping, read_id_keep) {
  signal_mapping[read_id %in% read_id_keep, ]
}





add_reference <- function(chunk, sequence) {
  chunk[
      strand == "+", base := toupper(sequence[pos_ref])
    ][
      strand == "-", base := get_complement_sequence(sequence[pos_ref])
    ]
}
