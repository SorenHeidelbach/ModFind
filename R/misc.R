#' Return list with vector and indices
#' 
#' Calculates index of one vector based on increment vector and returns list with both
#' 
#' @param increment_pos Vector of positions of index increment
#' @param vec Vector to add indices too
#' @return data.table with read mappings
#' @export
add_index_to_vector  <- function(increment_pos, vec) {
    assert::assert(
        is.numeric(increment_pos),
        msg = "Ensure increment_pos is numeric"
    )
    
    assert::assert(
        length(increment_pos) >= 2,
        msg = "increment_pos should contain a minimum of 2 values (start, end)"
    )

    end <- increment_pos[-1] - 1
    start <- head(increment_pos, -1)
    repeats <- end - start + 1

    assert::assert(
        max(end) <= length(vec),
        min(start) >= 1,
        msg = "Increment_pos out of bounds"
    )
    list(
      ind = rep(seq_along(start), repeats),
      vec = vec[min(start):max(end)]
    )
}


#' Replace all NA's in DT
#' 
#' Replaces all the NA's in all columns in data.table with the specified value 
#' 
#' @param DT data.table
#' @param replacement value to replace with
#' @return data.table with filled NA's
#' @export
replace_na_dt = function(DT, replacement = 0L) {
  for (i in seq_len(ncol(DT))) {
    set(DT, which(is.na(DT[[i]])), i, replacement)
  }
}

invert_range <- function(vec) {
  2*(median(range(vec))) - vec
}