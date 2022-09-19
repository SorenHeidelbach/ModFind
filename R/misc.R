#' Return list with vector and indices
#' 
#' 
#' 
#' @param increment_pos Positions where index increment by one, the range of the returned vectors is min:max of this vector
#' @param vec Vector to add indices too
#' @return data.table with read mappings
#' @export
add_index_to_vector  <- function(increment_pos, vec) {
    assert::assert(
        is.numeric(increment_pos),
        msg = "Ensure increment_pos is numeric"
    )
    assert::assert(
        max(increment_pos) <= length(vec),
        msg = "Increment_pos out of bounds"
    )
    assert::assert(
        length(increment_pos) >= 2,
        msg = "increment_pos should contain a minimum of 2 values (start, end)"
    )

    end <- increment_pos[-1] - 1
    start <- head(increment_pos, -1)
    repeats <- end - start + 1

    list(
      ind = rep(seq_along(start), repeats),
      vec = vec[min(start):max(end)]
    )
}
