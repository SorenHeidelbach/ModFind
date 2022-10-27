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

    assert::assert(
        max(end) <= length(vec),
        min(start) >= 1,
        msg = "Increment_pos out of bounds"
    )
    mapply(
      function(start, end) {
        vec[start:end]
      },
      start,
      end
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


#' Paste string to path
#' 
#' Ensures that leading or trailing / are correctly handled
#' 
#' @param ... Strings to be concatenated
#' @return Character
#' @export
paste_path <- function(...) {
  x <- list(...)

  # Check inputs
  position_of_vector <- lengths(x) > 1
  vec <- x[position_of_vector]
  checkmate::assert_true(sum(position_of_vector) <= 1)
  
  # Repeat non multi character inputs
  x <- lapply(
    x[!position_of_vector],
    function(chr) rep(chr[1], max(lengths(x)))
  )
  x[!position_of_vector] <- x
  x[position_of_vector] <- vec
  x <- transpose(x)

  # Paste the strings
  x <- lapply(
    x,
    function(paths) {
    paths <- unlist(paths)
    path  <- ""
    for (i in seq_along(paths)) {
      path <- stringr::str_remove(path, "\\/$")
      string <- stringr::str_remove(paths[i], "^\\/")
      if (i > 1) {
        path <- paste(path, string, sep = "/")
      } else {
        path <- paste0(path, string)
      }
    }
    return(path)
  }
  )
  if (sum(position_of_vector) == 0) {
    x <- unlist(x)
  }
  return(x)
}

# Transfer attributes from one dt to another
transfer_attributes <- function(dt_from, dt_to, override = FALSE){
  if (!override) {
    old_attr <- attributes(dt_to) %>% 
      names()
  } else {
    old_attr <- data.table() %>% 
      attributes() %>% 
      names()
  }
  
  new_attr <- attributes(dt_from) %>% 
    `[`(!names(.) %in% old_attr)
  
  invisible(
    lapply(names(new_attr), function(x) setattr(dt_to, x, new_attr[[x]]))
  )
}