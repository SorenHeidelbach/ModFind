#' Load sorted read mappings
#'
#' Load read mapping required for correct assignment of signal mapping reference. 
#'
#' @param path_mapping Preprocessed signal mappings
#' @return read mapping data.table
#' @import data.table
#' @export
load_mapping = function(path_mapping) {
  if (!file.exists(glue::glue("{path_mapping}.bai"))) {
    Rsamtools::indexBam(path_mapping)
  }
  what = c("qname", "rname", "pos", "strand", "qwidth")
  param = Rsamtools::ScanBamParam(what = what)
  read_mapping =
    Rsamtools::scanBam(
      path_mapping,
      param = param,
      index = glue::glue("{path_mapping}.bai")
    ) %>%
    unlist(recursive = FALSE)
  setDT(read_mapping)
}



#' Downsample either nat or pcr to least abundant type
#'
#' Removes read of the most bundant mapping type, such there are almost equal mapped bases of each type to each chunk.
#' 
#'
#' @param read_mapping Preprocessed signal mappings
#' @param min_cov Minimum coverage of type
#' @param chunk_size size of chunks
#' @return read mapping data.table
#' @import data.table
#' @export
downsample <- function(read_mapping, chunk_size, min_cov = 20){
  min_chunk_cov <- min_cov * chunk_size
  read_mapping <- read_mapping[sample(1:nrow(read_mapping))]
  read_mapping[
      , read_chunk_cov := fcase(
          pos + qwidth > (chunk + 1) * chunk_size,  as.integer((chunk + 1) * chunk_size - pos),
          pos < chunk * chunk_size, as.integer(pos + qwidth - chunk * chunk_size),
          rep(TRUE, length(pos)), qwidth
        )
    ][
      , coverage := cumsum(as.numeric(read_chunk_cov)), by = .(type, chunk, strand)
    ][
      , max_type_coverage_by_chunk := max(coverage), by = .(type, chunk, strand)
    ][
      , max_allowed_chunk_cov := min(max_type_coverage_by_chunk), by = .(chunk, strand)
    ][
      coverage < max_allowed_chunk_cov | coverage < min_chunk_cov,
    ][
      , `:=`(max_type_coverage_by_chunk = NULL, max_allowed_chunk_cov = NULL)
    ]
}


unnest_dt = function(dt, list_col, keep_col) {
  stopifnot(data.table::is.data.table(dt))
  dt_unnested <- dt[, get(list_col), by = mget(keep_col)]
  setnames(dt_unnested, "V1", list_col)
  return(dt_unnested)
}

#' Get read mappings that overextent their chunk
#'
#' Gets the read that overextent their assigned chunk and return dt with new chunk assignments
#' 
#' 
#' @param read_mapping Read mapping dt
#' @param chunk_size Size of chunks
#' @return read mapping data.table
#' @import data.table
#' @export
get_overextending_reads <- function(read_mapping, chunk_size){
  cnames <- colnames(read_mapping)[!(colnames(read_mapping) == "chunk")]
  max_chunk <- max(read_mapping$chunk)
  corrected_read_mappings <- read_mapping[
      pos + qwidth > (1 + chunk) * chunk_size,
    ][
      , over_extension := pos + qwidth - (1 + chunk) * chunk_size
    ][
      , .(chunk = list(1:(1 + over_extension %/% chunk_size) + chunk)),
      by = mget(cnames)
    ][
      chunk <= max_chunk
    ]
  unnest_dt(corrected_read_mappings, "chunk", cnames)
}
