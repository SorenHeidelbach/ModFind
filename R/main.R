#' Run preprocess 
#' 
#' Function to process signal and read mappings into preprocessed .tsv files, ready for current difference calculation
#' 
#' @param nat_mapping Path to NAT sorted bam file
#' @param pcr_mapping Path to PCR sorted bam file
#' @param nat_signal HDF5 object of NAT signal mappings
#' @param pcr_signal HDF5 object of PCR signal mappings
#' @param out Output path for preprocessed signal mappings
#' @param chunk_size Size to chunk reference into (default 1e4 bases)
#' @param threads Number of jobs to run parallel (default 1; >1 not supported on windows)
#' @return preprocessed chunk saved to 'out'
#' @export
process_chunks <- function(metainfo_list) {

  chunks <- names(metainfo_list)

  if (contigs != "all"){
    chunks <- chunks[contigs %>%
      lapply(function(contig) grepl(contig, chunks)) %>%
      transpose() %>%
      lapply(any) %>%
      unlist()]
    logger::log_info(paste0("Only processing chunks: ", paste0(chunks, collapse = ", ")))
  }

  lapply(
    chunks,
    function(chunk) {
      logger::log_debug(glue::glue("Processing {chunk}"))
      batches <- names(metainfo_list[[chunk]])
      successfully_read <- lapply(
        batches,
        function(batch) {
          logger::log_debug(glue::glue("\tProcessing {batch}"))
          types <- names(metainfo_list[[chunk]][[batch]])
          lapply(
            types,
            function(type) {
              logger::log_debug(glue::glue("\t\t{type}"))
              dacs <- hdf5[[type]][[glue::glue('/Batches/{batch}/Dacs')]][]
              ref_to_signal <-  hdf5[[type]][[glue::glue('/Batches/{batch}/Ref_to_signal')]][]
              tryCatch(
                add_signal(
                  metainfo_list[[chunk]][[batch]][[type]],
                  dacs,
                  ref_to_signal
                ),
                error = function(e) FALSE
              )
              rm("dacs", "ref_to_signal")
              gc()
            } # type
          )
        } # batch
      )
      signal <- rbindlist(unlist(metainfo_list[[chunk]], recursive = FALSE))
      metainfo_list[chunk] <- NULL
      signal <- get_reference_context(signal)

      # Current difference
      logger::log_debug(glue::glue("Calculating statistics"))
      signal_stats <- calculate_statistics(signal)


      fwrite(
          signal_stats,
          file = glue::glue("{out}/{chunk}.tsv"),
          sep = "\t")
      rm(signal); gc()
    } # chunk
  )
}