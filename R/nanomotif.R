#' Nanomotif
#'
#' Perform complete motif identification pipeline.
#'
#' @param nat_mapping Path to NAT sorted bam file
#' @param pcr_mapping Path to PCR sorted bam file
#' @param nat_signal Path to NAT signal_mappings.hdf5
#' @param pcr_signal Path to PCR signal_mappings.hdf5
#' @param out Output folder
#' @param reference_path Path to reference signal was mapped too
#' @param chunk_size Size of chunks (default 1e5)
#' @param max_cov Maximum coverage to include in motif detection (default: 200)
#' @param minPval Minimum p-value of current differences to considered event  (default: 1e-15)
#' @param hdbscan_args Arguments to HDBSCAN clustering (see documentation for possibilities)
#' @param umap_args Arguments to UMAP embedding (see documentation for possibilities)
#' @param align_event_sequences Align event sequence to increase possibility of detecting consensus sequence (default: TRUE)
#' @param iterations Number of clustering iterations (default 1)
#' @param chunks_per_contig Number of chunks to process pr. contig (deafault: 0 (All chunks))
#' @param min_entropy Minimum entropy of cluster alignment to consider a postion conserved (default: 1)
#' @return 
#' @import data.table
#' @export
nanomotif <- function(
  nat_mapping,
  pcr_mapping,
  nat_signal,
  pcr_signal,
  out,
  reference_path,
  chunk_size = 1e5,
  max_cov = 200,
  minPval = 1e-15,
  hdbscan_args = list(minPts = 10),
  umap_args = list(
      n_components = 2,
      min_dist = 0.05,
      n_neighbors = 30
    ),
  align_event_sequences = TRUE,
  iterations = 1,
  min_entropy = 1,
  chunks_per_contig = 0, # 0 = keep all
  select_contigs = "all",
  rolling_mean_k = 5
) {
  # Setting up logger
  logger::log_threshold(logger::TRACE)
  log_appender(appender_file(file.path(out, "log.txt")))

  # Load read metainfo
  metainfo <- prepare_metainfo(
    nat_mapping = nat_mapping,
    pcr_mapping = pcr_mapping,
    nat_hdf5 = nat_signal,
    pcr_hdf5 = pcr_signal,
    chunk_size = chunk_size,
    max_cov = max_cov
  )

  if (chunks_per_contig > 0) {
    # Select specified number of chunks per contig
    metainfo <- metainfo[
      , chunks_keep := (chunk+1) %in% sample(seq_len(max(chunk+1)), min(chunks_per_contig, max(chunk+1))), by = .(reference)
    ][
      chunks_keep == TRUE,
    ][
      , chunks_keep := NULL
    ]
  }
  if (select_contigs != "all") {
    # Select only some contigs
    metainfo <- metainfo[
      reference %in% select_contigs
    ]
  }
  metainfo_chunk <- split(
    metainfo,
    by = c("chunk_ref", "batch", "type"),
    flatten = FALSE
  )

  # Process chunked reads
  preprocess_all_chunks(
    metainfo = metainfo_chunk,
    nat_hdf5 = nat_signal,
    pcr_hdf5 = pcr_signal,
    out = out,
    chunk_size = chunk_size
  )
  find_motifs(
    path_chunk_stats = out,
    path_ref = reference_path,
    out = out,
    align_event_sequences = TRUE,
    hdbscan_args = hdbscan_args,
    umap_args = umap_args,
    event_sequence_frame = 7,
    event_pval_limit = minPval,
    entropy_threshold_motif = min_entropy,
    iteration_approach = "remove_noise",
    iterations = iterations
  )
  for (i in (seq_len(iterations))) {
    plot_motifs(
      cluster_path = file.path(out, "/", i, "/clusters.tsv"),
      plot_out  = file.path(out, "/", i),
      chunks_path = out,
      path_ref = reference_path
    )
  }
}


