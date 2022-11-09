
#' Embed joined chunks
#'
#' Wrapper function for umbedding and clustering. Take the data.table containing
#' the output the preprocess_chunk function. Precessed chunks can be loaded using
#' load_processed_chunks
#'
#' @param chunks preprocessed chunks generated with preprocess_all_chunks
#' @param ref_path path to the reference passed to megalodon
#' @param out character string with path to output results
#' @param event_pval_limit Numeric, max p-value for events
#' @param event_sequence_frame integer, number of nucleotide to include to each side
#' @param umap_args list of arguments passed directly to umap (See ?umap::umap)
#' https://cran.r-project.org/web/packages/umap/vignettes/umap.html#configuration-objects
#' @param hdbscan_Args list of arguments passed directly to hdbscan (See ?dbscan::hdbscan)
#' @return embedded and clustered positions
#' @export
embed_and_cluster <- function(
    chunks,
    ref_path,
    out,
    event_pval_limit = 1e-30,
    event_sequence_frame = 5,
    umap_args = list(),
    hdbscan_args = list()
  ) {
  checkmate::assert_data_table(chunks)

  # Prepare features for embedding
  features <- get_event_features(
    chunks,
    event_pval_limit
  )
  logger::log_debug(glue::glue("Events: {nrow(features)}"))
  if (nrow(features) < 50) {
    stop("Number of events too low")
  }
  features_metainfo <- features[
    , .SD, .SDcols = !(grepl("diff", names(features)))
  ]

  # Add reference sequence to events
  reference <- seqinr::read.fasta(ref_path)
  add_event_frame_sequence(
    features_metainfo,
    reference,
    width = event_sequence_frame
  )

  # Remove non feature values
  features[
    , `:=`(pos_ref = NULL, chunk_ref = NULL)
  ]

  replace_na_dt(features, 0)

  # Embedding
  logger::log_debug("Embedding")
  embedding <- do.call(embed, list(features, umap_args = umap_args))

  # Clustering
  logger::log_debug(glue::glue("Clustering"))
  do.call(cluster_embedding, list(embedding, hdbscan_args = hdbscan_args))
  embedding <- cbind(
    embedding,
    features_metainfo
  )

  # Save results
  fwrite(embedding, paste_path(out, "embedding.tsv"))
  return(embedding)
}

################################################################################

get_event_features <- function(signal, event_pval_limit = 1e-6) {
  signal[
    , event := min(p_val) < event_pval_limit, by = .(pos_ref)
    ]
  events <- gather_plus_and_minus_strand(signal)
  return(events)
}

################################################################################

gather_plus_and_minus_strand <- function(chunk) {
    chunk_merged <- merge(
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

################################################################################

load_processed_chunks <- function(
  preprocess_out
) {
  chunks <- rbindlist(
    lapply(
    paste_path(list.files(paste_path(preprocess_out, "/chunks"), full.names = TRUE, include.dirs = TRUE), "chunk.tsv"),
    function(name) {
      fread(name)
    })
  )
  replace_na_dt(chunks, 0)
  return(chunks)
}

################################################################################

embed <- function(
  features,
  umap_args = list()
) {
  checkmate::assert_data_table(
    features,
    types = c("numeric"),
    any.missing = FALSE
  )
  # Add data to UMAP arguments
  umap_args$d <- features
  embedding <- do.call(umap::umap, umap_args)
  dt <- data.table(embedding$layout)
  setnames(
    dt,
    old = names(dt),
    new = paste0("UMAP", seq_len(ncol(dt)))
  )
  return(dt)
}

################################################################################

cluster_embedding <- function(
  embedding,
  hdbscan_args = list()
) {
  # Argument check
  checkmate::assert_data_table(
    embedding,
    types = "numeric",
    any.missing = FALSE
  )
  # add embedding to cluster arguments
  hdbscan_args$x <- embedding
  clustering <- do.call(dbscan::hdbscan, hdbscan_args)
  logger::log_debug(paste0("Clusters found:", length(unique(clustering$cluster)) - 1))
  embedding[
      , cluster := clustering$cluster
    ][
      , membership_prob := clustering$membership_prob
    ]
  embedding[
      , cluster := fifelse(cluster == 0, NA_integer_, cluster)
    ][
      , cluster := as.factor(cluster)
    ]
}

################################################################################

add_event_frame_sequence <- function(
  feature,
  reference,
  width = 3L
) {
  checkmate::assert_subset(
    c("chunk_ref", "pos_ref"),
    names(feature)
  )
  checkmate::assert_int(width)
  checkmate::assert_list(reference)

  feature[
      , ref := stringr::str_remove(chunk_ref, "^\\d+_")
    ][
      , seq := paste(reference[[ref]][(pos_ref - width):(pos_ref + width)], collapse = ""), by = .(pos_ref, chunk_ref)
    ]
}

################################################################################

calculate_bit_score <- function(
  sequences,
  nucleotides = c("G", "C", "A", "T"),
  min_entropy = 0.4
) {
    # Argument checks
    checkmate::assert_vector(
        sequences,
        strict = TRUE,
        min.len = 2,
        any.missing = FALSE
    )
    checkmate::assert_character(sequences)
    checkmate::assert_true(min(nchar(sequences)) == max(nchar(sequences)))

    # Constants
    uniform_entropy <- 4 * -0.25 * log2(0.25)
    n_sequences <- length(sequences)
    sequence_length <- nchar(sequences)[1]
    small_n_correction <- function(n) {
      (1/log(2)) * 3 / (2 * n)
    }

    # Nucleotides count at each position
    count_positional_character  <- function(seq, character_regex) {
        sequence_length <- nchar(seq)[1]
        lapply(1:sequence_length,
            function(position) {
                paste0("^", paste(rep(".", position - 1), collapse = ""), character_regex) %>%
                    grepl(seq) %>%
                    sum()
            }
        ) %>%
        unlist()
    }
    # Identify number of positional matches of nucleotides
    dt <- data.table(
        do.call(
            cbind,
            lapply(
                paste0("[", nucleotides, "]"),
                function(x) count_positional_character(sequences, x)
            )
        )
    )

    dt[
      , position := seq_len(.N)
    ]
    setnames(
        dt,
        old = paste0("V", seq_along(nucleotides)),
        new = nucleotides
    )
    dt <- melt(
      dt,
      measure.vars = paste0(nucleotides)
    )
    setnames(
        dt,
        old = c("variable", "value"),
        new = c("nucleotide", "n")
    )

    # Positionwise frequency of nucleotides
    dt[
        , freq := n / n_sequences, by = position
      ][
        , entropy_loss := - freq * log2(freq)
      ]
    replace_na_dt(dt)

    # Entropy at each position
    dt[
        , entropy := uniform_entropy - sum(entropy_loss) - small_n_correction(n_sequences)
        , by = position
      ][
        , height := entropy * freq
      ]
    motif <- dt[
        , .({
          if (all(entropy > min_entropy)) {
            ind <- order(freq)
            nucleotide[ind][cumsum(freq[ind]) > 0.5]
          } else {
            factor(".")
          }
        }, {
          if (all(entropy > min_entropy)) {
            ind <- order(freq)
            freq[ind][cumsum(freq[ind]) > 0.5]
          } else {
            NA_real_
          }
        }),  by = position
      ]
    setnames(
      motif,
      old = c("V1", "V2"),
      new = c("nucleotide", "freq")
    )
    setorder(motif, position, freq)
    consensus <- paste0(
      motif[
        , .(paste(nucleotide, collapse = "|"))
        , by = position
      ]$V1,
      collapse = "  ")

    return(
      list(
        dt = dt,
        motif = motif,
        motif_consensus = consensus
      )
    )
}

################################################################################

visualise_clusters <- function(
  dt
) {
  # Check inputs
  checkmate::assert_data_table(dt)
  checkmate::assert_subset(
    c("cluster", "UMAP1", "UMAP2", "seq"),
    names(dt)
  )

  # Get list of event sequences by cluster
  cluster_sequences <- lapply(
    split(dt[!is.na(cluster)], by = "cluster"),
    function(x) {
      toupper(x$seq)
    }
  )
  # Meaningfull name
  names(cluster_sequences) <- paste0(
    "Cluster ", names(cluster_sequences), ", n = ",
    lengths(cluster_sequences)
    )
  plot_logo <- ggseqlogo::ggseqlogo(
    cluster_sequences,
    seq_type = "dna") +
    default_theme_SH() +
    theme(
      panel.grid.major.x = element_blank()
    )

  plot_scatter <- ggplot(dt) +
    aes(x = UMAP1, y = UMAP2, color = cluster) +
    geom_point(alpha = 0.3, size = 0.3) +
    default_theme_SH() +
    ggplot2::scale_colour_viridis_d(na.value = "#cccccc88")

  plot_scatter_motif <- ggplot(dt) +
    aes(x = UMAP1, y = UMAP2, color = motif) +
    geom_point(alpha = 0.3, size = 0.3) +
    default_theme_SH() +
    ggplot2::scale_colour_viridis_d(na.value = "#cccccc88")

  return(
    list(
      plot_scatter,
      plot_logo,
      plot_scatter_motif
    )
  )
}

################################################################################

#' Visualise motif profiles
#'
#' Visualise the identified clusters
#'
#' @param cluster_path path to motif_find clsuter.tsv
#' @param plot_out path to plot output
#' @param chunks_path path to preprocess chunks
#' @param n_extra_positions number of adjacent positions to include in plot
#' @param motifs_evaluated evaluation type (all, clustered, or vector of motifs)
#' @return plots
#' @export
plot_motifs  <- function(
  cluster_path,
  plot_out,
  chunks_path,
  path_ref,
  n_extra_positions = 7,
  motifs_evaluated = "all"
) {
  chunks <- load_processed_chunks(chunks_path)
  clusters <- fread(cluster_path)
  motifs <- unique(clusters$motif)
  motifs <- motifs[nchar(motifs) > 1]
  ref  <- seqinr::read.fasta(path_ref, as.string = TRUE)

  # Clean motif consensus
  clusters[
    , motif := motif %>%
        stringr::str_remove_all("^[. ]*") %>%
        stringr::str_remove_all("[. ]*$") %>%
        stringr::str_remove_all("[ ]")  %>%
        stringr::str_remove_all("^.\\.*.$") %>%
        stringr::str_remove_all("^.$")
  ][
    , cluster := as.factor(cluster)
  ][
    , small_cluster := .N < nrow(clusters) * 0.001, by = motif
  ][
    small_cluster == TRUE, motif := "", by = motif
  ]

  # Visualise embedd + clustering
  viz <- visualise_clusters(clusters)
  ggplot2::ggsave(
    file.path(plot_out, "cluster_scatter.png"),
    width = 6,
    height = 5,
    viz[[1]]
  )
  ggplot2::ggsave(
    file.path(plot_out, "cluster_logo.png"),
    width = 20,
    height = 16,
    viz[[2]]
  )
  ggplot2::ggsave(
    file.path(plot_out, "cluster_scatter_motif.png"),
    width = 6,
    height = 5,
    viz[[3]]
  )

  motif_position <- if ("clustered" %in% motifs_evaluated) {
    # Get only clustered motif positions
    clusters[
        nchar(motif) > 1
      ][
        , start := stringr::str_locate(toupper(seq), motif)[[1]], by = .(pos_ref, chunk_ref)
      ][
        , start := start - (nchar(clusters$seq[[1]]) - 1) / 2 - 1 + pos_ref
      ]
  } else if ("all" %in% motifs_evaluated) {
    # Get all motif positions in reference
    motif_position <- lapply(
      motifs,
      function(motif) {
        as.data.table(
          stringr::str_locate_all(toupper(ref), motif)
        )[
          , motif := motif
        ]
      }
    ) %>% rbindlist()
  } else {
    motif_position <- lapply(
      motifs_evaluated,
      function(motif) {
        as.data.table(
          stringr::str_locate_all(toupper(ref), motif)
        )[
          , motif := motif
        ]
      }
    ) %>% rbindlist()
  }

  # Expand relative to motif start
  motif_position <- motif_position[
      , .(N = c(
          rep(".", n_extra_positions),
          unlist(strsplit(motif, "")),
          rep(".", n_extra_positions)
        ))
      , by = .(start, motif)
    ][
      , rel_pos := seq_len(.N) - n_extra_positions - 1, by = .(start, motif)
    ][
      , pos_ref := start + rel_pos
    ][
      , pos_label := paste0(rel_pos, "\n", N)
    ]

  # Add statistics and strand direction
  motif_position <- merge(
      motif_position,
      chunks[, .SD, .SDcols = c("pos_ref", "strand", "mean_diff", "p_val")],
      all = TRUE
    )[
      !is.na(motif) & !is.na(mean_diff)
    ]

  motif_position <- split(motif_position, by = "motif")
  aggregated_motif_boxplot  <- function(motif_position_filt) {
    ggplot(motif_position_filt) +
      aes(y = mean_diff, x = reorder(pos_label, rel_pos)) +
      geom_jitter(size = 0.3, alpha = 0.1) +
      geom_boxplot(outlier.alpha = 0) +
      facet_wrap(~strand) +
      labs(
        x = "Relative motif position (from motif start)",
        y = "Mean difference",
        title = "Aggregation of mean difference for all motifs"
      ) +
      default_theme_SH()
  }

  lapply(
    names(motif_position),
    function(motif) {
      # Box plot of all motifs
      dir.create(file.path(plot_out, motif))
      p_boxplot  <- aggregated_motif_boxplot(motif_position[[motif]])
      ggsave(
        paste0(file.path(plot_out, motif, "boxplot"), ".pdf"),
        p_boxplot,
        width = 10,
        height = 5,
        device = "pdf"
      )
      ggsave(
        paste0(file.path(plot_out, motif, "boxplot"), ".png"),
        p_boxplot,
        width = 10,
        height = 5,
        device = "png"
      )


      motif_position_event_pval_limit  <- lapply(
          10^-(c(0, 2, 4, 6, 8, seq(10, 50, by = 10))),
          function(max_p_val) {
            motif_position_filt <- motif_position[[motif]][p_val < max_p_val]
            if(nrow(motif_position_filt > 0)) {
              # Box plot of motifs
              p_boxplot_filtered  <- aggregated_motif_boxplot(motif_position_filt)
              dir.create(file.path(plot_out, motif, "boxplot_filt"), showWarnings = FALSE)
              ggsave(
                paste0(file.path(plot_out, motif, "boxplot_filt"), "/", max_p_val, ".pdf"),
                p_boxplot_filtered,
                width = 10,
                height = 5,
                device = "pdf"
              )
              ggsave(
                paste0(file.path(plot_out, motif, "boxplot_filt"), "/", max_p_val, ".png"),
                p_boxplot_filtered,
                width = 10,
                height = 5,
                device = "png"
              )
            } else{
              return()
            }

            # Number of events at each relaitve motif position
            motif_position_filt[
                , .(n = .N), by = .(N, rel_pos, strand)
              ][
                , event_pval_limit := max_p_val
              ][
                , event_pval_limit := reorder(as.character(event_pval_limit), event_pval_limit, decreasing = TRUE)
              ][
                , labs := reorder(paste0(rel_pos, "\n", N), rel_pos)
              ]
          }
        )  %>%
        rbindlist()

      # Plot number of events remaining at different thresholds
      p_number_of_event <-  ggplot(motif_position_event_pval_limit) +
          aes(x = rel_pos, y = n, fill = event_pval_limit, text = event_pval_limit) +
          geom_area(position = "identity") +
          facet_wrap(~strand) +
          default_theme_SH() +
          ggplot2::scale_fill_viridis_d() +
          ggplot2::scale_x_continuous(
            labels = unique(motif_position_event_pval_limit$labs),
            breaks = unique(motif_position_event_pval_limit$rel_pos),
            expand = c(0, 0)
          ) +
          scale_y_continuous(
            expand = c(0, 0)
          )

      ggsave(
        file.path(plot_out, motif, "event_threshold.pdf"),
        p_number_of_event,
        width = 11,
        height = 5,
        device = "pdf"
      )
      ggsave(
        file.path(plot_out, motif, "event_threshold.png"),
        p_number_of_event,
        width = 11,
        height = 5,
        device = "png"
      )
    }
  )
}

################################################################################

#'
#'
#' Visualise the identified clusters
#'
#' @param out output folder
#' @param path_ref path to reference
#' @param path_chunk_stats path to preprocess chunks
#' @param n_extra_positions number of adjacent positions to include in plot
#' @param motifs_evaluated evaluation type (all, clustered, or vector of motifs)
#' @param entropy_threshold_motif Minimum entropy (how stringent motif detection is)
#' @return plots
#' @import msa
#' @export
find_motifs <- function(
  path_chunk_stats,
  path_ref,
  out,
  umap_args = list(
    n_components = 2,
    min_dist = 0.05,
    n_neighbors = 30
    ),
  hdbscan_args =  list(minPts = 15),
  event_pval_limit = 1e-20,
  event_sequence_frame = 8,
  align_event_sequences = TRUE,
  iterations = 3,
  iteration_approach = "remove_noise",
  entropy_threshold_motif = 1
) {
  # Load processed chunks
  chunks <- load_processed_chunks(path_chunk_stats)
  continue <- TRUE
  iter <- 1
  while (continue) {
    dir.create(file.path(out, iter), recursive = TRUE, showWarnings = FALSE)

    # Embed and cluster processed chunks
    clusters <- embed_and_cluster(
      chunks,
      ref_path = path_ref,
      out = file.path(out, iter),
      event_pval_limit = event_pval_limit,
      event_sequence_frame = event_sequence_frame,
      umap_args = umap_args,
      hdbscan_args = hdbscan_args
    )
    # Get cluster events sequences
    setorder(clusters, cluster)

    cluster_sequences <- lapply(
        split(clusters[!is.na(cluster)], by = "cluster"),
        function(x) {
          toupper(x$seq)
        }
      )
    if (length(cluster_sequences) == 0) {
      logger::log_debug("No clusters")
      clusters[
        , motif := ""
      ]
      break()
    }
    if (align_event_sequences) {
      cluster_sequences_aligned <- lapply(
        cluster_sequences,
        function(sequences) {
          alignment <- msa::msa(sequences, type = "dna", gapOpen = 1000) %>%
            msa::msaConvert(type = "seqinr::alignment")
          alignment$seq
        })
      names(cluster_sequences_aligned) <- names(cluster_sequences)
      cluster_sequences <- cluster_sequences_aligned
    }
    # Calculate entropy and select motifs
    add_motifs <- function(cluster_sequences) {
      cluster_entropy <- lapply(
          cluster_sequences,
          function(x) calculate_bit_score(x, min_entropy = entropy_threshold_motif)
        )
      cluster_motifs <- lapply(
        names(cluster_entropy),
        function(clust) {
          data.table(
            cluster = as.factor(clust),
            motif = cluster_entropy[[clust]]$motif_consensus
      )}) %>%
        rbindlist()
      return(cluster_motifs)
    }
    # Join cluster motif
    clusters[
      add_motifs(cluster_sequences), on = "cluster", motif := i.motif
    ]
    # Save clustering
    fwrite(clusters, file.path(out, iter, "clusters.tsv"), sep = "\t")

    switch(iteration_approach,
      "remove_biggest_cluster" = {
        # Identify biggest cluster
        clusters_to_keep <- clusters$cluster %>%
          table() %>%
          sort(decreasing = TRUE) %>%
          `[`(-1) %>%
          names() %>%
          as.numeric()
        # Remove biggest cluster events from `chunks`
        clusters_kept <- clusters[
            cluster %in% c(clusters_to_keep, NA), .SD, .SDcols = c("pos_ref", "chunk_ref", "cluster")
          ][
            , keep := TRUE
          ]
        chunks[, keep := FALSE]
        chunks <- chunks[
            clusters_kept, on = c("pos_ref", "chunk_ref"), keep := i.keep
          ][
            keep == TRUE
          ]
      },
      "remove_noise" = {
        # Remove noise events from `chunks`
        clusters_kept <- clusters[
            !is.na(cluster), .SD, .SDcols = c("pos_ref", "chunk_ref", "cluster")
          ][
            , keep := TRUE
          ]
        chunks[, keep := FALSE]
        chunks <- chunks[
            clusters_kept, on = c("pos_ref", "chunk_ref"), keep := i.keep
          ][
            keep == TRUE
          ]
      },
    stop("Invalid `iteration_approach` value ('remove_biggest_cluster' or 
          'remove_noise')")
    )
    iter <- iter + 1
    if (iter > iterations) {
      continue <- FALSE
    }
  }
}
