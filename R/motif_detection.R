
#' Embed joined chunks
#' 
#' Wrapper function for umbedding and clustering. Take the data.table containing
#' the output the preprocess_chunk function. Precessed chunks can be loaded using 
#' load_processed_chunks
#' 
#' 
#' @param chunks preprocessed chunks generated with preprocess_all_chunks
#' @param ref_path path to the reference passed to megalodon
#' @param out character string with path to output results
#' @param p_value_threshold Numeric, max p-value for events
#' @param event_sequence_frame integer, number of nucleotide to include to each side
#' @param umap_args list of arguments passed directly to umap (See ?umap::umap)
#' https://cran.r-project.org/web/packages/umap/vignettes/umap.html
#' @param hdbscan_Args list of arguments passed directly to hdbscan (See ?dbscan::hdbscan)
#' @return embedded and clustered positions
#' @export
embed_and_cluster <- function(
    chunks,
    ref_path,
    out,
    p_value_threshold = 1e-6,
    event_sequence_frame = 5,
    umap_args = list(),
    hdbscan_args = list()
  ) {
  checkmate::assert_data_table(chunks)

  # Prepare features for embedding
  features <- get_event_features(
    chunks, 
    p_value_threshold
  )

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
  logger::log_debug(glue::glue("Events: {nrow(features)}"))

  # Embedding
  logger::log_debug("Embedding")
  emb_umap <- do.call(embed_umap, list(features, umap_args = umap_args))
  fwrite(emb_umap, paste0(out, "umap.tsv"))

  logger::log_debug(glue::glue("Clustering"))
  do.call(cluster_emb, list(emb_umap, hdbscan_args = hdbscan_args))
  emb_umap <- cbind(
    emb_umap,
    features_metainfo
  )
  return(emb_umap)
}

#' Load preprocessed chunks
#' 
#' Loads the output of the preprocess function
#' 
#' @param preprocess_out path specified by 'out' in the preprocess function
#' @return joined chunks
#' @export
load_processed_chunks <- function(preprocess_out) {
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

embed_umap <- function(
    features,
    umap_args = list()
  ) {
  checkmate::assert_data_table(
    features,
    types = c("numeric"),
    any.missing = FALSE
  )
  umap_args$d <- features
  embedding <- do.call(umap::umap, umap_args)
  dt <- data.table(embedding$layout)
  setnames(
    dt,
    old = names(dt),
    new = paste0("UMAP", 1:ncol(dt))
  )
  return(dt)
}



cluster_emb <- function(
    embedding,
    hdbscan_args = list()
  ){
  # Argument check
  checkmate::assert_data_table(
    embedding,
    types = "numeric",
    any.missing = FALSE
  )

  hdbscan_args$x <- embedding
  clusters_hdbscan <- do.call(dbscan::hdbscan, hdbscan_args)
  embedding[
      , HDBSCAN := clusters_hdbscan$cluster
    ][
      , cluster_prob := clusters_hdbscan$membership_prob
    ]
  embedding[
      , HDBSCAN := fifelse(HDBSCAN == 0, NA_integer_, HDBSCAN)
    ][
      , HDBSCAN := as.factor(HDBSCAN)
    ]
}


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

#' Calculate bit score and extract motifs
#' 
#' The entropy at each postion in supplied sequences. The entropy
#' calculated based on
#' https://en.wikipedia.org/wiki/Sequence_logo#Logo_creation.
#' 
#' @param sequences list of sequences of the same length
#' @return list(entropy_dt, motif_dt, motif_vec)
#' @export
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
    small_n_correction <- (1/log(2)) * 3/(2*n_sequences)

    # Nucleotides count at each position
    count_positional_character  <- function(seq, character_regex){
        sequence_length = nchar(seq)[1]
        lapply(1:sequence_length,
            function(position) {
                paste0("^", paste(rep(".", position - 1), collapse = ""), character_regex) %>%
                    grepl(seq) %>%
                    sum()
            }
        ) %>% unlist()
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
      , position := 1:.N
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
        , freq := n / n_sequences
      ][
        , entropy_loss := - freq * log2(freq)
      ]
    replace_na_dt(dt)

    # Entropy at each position
    dt[
        , entropy := uniform_entropy - sum(entropy_loss) - small_n_correction
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
        , .(paste(nucleotide, collapse = ">"))
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


#' Visualise clusters
#' 
#' Visualise the identified clusters 
#' 
#' @param dt data.table outputted from embed_and_cluster
#' @return list(scatter_plot, logo_plot)
#' @export
visualise_clusters <- function(
    dt
  ) {
    # Check inputs
    checkmate::assert_data_table(dt)
    checkmate::assert_subset(
      c("HDBSCAN", "UMAP1", "UMAP2", "seq"),
      names(dt)
    )

    # Get list of event sequences by cluster
    cluster_sequences <- lapply(
      split(dt[!is.na(HDBSCAN)], by = "HDBSCAN"),
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
      aes(x = UMAP1, y = UMAP2, color = HDBSCAN) +
      geom_point() +
      default_theme_SH() +
      ggplot2::scale_colour_viridis_d(na.value = "gray80")

    plot_scatter_motif <- ggplot(dt) +
      aes(x = UMAP1, y = UMAP2, color = motif) +
      geom_point() +
      default_theme_SH() +
      ggplot2::scale_colour_viridis_d(na.value = "gray80")

    return(
      list(
        plot_scatter,
        plot_logo,
        plot_scatter_motif
      )
    )
  }



