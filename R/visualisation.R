range_01 <- function(vec) {
  return((vec - min(vec))/max(vec - min(vec)))
}

add_rel_pos_current_pos <- function(chunk) {
  chunk[
      ,
      current_pos_decimal := round(
        pos_ref + range_01(pos_signal),
        4),
      by = .(read_id, pos_ref)
    ]
}

get_chunk_zoom <- function(chunk, start_pos, frame = 10) {
  chunk[
    (pos_ref > start_pos) & (pos_ref < (start_pos + frame)),
  ][
    , mean_current := mean(signal), by = .(type, pos_ref, strand)
  ][
    , strand := fifelse(strand == "+", "Plus strand", "Minus strand")
  ]
}

get_n_reads_in_zoom <- function(chunk_zoom) {
  chunk_zoom[
    , min_y := min(signal)
  ][
    , .(n_reads = length(unique(read_id)),
        min_y = min_y[1]),
    by = .(strand, type)
  ][
    , .(label = paste0("N reads: ",
                        paste0(paste0(type, " ", n_reads), collapse = ": ")),
          min_y = min_y[1]),
    by = strand
  ]
}

add_n_mapped_label <- function(chunk_zoom) {
      geom_label(
        data = get_n_reads_in_zoom(chunk_zoom),
        aes(label = label, group = strand, y = min_y),
        x = 1,
        hjust = 0,
        vjust = -0.1,
        color = "#757575"
    )
}

plot_mapped_current_boxplot <- function(chunk_zoom, start_plot, 
                                        col = c("#1f1da5", "#c08731"),
                                        fill = c("#1f1da580", "#c0873180")) {
  chunk_zoom %>%
    ggplot() +
    #default_theme_SH() +
    add_n_mapped_label(chunk_zoom) +
    aes(x = as.factor(pos_ref - start_plot), y = signal) +
    geom_boxplot(aes(fill = type)) +
    geom_line(aes(group = read_id, y = mean_current, col = type), size = 2) +
    labs(
      x = "Reference Position",
      y = "Normalised current values",
      title = glue("Mapped current at position {start_plot}")
    ) +
    scale_color_manual(values = col) +
    scale_fill_manual(values = fill) +
    facet_wrap(~strand) 
}
plot_mapped_current_lines <- function(chunk_zoom, start_plot, frame,
                                        col = c("#1f1da5", "#c08731"),
                                        fill = c("#1f1da580", "#c0873180")) {
  chunk_zoom %>%
    ggplot() +
    default_theme_SH() +
    add_n_mapped_label(chunk_zoom) +
    aes(x = current_pos_decimal - start_plot, col = type) +
    geom_line(aes(group = read_id, y = signal), size = 0.2, alpha = 0.7) +
    geom_line(aes(group = read_id, y = mean_current), size = 1) +
    scale_x_continuous(breaks = seq(0, frame, 1)) +
    theme(
      panel.grid.major.x = element_line(color = "#616161")
    ) +
    labs(
      x = "Reference Position",
      y = "Normalised current values",
      title = glue("Mapped current at position {start_plot}")
    ) +
    scale_color_manual(values = col) +
    scale_fill_manual(values = fill) +
    facet_wrap(~strand) +
    geom_text(
      aes(x = pos_ref - start_plot, label = base),
      y = -Inf, vjust = 0, col = "#757575"
    )
}



get_chunk_coverage <- function(chunk, batch_size) {
  chunk[
      , .(N = .N), by = .(read_id, type, pos_ref, strand)
    ][
      , .(N = .N), by = .(type, pos_ref, strand)
    ][
      , batch := pos_ref %/% batch_size
    ][
      , .(N = mean(N)), by = .(type, batch, strand)
    ] 
}

plot_chunk_coverage <- function(chunk, batch_size = 100) {
  get_chunk_coverage(chunk, batch_size) %>% 
  ggplot(aes(x = batch * batch_size, y = N, col = type)) +
  geom_line(size = 1.5) +
  scale_color_manual(values = c("#1f1da5", "#c08731")) +
  labs(
    x = "Contig position",
    y = "N mapped reads"
  ) +
  facet_wrap(~strand)
}


plot_chunk_current <- function(chunk,
                               path_out,
                               start_plot = "random",
                               frame = 10,
                               batch_size = 100) {
  if (start_plot == "random") {
    start_plot <- sample(chunk$pos_ref, 1)
  }

  chunk_zoom <- get_chunk_zoom(chunk, start_plot, frame = frame)
  add_rel_pos_current_pos(chunk_zoom)

  dir.create(path_out, showWarnings = FALSE, recursive = TRUE)
  ggsave(
    glue("{path_out}/reference_mapped_curents.png"),
    width = 10,
    height = 5,
    plot_mapped_current_lines(
      chunk_zoom,
      start_plot = start_plot,
      frame = frame
    )

  )
  ggsave(
    glue("{path_out}/reference_mapped_curents_boxplot.png"),
    width = 10,
    height = 5,
    plot_mapped_current_boxplot(
      chunk_zoom,
      start_plot
    )
  )

  ggsave(
    glue("{path_out}/reference_n_mapped_reads.png"),
    width = 10,
    height = 5,
    plot_chunk_coverage(
      chunk,
      batch_size
    )
  )
}

