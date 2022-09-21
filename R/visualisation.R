
range_01 <- function(vec) {
  return((vec - min(vec))/max(vec - min(vec)))
}

add_rel_pos_current_pos <- function(chunk) {
  chunk[
      ,
      current_pos_decimal := round(
        pos_ref + range_01(pos_signal),
        4),
      by = .(read_id, pos_ref, strand, type)
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

plot_mapped_current_boxplot <- function(chunk, start_plot, 
                                        col = c("#1f1da5", "#c08731"),
                                        fill = c("#1f1da580", "#c0873180")) {
  chunk_zoom <- get_chunk_zoom(chunk, start_plot)
  setorder(chunk_zoom, pos_ref)
  chunk_zoom %>%
    ggplot() +
    default_theme_SH() +
    add_n_mapped_label(chunk_zoom) +
    aes(x = as.factor(pos_ref - start_plot), y = signal) +
    geom_boxplot(aes(fill = type)) +
    geom_line(aes(group = type, y = mean_current, col = type), size = 2) +
    labs(
      x = "Relative Reference Position",
      y = "Normalised current values",
      title = glue::glue("Mapped current at position {start_plot}")
    ) +
    scale_color_manual(values = col) +
    scale_fill_manual(values = fill) +
    facet_wrap(~strand) 
}

plot_mapped_current_lines <- function(chunk, start_plot,
                                        col = c("#1f1da5", "#c08731"),
                                        fill = c("#1f1da580", "#c0873180")) {
  chunk_zoom <- get_chunk_zoom(chunk, start_plot)
  add_rel_pos_current_pos(chunk_zoom)
  setorder(chunk_zoom, pos_ref)
  chunk_zoom %>%
    ggplot() +
    default_theme_SH() +
    add_n_mapped_label(chunk_zoom) +
    aes(x = current_pos_decimal - start_plot, col = type) +
    geom_line(aes(group = read_id, y = signal), size = 0.2, alpha = 0.7) +
    geom_line(aes(group = type, y = mean_current), size = 1) +
    theme(
      panel.grid.major.x = element_line(color = "#616161")
    ) +
    labs(
      x = "Reference Position",
      y = "Normalised current values",
      title = glue::glue("Mapped current at position {start_plot}")
    ) +
    scale_color_manual(values = col) +
    scale_fill_manual(values = fill) +
    facet_wrap(~strand)
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
  default_theme_SH() +
  labs(
    x = "Contig position",
    y = "N mapped reads"
  ) +
  facet_wrap(~strand)
}

#' Plot signal of chunk
#' 
#' Plots the signal and estimated coverage of a chunk.
#' 
#' @param chunk data.table, chunk with signal and pos_signal
#' @param path_out character, output path of plots
#' @param start_plot integer, start position of zoom plots (default: "random")
#' @param frame ineteger, Number of positions to include in zoom plots (default: 10)
#' @param batch_size integer, size of batching when estimating coverage (default: 100)
#' @return nothing, save plot to path_out
#' @import ggplot2
#' @export
plot_chunk_current <- function(chunk,
                               path_out,
                               start_plot = "random",
                               frame = 10,
                               batch_size = 100) {
  if (start_plot == "random") {
    start_plot <- sample(chunk$pos_ref, 1)
  }


  dir.create(path_out, showWarnings = FALSE, recursive = TRUE)
  ggsave(
    glue::glue("{path_out}/reference_mapped_curents.png"),
    width = 10,
    height = 5,
    plot_mapped_current_lines(
      chunk,
      start_plot = start_plot
    )

  )
  ggsave(
    glue::glue("{path_out}/reference_mapped_curents_boxplot.png"),
    width = 10,
    height = 5,
    plot_mapped_current_boxplot(
      chunk,
      start_plot
    )
  )

  ggsave(
    glue::glue("{path_out}/reference_n_mapped_reads.png"),
    width = 10,
    height = 5,
    plot_chunk_coverage(
      chunk,
      batch_size
    )
  )
}



default_theme_SH = function(){
  theme_grey(base_size = 14, base_family = "sans") %+replace% 
    theme(
      # plot margin
      plot.margin = unit(rep(0.5, 4), "cm"),
      # plot background and border
      plot.background = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank(),
      # grid lines
      panel.grid.major.x = element_line(size = 0.5, color = "#cbcbcb00"),
      panel.grid.major.y = element_line(size = 0.5, color = "#cbcbcb"),
      panel.grid.minor = element_blank(),
      # axis ticks and lines
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      # title, subtitle and caption
      plot.title = element_text(size = 25, face = "bold", colour = "#757575",
                                hjust = 0, margin = margin(9, 20, 9, 0)),
      plot.subtitle = element_text(size = 24, colour = "#757575", hjust = 0,
                                   margin = margin(9, 0, 18, 9)),
      plot.caption = element_text(size = 15, color = "grey50", hjust = 1,
                                  margin = margin(t = 15)),
      # axes titles
      axis.title = element_text(size = 15, colour = "#757575", hjust = 1),
      axis.text.x = element_text(size = 10, margin = margin(b = 14)),
      axis.text.y = element_text(size = 10, margin = margin(l = 14)),
      # legend
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.title = element_text(size = 12, colour = "#757575"),
      legend.text.align = 0,
      legend.text = element_text(size = 14, colour = "#757575"),
      # facetting
      strip.background = element_rect(fill = "transparent", colour = NA),
      strip.text = element_text(size = 12, face = "bold", colour = "#757575",
                                hjust = 0)
    )
}



