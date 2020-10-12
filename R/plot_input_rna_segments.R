#' Title
#'
#' @param x
#' @param normalised
#'
#' @return
#' @export
#'
#' @examples
#'
#' input = get_counts_matrix(x)
plot_counts_rna_segments = function(x, normalised = TRUE, z_score = FALSE)
{
  # Get segments_input
  segments_input = x$data$counts

  # Input segmentation - get size in Megabases (Mb)
  input_segments = get_input_segmentation(x) %>%
    dplyr::mutate(size = ceiling((to - from) / 10 ^ 6),
                  label_chr = paste0(chr, " (", size, 'Mb)')) %>%
    idify()

  # RNA data
  RNA = get_counts(x, normalise = normalised, z_score = z_score) %>%
    idify()

  RNA = dplyr::left_join(RNA, input_segments %>% dplyr::select(segment_id, label_chr), by = "segment_id") %>%
    dplyr::rename(segment = segment_id)

  # summary stats
  ngenes = sum(get_input_segmentation(x)$mu)
  nsegments = nrow(get_input_segmentation(x))
  MB_covered = round(sum(get_input_segmentation(x)$size) / 10 ^ 6)

  # Cluster assignments
  clustering = get_clusters(x) %>%
    arrange(desc(cluster))

  # Clustering assignments plot
  clusters_colors = get_clusters_colors(clustering$cluster)

  clustering_sideplot = ggplot(clustering) +
    geom_tile(aes(
      x = 'Cluster',
      y = factor(cell, levels = cell),
      fill = cluster
    )) +
    CNAqc:::my_ggplot_theme() +
    scale_fill_manual(values = clusters_colors) +
    labs(y = 'Cell', x = 'Segment') +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1)) +
    guides(fill = guide_legend(NULL, ncol = 1)) +
    labs(y = NULL, x = NULL)

  # RNA plot
  rna_plot = ggplot(RNA) +
    geom_tile(aes(
      x = factor(label_chr, levels = gtools::mixedsort(label_chr) %>% unique()),
      y = factor(cell, levels = clustering$cell),
      fill = n
    )) +
    CNAqc:::my_ggplot_theme() +
    scale_fill_distiller(palette = 'GnBu', direction = 1) +
    labs(y = "Cell", x = "Segment") +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1)) +
    guides(fill = guide_colorbar("Counts per segment", barwidth = unit(3, 'cm'))) +
    labs(
      title = paste0("Raw CONGAS dataset (k =", get_k(x), ')'),
      subtitle = paste0(ngenes, " genes mapped to ", nsegments, " segments."),
      caption = paste0(
        "Input segments span ",
        MB_covered,
        " Mb. Counts are ",
        ifelse(normalised, "noramalised.", "not normalised.")
      )
    )
  # +
  #   scale_x_discrete(labels = RNA$label_chr)

  cowplot::plot_grid(
    rna_plot,
    clustering_sideplot,
    nrow = 1,
    align = 'h',
    rel_widths = c(1, .1)
  )


}