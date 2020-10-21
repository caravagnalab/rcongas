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
plot_counts_rna_segments = function(x,
                                    normalised = TRUE,
                                    z_score = FALSE,
                                    cutoff_p = 0.01)
{
  # Get segments_input
  segments_input = x$data$counts

  # Input segmentation - get size in Megabases (Mb)
  input_segments = Rcongas::get_input_segmentation(x) %>%
    dplyr::mutate(size = ceiling((to - from) / 10 ^ 6),
                  label_chr = paste0(chr, " (", size, 'Mb)')) %>%
    Rcongas:::idify()

  # RNA data
  RNA = Rcongas::get_counts(x, normalise = normalised, z_score = z_score) %>%
    Rcongas:::idify()

  RNA = dplyr::left_join(RNA,
                         input_segments %>% dplyr::select(segment_id, label_chr),
                         by = "segment_id") %>%
    dplyr::rename(segment = segment_id)

  # summary stats
  ngenes = sum(Rcongas::get_input_segmentation(x)$mu)
  nsegments = nrow(Rcongas::get_input_segmentation(x))
  MB_covered = round(sum(Rcongas::get_input_segmentation(x)$size) / 10 ^ 6)

  # prepare plot caption
  # caption = paste0("RNA counts are ",
  #                  ifelse(normalised, "normalised.", "not normalised."))
  caption = ""

  # Cluster assignments
  clustering = Rcongas::get_clusters(x) %>%
    dplyr::arrange(desc(cluster))

  # Clustering assignments plot
  clusters_colors = Rcongas:::get_clusters_colors(clustering$cluster)

  # clustering_sideplot = ggplot(clustering) +
  #   geom_tile(aes(
  #     x = 'Cluster',
  #     y = factor(cell, levels = cell),
  #     fill = cluster
  #   )) +
  #   CNAqc:::my_ggplot_theme() +
  #   scale_fill_manual(values = clusters_colors) +
  #   labs(y = 'Cell', x = 'Segment') +
  #   theme(axis.text.y = element_blank(),
  #         axis.text.x = element_text(angle = 90, hjust = 1)) +
  #   guides(fill = guide_legend(NULL, ncol = 1)) +
  #   labs(y = NULL, x = NULL)

  # RNA plot
  rna_plot = ggplot(RNA) +
    geom_tile(aes(
      x = factor(label_chr, levels = gtools::mixedsort(label_chr) %>% unique()),
      y = factor(cell, levels = clustering$cell),
      fill = n
    )) +
    CNAqc:::my_ggplot_theme() +
    labs(y = "Cell", x = "Segment") +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1)) +
    guides(fill = guide_colorbar("Counts", barwidth = unit(3, 'cm'))) +
    labs(
      # titl = paste0("Raw CONGAS dataset (k =", Rcongas::get_k(x), ')'),
      title = paste0(
        ngenes,
        " genes (",
        nsegments,
        " segments, ",
        MB_covered,
        " Mb)"
      )
    )

  if (!z_score)
    rna_plot = rna_plot + scale_fill_distiller(palette = 'GnBu', direction = 1)
  else
  {
    # Change the legend and palettes for the z-score
    rna_plot = rna_plot +
      scale_fill_gradient2(low = "steelblue", high = 'indianred3') +
      guides(fill = guide_colorbar("Z-score", barwidth = unit(3, 'cm')))
  }

  # Annotate special areas of the tile - first split the clusters
  nclusters = Rcongas::get_k(x)
  cluster_size = Rcongas::get_clusters_size(x)

  if (nclusters > 1)
    rna_plot = rna_plot +
    geom_hline(yintercept = cumsum(cluster_size[-1]),
               size = 0.3,
               linetype = 'dashed')

  # .. then put some squares around certain segments
  chrs_to_annotate = Rcongas::get_segment_test_counts(x,
                                                      group1 = 1,
                                                      group2 = 2,
                                                      cutoff_p = cutoff_p) %>%
    dplyr::filter(sign) %>%
    Rcongas:::idify() %>%
    dplyr::left_join(RNA %>% Rcongas:::idify(), by = 'segment_id') %>%
    dplyr::select(segment_id, label_chr) %>%
    dplyr::pull(label_chr) %>%
    unique()

  # Annotation of squares uses the x-axis ordering to determine
  # the perimeter of rectangles that will highlight the data
  if (length(chrs_to_annotate) >= 1)
  {
    annotation_color = 'mediumseagreen'
    order_x_axis = gtools::mixedsort(RNA$label_chr) %>% unique()

    chrs_to_annotate = Reduce(dplyr::bind_rows,
                              lapply(chrs_to_annotate, function(y) {
                                data.frame(
                                  xmin = which(y == order_x_axis) - 0.5,
                                  xmax = which(y == order_x_axis) + 0.5,
                                  ymin = 0,
                                  ymax = Inf
                                )
                              }))

    rna_plot = rna_plot +
      geom_rect(
        data = chrs_to_annotate,
        aes(
          xmin = xmin,
          xmax = xmax,
          ymin = ymin,
          ymax = ymax
        ),
        fill = NA,
        color = annotation_color,
        alpha = .2,
        size = 0.1
      )

    # edit caption
    caption = paste0(caption,
                     "Shaded area: Poisson test significant at level ",
                     cutoff_p)
  }

  # Add caption
  rna_plot = rna_plot + labs(caption = caption)

  # cowplot::plot_grid(
  #   rna_plot,
  #   clustering_sideplot,
  #   nrow = 1,
  #   align = 'h',
  #   rel_widths = c(1, .1)
  # )


  rna_plot = rna_plot +
    geom_point(data = clustering,
               aes(x = 0, y = cell, color = cluster),
               size = 3) +
    scale_color_manual(values = clusters_colors)


}