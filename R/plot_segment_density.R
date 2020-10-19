#' Title
#'
#' @param x
#' @param segment_id
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plot_segment_density = function(x, segments_ids, group1 = 1, group2 = 2, ...)
{
  counts_data = Rcongas:::get_counts(x) %>% Rcongas:::idify() %>%
    dplyr::filter(segment_id %in% !!segments_ids)

  mean_data = counts_data %>%
    dplyr::group_by(segment_id, cluster) %>%
    dplyr::summarise(mean = mean(n), .groups = 'keep')

  # Poisson p-value
  test_pvalue = Rcongas:::get_segment_test_counts(x, group1 = group1, group2 = group2, ...) %>%
    Rcongas:::idify() %>%
    dplyr::filter(segment_id %in% !!segments_ids)


  # Poisson parameters - plot density
  bm = Rcongas:::get_best_model(x)


  lighten <- function(color, factor=0.5){
    if ((factor > 1) | (factor < 0)) stop("factor needs to be within [0,1]")
    col <- col2rgb(color)
    col <- col + (255 - col)*factor
    col <- rgb(t(as.matrix(apply(col, 1, function(x) if (x > 255) 255 else x))), maxColorValue=255)
    col
  }

  clusters_colors = Rcongas:::get_clusters_colors(counts_data$cluster)


  counts_data %>%
    ggplot(aes(n, fill = cluster)) +
    facet_wrap(~segment_id, ncol = 2, scales = 'free') +
    geom_histogram(bins = 100) +
    CNAqc:::my_ggplot_theme() +
    scale_fill_manual(values = sapply(clusters_colors, lighten, factor = .3)) +
    geom_vline(
      data = mean_data,
      aes(xintercept = mean, color = cluster),
      size = .5,
      color = 'black',
      linetype = 'dashed'
    ) +
    scale_color_manual(values = Rcongas:::get_clusters_colors(counts_data$cluster)) +
    labs(title = "Poisson counts") +
    labs(x = "Counts", y = "Observations") +
   geom_text(
     data = test_pvalue,
     aes(x = 0, y = Inf, label = paste0('p ', format.pval(p, eps = 0.001))),
     inherit.aes = FALSE,
     size = 2.5,
     hjust = 0,
     vjust = 1.3
   ) +
  coord_cartesian(clip = 'off')
}
