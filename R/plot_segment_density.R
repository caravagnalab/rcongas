#' Plots the data and fit density for a set of segments.
#'
#' @description This function works for a clustering result, whatever is the
#' adopted model. It returns a list of plots for each one of the input segments.
#' Each plot shows the data, plus the fit density per cluster. It can be
#' used to assess the quality of the data and fits.
#'
#' @param x Input object with clusters.
#' @param segment_ids Vector of segment ids in the CONGAS format (e.g, \code{"chr12:40515401:133219374"}).
#' By default all are used via function \code{get_segment_ids}.
#' @param ...
#'
#' @return A list of ggplot objects.
#' @export
#'
#' @examples
#'
#' x = Rcongas::congas_example
#'
#' # Plot the first segment
#' plot_segment_density(x, segments_ids = get_segment_ids(x)[1])
#'
#' # Plot the first one that is highlighted
#' plot_segment_density(x, segments_ids = get_segment_ids(x, highlight = TRUE)[1])
#'
#' # Plot the top 4 highlighted, in 2x2 format, specifying alpha value
#'
#' what_to_plot = get_segment_ids(x, highlight = TRUE, alpha = 0.05)[1:4]
#'
#' ggpubr::ggarrange(
#'  plotlist = plot_segment_density(x, segments_ids = what_to_plot),
#'  nrow = 2,
#'  ncol = 2
#')
plot_segment_density = function(x,
                                segments_ids = get_segment_ids(x),
                                sum_denominator = TRUE)
{
  if (!has_inference(x))
    stop("Cannot plot this density without clusters, or the input is not a CONGAS object.")
  
  if (length(segments_ids) == 0)
  {
    cli::cli_alert_warning("No input segments for plotting.")
    return(ggplot())
  }
  
  if (is_gaussian(x))
    plots = lapply(
      segments_ids,
      plot_single_segment_gaussian,
      x = x,
      sum_denominator = sum_denominator
    )
  else
    plots = lapply(
      segments_ids,
      plot_single_segment_poisson,
      x = x,
      sum_denominator = sum_denominator
    )
  
  names(plots) = segments_ids
  
  return(plots)
}

# Visualisation for a Poisson-based input
plot_single_segment_poisson = function(x, segment, sum_denominator)
{
  # Counts data
  counts_data = get_counts(x, normalise = TRUE,  sum_denominator = sum_denominator) %>%
    idify() %>%
    dplyr::filter(segment_id == segment)
  
  # Coloring
  clusters_colors = get_clusters_colors(counts_data$cluster)
  clusters = get_clusters_size(x) %>% names()
  
  # This makes a progress bar
  density_points = easypar::run(
    FUN = function(cluster) {
      get_poisson_density_values(x = x,
                                 segment_id = segment,
                                 cluster = cluster)
    },
    lapply(clusters, list),
    parallel = FALSE
  )
  
  density_points = Reduce(dplyr::bind_rows, density_points)
  
  # 2 plots
  nbins = 100
  
  density_plot = density_points %>%
    ggplot(aes(x = x, y = y, color = cluster)) +
    facet_wrap( ~ paste0(segment_id, ' (Poisson)')) +
    geom_point(size = .8) +
    geom_line(size = .3) +
    # CNAqc:::my_ggplot_theme() +
    theme_linedraw(base_size = 10) +
    scale_color_manual(values = clusters_colors) +
    labs(x = NULL, y = "Density") +
    coord_cartesian(clip = 'off') +
    guides(color = FALSE) +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
  
  data_plot = counts_data %>%
    ggplot(aes(
      n,
      y = ..count..,
      fill = cluster,
      group = cluster
    )) +
    geom_histogram(bins = nbins) +
    # facet_wrap(~ segment_id) +
    CNAqc:::my_ggplot_theme() +
    scale_fill_manual(values = sapply(clusters_colors, Rcongas:::lighten, factor = .1)) +
    labs(x = "RNA counts", y = "Observations") +
    coord_cartesian(clip = 'off') +
    guides(fill = guide_legend("Cluster", ncol = 1)) +
    theme(legend.position = c(.9, 0.85), legend.background = element_blank())
  
  cowplot::plot_grid(
    density_plot,
    data_plot,
    ncol = 1,
    align = 'v',
    axis = 'lr',
    rel_heights = c(1, 2)
  )
}

# Visualisation for a Gaussian input
plot_single_segment_gaussian = function(x, segment, sum_denominator)
{
  # Counts data
  counts_data = Rcongas:::get_counts(x, normalise = TRUE,  sum_denominator = sum_denominator) %>%
    Rcongas:::idify() %>%
    dplyr::filter(segment_id == segment)
  
  # Coloring
  clusters_colors = get_clusters_colors(counts_data$cluster)
  clusters = Rcongas::get_clusters_size(x) %>% names()
  
  # We need to make some adjustments to get the right binning, scaling etc
  nbins = 30
  
  binsize = counts_data %>%
    dplyr::group_by(segment_id) %>%
    dplyr::summarise(min = min(n),
                     max = max(n),
                     .groups = 'keep') %>%
    dplyr::mutate(binsize = (max - min) / nbins) %>%
    dplyr::pull(binsize)
  
  # This makes a progress bar
  density_points = easypar::run(
    FUN = function(cluster) {
      get_gaussian_density_values(x = x,
                                  segment_id = segment,
                                  cluster = cluster)
    },
    lapply(clusters, list),
    parallel = FALSE
  )
  
  density_points = Reduce(dplyr::bind_rows, density_points)
  
  # 2 plots
  
  density_plot = density_points %>%
    ggplot(aes(x = x, y = y, color = cluster)) +
    facet_wrap( ~ segment_id) +
    geom_point(size = .8) +
    geom_line(size = .3) +
    CNAqc:::my_ggplot_theme() +
    scale_color_manual(values = clusters_colors) +
    labs(title = "Gaussian density and counts") +
    labs(x = NULL, y = "Density") +
    coord_cartesian(clip = 'off') +
    guides(color = FALSE) +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
  
  data_plot = counts_data %>%
    ggplot(aes(
      n,
      y = ..count..,
      fill = cluster,
      group = cluster
    )) +
    geom_histogram(bins = nbins) +
    CNAqc:::my_ggplot_theme() +
    scale_fill_manual(values = sapply(clusters_colors, Rcongas:::lighten, factor = .1)) +
    labs(x = "RNA counts", y = "Observations") +
    coord_cartesian(clip = 'off') +
    guides(fill = guide_legend("Cluster", nrow = 1))
  
  cowplot::plot_grid(
    density_plot,
    data_plot,
    ncol = 1,
    align = 'v',
    axis = 'lr',
    rel_heights = c(1, 2)
  )
}



# Density functions for the fit plots - Poisson and Gaussian models
get_poisson_density_values = function(x,
                                      segment_id,
                                      cluster)
{
  # Poisson parameters - for the density
  poisson_params = get_poisson_parameters(x) %>%
    Rcongas:::idify() %>%
    dplyr::filter(segment_id == !!segment_id, cluster == !!cluster)
  
  if (nrow(poisson_params) == 0) {
    cli::cli_alert_warning(
      "Required density forcluster {.field {cluster}}, with no Poisson parameters associated!"
    )
    return(NULL)
  }
  
  # Counts data - to determine density range
  counts_data = Rcongas:::get_counts(x, normalise = TRUE) %>%
    Rcongas:::idify() %>%
    dplyr::filter(segment_id == !!segment_id, cluster == !!cluster)
  
  # Ranges
  min_x = counts_data %>% dplyr::pull(n) %>% min %>% floor()
  max_x = counts_data %>% dplyr::pull(n) %>% max %>% ceiling()
  
  min_x = ifelse(min_x > 10, min_x - 10, 0)
  max_x = max_x + 10
  
  # Mixture component likelihood
  lambda  = poisson_params$lambda
  pi = Rcongas::get_clusters_size(x, normalised = TRUE)[cluster]
  
  data.frame(
    segment_id = segment_id,
    cluster = paste(cluster),
    x = seq(min_x, max_x, by = (max_x - min_x) / 100) %>% round,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(y = dpois(x, lambda) * pi) %>%
    as_tibble()
}



get_gaussian_density_values = function(x,
                                       segment_id,
                                       cluster)
{
  # Gaussian parameters - for the density
  gaussian_params = Rcongas:::get_gaussian_parameters(x) %>%
    Rcongas:::idify() %>%
    dplyr::filter(segment_id == !!segment_id, cluster == !!cluster)
  
  # Counts data - to determine density range
  counts_data = Rcongas:::get_counts(x, normalise = TRUE, sum_denominator = FALSE) %>%
    Rcongas:::idify() %>%
    dplyr::filter(segment_id == !!segment_id, cluster == !!cluster)
  
  # Ranges
  min_x = counts_data %>% dplyr::pull(n) %>% min
  max_x = counts_data %>% dplyr::pull(n) %>% max
  
  min_x = min_x - 0.2
  max_x = max_x + 0.2
  
  # Mixture component likelihood
  mean  = gaussian_params$mean
  sd = gaussian_params$sd
  pi = Rcongas::get_clusters_size(x, normalised = TRUE)[cluster]
  
  data.frame(
    segment_id = segment_id,
    cluster = paste(cluster),
    x = seq(min_x, max_x, by = (max_x - min_x) / 100),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(y = dnorm(x, mean, sd) * pi) %>%
    as_tibble()
}


# OLD
#
# plot_single_segment = function(x, segment, sum_denominator)
# {
#   # Counts data
#   counts_data = Rcongas:::get_counts(x, normalise = TRUE,  sum_denominator = sum_denominator) %>%
#     Rcongas:::idify() %>%
#     dplyr::filter(segment_id == segment)
#
#   # We need to make some adjustments to get the right binning, scaling etc
#   if (is_gaussian(x)) {
#     nbins = 30
#     binsize = counts_data %>%
#       dplyr::group_by(segment_id) %>%
#       dplyr::summarise(min = min(n),
#                        max = max(n),
#                        .groups = 'keep') %>%
#       dplyr::mutate(binsize = (max - min) / nbins) %>%
#       dplyr::pull(binsize)
#   }
#   else
#   {
#     nbins = 100
#     binsize = 1
#   }
#
#   if (!is_gaussian(x)) {
#     # Poisson likelihood
#
#     # Poisson parameters
#     clusters = Rcongas::get_clusters_size(x) %>% names()
#
#     # density_points = lapply(clusters,
#     #                         get_poisson_density_values,
#     #                         x = x,
#     #                         segment_id = segment)
#
#     # This makes a progress bar
#     density_points = easypar::run(
#       FUN = function(cluster) {
#         get_poisson_density_values(x = x,
#                                    segment_id = segment,
#                                    cluster = cluster)
#       },
#       lapply(clusters, list),
#       parallel = FALSE
#     )
#
#     density_points = Reduce(dplyr::bind_rows, density_points)
#
#     # Scale by binsize
#     density_points$y = density_points$y * binsize
#
#   } else {
#     # Gaussian parameters
#     clusters = Rcongas::get_clusters_size(x) %>% names()
#
#     density_points = lapply(clusters,
#                             get_gaussian_density_values,
#                             x = x,
#                             segment_id = segment)
#
#     density_points = Reduce(dplyr::bind_rows, density_points)
#   }
#
#   # Coloring
#   clusters_colors = get_clusters_colors(counts_data$cluster)
#
#   # Plot
#   counts_data %>%
#     ggplot(aes(
#       n,
#       y = ..count.. / sum(..count..) * binsize,
#       fill = cluster,
#       group = cluster
#     )) +
#     # geom_histogram(bins = nbins) +
#     facet_wrap( ~ segment_id, ncol = 2, scales = 'free') +
#     geom_point(data = density_points,
#                aes(x = x, y = y, color = cluster),
#                inherit.aes = FALSE) +
#     geom_line(data = density_points,
#               aes(x = x, y = y, color = cluster),
#               inherit.aes = FALSE) +
#     CNAqc:::my_ggplot_theme() +
#     scale_fill_manual(values = sapply(clusters_colors, Rcongas:::lighten, factor = .5)) +
#     scale_color_manual(values = clusters_colors) +
#     labs(title = "Poisson counts") +
#     labs(x = "RNA counts", y = "Density") +
#     coord_cartesian(clip = 'off')
#
#   counts_data %>%
#     ggplot(aes(
#       n,
#       y = ..count.. / sum(..count..) / binsize,
#       fill = cluster,
#       group = cluster
#     )) +
#     geom_histogram(bins = nbins) +
#     facet_wrap( ~ segment_id, ncol = 2, scales = 'free') +
#     geom_point(data = density_points,
#                aes(x = x, y = y, color = cluster),
#                inherit.aes = FALSE) +
#     geom_line(data = density_points,
#               aes(x = x, y = y, color = cluster),
#               inherit.aes = FALSE) +
#     CNAqc:::my_ggplot_theme() +
#     scale_fill_manual(values = sapply(clusters_colors, Rcongas:::lighten, factor = .5)) +
#     scale_color_manual(values = clusters_colors) +
#     labs(title = "Poisson counts") +
#     labs(x = "RNA counts", y = "Density") +
#     coord_cartesian(clip = 'off')
#
# }

# x = Rcongas::congas_example
# plot_segment_density(x, get_segment_ids(x, highlight = TRUE, alpha = 0.05)[1])
