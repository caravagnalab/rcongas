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
plot_segment_density = function(x,
                                segments_ids,
                                group1 = 1,
                                group2 = 2,
                                ...)
{
  plot_single_segment = function(x, segment)
  {
    # Counts data
    counts_data = Rcongas:::get_counts(x, normalise = TRUE) %>%
      Rcongas:::idify() %>%
      dplyr::filter(segment_id == segment)

    # Poisson p-value from the usual test
    test_pvalue = Rcongas:::get_segment_test_counts(x, group1 = group1, group2 = group2, ...) %>%
      Rcongas:::idify() %>%
      dplyr::filter(segment_id == segment)

    # Poisson parameters - used to plot the density
    poisson_params = Rcongas:::get_poisson_parameters(x)

    all_densities = Reduce(
      dplyr::bind_rows,
      lapply(
        segment,
        get_poisson_density_values,
        counts_data = counts_data,
        poisson_params = poisson_params
      )
    ) %>% as_tibble()

    # ggplot(all_densities, aes(x = x, y = y, color = cluster)) + geom_line()

    # We need to make some adjustments to get the right binning scaling etc
    nbins = 100

    range_data = counts_data %>%
      dplyr::group_by(segment_id) %>%
      dplyr::summarise(min = min(n), max = max(n), .groups = 'keep') %>%
      dplyr::mutate(binsize = (max - min)/nbins)

    all_densities = all_densities %>%
      dplyr::left_join(range_data, by = 'segment_id') %>%
      dplyr::mutate(y = y * binsize)

    # Coloring
    clusters_colors = Rcongas:::get_clusters_colors(counts_data$cluster)

    # Plot
    counts_data %>%
      # filter(segment_id == "chr15:1:102600000" ) %>%
      # ggplot(aes(n, y = ..count.., fill = cluster)) +
      ggplot(aes(n, y = ..count../sum(..count..), fill = cluster, group = cluster)) +
      geom_histogram(bins = nbins) +
      facet_wrap(~ segment_id, ncol = 2, scales = 'free') +
      geom_point(data = all_densities,
                 aes(x = x, y = y, color = cluster),
                 # pch = 21,
                 inherit.aes = FALSE) +
      geom_line(
        data = all_densities,
        aes(x = x, y = y, color = cluster),
        inherit.aes = FALSE
      ) +
      CNAqc:::my_ggplot_theme() +
      scale_fill_manual(values = sapply(clusters_colors, lighten, factor = .5)) +
      scale_color_manual(values = Rcongas:::get_clusters_colors(counts_data$cluster)) +
      labs(title = "Poisson counts") +
      labs(x = "RNA counts", y = "Density") +
      geom_text(
        data = test_pvalue,
        aes(
          x = 0,
          y = Inf,
          label = p_value_format(p)
        ),
        inherit.aes = FALSE,
        size = 2.5,
        hjust = 0,
        vjust = 1.3
      ) +
      coord_cartesian(clip = 'off')

  }


  plots = lapply(segments_ids, plot_single_segment, x = x)
  names(plots) = segments_ids

  return(plots)
}

get_poisson_density_values = function(counts_data,
                                      segment_id,
                                      poisson_params)
{
  min_x = counts_data %>% dplyr::filter(segment_id == !!segment_id) %>% dplyr::pull(n) %>% min
  max_x = counts_data %>% dplyr::filter(segment_id == !!segment_id) %>% dplyr::pull(n) %>% max

  min_x = min_x - 10
  max_x = max_x + 10

  segment_params = poisson_params %>% dplyr::filter(segment_id == !!segment_id)

  all_density = NULL

  # Mixture component likelihood
  for (i in 1:nrow(segment_params))
  {
    lambda_i  = segment_params$lambda[i]
    pi_i = Rcongas::get_clusters_size(x, normalised = TRUE)[segment_params$cluster[i]]

    df_i = data.frame(
      segment_id = segment_params$segment_id[i],
      cluster = segment_params$cluster[i] %>% paste,
      x = seq(min_x, max_x, by = (max_x - min_x) / 100) %>% round,
      stringsAsFactors = FALSE
    ) %>%
      dplyr::mutate(y = dpois(x, lambda_i) * pi_i)

    all_density = dplyr::bind_rows(all_density, df_i)
  }

  all_density
}


