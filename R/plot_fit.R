#' Plot fits.
#'
#' @description This function is like \code{\link{plot_data}}, but shows information about a fit model.
#'
#' As \code{\link{plot_data}}, this function uses a \code{what} parameter
#' to dispatch visualization to a number of internal functions. For the fits one can visualize:
#'
#' * (\code{what = "CNA"}) A genome-wide plot of Copy Number Alteration profiles (inferred) per cluster;
#' * (\code{what = "mixing_proportions"}) The normalized size of each cluster, per modality;
#' * (\code{what = "density"}) The density per cluster and segment, split by modality. By default this is shown
#' only for segments that differ for a segment CNA value among one of the inferred clusters;
#' * (\code{what = "heatmap"}) The same heatmap plot from \code{\link{plot_data}} with
#' \code{what = "heatmap"}, where rows are sorted by cluster and cluster annotations reported;
#' * (\code{what = "scores"}) The scores used for model selection;
#' * (\code{what = "posterior_CNA"}) The posterior probability for CNA values;
#'
#' This function has the same logic of \code{\link{plot_data}} with respect to the ellipsis and
#' input parameters.
#'
#' @param x An object of class \code{rcongasplus}.
#' @param what Any of \code{"CNA"},  \code{"density"}, \code{"mixing_proportions"}, \code{"heatmap"}
#' or \code{"scores"}.
#' @param ... Parameters forwarded to the internal plotting functions.
#'
#' @return A \code{ggplot} plot, or a more complex \code{cowplot} figure.
#'
#' @export
#'
#' @examples
#' data("example_object")
#'
#' # Genome-wide segments plo
#' plot_fit(example_object, what = 'CNA')
#'
#' # Density plot
#' plot_fit(example_object, what = 'density')
#'
#' # Mixing proportions
#' plot_fit(example_object, what = 'mixing_proportions')
#'
#' # Fit heatmap
#' plot_fit(example_object, what = 'heatmap')
#'
#' # Scores for model selection
#' plot_fit(example_object, what = 'scores')
#'
#' # Posterior for CNAs
#' plot_fit(example_object, what = 'posterior_CNA')
plot_fit = function(x, what = 'CNA', ...)
{
  x %>% sanitize_obj()

  if(!('best_fit' %in% (x %>% names))) {
    cli::cli_alert_danger("No fits available, returning empty plot.")

    return(ggplot())
  }

  if (what == 'CNA')
    return(x %>% plot_fit_CNA())

  if (what == 'density')
    return(x %>% plot_fit_density(...))

  if (what == 'mixing_proportions')
    return(x %>% plot_mixing_proportions(...))

  if (what == 'heatmap')
    return(x %>% plot_fit_heatmap(...))

  if (what == 'scores')
    return(x %>% plot_fit_scores())

  if (what == 'posterior_CNA')
    return(x %>% plot_fit_posterior_CNA())

  what_supported = c("CNA", "density", "mixing_proportions", "heatmap", 'scores', "posterior_CNA")

  stop(paste0("Unrecognised 'what': use any of ", what_supported %>% paste(collapse = ', '), '.'))
}

plot_fit_CNA = function(x)
{
  # Fit CNAs in absolute coordinates
  CNAs = get_fit(x, 'CNA') %>%
    deidify() %>%
    CNAqc:::relative_to_absolute_coordinates(x = list(reference_genome = x$reference_genome))

  # Shift equal values
  CNAs = CNAs %>%
    group_by(chr, from, to, value) %>%
    mutate(row_id = row_number(),
           grp_size = n()) %>%
    mutate(value = ifelse(grp_size > 1,
                          value + (row_id - 1) * 0.05,
                          value)) %>%
    ungroup()

  # Segments to highlight
  specials = CNAs %>%
    filter(grp_size != CNAs$cluster %>% unique %>% length) %>%
    distinct(chr, from, to)

  # Determine plot bounds
  ymin = CNAs$value %>% min
  ymax = CNAs$value %>% max
  ymin = ymin - ymin * .15
  ymax = ymax + ymax * .15

  if (ymax < 4)
    ymax = 4

  # Plot template
  segments_plot = blank_genome(x$reference_genome) + if(nrow(specials > 0)){
    geom_rect(
      data = specials ,
      aes(
        xmin = from,
        xmax = to,
        ymin = ymin,
        ymax = ymax
      ),
      fill = "indianred3",
      alpha = .2
    )} else {
      geom_blank()
    }

  # Add segments
  segments_plot = segments_plot + geom_segment(
    data = CNAs,
    aes(
      x = from,
      xend = to,
      y = value ,
      yend = value ,
      colour = cluster
    ),
    size = 1
  ) +
    theme_linedraw(base_size = 9) +
    scale_color_brewer(palette = 'Set1') +
    theme(legend.position = 'bottom') +
    ylim(ymin, ymax) +
    labs(title = x$description, y = "Ploidy")

  return(segments_plot)
}

plot_fit_density = function(x, highlights = TRUE)
{
  # CNAs where these are different
  CNAs = get_fit(x, 'CNA')

  if (highlights)
  {
    nclusters = CNAs$cluster %>% unique %>% length()

    CNAs = CNAs %>%
      group_by(segment_id, value) %>%
      mutate(grp_size = n()) %>%
      filter(grp_size != nclusters) %>%
      pull(segment_id) %>%
      unique

    cli::cli_alert("Plotting segments where different CNAs are present: {.field {CNAs}}.")
  }
  else
  {
    CNAs = CNAs %>%
      pull(segment_id) %>%
      unique

    cli::cli_alert("Showing all segments (this plot can be large).")
  }

  #plots <- plot_data_histogram(x, CNAs)


  plots <- plot_data_histogram(x, segment_id)

  # Per cell clustering assignments
  clustering_assignments = get_fit(x, what = 'cluster_assignments') %>%
    select(-modality)

  plots$data$modality <- clustering_assignments$

  what = what %>% left_join(clustering_assignments, by = 'cell')

  # Call
  what %>%
    ggplot(aes(value, fill = cluster)) +
    geom_histogram(bins = 50) +
    facet_grid(segment_id ~ modality, scales = 'free_x') +
    labs(title = x$description
         ) +
    guides(fill = FALSE) +
    theme_linedraw(base_size = 9) +
    scale_fill_brewer(palette = 'Set1') +
    labs(x = "Input",
         y = 'Observations')
}

plot_mixing_proportions = function(x)
{
  clustering_assignments =
    rna_clustering_assignments =
    atac_clustering_assignments = NULL

  if(x %>% has_rna)
    rna_clustering_assignments = get_fit(x, what = 'cluster_assignments') %>%
      filter(modality == 'RNA') %>%
      group_by(cluster) %>%
      summarise(n = n()) %>%
      mutate(prop = (100 * n/sum(n)) %>% round(2)) %>%
      arrange(desc(cluster)) %>%
      mutate(
        lab.ypos = cumsum(prop) - 0.5 * prop,
        modality = 'RNA'
        )

  if(x %>% has_atac)
    atac_clustering_assignments = get_fit(x, what = 'cluster_assignments') %>%
      filter(modality == 'ATAC') %>%
      group_by(cluster) %>%
      summarise(n = n()) %>%
      mutate(prop = (100 * n/sum(n)) %>% round(2)) %>%
      arrange(desc(cluster)) %>%
      mutate(
        lab.ypos = cumsum(prop) - 0.5 * prop,
        modality = 'ATAC'
        )

  clustering_assignments = bind_rows(
    rna_clustering_assignments,
    atac_clustering_assignments
  )

  ggplot(clustering_assignments,
         aes(x = "", y = prop, fill = cluster)) +
    geom_bar(width = 1,
             stat = "identity",
             color = "white") +
    coord_polar("y", start = 0) +
    scale_fill_brewer(palette = 'Set1') +
    geom_text(aes(
      y = lab.ypos,
      label = paste0(cluster, ' (', prop, '%\nn = ', n, ')'),
    ),       size = 2, color = "white") +
    theme_linedraw(base_size = 9) +
    labs(x = NULL, y = NULL) +
    theme(
      axis.text.x = element_blank(),
      panel.border = element_blank(),
      legend.position = 'right'
    ) +
    labs(
      title = "Mixing proportions"
    ) +
    guides(fill = guide_legend('')) +
    facet_wrap(~modality, nrow = 1)
}

plot_fit_density = function(x, highlights = TRUE)
{
  # CNAs where these are different
  CNAs = get_fit(x, 'CNA')

  if (highlights)
  {
    nclusters = CNAs$cluster %>% unique %>% length()

    CNAs = CNAs %>%
      group_by(segment_id, value) %>%
      mutate(grp_size = n()) %>%
      filter(grp_size != nclusters) %>%
      pull(segment_id) %>%
      unique

    cli::cli_alert("Plotting segments where different CNAs are present: {.field {CNAs}}.")
  }
  else
  {
    CNAs = CNAs %>%
      pull(segment_id) %>%
      unique

    cli::cli_alert("Showing all segments (this plot can be large).")
  }

  #CNAs <-  segment_id

  if(length(CNAs) == 0){
    return(ggplot()+ geom_blank())
  }

  plots <- plot_data_histogram(x, CNAs)

  data_hist <-  plots$data


  # Per cell clustering assignments
  clustering_assignments = get_fit(x, what = 'cluster_assignments') %>%
    select(-modality)

  data_hist$modality <- sapply(data_hist$modality %>%  strsplit(" "), function(y) y[1])

  data_hist <-  dplyr::left_join(data_hist, clustering_assignments)

  mixing_props <- get_mixing_proportions(x)

  densities <- assemble_likelihood_tibble(x, CNAs)
  colnames(densities)[c(3,5)] <-  c("cluster", "segment_id")
  densities <- dplyr::left_join(densities,mixing_props) %>%  mutate(value = value * (mixing + 1e-3))

  CNA <-  get_CNA(x)

  colnames(CNA)[3] <-  "CNA"

  data_hist <-  dplyr::left_join(data_hist, CNA)

  ret <- lapply(CNAs, function(s) plot_fit_density_aux(data_hist, densities,s, x))

  return(ret)


}


plot_fit_density_aux <-  function(df, densities, segment, x){

  df <-  df %>%  filter(segment_id == segment)
  densities <-  densities %>%  filter(segment_id == segment)


  cols <- RColorBrewer::brewer.pal(9, "Set1")
  CNA <- df %>% arrange(cluster) %>%  select(cluster,CNA) %>%  unique() %>%  pull(CNA)
  clts <- sort(df$cluster %>%  unique())


  p1 <- ggplot( ) +
    geom_histogram(aes(x = value, fill = factor(cluster, levels = clts)
    ),bins = 50,  data = df, color = "black", alpha = 0.4) +
    facet_wrap(segment_id ~ modality, scales = 'free') +
    guides() +
    theme_linedraw(base_size = 9) +
    scale_fill_manual("CN value", values = cols, labels = CNA, drop=FALSE) +
    labs(x = "Input",
         y = 'Observations') +
    theme(strip.text.y.right = element_text(angle = 0), , axis.text.y = element_blank(), axis.ticks.y = element_blank())

  p2 <-ggplot() +
    geom_line(aes(x = X, y = value, color = factor(cluster, levels = clts)), data = densities, size = 0.8) +
    geom_point(aes(x = X, y = value, color = factor(cluster, levels = clts)), data = densities, size = 0.3) +
    facet_wrap(segment_id ~ modality, scales = 'free') +
    labs(title = x$description, subtitle = "Cell cluster assignments and densities") +
    guides(fill = FALSE) +
    theme_linedraw(base_size = 9) +
    scale_color_manual("Clusters",values = cols, drop=FALSE) +
    labs(x = "Input",
         y = 'Density') +
    theme(strip.text.y.right = element_text(angle = 0), axis.text.y = element_blank(), axis.ticks.y = element_blank())

  cowplot::plot_grid(p2,p1, align = "hv", axis = "lrtb", ncol = 1)



}


plot_fit_heatmap = function(x, segments = get_input(x, what = 'segmentation') %>% pull(segment_id))
{
  # devtools::load_all()

  # Create plots with internal function
  data_plot = plot_data_heatmap(x, segments = segments)

  # Get assignments and upgrade plot row ordering
  rna_clustering_assignments = atac_clustering_assignments = NULL
  rna_plot = atac_plot = NULL

  if(x %>% has_rna)
  {
    rna_clustering_assignments = get_fit(x, what = 'cluster_assignments') %>%
      filter(modality == 'RNA') %>%
      arrange(cluster)

    c_size = rna_clustering_assignments %>% group_by(cluster) %>% summarise(n = n())

    rna_plot = data_plot$plots[[1]] +
      scale_y_discrete(limits = rna_clustering_assignments$cell) +
      geom_hline(
        yintercept = cumsum(c_size$n)[-nrow(c_size)],
        size = .6,
        linetype = 'dashed'
      )  +
      geom_point(data = rna_clustering_assignments,
               aes(x = 0, y = cell, color = cluster),
               size = 2) +
      scale_color_brewer(palette = "Set1")
  }

  if(x %>% has_atac)
  {
    atac_clustering_assignments = get_fit(x, what = 'cluster_assignments') %>%
      filter(modality == 'ATAC') %>%
      arrange(cluster)

    c_size = atac_clustering_assignments %>% group_by(cluster) %>% summarise(n = n())

    atac_plot = data_plot$plots[[2]] +
      scale_y_discrete(limits = atac_clustering_assignments$cell) +
      geom_hline(
        yintercept = cumsum(c_size$n)[-nrow(c_size)],
        size = .6,
        linetype = 'dashed'
      )  +
      geom_point(data = atac_clustering_assignments,
                 aes(x = 0, y = cell, color = cluster),
                 size = 2) +
      scale_color_brewer(palette = "Set1")
  }

  stats_data = stat(x)

  if(length(stats_data$modalities) > 1)
  {
    # Proportional to number of cells per assay
    rel_height = c(stats_data$ncells_ATAC, stats_data$ncells_RNA) /
      (stats_data$ncells_RNA + stats_data$ncells_ATAC)

    f = cowplot::plot_grid(
      atac_plot + theme(axis.text.x = element_blank()) + labs(x = NULL),
      rna_plot,
      rel_heights = rel_height,
      ncol = 1,
      align = 'v',
      axis = 'lr'
    )

    return(f)
  }

  if(x %>% has_rna) return(rna_plot)
  if(x %>% has_atac) return(atac_plot)
}


plot_fit_scores = function(x)
{
scores = reshape2::melt(
  x$model_selection %>%
    select(-n_observations),
  id = c('K', 'lambda')
)
scores$K = as.factor(scores$K)
IC_best = scores %>%
  group_by(variable) %>%
  filter(value == min(value)) %>%
  filter(variable != 'entropy')

H_best = scores %>%
  group_by(variable) %>%
  filter(value == max(value)) %>%
  filter(variable == 'entropy')
scores %>%
  ggplot(aes(x = lambda, y = value, color = K)) +
  geom_line() +
  geom_point() +
  facet_wrap(~variable, scales = 'free') +
  geom_point(
    data = H_best,
    color = 'red'
  ) +
  geom_point(
    data = IC_best,
    color = 'red'
  ) +
  labs(title = x$description) +
  theme_linedraw(base_size = 9) +
  scale_color_brewer(palette="Dark2") +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE, digits = 3))

}

plot_fit_posterior_CNA = function(x)
{
  x = get_fit(x, what = 'posterior_CNA')

  x %>%
    ggplot(aes(y = segment_id, x = value, fill = probability)) +
    geom_tile() +
    scale_fill_distiller(palette = 'Spectral') +
    theme_linedraw(base_size = 9) +
    theme(legend.position = 'bottom') +
    labs(
      x = "Copy Number",
      y = "Segment"
    ) +
    guides(fill = guide_colorbar("Posterior", barheight = .5)) +
    scale_y_discrete(limits = gtools::mixedsort(x$segment_id %>% unique) %>% rev) +
    facet_wrap(~cluster)

}


