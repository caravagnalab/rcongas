#' Plot fits.
#'
#' @description This function is like \code{\link{plot_data}}, but shows information about a fit model.
#'
#' As \code{\link{plot_data}}, this function uses a \code{what} parameter
#' to dispatch visualization to a number of internal functions. For the fits one can visualize:
#'
#' * (\code{what = "CNA"}) a genome-wide plot of Copy Number Alteration profiles (inferred) per cluster;
#' * (\code{what = "mixing_proportions"}) The normalized size of each cluster, per modality;
#' * (\code{what = "density"}) The density per cluster and segment, split by modality. By default this is shown
#' only for segments that differ for a segment CNA value among one of the inferred clusters.
#'
#' This function has the same logic of \code{\link{plot_data}} (with the ellipsis and parameters).
#'
#' @param x An object of class \code{rcongasplus}.
#' @param what Any of \code{"CNA"},  \code{"density"},  or \code{"mixing_proportions"}.
#' @param ... Parameters forwarded to the internal plotting functions.
#'
#' @return A \code{ggplot} plot.
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
plot_fit = function(x, what = 'CNA', ...)
{
  x %>% sanitize_obj()

  if(!('best_model' %in% (x %>% names))) {
    cli::cli_alert_danger("No fits available, returning empty plot.")

    return(ggplot())
  }

  if (what == 'CNA')
    return(x %>% plot_fit_CNA())

  if (what == 'density')
    return(x %>% plot_fit_density(...))

  if (what == 'mixing_proportions')
    return(x %>% plot_mixing_proportions(...))

  stop("Unrecognised 'what': use any of 'CNA', 'density' or 'plot_mixing_proportions'.")
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


plot_fit_density_aux <- function(x, segment_id){


  dplot <- plot_data_histogram(x, segments = segment_id)
  dplot()



}
