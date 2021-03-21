#' Plot fits
#'
#' @description General plot for fits, showing them as .....
#'
#' @param x
#' @param segments Can subset the plot only to certain segments, which are indexed
#' by the \code{segment_id} format in the tibble returned by function \code{get_input}
#' with \code{what = 'segmentation'} as parameter.
#' @param what Any of \code{"histogram"},  \code{"lineplot"},  or \code{"heatmap"}.
#'
#' @return A \code{ggplot} or \code{cowplot} figure, depending on the number of
#' modalities and plot required.
#'
#' @export
#'
#' @examples
plot_fit = function(x, ...)
{
  x %>% sanitize_obj()
  
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
  ymin = ymin + ymin * .15
  ymax = ymax + ymax * .15
  
  if (ymax < 4)
    ymax = 6
  
  # Plot template
  segments_plot = blank_genome(x$reference_genome) +
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
    )
  
  # Add segments
  segments_plot = segments_plot + geom_segment(
    data = CNAs,
    aes(
      x = from,
      xend = to,
      y = value,
      yend = value,
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
  
  # Normalisation factors
  norm_factors = get_input(x, what = 'normalisation')
  
  # RNA_data values
  what_RNA = get_input(x, what = 'data') %>%
    filter(modality == "RNA")
  
  if (nrow(what_RNA) > 0  && stats_data$rna_dtype != "G")
    what_RNA = normalise_modality(what_RNA, norm_factors %>% filter(modality == "RNA"))
  
  # ATAC_data values
  what_ATAC = get_input(x, what = 'data') %>%
    filter(modality == "ATAC")
  
  if (nrow(what_ATAC) > 0 && stats_data$atac_dtype != "G")
    what_ATAC = normalise_modality(what_ATAC, norm_factors %>% filter(modality == "ATAC"))
  
  what = bind_rows(what_RNA, what_ATAC) %>%
    filter(segment_id %in% CNAs)
  
  # Per cell clustering assignments
  clustering_assignments = get_fit(x, what = 'cluster_assignments') %>%
    select(-modality)
  
  what = what %>% left_join(clustering_assignments, by = 'cell')
  
  # Call
  what %>%
    ggplot(aes(value, fill = cluster)) +
    geom_histogram(bins = 50) +
    facet_grid(segment_id ~ modality, scales = 'free_x') +
    labs(title = x$description,
         subtitle = subtitle) +
    guides(fill = FALSE) +
    theme_linedraw(base_size = 9) +
    scale_fill_brewer(palette = 'Set1') +
    labs(x = "Input",
         y = 'Observations')
}

plot_mixing_proportions = function(x)
{
  clustering_assignments = get_fit(x, what = 'cluster_assignments') %>% 
    group_by(cluster) %>% 
    summarise(n = n()) %>% 
    mutate(prop = (100 * n/sum(n)) %>% round(2)) %>%
    arrange(desc(cluster)) %>%
    mutate(lab.ypos = cumsum(prop) - 0.5 * prop)
  
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
    guides(fill = guide_legend(''))
}
