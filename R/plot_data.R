#' Plot data.
#'
#' @description General plot for data, showing data as a histogram, a lineplot
#' or a heatmap. Discrete count data are normalised by input factors.
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
plot_data = function(x,
                     segments = get_input(x, what = 'segmentation') %>% pull(segment_id),
                     what = 'histogram',
                     ...)
{
  x %>% sanitize_obj()
  
  if (what == 'histogram')
    return(x %>% plot_data_histogram(segments = segments))
  
  if (what == 'lineplot')
    return(x %>% plot_data_lineplot(segments = segments, ...))
  
  if (what == 'heatmap')
    return(x %>% plot_data_heatmap(segments = segments))
  
  stop("Unrecognised 'what': use any of 'histogram', 'lineplot' or 'heatmap'.")
}


plot_data_histogram = function(x, segments)
{
  stats_data = stat(x, what = 'data')
  
  # Data stats
  subtitle = paste0(
    "RNA (",
    stats_data$ncells_RNA,
    " cells, ",
    stats_data$rna_genes,
    " genes) ATAC (",
    stats_data$ncells_ATAC,
    " cells, ",
    stats_data$atac_peaks,
    " peaks)"
  )
  
  # What to print, enhanced facets
  what_rna_lik = case_when(
    stats_data$rna_dtype == "NB" ~ "Negative Binomial",
    stats_data$rna_dtype == "P" ~ "Poisson",
    stats_data$rna_dtype == "G" ~ "Gaussian"
  )
  
  what_atac_lik = case_when(
    stats_data$atac_dtype == "NB" ~ "Negative Binomial",
    stats_data$atac_dtype == "P" ~ "Poisson",
    stats_data$atac_dtype == "G" ~ "Gaussian"
  )
  
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
    filter(segment_id %in% segments) %>%
    mutate(modality = case_when(
      modality == "RNA" ~ paste(modality, '(', what_rna_lik, ')'),
      modality == "ATAC" ~ paste(modality, '(', what_atac_lik, ')')
    ))
  
  # Call
  what %>%
    ggplot(aes(value, fill = modality)) +
    geom_histogram(bins = 50) +
    facet_grid(segment_id ~ modality, scales = 'free_x') +
    labs(title = x$description,
         subtitle = subtitle) +
    guides(fill = FALSE) +
    theme_linedraw(base_size = 9) +
    scale_fill_manual(values = modality_colors(what$modality %>% unique)) +
    labs(x = "Input",
         y = 'Observations')
}

plot_data_lineplot = function(x, segments, alpha)
{
  stats_data = stat(x, what = 'data')
  
  # Data stats
  subtitle = paste0(
    "RNA (",
    stats_data$ncells_RNA,
    " cells, ",
    stats_data$rna_genes,
    " genes) ATAC (",
    stats_data$ncells_ATAC,
    " cells, ",
    stats_data$atac_peaks,
    " peaks)"
  )
  
  # What to print, enhanced facets
  what_rna_lik = case_when(
    stats_data$rna_dtype == "NB" ~ "Negative Binomial",
    stats_data$rna_dtype == "P" ~ "Poisson",
    stats_data$rna_dtype == "G" ~ "Gaussian"
  )
  
  what_atac_lik = case_when(
    stats_data$atac_dtype == "NB" ~ "Negative Binomial",
    stats_data$atac_dtype == "P" ~ "Poisson",
    stats_data$atac_dtype == "G" ~ "Gaussian"
  )
  
  # Normalisation factors
  norm_factors = get_input(x, what = 'normalisation')
  
  # RNA_data values
  what_RNA = get_input(x, what = 'data') %>%
    filter(modality == "RNA") %>% 
    mutate(
      alpha = alpha,
      size = .5
    )
  
  if (nrow(what_RNA) > 0 && stats_data$rna_dtype != "G")
    what_RNA = normalise_modality(what_RNA, norm_factors %>% filter(modality == "RNA"))
  
  # ATAC_data values
  what_ATAC = get_input(x, what = 'data') %>%
    filter(modality == "ATAC") %>% 
    mutate(
      alpha = alpha,
      size = .5
    )
  
  if (nrow(what_ATAC) > 0 && stats_data$atac_dtype != "G")
    what_ATAC = normalise_modality(what_ATAC, norm_factors %>% filter(modality == "ATAC"))
  
  # Normalise observed values by number of events per segment to make them comparable
  if (nrow(what_RNA) > 0)
  {
    cli::cli_alert("Scaling RNA observed values by number of genes mapped per segment.")
    
    what_RNA = what_RNA %>% 
      left_join(
        x %>% get_input(what = 'segmentation') %>% select(segment_id, RNA_genes),
        by = 'segment_id'
      ) %>% 
      mutate(value = value/RNA_genes)
  }
  
  if (nrow(what_ATAC) > 0)
  {
    cli::cli_alert("Scaling ATAC observed values by number of peaks mapped per segment.")
    
    what_ATAC = what_ATAC %>% 
      left_join(
        x %>% get_input(what = 'segmentation') %>% select(segment_id, ATAC_peaks),
        by = 'segment_id'
      ) %>% 
      mutate(value = value/ATAC_peaks)
  }
  
  # Pool all modalities plus the segmentation, apply filters
  what = bind_rows(
    what_RNA %>% deidify(),
    what_ATAC %>% deidify(),
    get_input(x, what = 'segmentation') %>%
      mutate(modality = " Input segmentation") %>% # add space for factors ordering
      rename(value = copies) %>% 
      mutate(
        alpha = 1,
        size = 2
      )
  ) %>%
    idify() %>%
    filter(segment_id %in% segments) %>%
    mutate(modality = case_when(
      modality == "RNA" ~ paste(modality, '(', what_rna_lik, ')'),
      modality == "ATAC" ~ paste(modality, '(', what_atac_lik, ')'),
      TRUE ~ modality
    ))
  
  # Shift CNA coordinates to absolute
  what = CNAqc:::relative_to_absolute_coordinates(x,
                                                  what)
  
  # Genome plot call
  blank_genome(ref = x$reference_genome,
                       chromosomes = what$chr %>% unique()) +
    geom_segment(
      data = what,
      aes(
        x = from,
        xend = to,
        y = value,
        yend =  value,
        color = modality,
        alpha = alpha,
        size = size
      ),
      inherit.aes = FALSE
    ) +
    facet_wrap( ~ modality, ncol = 1, scales = 'free_y') +
    guides(color = FALSE) +
    scale_color_manual(values = modality_colors(what$modality %>% unique)) +
    labs(title = x$description,
         subtitle = subtitle) +
    theme_linedraw(base_size = 9) +
    labs(x = "Chromosomes",
         y = 'Input',
         caption = "Per segment values are normalised by number of mapped RNA genes/ATAC peaks.") +
    guides(
      alpha = FALSE,
      size = FALSE
    )
}

plot_data_heatmap = function(x, segments)
{
  stats_data = stat(x, what = 'data')
  
  # Data stats
  subtitle_RNA = paste0("RNA (",
                        stats_data$ncells_RNA,
                        " cells, ",
                        stats_data$rna_genes,
                        " genes)")
  
  subtitle_ATAC = paste0("ATAC (",
                         stats_data$ncells_ATAC,
                         " cells, ",
                         stats_data$atac_peaks,
                         " peaks)")
  
  # What to print, enhanced facets
  what_rna_lik = case_when(
    stats_data$rna_dtype == "NB" ~ "Negative Binomial",
    stats_data$rna_dtype == "P" ~ "Poisson",
    stats_data$rna_dtype == "G" ~ "Gaussian"
  )
  
  what_atac_lik = case_when(
    stats_data$atac_dtype == "NB" ~ "Negative Binomial",
    stats_data$atac_dtype == "P" ~ "Poisson",
    stats_data$atac_dtype == "G" ~ "Gaussian"
  )
  
  # Normalisation factors
  norm_factors = get_input(x, what = 'normalisation')
  
  # RNA_data values
  what_RNA = get_input(x, what = 'data') %>%
    filter(modality == "RNA")
  
  if (nrow(what_RNA) > 0 && stats_data$rna_dtype != "G")
    what_RNA = normalise_modality(what_RNA, norm_factors %>% filter(modality == "RNA"))
  
  # ATAC_data values
  what_ATAC = get_input(x, what = 'data') %>%
    filter(modality == "ATAC")
  
  if (nrow(what_ATAC) > 0 && stats_data$atac_dtype != "G")
    what_ATAC = normalise_modality(what_ATAC, norm_factors %>% filter(modality == "ATAC"))
  
  # Pool all modalities plus the segmentation, apply filters
  what = bind_rows(what_RNA %>% deidify(),
                   what_ATAC %>% deidify()) %>%
    idify() %>%
    filter(segment_id %in% segments) %>%
    mutate(modality = case_when(
      modality == "RNA" ~ paste(modality, '(', what_rna_lik, ')'),
      modality == "ATAC" ~ paste(modality, '(', what_atac_lik, ')'),
      TRUE ~ modality
    ))
  
  # Genome plot call
  RNA_plot = ATAC_plot = NULL
  
  if (nrow(what_RNA) > 0)
    # RNA
  {
    RNA_plot = ggplot(what_RNA %>% filter(segment_id %in% segments)) +
      geom_tile(aes(x = segment_id, y = cell, fill = value)) +
      theme_linedraw(base_size = 9) +
      labs(x = "Segments",
           y = subtitle_RNA) +
      theme(axis.text.y = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
      ) +
      guides(fill = guide_colorbar(paste('RNA (', what_rna_lik, ')'),
                                   barheight = unit(3, 'cm')))
    
    if (stats_data$rna_dtype == "G")
      RNA_plot = RNA_plot +
        scale_fill_gradient2(low = "steelblue", high = 'indianred3')
    else
      RNA_plot = RNA_plot +
        scale_fill_distiller(palette = 'GnBu', direction = 1)
  }
  
  if (nrow(what_ATAC) > 0)
    # RNA
  {
    ATAC_plot = ggplot(what_ATAC %>% filter(segment_id %in% segments)) +
      geom_tile(aes(x = segment_id, y = cell, fill = value)) +
      theme_linedraw(base_size = 9) +
      labs(x = "Segments",
           y = subtitle_ATAC) +
      theme(axis.text.y = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
      ) +
      guides(fill = guide_colorbar(paste('ATAC (', what_atac_lik, ')'),
                                   barheight = unit(3, 'cm')))
    
    if (stats_data$atac_dtype == "G")
      ATAC_plot = ATAC_plot +
        scale_fill_gradient2(low = "steelblue", high = 'indianred3')
    else
      ATAC_plot = ATAC_plot + scale_fill_distiller(palette = 'GnBu', direction = 1)
  }
  
  # Figure assembly
  if (stats_data$nmodalities == 1)
  {
    if ("RNA" %in% stats_data$modalities)
      return(RNA_plot)
    if ("ATAC" %in% stats_data$modalities)
      return(ATAC_plot)
  }
  else
  {
    # Proportional to number of cells per assay
    rel_height = c(stats_data$ncells_ATAC, stats_data$ncells_RNA) /
      (stats_data$ncells_RNA + stats_data$ncells_ATAC)
    
    # Place ATAC on top of RNA, remove ATAC segment ids, move labels
    figure =
      cowplot::plot_grid(
        ATAC_plot + theme(axis.text.x = element_blank()) + labs(x = NULL),
        RNA_plot,
        rel_heights = rel_height,
        ncol = 1,
        align = 'v',
        axis = 'lr'
      )
    
    return(figure)
  }
}


blank_genome = function(ref = "GRCh38",
                        chromosomes = paste0("chr", c(1:22, "X", "Y")),
                        label_chr = -0.5,
                        cex = 1)
{
  reference_coordinates = CNAqc:::get_reference(ref) %>%
    filter(chr %in% chromosomes)
  
  low = min(reference_coordinates$from)
  upp = max(reference_coordinates$to)
  
  pl = ggplot(reference_coordinates) +
    geom_rect(
      aes(
        xmin = centromerStart,
        xmax = centromerEnd,
        ymin = 0,
        ymax = Inf
      ),
      alpha = 0.3,
      colour = "gainsboro"
    ) +
    geom_segment(
      data = reference_coordinates,
      aes(
        x = from,
        xend = from,
        y = 0,
        yend = Inf
      ),
      size = 0.1,
      color = "black",
      linetype = 8
    ) + 
    labs(x = "Chromosome",
         y = "Value") + ggpubr::rotate_y_text() +
    scale_x_continuous(
      breaks = c(0, reference_coordinates$from,
                 upp),
      labels = c(
        "",
        gsub(
          pattern = "chr",
          replacement = "",
          reference_coordinates$chr
        ),
        ""
      )
    )
  return(pl)
}
