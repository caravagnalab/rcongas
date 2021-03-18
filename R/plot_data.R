#' Plot data.
#'
#' @description General plot for data, accessing getter function to return
#' segment_factors, CNAs (point estimates), posterior CNA (distribution),
#' mixing proportions, clustering assignments or the latent responsibilities.
#'
#' @param x
#' @param segments Can subset the plot only to certain segments, which are indexed
#' by the \code{segment_id} format in the tibble returned by function \code{get_input}
#' with \code{what = 'segmentation'} as parameter.
#' @param what Any of \code{"segment_factors"}, \code{"CNA"},  \code{"posterior_CNA"},
#'  \code{"mixing_proportions"},  \code{"cluster_assignments"},  or \code{"z_nk"}.
#'
#' @return A tibble.
#' @export
#'
#' @examples
plot_data = function(x, 
                     segments = get_input(x, what = 'segmentation') %>% pull(segment_id), 
                     what = 'histogram')
{
  x %>% sanitize_obj()
  
  if (what == 'histogram')
    return(x %>% plot_data_histogram(segments = segments))

  if (what == 'lineplot')
    return(x %>% plot_data_lineplot(segments = segments))
  
  stop("Unrecognised 'what': use any of 'histogram', 'lineplot'")
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

  if(nrow(what_RNA) > 0 & stats_data$rna_dtype != "G")  
    what_RNA = normalise_modality(what_RNA, norm_factors %>% filter(modality == "RNA"))
    
  # ATAC_data values
  what_ATAC = get_input(x, what = 'data') %>% 
    filter(modality == "ATAC")
  
  if(nrow(what_ATAC) > 0 & stats_data$atac_dtype != "G")  
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

plot_data_lineplot = function(x, segments)
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
  
  # RNA_data values
  what_RNA = get_input(x, what = 'data') %>% 
    filter(modality == "RNA")
  
  if(nrow(what_RNA) > 0 & stats_data$rna_dtype != "G")  
    what_RNA = normalise_modality(what_RNA, norm_factors %>% filter(modality == "RNA"))
  
  # ATAC_data values
  what_ATAC = get_input(x, what = 'data') %>% 
    filter(modality == "ATAC")
  
  if(nrow(what_ATAC) > 0 & stats_data$atac_dtype != "G")  
    what_ATAC = normalise_modality(what_ATAC, norm_factors %>% filter(modality == "ATAC"))
  
  # Pool all modalities plus the segmentation, apply filters
  what = bind_rows(
    what_RNA %>% deidify(),
    what_ATAC %>% deidify(),
    get_input(x, what = 'segmentation') %>%
      mutate(modality = " Input segmentation") %>% # add space for factors ordering
      rename(value = copies) 
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
  CNAqc:::blank_genome(
    ref = x$reference_genome,
    chromosomes = what$chr %>% unique()
  ) +
    geom_segment(data = what,
                 aes(
                   x = from,
                   xend = to,
                   y = value,
                   yend =  value,
                   color = modality
                 ),
                 inherit.aes = FALSE) +
    facet_wrap( ~ modality, ncol = 1, scales = 'free_y') +
    guides(color = FALSE) +
    scale_color_manual(values = modality_colors(what$modality %>% unique)) +
    labs(title = x$description,
         subtitle = subtitle) +
    theme_linedraw(base_size = 9) +
    labs(x = "Chromosomes",
         y = 'Input')
}
