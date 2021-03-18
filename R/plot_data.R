plot_data = function(x,
                     segments = get_input(x, what = 'segmentation') %>% pull(segment_id))
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
  
  what = get_input(x, what = 'data') %>%
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
    ggsci::scale_color_nejm() +
    labs(x = "Input",
         y = 'Observations')
}

# plot_data(x,
#           c("chr15:0:102531392", "chr16:0:90354753",  "chr17:0:81195210"))
