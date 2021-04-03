







#' Filter cells by upper and lower quantile of their normalisation factors.
#'
#' @param x A tibble with a cell and normalisation_factor column.
#' @param lower_quantile The minimum quantile to determine cuts. 
#' @param upper_quantile The maximum quantile to determine cuts. 
#'
#' @return The input data with normalisation factors in between the quantiles range.
#' @export
#'
#' @examples
#' data('example_input')
filter_cells_by_quantile = function(x, lower_quantile = 0.05, upper_quantile = .95)
{
  m_quantile = quantile(x$normalisation_factor, upper_quantile)
  l_quantile = quantile(x$normalisation_factor, lower_quantile)
  
  x = x %>% 
    mutate(del = normalisation_factor > m_quantile | normalisation_factor < l_quantile)
  
  r_cap = x %>% filter(del)
  
  if(nrow(r_cap) > 0)
  {
    cli::cli_h3("Upper quantile {.field {upper_quantile}}")
    
    cli::cli_alert_info("{.field n = {nrow(r_cap)}} entries to remove")
    
    cat("\n")
    r_cap  %>% print
    
    x = x %>% 
      filter(!del)
  }
  
  return(x %>% select(-del))
}






outliers_inspector = function(x, x_original)
{
  all_plots = lapply(x$input$segmentation$segment_id, function(s_id)
  {
    # x
    what_x_rna = x$input$dataset %>% filter(segment_id == s_id, modality == 'RNA')
    what_x_atac = x$input$dataset %>% filter(segment_id == s_id, modality == 'ATAC')
    
    if((x %>% has_rna) && what_x_rna$value_type[1] == 'NB')
      what_x_rna = normalise_modality(what_x_rna, x$input$normalisation %>% filter(modality == 'RNA'))
    
    if((x %>% has_atac) && what_x_atac$value_type[1] == 'NB')
      what_x_atac = normalise_modality(what_x_atac, x$input$normalisation %>% filter(modality == 'ATAC'))
    
    what_x = bind_rows(what_x_rna, what_x_atac)
    what_x_cellids = what_x %>% filter(outlier_cap) %>% pull(cell)
    
    # x original
    what_x_rna2 = x_original$input$dataset %>% filter(segment_id == s_id, modality == 'RNA')
    what_x_atac2 = x_original$input$dataset %>% filter(segment_id == s_id, modality == 'ATAC')
    
    if((x_original %>% has_rna) && what_x_rna2$value_type[1] == 'NB')
      what_x_rna2 = normalise_modality(what_x_rna2, x_original$input$normalisation %>% filter(modality == 'RNA'))
    
    if((x_original %>% has_atac) && what_x_atac2$value_type[1] == 'NB')
      what_x_atac2 = normalise_modality(what_x_atac2, x_original$input$normalisation %>% filter(modality == 'ATAC'))
    
    what_x_original = bind_rows(what_x_rna2, what_x_atac2)
    what_x_original$outlier = what_x_original$cell %in% what_x_cellids
    
    
    cowplot::plot_grid(
      ggplot(what_x) +
        geom_histogram(aes(value, fill = outlier_cap), bins = 70) +
        facet_wrap(segment_id~modality, scales = 'free') +
        theme_linedraw(base_size = 9) +
        labs(title = "After filtering (red: capped)") +
        scale_fill_manual(values = c(`TRUE` = 'indianred3', `FALSE` = 'gray'))  +
        guides(fill = FALSE),
      ggplot(what_x_original) +
        geom_histogram(aes(value, fill = outlier), bins = 70) +
        facet_wrap(segment_id~modality, scales = 'free') +
        theme_linedraw(base_size = 9) +
        labs(title = "Before filtering (red: outliers)") +
        scale_fill_manual(values = c(`TRUE` = 'indianred3', `FALSE` = 'gray'))  +
        guides(fill = FALSE),
      ncol = 1, 
      nrow = 2
    )
  })
  
  return(all_plots)
}





outliers_inspector = function(x, x_original)
{
  all_plots = lapply(x$input$segmentation$segment_id, function(s_id)
  {
    # x
    what_x_rna = x$input$dataset %>% filter(segment_id == s_id, modality == 'RNA')
    what_x_atac = x$input$dataset %>% filter(segment_id == s_id, modality == 'ATAC')
    
    if((x %>% has_rna) && what_x_rna$value_type[1] == 'NB')
      what_x_rna = normalise_modality(what_x_rna, x$input$normalisation %>% filter(modality == 'RNA'))
    
    if((x %>% has_atac) && what_x_atac$value_type[1] == 'NB')
      what_x_atac = normalise_modality(what_x_atac, x$input$normalisation %>% filter(modality == 'ATAC'))
    
    what_x = bind_rows(what_x_rna, what_x_atac)
    what_x_cellids = what_x %>% filter(outlier_cap) %>% pull(cell)
    
    # x original
    what_x_rna2 = x_original$input$dataset %>% filter(segment_id == s_id, modality == 'RNA')
    what_x_atac2 = x_original$input$dataset %>% filter(segment_id == s_id, modality == 'ATAC')
    
    if((x_original %>% has_rna) && what_x_rna2$value_type[1] == 'NB')
      what_x_rna2 = normalise_modality(what_x_rna2, x_original$input$normalisation %>% filter(modality == 'RNA'))
    
    if((x_original %>% has_atac) && what_x_atac2$value_type[1] == 'NB')
      what_x_atac2 = normalise_modality(what_x_atac2, x_original$input$normalisation %>% filter(modality == 'ATAC'))
    
    what_x_original = bind_rows(what_x_rna2, what_x_atac2)
    what_x_original$outlier = what_x_original$cell %in% what_x_cellids
    
    
    cowplot::plot_grid(
      ggplot(what_x) +
        geom_histogram(aes(value, fill = outlier_cap), bins = 70) +
        facet_wrap(segment_id~modality, scales = 'free') +
        theme_linedraw(base_size = 9) +
        labs(title = "After filtering (red: capped)") +
        scale_fill_manual(values = c(`TRUE` = 'indianred3', `FALSE` = 'gray'))  +
        guides(fill = FALSE),
      ggplot(what_x_original) +
        geom_histogram(aes(value, fill = outlier), bins = 70) +
        facet_wrap(segment_id~modality, scales = 'free') +
        theme_linedraw(base_size = 9) +
        labs(title = "Before filtering (red: outliers)") +
        scale_fill_manual(values = c(`TRUE` = 'indianred3', `FALSE` = 'gray'))  +
        guides(fill = FALSE),
      ncol = 1, 
      nrow = 2
    )
  })
  
  return(all_plots)
}
