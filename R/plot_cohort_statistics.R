#' Title
#'
#' @param x
#' @param cutoff_zeroes
#' @param cutoff_counts
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plot_cohort_statistics = function(x,
                                cutoff = 0.9,
                                assembly = TRUE,
                                ...)
{
  stopifnot(inherits(x, 'rcongas'))
  
  
  if(!has_inference(x))
  {
    input_rna = get_input_raw_data(x, ...)
    input_rna$cluster = 'Not Available'
  }
  else
    input_rna = get_input_raw_data(x, add_clusters = TRUE, ...)
  
  p1 = aux_plot_cells_zerocounts(x, input_rna, cutoff)
  p2 = aux_plot_counts_distribution(x, input_rna, cutoff)
  p3 = aux_plot_genes_zerocounts(x, input_rna, cutoff)
    
  
  if(assembly)
    return(cowplot::plot_grid(p1, p2, p3, align = 'h', axis = 'tb'))
  else
    return(list(zeroes = p1, counts = p2, genes = p3))
}

aux_plot_cells_zerocounts = function(x, input_rna, cutoff)
{
  # Counts per cell
  # zcounts = apply(input_rna,
  #                 2,
  #                 function(x) {
  #                   sum(x == 0)
  #                 })
  
  norm_factor = get_dataset_stats(x)$ngenes  

  # Already filtered (these are >0, so 1 - % is what i want)
  foo = input_rna %>%
    group_by(cell, cluster) %>%
    summarise(n = n()) %>% 
    mutate(p = 1 - n/norm_factor) %>% 
    arrange(p) %>% 
    ungroup() %>% 
    mutate(x = row_number())
  
  # zpcounts = zcounts / nrow(input_rna)
  
  # zpcounts_tb = zpcounts %>% as_tibble()
  # zpcounts_tb$cell = names(zcounts)
  # 
  # zpcounts_tb = zpcounts_tb %>%
  #   arrange(value)
  # 
  # if (has_inference(x))
  #   zpcounts_tb = zpcounts_tb %>%
  #   left_join(get_clusters(x),
  #             by = 'cell')
  
  p1 = ggplot(foo) +
    geom_bar(aes(x = x, y = p), stat = 'identity') +
    CNAqc:::my_ggplot_theme() +
    theme(
      axis.text.x = element_blank(),
      legend.position = 'right',
      axis.ticks.x = element_blank()
    ) +
    ylim(0, 1) +
    geom_hline(yintercept = cutoff,
               color = 'indianred3',
               linetype = 'dashed') +
    labs(x = "Cell",
         y = "% genes with 0 reads",
         title = "Zero-counts across cells")
  
  if (has_inference(x))
    p1 = p1 + facet_wrap(~ cluster)
  
  p1
}

aux_plot_counts_distribution = function(x, input_rna, cutoff)
{
  ggplot(input_rna %>% filter(n > cutoff)) +
    geom_histogram(aes(n), bins = 100) +
    labs(x = "Counts" ,
         y = "Observations (log)",
         title = "Total counts per cell") +
    scale_y_log10() +
    CNAqc:::my_ggplot_theme() +
    geom_vline(xintercept = cutoff,               
               color = 'indianred3',
               linetype = 'dashed') 
}

aux_plot_genes_zerocounts = function(x, input_rna, cutoff)
{
  norm_factor = get_dataset_stats(x)$ncells  
  
  # Already filtered (these are >0, so 1 - % is what i want)
  foo = input_rna %>%
    group_by(gene, cluster) %>%
    summarise(n = n()) %>% 
    mutate(p = 1 - n/norm_factor) %>% 
    arrange(p) %>% 
    ungroup() %>% 
    mutate(x = row_number())
  
  
  p1 = ggplot(foo) +
    geom_bar(aes(x = x, y = p), stat = 'identity') +
    CNAqc:::my_ggplot_theme() +
    theme(
      axis.text.x = element_blank(),
      legend.position = 'right',
      axis.ticks.x = element_blank()
    ) +
    ylim(0, 1) +
    geom_hline(yintercept = cutoff,
               color = 'indianred3',
               linetype = 'dashed') +
    labs(x = "Gene",
         y = "% cells with 0 reads",
         title = "Zero-counts across cells")

  if (has_inference(x))
    p1 = p1 + facet_wrap(~ cluster)
  
  p1
  
}





