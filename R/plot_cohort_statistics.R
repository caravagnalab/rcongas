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
  
  input_rna = get_input_raw_data(x, ...)
  
 
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
  zcounts = apply(input_rna,
                  2,
                  function(x) {
                    sum(x == 0)
                  })
  
  zpcounts = zcounts / nrow(input_rna)
  
  zpcounts_tb = zpcounts %>% as_tibble()
  zpcounts_tb$cell = names(zcounts)
  
  zpcounts_tb = zpcounts_tb %>%
    arrange(value)
  
  if (has_inference(x))
    zpcounts_tb = zpcounts_tb %>%
    left_join(get_clusters(x),
              by = 'cell')
  
  p1 = ggplot(zpcounts_tb) +
    geom_bar(aes(x = cell, y = value), stat = 'identity') +
    theme(
      axis.text.x = element_blank(),
      legend.position = 'right',
      axis.ticks.x = element_blank()
    ) +
    ylim(0, 1) +
    geom_hline(yintercept = cutoff,
               color = 'indianred3',
               linetype = 'dashed') +
    scale_x_discrete(limits = zpcounts_tb$cell) +
    labs(x = "Cell",
         y = "% genes with 0 reads",
         title = "Cell counts")
  
  if (has_inference(x))
    p1 = p1 + facet_wrap(~ cluster)
  
  p1
}

aux_plot_counts_distribution = function(x, input_rna, cutoff)
{
  counts_total = input_rna %>% as.vector()
  counts_total = counts_total[counts_total > cutoff]
  counts_total = as_tibble(counts_total)
  
  ggplot(counts_total) +
    geom_histogram(aes(value), bins = 100) +
    labs(x = "Counts" ,
         y = "Observations (log)",
         title = "Counts distribution") +
    scale_y_log10() +
    CNAqc:::my_ggplot_theme() +
    geom_vline(xintercept = cutoff,               
               color = 'indianred3',
               linetype = 'dashed') 
}

aux_plot_genes_zerocounts = function(x, input_rna, cutoff)
{
  # Counts per cell
  gcounts = apply(input_rna,
                  1,
                  function(x) {
                    sum(x == 0)
                  })
  
  gpcounts = gcounts / ncol(input_rna)
  
  gpcounts_tb = gpcounts %>% as_tibble()
  gpcounts_tb$gene = names(gcounts)
  
  gpcounts_tb = gpcounts_tb %>%
    arrange(value)
  
  ggplot(gpcounts_tb) +
    geom_bar(aes(x = gene, y = value), stat = 'identity') +
    theme(
      axis.text.x = element_blank(),
      legend.position = 'right',
      axis.ticks.x = element_blank()
    ) +
    ylim(0, 1) +
    geom_hline(yintercept = cutoff,
               color = 'indianred3',
               linetype = 'dashed') +
    scale_x_discrete(limits = gpcounts_tb$gene) +
    labs(x = "Gene",
         y = "% cells with 0 reads",
         title = "Gene counts")
}





