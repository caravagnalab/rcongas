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
plot_cell_statistics = function(x,
                                cutoff_zeroes = 0.9,
                                cutoff_counts = 0,
                                assembly = TRUE,
                                ...)
{
  stopifnot(inherits(x, 'rcongas'))
  
  input_rna = get_input_raw_data(x, ...)
  
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
    geom_hline(yintercept = cutoff_zeroes,
               color = 'indianred3',
               linetype = 'dashed') +
    scale_x_discrete(limits = zpcounts_tb$cell) +
    labs(x = "Cell",
         y = "%",
         title = "Zeroes proportions")
  
  if (has_inference(x))
    p1 = p1 + facet_wrap(~ cluster)
  
  
  counts_total = input_rna %>% as.vector()
  counts_total = counts_total[counts_total > cutoff_counts]
  counts_total = as_tibble(counts_total)
  
  p2 = ggplot(counts_total) +
    geom_histogram(aes(value), bins = 100) +
    labs(x = "Counts",
         y = "Observations (log)",
         title = "Counts distribution") +
    scale_y_log10() +
    CNAqc:::my_ggplot_theme()
  
  if(assembly)
    return(cowplot::plot_grid(p1, p2, align = 'h', axis = 'tb'))
  else
    return(list(zeroes = p1, counts = p2))
}
