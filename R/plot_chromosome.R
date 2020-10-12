#' Title
#'
#' @param x
#' @param input
#' @param chr
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plot_chromosome = function(x, chr = 'chr1', ...)
{
  # Get all counts
  all_counts = get_counts(x, chromosomes = chr, normalise = TRUE, ...)

  if (nrow(all_counts) == 0)
    return(CNAqc:::eplot())

  # Segments plot with Salvatore' trick
  range = get_clones_ploidy(x)
  range_low = min(range$CN)
  range_up = max(range$CN)

  segments_plot = plot_gw_cna_profiles(x, chromosomes = chr) + ylim(range_low, range_up) +
    labs(title = paste0("Chromosome: ", chr))

  # Counts plot
  counts_plot = all_counts %>%
    ggplot() +
    geom_histogram(aes(n, fill = cluster),
                   bins = 100) +
    facet_grid(from ~ chr) +
    CNAqc:::my_ggplot_theme() +
    labs(x = NULL, y = NULL, title = 'Normalised counts') +
    scale_fill_manual(values = get_clusters_colors(all_counts$cluster))

  # Volcano - if DE is available
  volc_plot = CNAqc:::eplot()
  if(has_DE(x))
    volc_plot = plot_DE_volcano(x, chromosomes = chr)

  # Assembly
  cowplot::plot_grid(
    segments_plot,
    counts_plot,
    volc_plot,
    nrow = 1,
    align = 'h',
    axis = 'tb'
  )

}
