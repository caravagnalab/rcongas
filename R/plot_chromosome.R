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


plot_chromosome = function(x, chr = 'chr1', plot = c("segments", "counts", "DE"), ...)
{
  # Get all counts
  all_counts = get_counts(x, chromosomes = chr, normalise = TRUE, ...)

  if (nrow(all_counts) == 0)
    return(CNAqc:::eplot())

  # Segments plot with Salvatore' trick
  range = get_clusters_ploidy(x)
  range_low = min(range$CN)
  range_up = max(range$CN)

  segments_plot = NULL

  if("segments" %in% plot)
    segments_plot = plot_gw_cna_profiles(x, chromosomes = chr) + ylim(range_low, range_up) +
      labs(title = paste0("Chromosome: ", chr))

  # Counts plot
  counts_plot = NULL

  if("counts" %in% plot)
    counts_plot = all_counts %>%
      ggplot() +
      geom_histogram(aes(n, fill = cluster),
                     bins = 100) +
      facet_grid(from ~ chr) +
      CNAqc:::my_ggplot_theme() +
      labs(x = NULL, y = NULL, title = 'Normalised counts') +
      scale_fill_manual(values = Rcongas:::get_clusters_colors(all_counts$cluster))

  # Volcano - if DE is available
  volc_plot = NULL

  if("DE" %in% plot)
  {
    volc_plot = CNAqc:::eplot()
    if(has_DE(x))
      volc_plot = plot_DE_volcano(x, chromosomes = chr)
  }

  # Assembly - depends on plot = c("segment", "counts", "DE")
  all_plots = list(
    segments_plot,
    counts_plot,
    volc_plot
  )

  all_plots = all_plots[!sapply(all_plots, is.null)]

  if(length(all_plots) == 1) return(all_plots[[1]])

  cowplot::plot_grid(
    plotlist = all_plots,
    nrow = 1,
    align = 'h',
    axis = 'tb'
  )

}
