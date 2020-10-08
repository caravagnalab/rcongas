plot_chromosome = function(x, chr = 'chr1', ...)
{
  # Get all counts
  all_counts = get_counts(x, input, chromosomes = chr, ...)

  if (nrow(all_counts) == 0)
    return(CNAqc:::eplot())

  # Segments plot with Salvatore' trick
  range = get_clones_ploidy(x)
  range_low = min(range$CN)
  range_up = max(range$CN)

  segments_plot = plot_gw_cna_profiles(x, chromosomes = chr) + ylim(range_low, range_up)

  # Counts plot
  counts_plot = all_counts %>%
    ggplot() +
    geom_histogram(aes(n, fill = cluster),
                   bins = 100) +
    facet_grid(from ~ chr) +
    CNAqc:::my_ggplot_theme() +
    # theme(
    #   axis.text.x = element_text(angle = 90, hjust = 1),
    #   axis.text.y = element_blank(),
    #   axis.ticks.y = element_blank()
    # ) +
    labs(x = NULL, y = NULL) +
    scale_fill_manual(values = get_clusters_colors(all_counts$cluster))


  # Assembly
  ggpubr::ggarrange(
    segments_plot,
    counts_plot,
    nrow = 1,
    common.legend = TRUE,
    legend = 'bottom'
  )

}
