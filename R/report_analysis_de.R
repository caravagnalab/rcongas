#' Report for a template set of results after CONGAS
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
#' 
#' x = Rcongas::congas_example
#'
#' print(x)
#' 
report_analysis_de = function(x)
{
  require(tidyverse)

  # Top plot
  DE_gw = Rcongas::plot_DE_gw(x)
  DE_gw = ggpubr::annotate_figure(DE_gw,
                                  top = ggpubr::text_grob(label = "   Clones inferred, and DE analysis",  x = 0, hjust = 0))

  counts_plot = Rcongas::plot_counts_rna_segments(x, normalised = TRUE, z_score = TRUE)

  DE_volcano = Rcongas::plot_DE_volcano(x)

  mid_strip = ggpubr::ggarrange(
    counts_plot,
    DE_volcano,
    nrow = 1,
    widths = c(1, .5),
    labels = c("B", "C")
  )

  # Counts data for the chromosomes
  M = Rcongas::get_counts(x) %>% Rcongas:::idify()

  bottom_strip = ggplot(M) +
    geom_histogram(aes(n, fill = cluster), bins = 100) +
    CNAqc:::my_ggplot_theme() +
    facet_wrap( ~ segment_id, scales = 'free') +
    scale_fill_manual(values = Rcongas:::get_clusters_colors(M$cluster)) +
    labs(title = "Counts per chromosome")

  figure = ggpubr::ggarrange(
    DE_gw,
    mid_strip,
    bottom_strip,
    ncol = 1,
    heights = c(.5, .5, 1),
    labels = c("A", "", "D")
  )

  # ggsave(figure, filename = "breast.pdf", height = 18, width = 13)
  return(figure)
}
