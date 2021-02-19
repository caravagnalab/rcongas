#' Title
#'
#' @param x
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
report_data = function(x, ...)
{
  curate = function(p) {
    p +
      theme(
        plot.title = element_text(
          size = 11 * cex,
          color = "indianred3",
          face = "bold",
          family = "Roboto Black"
        ),
        plot.subtitle = element_text(color = "#666666", size = 7 * cex),
        plot.caption = element_text(color = "#AAAAAA", size = 5 * cex)
      )
  }
  
  # Mid and bottom panels: what we show
  to_show = get_input_segmentation(x) %>% pull(chr) %>% unique
  
  if (has_inference(x))
    to_show = highlights(x, ...) %>%
    filter(highlight) %>%
    pull(chr) %>%
    unique
  
  # Top panel
  prna = plot_counts_rna_segments(x, z_score = TRUE, chromosomes = to_show, ...) +
    theme(axis.text.x = element_blank()) +
    curate()
  
  pcs = plot_cohort_statistics(x, assembly = FALSE, ...)
  
  tp_strip = cowplot::plot_grid(plotlist = lapply(pcs,  curate),
                                ncol = 3,
                                nrow = 1)
  
  
  # tp_panel = cowplot::plot_grid(
  #   plotlist = append(list(prna), pcs),w
  #   ncol = 4,
  #   nrow = 1,
  #   rel_widths = c(1, .5, .5, .5),
  #   axis = 'tb',
  #   align = 'h'
  # )
  
  # return(tp_panel)
  
  pcps = plot_counts_per_segment(x, chromosomes = to_show, ...) +
    curate()
  
  tp_panel = cowplot::plot_grid(
    prna,
    pcps,
    ncol = 2,
    nrow = 1,
    # rel_widths = c(1, .5, .5, .5),
    axis = 'tb',
    align = 'h'
  )
  
  # tp_panel = ggpubr::ggarrange(
  #   tp_strip,
  #   tp_panel,
  #   ncol = 1,
  #   nrow = 2,
  #   heights = c(1, 2)
  # )
  #
  # tp_panel = cowplot::plot_grid(
  #   tp_panel,
  #   pcps,
  #   ncol = 2,
  #   nrow = 1,
  #   axis = 'tb',
  #   align = 'h'
  # )
  
  pgcg = plot_gene_counts_on_genome(x, chromosomes = to_show, assembly = FALSE, ...)
  
  pgcg[[1]] = pgcg[[1]] +
    labs(title ="Counts across the genome")
    curate()
  
  pgcg = ggpubr::ggarrange(
    plotlist = pgcg,
    ncol = 4,
    nrow = (length(to_show) / 4) %>% ceiling
  )
  
  pgst = plot_gene_rank(x,
                        chromosomes = to_show,
                        assembly = FALSE,
                        top = 100,
                        ...)
  
  pgst[[1]] = pgst[[1]] +
    labs(title ="Gene rank")
  curate()
  
  
  pgst = ggpubr::ggarrange(
    plotlist = pgst,
    ncol = 4,
    nrow = (length(to_show) / 4) %>% ceiling
  )
  
  # bo_mi_panel = cowplot::plot_grid(
  #   pgcg,
  #   pgst,
  #   ncol = 1,
  #   nrow = 2,
  #   rel_heights = c(2, 1)
  # )
  
  
  # cowplot::plot_grid(
  #   tp_panel,
  #   bo_mi_panel,
  #   ncol = 1,
  #   nrow = 2,
  #   rel_heights = c(1, 2)
  # )
  
  ggpubr::ggarrange(
    tp_strip,
    cowplot::plot_grid(
      tp_panel,
      pgcg,
      pgst,
      ncol = 1,
      nrow = 3,
      rel_heights = c(1.7, 2, 1)
    ),
    ncol = 1,
    nrow = 2,
    heights = c(1, 7)
  )
  
}
