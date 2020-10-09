# x = calculate_DE(x, mat_pre2 %>% t, 1, 2)

#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
plot_gw_DE = function(x,
                      chromosomes = paste0("chr", c(1:22, "X", "Y")),
                      ...
                      )
{
  # Special case - analysis not available
  if(!has_DE(x)) return(CNAqc:::eplot())

  # Load DE results - forward params
  DE_table = get_DE_table(x, chromosomes = chromosomes, ...)

  # Get segments plot - gw
  segments_plot = plot_gw_cna_profiles(x, whole_genome = TRUE, chromosomes = chromosomes)

  # Prepare a blank WG plot
  blank_wg = get_plain_chrplot(x$reference_genome, chromosomes = chromosomes) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(angle = 0))

  gw_n_plot = blank_wg +
    labs(y = 'Genes', x = NULL) +
    geom_rect(data = y,
              aes(
                xmin = from,
                xmax = to,
                ymin = 0,
                ymax = n,
                fill = n
              )) +
    scale_y_continuous(breaks = c(0, max(y$n))) +
    scale_fill_distiller(palette = 'Greys', direction = 1) +
    guides(fill = FALSE)

  gw_p_plot = blank_wg +
    labs(y = "DEG", x = NULL) +
    geom_rect(data = y,
              aes(
                xmin = from,
                xmax = to,
                ymin = 0,
                ymax = p,
                fill = p
              )) +
    scale_y_continuous(breaks = c(0, 1)) +
    scale_fill_distiller(palette = 'Purples', direction = 1) +
    guides(fill = FALSE)



  cowplot::plot_grid(
    gw_p_plot,
    gw_n_plot,
    segments_plot,
    rel_heights = c(.2, .2, 1),
    ncol = 1,
    align = 'v'
  )
}
