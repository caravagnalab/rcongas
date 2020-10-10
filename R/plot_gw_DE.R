# x = calculate_DE(x, mat_pre2 %>% t, 1, 2)
# plot_gw_DE(x)

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
                      annotate_top_DE = 5,
                      ...
                      )
{
  # Special case - analysis not available
  if(!has_DE(x)) return(CNAqc:::eplot())

  # Load DE results - forward params
  DE_table = get_DE_table(x, chromosomes = chromosomes, ...)

  # Get gene locations (this and the next steps shoud be
  # done when we compute DE!)
  gene_locations = x$data$gene_locations %>%
    dplyr::filter(gene %in% DE_table$gene) %>%
    dplyr::mutate(from = as.numeric(from), to = as.numeric(to))

  # Bind genes to locations, and count NUM of DEs
  n_sign = DE_table %>%
    dplyr::select(-chr, -from, -to) %>%
    dplyr::left_join(gene_locations, by = 'gene', suffix = c('.gene', '.segment')) %>%
    dplyr::group_by(chr, from, to) %>%
    dplyr::summarise(n_signif = n())

  n_sign = n_sign[complete.cases(n_sign), ]

  n_sign = CNAqc:::relative_to_absolute_coordinates(list(reference_genome = x$reference_genome), n_sign)

  # Bind genes to locations, and count NUM of genes
  n_tot = x$DE$table %>%
    dplyr::select(-chr, -from, -to) %>%
    dplyr::left_join(gene_locations, by = 'gene', suffix = c('.gene', '.segment')) %>%
    dplyr::group_by(chr, from, to) %>%
    dplyr::summarise(n_tot = n())

  n_tot = n_tot[complete.cases(n_tot), ]

  n_tot = CNAqc:::relative_to_absolute_coordinates(list(reference_genome = x$reference_genome), n_tot)


  # top_ten_DE = DE_table %>%
  #   dplyr::arrange(p_val_adj) %>%
  #   dplyr::filter(row_number() <= annotate_top_DE)
  #
  # top_ten_DE = CNAqc:::relative_to_absolute_coordinates(list(reference_genome = x$reference_genome), top_ten_DE)

  # Get segments plot - gw
  segments_plot = plot_gw_cna_profiles(x, whole_genome = TRUE, chromosomes = chromosomes)

  # Annotate a dashed line
  # segments_plot = segments_plot +
  #   geom_vline(
  #     data = top_ten_DE,
  #     aes(xintercept = from),
  #     color = 'black',
  #     linetype = "dashed",
  #     size = .25
  #   )

  # segments_plot = segments_plot +
  # geom_text(
  #   data = top_ten_DE,
  #   aes(x = from, y = Inf, label = gene),
  #   # nudge_y      = 0.0,
  #   # direction    = "x",
  #   angle        = 90,
  #   vjust        = 0,
  #   hjust = 1,
  #   size = 1.5
  # ) +
  #   coord_cartesian(clip = 'off')

  # segments_plot = segments_plot +
  #   ggrepel::geom_label_repel(
  #     data = top_ten_DE,
  #     aes(x = from, y = Inf, label = gene),
  #     size = 1.5,
  #     color = 'white',
  #     fill = 'black'
  #   ) +
  #   coord_cartesian(clip = 'off')

  # Prepare a blank WG plot
  blank_wg = get_plain_chrplot(x$reference_genome, chromosomes = chromosomes) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(angle = 0))

  gw_n_plot = blank_wg +
    labs(y = 'Genes', x = NULL) +
    geom_rect(data = n_sign,
              aes(
                xmin = as.numeric(from),
                xmax = as.numeric(to),
                ymin = 0,
                ymax = n_signif,
                fill = n_signif
              )) +
    scale_y_continuous(breaks = c(0, ceiling(max(n_sign$n_signif)))) +
    scale_fill_distiller(palette = 'Greys', direction = 1) +
    guides(fill = FALSE)

  gw_t_plot = blank_wg +
    labs(y = "DEG", x = NULL) +
    geom_rect(data = n_tot,
              aes(
                xmin = from,
                xmax = to,
                ymin = 0,
                ymax = n_tot,
                fill = n_tot
              )) +
    scale_y_continuous(breaks = c(0, ceiling(max(n_tot$n_tot)))) +
    scale_fill_distiller(palette = 'Purples', direction = 1) +
    guides(fill = FALSE)



  cowplot::plot_grid(
    gw_t_plot,
    gw_n_plot,
    segments_plot,
    rel_heights = c(.2, .2, 1),
    ncol = 1,
    align = 'v'
  )
}
