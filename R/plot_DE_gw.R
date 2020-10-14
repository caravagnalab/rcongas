# x = calculate_DE(x, mat_pre2 %>% t, 1, 2)
# plot_gw_DE(x)
# plot_gw_DE(x, color_DEG = 'Oranges', color_secondary = 'darkgreen')

#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
plot_DE_gw = function(x,
                      chromosomes = paste0("chr", c(1:22, "X", "Y")),
                      cut_pvalue = 0.001,
                      cut_lfc = 0.25
                      )
{
  # Special case - analysis not available
  if (!has_DE(x))
    return(CNAqc:::eplot())

  # Plots colours
  color_DEG = "Greys"
  color_secondary = 'orange'

  # Load DE results - forward params
  DE_full_table = get_DE_table(x,
                               chromosomes = chromosomes,
                               cut_pvalue = 1,
                               cut_lfc = 0)
  DE_table = get_DE_table(x, chromosomes = chromosomes, cut_pvalue = cut_pvalue, cut_lfc = cut_lfc)

  # Mapping function for all genes (gives us the total genes per segment)
  mapping_DE_total = map_de_to_segments(DE_full_table, x)

  mapping_DE_total_mapped = mapping_DE_total$mapping
  mapping_DE_total_mapped_stats_total = mapping_DE_total_mapped %>%
    dplyr::group_by(segment_id) %>%
    dplyr::summarise(n_mappable = n(), .groups = 'keep') %>%
    dplyr::ungroup()

  # Mapping function for DE genes (gives us the significant genes per segment)
  mapping_DE_sign = map_de_to_segments(DE_table, x)

  mapping_DE_sign = mapping_DE_sign$mapping
  mapping_DE_sign_total = mapping_DE_sign %>%
    dplyr::group_by(segment_id) %>%
    dplyr::summarise(n_mappable_sign = n(), .groups = 'keep') %>%
    dplyr::ungroup()

  joint_counts = dplyr::full_join(mapping_DE_total_mapped_stats_total,
                                  mapping_DE_sign_total,
                                  by = "segment_id") %>%
    dplyr::mutate(prop = n_mappable_sign / n_mappable) %>%
    deidify()

  joint_counts = CNAqc:::relative_to_absolute_coordinates(list(reference_genome = x$reference_genome),
                                                          joint_counts)

  max_joint_counts = max(joint_counts$n_mappable_sign, na.rm = TRUE)

  # Get segments plot - gw, and add some caption
  segments_plot = plot_gw_cna_profiles(x, whole_genome = TRUE, chromosomes = chromosomes) +
    labs(caption = paste("Cutoffs: p-value p <", cut_pvalue, "LFC < ",
                         cut_lfc, '.'))

  # Prepare a blank WG plot
  blank_wg = get_plain_chrplot(x$reference_genome, chromosomes = chromosomes) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(angle = 0))

  # color_secondary = 'orange'

  # Plot the
  gw_n_plot =
    blank_wg +
    labs(y = 'DEG', x = NULL) +
    geom_rect(
      data = joint_counts,
      aes(
        xmin = as.numeric(from),
        xmax = as.numeric(to),
        ymin = 0,
        ymax = n_mappable_sign,
        fill = n_mappable_sign
      )
    ) +
    # scale_y_continuous(breaks = c(0, ceiling(max(
    #   joint_counts$n_mappable_sign
    # )))) +
    scale_fill_distiller(palette = color_DEG, direction = 1) +
    guides(fill = FALSE) +
    geom_line(
      data = joint_counts,
      mapping = aes(x = from + (to - from) / 2, y = prop * max_joint_counts),
      size = .45,
      color = color_secondary
    ) +
    geom_point(
      data = joint_counts %>% dplyr::mutate(prop = n_mappable_sign / n_mappable),
      mapping = aes(x = from + (to - from) / 2, y = prop * max_joint_counts),
      size = .7,
      color = color_secondary
    ) +
    scale_y_continuous(
      breaks = c(0, ceiling(
        max(joint_counts$n_mappable_sign, na.rm = TRUE)
      )),
      sec.axis = sec_axis(
        ~ . / max_joint_counts,
        name = "% DEG",
        labels = function(b) {
          paste0(round(b * 100, 0), "%")
        }
      )
    )

  # gw_t_plot = blank_wg +
  #   labs(y = "Genes", x = NULL) +
  #   geom_rect(data = joint_counts,
  #             aes(
  #               xmin = from,
  #               xmax = to,
  #               ymin = 0,
  #               ymax = n_mappable,
  #               fill = n_mappable
  #             )) +
  #   scale_y_continuous(breaks = c(0, ceiling(max(
  #     joint_counts$n_mappable
  #   )))) +
  #   scale_fill_distiller(palette = 'Purples', direction = 1) +
  #   guides(fill = FALSE)

  # Off-target top strip
  mapping_DE_sign_offt = map_de_to_segments(DE_table, x)$off_target
  mapping_DE_sign_offt = CNAqc:::relative_to_absolute_coordinates(list(reference_genome = x$reference_genome),
                                                                  mapping_DE_sign_offt)

  offtarget_plot =
    blank_wg +
    geom_rect(
      data = joint_counts,
      aes(
        xmin = as.numeric(from),
        xmax = as.numeric(to),
        ymin = 0,
        ymax = n_mappable_sign
      ),
      fill = NA
    ) +
    geom_vline(data = mapping_DE_sign_offt,
               aes(xintercept = from),
               size = .5) +
    labs(y = "Off", x = NULL) +
    theme(axis.text.y = element_blank())

  # Figure assembly
  cowplot::plot_grid(
    offtarget_plot,
    gw_n_plot,
    segments_plot,
    rel_heights = c(.075, .2, 1),
    ncol = 1,
    align = 'v'
  )
}


# m - DE table
# x - rcongas obj
map_de_to_segments = function(m, x)
{
  s = get_input_segmentation(x)

  mapping = NULL

  for (i in 1:nrow(s))
  {
    matched =  m %>%
      dplyr::filter(chr == s$chr[i],
                    from >= s$from[i],
                    to <= s$to[i])

    if (nrow(matched) > 0)
    {
      matched = matched %>%
        dplyr::mutate(segment_id = paste(s$chr[i], s$from[i], s$to[i], sep =
                                           ':'))

      mapping = dplyr::bind_rows(mapping, matched)
    }
  }

  # Those off the segments
  off_target = m %>%
    dplyr::filter(!(gene %in% mapping$gene))

  stopifnot(nrow(m) == nrow(off_target) + nrow(mapping))

  return(list(mapping = mapping,
              off_target = off_target))

}
