# inp3 = input congas
# res = output congas
# mat_pre2 = ...
# input = inp3
# x = res
# head(get_counts(res, inp3))

#
# plot_gw_
#
plot_gw_cna_profiles = function(x, whole_genome = FALSE, chromosomes = paste0("chr", c(1:22, "X", "Y")))
{
  # Auxiliary plain plot function for whole_genome views
  get_plain_chrplot = function(
    reference = 'hg19',
    chromosomes = paste0("chr", c(1:22, "X", "Y"))
    )
  {
    reference_coordinates = CNAqc:::get_reference(reference) %>%
      dplyr::filter(chr %in% chromosomes)

    low = min(reference_coordinates$from)
    upp = max(reference_coordinates$to)

    ggplot(reference_coordinates) +
      CNAqc:::my_ggplot_theme() +
      geom_rect(
        aes(
          xmin = centromerStart,
          xmax = centromerEnd,
          ymin = -Inf,
          ymax = Inf
        ),
        alpha = 0.3,
        colour = "gainsboro"
      ) +
      labs(x = "Chromosome",
           y = "Normalised Copy Number") +
      ggpubr::rotate_y_text() +
      scale_x_continuous(
        breaks = c(0, reference_coordinates$from,
                   upp),
        labels = c(
          "",
          gsub(
            pattern = "chr",
            replacement = "",
            reference_coordinates$chr
          ),
          ""
        )
      )

  }

  # Returned plot objects
  segments_plot = NULL

  # Segments ploidy
  segments = get_clones_ploidy(x, chromosomes)

  # Two distinct view
  if (whole_genome)
  {
    # plain chr plot like in CNAqc
    plain_plot = NULL

    if (reference == "hg38")
      plain_plot = get_plain_chrplot(reference = 'GRCh38', chromosomes)
    else
      plain_plot = get_plain_chrplot(reference, chromosomes)

    # Adjustment for the view
    segments = CNAqc:::relative_to_absolute_coordinates(list(reference_genome = reference), segments)

    segments_plot = plain_plot +
      geom_segment(
        data = segments,
        aes(
          x = from,
          xend = to,
          y = CN,
          yend = CN,
          colour = cluster
        ),
        size = 1.5
      ) +
      scale_colour_manual(values = get_clusters_colors(segments$cluster))
  }

  if(!whole_genome)
  {
    # chr ordering
    levels_chr_ordering = paste0("chr", c(1:22, "X", "Y"))

    segments_plot = ggplot() +
      geom_segment(
        data = segments,
        aes(
          x = from,
          xend = to,
          y = CN,
          yend = CN,
          colour = cluster
        ),
        size = 1.5
      ) +
      scale_colour_manual(values = get_clusters_colors(segments$cluster)) +
      facet_wrap(~factor(chr, levels = levels_chr_ordering), nrow = 1) +
      CNAqc:::my_ggplot_theme() +
      labs(
        x = "Chromosome",
        y = "Normalised Copy Number"
      ) +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
  }

  # Final result
  return(segments_plot)
}


