#' Title
#'
#' @param x
#' @param whole_genome
#' @param chromosomes
#'
#' @return
#' @export
#'
#' @examples
#' x = Rcongas::congas_example
#' 
#' # Default view 
#' plot_gw_cna_profiles(x)
#' 
#' # Avoid facet (shows real chromosome width)
#' plot_gw_cna_profiles(x, whole_genome = TRUE)
#' 
#' # Subset to one chromosome
#' plot_gw_cna_profiles(x, chromosomes = c('chr1'))
#' 
#' # Forward alpha to get_clusters_ploidy (changes what is highlighted)
#' plot_gw_cna_profiles(x, whole_genome = TRUE, alpha = 0.001)

plot_gw_cna_profiles = function(x,
                                whole_genome = FALSE,
                                chromosomes = paste0("chr", c(1:22, "X", "Y")),
                                ...
                                )
{
  # Returned plot objects
  segments_plot = NULL

  # Segments ploidy
  segments = get_clusters_ploidy(x, chromosomes, ...)
  
  segments_h = CNAqc:::relative_to_absolute_coordinates(list(reference_genome = x$reference_genome), segments %>% dplyr::filter(highlight))
  segments_nh = CNAqc:::relative_to_absolute_coordinates(list(reference_genome = x$reference_genome), segments %>% dplyr::filter(!highlight))

  # Highlights
  
  l <- list(...)
  alpha = 0.05
  
  if('alpha' %in% names(l)) alpha = l['alpha'] %>% as.numeric()
  
  table_highlights = highlights(x, alpha) %>% filter(highlight)

  # Test for the difference
  # test_pvalue = Rcongas::get_segment_test_counts(x,
  #                                                group1 = 1,
  #                                                group2 = 2,
  #                                                cutoff_p = cutoff_p) %>%
  #   dplyr::filter(sign) %>%
  #   dplyr::filter(chr %in% chromosomes)

  ymin = segments$CN %>% min
  ymax = segments$CN %>% max
  ymin = ymin + ymin * .15
  ymax = ymax + ymax * .15

  shading_color = 'mediumseagreen'
  shading_color = 'darkgoldenrod2'
  shading_color = 'indianred3'
  

  # Two distinct view
  if (whole_genome)
  {
    # plain chr plot like in CNAqc
    plain_plot = get_plain_chrplot(x$reference_genome, chromosomes)

    plain_plot = plain_plot +
      geom_rect(
        data = segments_nh ,
        aes(
          xmin = as.numeric(from),
          xmax = as.numeric(to),
          ymin = ymin,
          ymax = ymax
        ),
        fill = "gainsboro",
        alpha = .2
      )

    # Add test_pvalue info
    if(nrow(segments_h) > 0)
    {
      # test_pvalue = CNAqc:::relative_to_absolute_coordinates(list(reference_genome = x$reference_genome), test_pvalue)

      plain_plot = plain_plot +
      geom_rect(
        data = segments_h ,
        aes(
          xmin = as.numeric(from),
          xmax = as.numeric(to),
          ymin = ymin,
          ymax = ymax
        ),
        fill = shading_color,
        alpha = .2
      ) #+
        # labs(caption = paste0("Shaded area: Poisson test significant at level ", cutoff_p))
    }

    # Adjustment for the view
    segments_sc = CNAqc:::relative_to_absolute_coordinates(list(reference_genome = x$reference_genome), segments)

    segments_plot = plain_plot +
      geom_segment(
        data = segments_sc,
        aes(
          x = from,
          xend = to,
          y = CN,
          yend = CN,
          colour = cluster
        ),
        size = .5
      ) 
    
    # Thicker
    segments_sc_highlight = NULL
    
    if(nrow(table_highlights) > 0)
    {
      for (i in 1:nrow(table_highlights))
        segments_sc_highlight = segments_sc_highlight %>%
        bind_rows(
          segments_sc %>%
            filter(table_highlights$segment_id[i] == segment_id) %>%
            filter(
              cluster %in% c(table_highlights$cluster[i], table_highlights$versus[i])
            )
        ) %>% 
        distinct(segment_id, cluster, .keep_all =T)
    
    
    segments_plot = segments_plot+
      geom_segment(
        data = segments_sc_highlight,
        aes(
          x = from,
          xend = to,
          y = CN,
          yend = CN,
          colour = cluster
        ),
        size = 1.5
      ) +
      scale_colour_manual(values = get_clusters_colors(segments_sc$cluster)) +
      guides(color = guide_legend(bquote("Cluster ("*alpha*' = '*.(alpha)*")")))
    }
  }

  if (!whole_genome)
  {
    # chr ordering
    levels_chr_ordering = paste0("chr", c(1:22, "X", "Y"))

    segments_plot = ggplot()

    if(nrow((segments %>% dplyr::filter(highlight)) > 0))
    {
      segments_plot = segments_plot +
        geom_rect(
          data = segments %>% dplyr::filter(highlight) ,
          aes(
            xmin = as.numeric(from),
            xmax = as.numeric(to),
            ymin = ymin,
            ymax = ymax
          ),
          fill = shading_color,
          alpha = .2
        ) #+
       # labs(caption = paste0("Shaded area: Poisson test significant at level ", cutoff_p))
    }

    segments_plot = segments_plot +
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
      facet_wrap( ~ factor(chr, levels = levels_chr_ordering), nrow = 1) +
      CNAqc:::my_ggplot_theme() +
      labs(x = "Chromosome",
           y = "Normalised CN") +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
  }

  # Add DE information on the plot
  # if(has_DE(x))
  # {
  #   DE_table = get_DE_table(x, chromosomes = chromosomes, cut_pvalue = pvalue_cut_DE, cut_lfc = lfc_cut_DE)
  #
  #   if (whole_genome)
  #     DE_table = CNAqc:::relative_to_absolute_coordinates(list(reference_genome = x$reference_genome), DE_table)
  #
  #   # TODO aggiornare messaggio stampa
  #   cli::cli_alert_info("Found DE analysis results, annotating n = {.value {nrow(DE_table)}} genes with adjusted p-value < {.field {p_cut_DE}}.")
  #
  #   if(nrow(DE_table) > 0)
  #   {
  #     segments_plot +
  #       geom_vline(
  #         data = DE_table,
  #         aes(xintercept = from),
  #         size = .1,
  #         linetype = 'dashed'
  #       )
  #
  #
  #   }
  #   else
  #   {
  #     cli::cli_alert_warning("No genes with significant DE with the requested parameters.")
  #   }
  #
  # }

  # Final result
  return(segments_plot)
}


# Auxiliary plain plot function for whole_genome views
get_plain_chrplot = function(reference = 'hg19',
                             chromosomes = paste0("chr", c(1:22, "X", "Y")))
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
