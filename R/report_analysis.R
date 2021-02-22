

#' Title
#'
#' @param x
#' @param ...
#' @param cex 
#'
#' @return
#' @export
#'
#' @importFrom ggpubr ggarrange
#' @importFrom CNAqc plot_segments
#' @importFrom patchwork plot_layout
#'
#' @examples
report_analysis <- function(x, cex = 1, ...)
{
  highlights <-
    x %>% Rcongas::get_clusters_ploidy() %>% filter(highlight == TRUE)
  
  chr <- unique(highlights$chr)
  
  gene_locations <- as.data.frame(x$data$gene_locations$segment_id)
  gene_locations <- distinct(gene_locations)
  
  tot_counts <- as.data.frame(rowSums(x$data$counts))
  tot_counts <- log(tot_counts)
  
  colnames(tot_counts) = "log_counts"
  
  clusters <- as.data.frame(Rcongas:::get_cluster_assignments(x))
  colnames(clusters) = "clusters"
  annotations <- cbind(tot_counts, clusters)
  
  curate = function(p) {
    p +
      # theme_minimal(base_size=10 * cex, base_family="Roboto") +
      theme(
        plot.title = element_text(
          size = 11 * cex,
          color = "indianred3",
          face = "bold",
          family = "Roboto Black"
        ),
        plot.subtitle = element_text(color = "#666666", size = 7 * cex),
        # plot.title = element_text(family="Roboto Black"),
        plot.caption = element_text(color = "#AAAAAA", size = 5 * cex)
        # legend.position = 'bottom'
        
      )
    
  }
  
  # p1 <- Rcongas:::plot_highlights(x) %>% curate
  # p1 = ggplot()
  
  stats_data = Rcongas::get_dataset_stats(x)
  
  highlights_segments = Rcongas:::get_segment_ids(x, highlight = TRUE) %>% unique
  
  ph = plot_highlights(x)
  
  # 2p
  p1 = lapply(plot_cohort_statistics(x, assembly = FALSE, ...), curate)
  
  # 1p
  p2 <- plot_mixing_proportions(x) %>% curate
  
  p12 = p1[[1]] | p1[[2]]
  
  # 3p
  bm = get_best_model(x)
  
  p3a = (plot_loss(x) +
           labs(
             title = paste("ELBO"),
             subtitle = paste0(
               stats_data$score_type,
               ' = ',
               stats_data$score %>% round(3),
               " (",
               bm$loss %>% length,
               ' steps)'
             )
           )) %>%
    curate
  
  p3b = (plot_normalization_factors(x)  +
           labs(
             title = paste("Library size"),
             subtitle = paste(
               "Theta (Gamma prior): shape = ",
               bm$run_information$input_hyper_params$theta_shape,
               '; rate  = ',
               bm$run_information$input_hyper_params$theta_rate
             )
           )) %>%
    curate
  
  p3c = plot_latent_variables(x) %>%
    curate
  
  # p3 = p3a | p3b | p3c
  
  p3 =  p3a | p3b
  # p12 = p3c | p1[[1]]
  p12 = p3c | ph
  
  p4b <-
    Rcongas:::plot_counts_rna_segments(x, normalised = T, z_score = TRUE) %>%
    curate
  
  p5 <-
    (
      Rcongas:::plot_gw_cna_profiles(x, whole_genome = TRUE) +
        labs(title = "Tumour CNA segments", subtitle = "Default highlight parameters.")
    ) %>%
    curate
  
  
  
  segments_ids = Rcongas:::get_segment_ids(x, highlight = TRUE) %>% unique
  ns = segments_ids %>% length %>% sqrt %>% ceiling
  
  # Density plots
  density_plots = lapply(plot_segment_density(x, segments_ids = segments_ids),
                         function(x)
                           x %>% curate)
  
  if (segments_ids %>% length <= 4) 
  {
    ns = 1
    
    bottom_panel = ggpubr::ggarrange(
      plotlist = density_plots,
      nrow = 1,
      ncol = segments_ids %>% length,
      common.legend = TRUE,
      legend = "bottom",
      labels = 'h'
    )
  }
  else {
    bottom_panel = ggpubr::ggarrange(
      plotlist = density_plots,
      nrow = ifelse(ns * (ns - 1) < (segments_ids %>% length), ns, ns - 1),
      ncol = ns,
      common.legend = TRUE,
      legend = "bottom",
      labels = 'h'
    )
  }
  
  strip_top =  p2 + p5  +  plot_layout(widths =  c(1, 4),
                                       ncol = 2,
                                       nrow = 1)
  
  strip_stat = ((p3) / p12)
  strip_stat = strip_stat | p4b
  
  ((strip_top / strip_stat) /
      bottom_panel)  +
    plot_layout(heights = c(.5, 1, ns),
                ncol = 1,
                nrow = 3)
  
  
  # ggpubr::ggarrange(top_panel,
  #                   mid_panel,
  #                   bottom_panel,
  #                   heights = c(1, 1, ns),
  #                   nrow = 3)
  
}


xy_cell_plot <- function(x,
                         f_x = min,
                         f_y = max,
                         f_x_label = "min",
                         f_y_label = "max",
                         quantiles = c(0.1, 0.9),
                         segments_ids = Rcongas:::get_segment_ids(x)) {
  counts_data = get_counts(
    x,
    normalise = FALSE,
    z_score = FALSE,
    sum_denominator = FALSE
  ) %>%
    Rcongas:::idify() %>%
    filter(segment_id %in% segments_ids)
  
  nsegments = segments_ids %>% length
  
  max_counts = counts_data %>%
    group_by(cell) %>%
    summarise(f_x = f_x(n))
  qmax_counts = quantile(max_counts$f_x, quantiles)
  
  max_counts = max_counts %>%
    mutate(off_qx = f_x < qmax_counts[1] | f_x > qmax_counts[2])
  
  min_counts = counts_data %>%
    group_by(cell) %>%
    summarise(f_y = f_y(n))
  qmin_counts = quantile(min_counts$f_y, quantiles)
  
  min_counts = min_counts %>%
    mutate(off_qy = f_y < qmin_counts[1] | f_y > qmin_counts[2])
  
  joined = full_join(max_counts, min_counts, by = "cell")
  
  # If there are clustering result inside
  if (Rcongas:::has_inference(x))
  {
    joined = joined %>%
      full_join(get_clusters(x) %>%
                  dplyr::select(cell, cluster),
                by = 'cell')
  }
  else{
    joined$cluster = "Not clustered"
  }
  
  n_true = sum(joined$off_qx | joined$off_qy)
  
  plot = ggplot(joined, aes(x = f_x, y = f_y, color = cluster)) +
    geom_point(size = 1, aes(shape = off_qx | off_qy)) +
    labs(title = paste0("XY statistics with ", nsegments, " segments")) +
    scale_color_brewer(palette = "Set1") +
    CNAqc:::my_ggplot_theme() +
    geom_abline()  +
    labs(x = f_x_label,
         y = f_y_label) +
    geom_vline(xintercept = qmax_counts,
               linetype = 'dashed',
               color = 'indianred3') +
    geom_hline(yintercept = qmin_counts,
               linetype = 'dashed',
               color = 'indianred3') +
    guides(shape = guide_legend(paste0(
      "Outside q={", paste(quantiles, collapse = ', '), '} n = ', n_true
    )))
  
  return(list(dataframe = joined, plot = plot))
}



# install.packages('patchwork')
#
# library(patchwork)
#
# cowplot::plot_grid(
#   p3,
#   p2,
#   p1,
#   rel_widths = c(1, .5, 1),
#   ncol = 3,
#   align = 'h',
#   axis = 'tb'
# )

#
#
# strip_stat = p2 + p1 +
#   plot_layout(widths = c(1, 2))
#
# strip_stat = ((p3)/strip_stat)
# strip_stat = strip_stat | p4b
#
#
# cowplot::plot_grid(
#   # p4,
#   p4b,
#   p5,
#   ncol = 2,
#   rel_widths = c(1, .8),
#   align = 'h',
#   axis = 'tb'
# )
