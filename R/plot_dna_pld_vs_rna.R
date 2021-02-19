#' Title
#'
#' @param x
#' @param library_size_norm
#'
#' @return
#' @export
#'
#' @examples
plot_dna_pld_vs_rna <- function(x, library_size_norm = FALSE) {
  counts <-
    get_counts(
      x,
      normalise = F,
      z_score = F,
      sum_denominator = F
    ) %>%  Rcongas:::idify()
  segments <- x$data$cnv[, c("ploidy_real", "segment_id", "mu")]
  
  
  plot_df <- merge(counts, segments, by = "segment_id")
  
  
  
  plot_df <- relative_to_absolute_coordinates(plot_df)
  
  plot_df$n <-  plot_df$n / plot_df$mu
  
  
  
  
  if (library_size_norm) {
    clb <- plot_df %>% group_by(cell) %>% summarize(n_lb = sum(n))
    plot_df <-  merge(plot_df, clb, by = "cell")
    plot_df$n <-  plot_df$n / plot_df$n_lb * 10
  }
  
  cm <-
    plot_df %>% group_by(segment_id) %>% summarize(n_m = mean(n))
  
  
  plot_df <-  merge(plot_df, cm, by = "segment_id")
  
  empty <- CNAqc:::blank_genome()
  
  # p1 <-
  #   empty + geom_rect(
  #     aes(
  #       xmin = from,
  #       xmax = to,
  #       ymin = ploidy_real - 0.2,
  #       ymax =  ploidy_real + 0.2,
  #       fill = paste(ploidy_real)
  #     ),
  #     data = plot_df
  #   ) +
  #   scale_fill_discrete("Ploidy") + CNAqc:::my_ggplot_theme()
  
  # p1 <-
  #   empty + geom_rect(
  #     aes(
  #       xmin = from,
  #       xmax = to,
  #       ymin = ploidy_real - 0.2,
  #       ymax =  ploidy_real + 0.2,
  #       fill = paste(ploidy_real)
  #     ),
  #     data = plot_df %>% distinct(segment_id, from, to, ploidy_real),
  #     color = 'indianred3'
  #   ) +
  #   # scale_fill_brewer(palette = "YlGnBu", "Ploidy") +
  #   # scale_fill_discrete("Ploidy") + 
  #   CNAqc:::my_ggplot_theme()  +
  #   labs(y='Ploidy')
  # 
  # if(max(plot_df$ploidy_real) < 3.5) p1 = p1 + ylim(0,5)
  # 
  # 
  # p2 <-
  #   empty + geom_rect(
  #     aes(
  #       xmin = from,
  #       xmax = to,
  #       ymin = n - 0.005,
  #       ymax =  n + 0.005
  #     ),
  #     fill = "black",
  #     alpha = 0.6,
  #     data = plot_df
  #   ) +
  #   geom_rect(aes(
  #     xmin = from,
  #     xmax = to,
  #     ymin = n_m - (0.05 * (max(n_m) - min(n_m))) ,
  #     ymax = n_m + (0.05 * (max(n_m) - min(n_m)))
  #   ),
  #   fill = "red",
  #   data = plot_df) + CNAqc:::my_ggplot_theme() +
  #   CNAqc:::my_ggplot_theme()   +
  #   labs(y='Counts per cell')
  # 
  #   
  # 
  # plt <-
  #   cowplot::plot_grid(
  #     p1,
  #     p2,
  #     align = "hv",
  #     axis = "tblr",
  #     nrow = 2,
  #     rel_heights = c(1, 1.2)
  #   )
  
  scaling = max(plot_df$ploidy_real)
  
  p1 <-
    empty + 
    geom_rect(
      aes(
        xmin = from,
        xmax = to,
        ymin = ploidy_real - 0.2 - scaling - 2,
        ymax =  ploidy_real + 0.2 - scaling - 2,
      ),
      data = plot_df %>% distinct(segment_id, from, to, ploidy_real),
      fill = 'steelblue'
    ) +
    CNAqc:::my_ggplot_theme()  +
    labs(y='Counts') +
    geom_rect(
      aes(
        xmin = from,
        xmax = to,
        ymin = n - 0.005,
        ymax =  n + 0.005
      ),
      fill = "black",
      alpha = 0.6,
      data = plot_df
    ) +
    geom_rect(aes(
      xmin = from,
      xmax = to,
      ymin = n_m - (0.05 * (max(n_m) - min(n_m))) ,
      ymax = n_m + (0.05 * (max(n_m) - min(n_m)))
    ),fill = "red",
    data = plot_df
    ) +
    labs(title = "Counts along the genome", subtitle = "Input segmentation: bottom line.")
  
  if(library_size_norm) +
    p1 = labs(caption = "Counts normalised by library size.")
  
  # if(max(plot_df$ploidy_real) < 3.5) p1 = p1 + ylim(0,5)
  
  
  return(p1)
}
