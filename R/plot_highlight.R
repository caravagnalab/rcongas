#' Plot highlighted segments
#'
#' @param x
#' @param alpha
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
#' # Default view
#' plot_highlights(x)
#'
#' # More stringent
#' plot_highlights(x, alpha = 0.01)
#'
#' # Even more stringent
#' plot_highlights(x, alpha = 0.001)
plot_highlights = function(x, alpha = 0.05)
{
  stopifnot(has_inference(x))
  
  tests_table = highlights(x, alpha) %>%
    mutate(label = paste(cluster , 'vs', versus))
  
  nh = tests_table %>% filter(highlight) %>% nrow
  
  if(nh == 0) return(CNAqc:::eplot())
  
  # Density points
  params = tests_table %>% distinct(cluster, versus, shape, rate)
  
  df_points = lapply(1:nrow(params), function(i) {
    x_ij = tests_table %>%
      filter(cluster == params$cluster[i], versus == params$versus[i])
    
    domain_x = seq(0, max(x_ij$diff) * 1.25, length.out = 30)
    image_y = dgamma(domain_x,
                     params$shape[i],
                     params$rate[i])
    
    data.frame(
      x = domain_x,
      y = image_y,
      cluster = params$cluster[i],
      versus = params$versus[i]
    )
  }) %>%
    Reduce(f = bind_rows) %>%
    as_tibble() %>%
    mutate(label = paste(cluster , 'vs', versus))
  
  
  ann = tests_table %>% 
    filter(highlight) %>% 
    group_by(label) %>% 
    summarise(slabel = paste(segment_id, collapse = '\n'), label = label) %>% 
    distinct(label, slabel)

  
  # Plot
  ggplot(tests_table,
         aes(diff, fill = highlight)) +
    geom_histogram(bins = 30) +
    facet_wrap(~ label) +
    CNAqc:::my_ggplot_theme() +
    scale_fill_manual(bquote("Above " * q[alpha] * ""),
                      values = c(`TRUE` = 'indianred3', `FALSE` = 'gray')) +
    labs(
      x = "Absolute CN difference",
      y = "Observations",
      title = bquote(.(nh)~"highlight(s) (" * q[alpha] * " = " * .(1 - alpha) *
                       ")")
    ) +
    geom_line(
      data = df_points,
      aes(x = x, y = y),
      color = 'black',
      inherit.aes = F,
      size = .5
    ) +
    ggrepel::geom_text_repel(
      data = ann,
      aes(x = 0, y = Inf, label = slabel),
      inherit.aes = F,
      size = 2,
      hjust = 0,
      color = 'indianred3'
    )
}