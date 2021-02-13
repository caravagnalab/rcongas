#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
#' x = Rcongas::congas_example
#'
#' plot_mixing_proportions(x)
plot_mixing_proportions = function(x)
{
  clusters_size = Rcongas::get_clusters_size(x)
  p_clusters_size = Rcongas::get_clusters_size(x, normalised = TRUE)
  
  clusters_data <- data.frame(
    cluster = names(clusters_size),
    n = clusters_size,
    prop = round(p_clusters_size * 100, 2),
    stringsAsFactors = FALSE
  ) %>%
    arrange(desc(cluster)) %>%
    mutate(lab.ypos = cumsum(prop) - 0.5 * prop)
  
  
  ggplot(clusters_data,
         aes(x = "", y = prop, fill = cluster)) +
    geom_bar(width = 1,
             stat = "identity",
             color = "white") +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = Rcongas:::get_clusters_colors(names(clusters_size))) +
    geom_text(aes(
      y = lab.ypos,
      label = paste0(cluster, '-', prop, '%\nn = ', n, '')
    ), color = "white") +
    CNAqc:::my_ggplot_theme() +
    labs(x = NULL, y = NULL) +
    theme(
      axis.text.x = element_blank(),
      panel.border = element_blank(),
      legend.position = 'right'
    ) +
    labs(
      title = "Mixing proportions",
      subtitle = paste0(
        sum(clusters_size),
        " cells and ",
        length(clusters_size),
        ' clusters'
      )
    ) +
    guides(fill = guide_legend(''))
  
  
}