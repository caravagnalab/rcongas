#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
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
    geom_text(aes(y = lab.ypos, label = paste0(prop, '%\nn = ', n, '')), color = "white") +
    CNAqc:::my_ggplot_theme() +
    labs(x = NULL, y = NULL) +
    theme(axis.text.x = element_blank(), panel.border = element_blank()) +
    labs(title = paste0(
      "n = ",
      sum(clusters_size),
      " cells, k = ",
      length(clusters_size),
      ' clusters'
    ))


}