#' Title
#'
#' @param x
#' @param cutoff_assignment
#'
#' @return
#' @export
#'
#' @examples
plot_latent_variables = function(x, cutoff_assignment = 0)
{
  assignments = Rcongas::get_clusters(x, cut_znk = cutoff_assignment) %>%
    dplyr::arrange(cluster, p_assignment)

  clusters_table = assignments %>%
    dplyr::select(-p_assignment,-cell)

  clusters_names = names(Rcongas::get_clusters_size(x))

  not_assign = is.na(assignments$cluster)
  n = sum(not_assign)
  p = round((n / nrow(assignments)) * 100)
  lv = reshape2::melt(
    assignments %>% dplyr::select(tidyselect::all_of(clusters_names)) %>%
      dplyr::mutate(pos = row_number()),
    id = "pos"
  )
  colnames(lv) = c("Point", "Cluster", "Value")

  lv <-  lv %>% dplyr::arrange(Cluster, Value)

  ggplot(lv, aes(x = Point, y = Cluster, fill = Value)) + geom_raster() +
    scale_fill_viridis_c(direction = -1) + mobster:::my_ggplot_theme() +
    guides(
      fill = guide_colorbar(bquote(z["nk"] ~ " "), barwidth = unit(3,"cm"))
      ) +
    labs(
      title = bquote("Latent variables"),
      subtitle = bquote(z["nk"] ~ " > " * .(cutoff_assignment) ~ ": n =" ~ .(n) ~ " NAs (" * .(p) * "%)"),
      y = paste0("Points (n =", nrow(assignments), ")")
    )

}
