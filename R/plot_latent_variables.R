#' Plot latent varibles.
#' 
#' @description Plots the clustering responsibilities (posterior over $z$)
#' for the available clusters. It can put \code{NA} to all entries where the hard 
#' clustering assignments - max of the posterior - is below a cutoff.
#'
#' @param x Input object with clusters
#' @param cutoff_assignment A value in $[0, 1)$ so that values below this are
#' set to \code{NA}
#'
#' @return A ggplot object
#' 
#' @export
#'
#' @examples
#' x = Rcongas::congas_example
#' 
#' # Default (cutoff 0)
#' plot_latent_variables(x) 
#' 
#' # >50% assignment probability
#' plot_latent_variables(x, cutoff_assignment = .5) 
#' 
#' # >90% assignment probability
#' plot_latent_variables(x, cutoff_assignment = .9) 
#' 
#' # >97.5% assignment probability
#' plot_latent_variables(x, cutoff_assignment = .975) 
plot_latent_variables = function(x, cutoff_assignment = 0)
{
  if(!has_inference(x)) stop("Cannot extract clustering information if not avaiable, or not a CONGAS object.")
  
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
