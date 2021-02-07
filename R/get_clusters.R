#' Extract cell clustering assignments.
#'
#' @description Extract a tibble for clustering assignments of the input cells,
#' possibly subsetting by cluster or leaving \code{NA} those cell assignments
#' that have probability ($z_{nk}$) below a cutoff (via latent variables).
#'
#' @param x Input object with clusters.
#' @param cluster_label Cluster id to subset the assignments.
#' @param cut_znk Cutoff for latent variables.
#'
#' @return A tibble with the cell id, cluster labels, latent variables values
#' for every clustering id, and \code{p_assignment} the hard-clustering assignment
#' probability for the cells (max among the latent variables).
#'
#' @export
#'
#' @examples
#'
#' x = Rcongas::congas_example
#'
#' # Default view
#' get_clusters(x)
#'
#' # Assignments to specific cluster
#' get_clusters(x, clusters = "c1")
#'
#' # Set to NA all clustering assignments with p_assignment below a 99%
#' get_clusters(x, cut_znk = .99)
get_clusters = function(x,
                        clusters = NULL,
                        cut_znk = 0)
{
  if(!has_inference(x)) stop("Cannot extract clustering information if not avaiable, or not a CONGAS object.")

  # TODO - take all assingments with z_nk > c, c >= 0
  best_model = Rcongas:::get_best_model(x)

  clusters_table = data.frame(
    cell = names(best_model$parameters$assignement),
    cluster = paste(best_model$parameters$assignement),
    stringsAsFactors = FALSE
  ) %>%
    as_tibble()

  # Latent variables
  z_nk = best_model$parameters$assignment_probs %>%  as.data.frame()
  colnames(z_nk) = Rcongas::get_clusters_size(x) %>% names

  z_nk$p_assignment = as.double(unlist(apply(z_nk, 1, function(x) {
    max(x, na.rm = TRUE)
  })))


  # This is where the error happens.
  z_nk$cell = names(best_model$parameters$assignement)
  
  # which(!(names(best_model$parameters$assignement) %in% clusters_table$cell))
  
  clusters_table = clusters_table %>%
    full_join(z_nk, by = "cell") %>%
    dplyr::mutate(cluster = ifelse(p_assignment > cut_znk, cluster, NA)) %>%
    as_tibble()

  if (!is.null(clusters))
    clusters_table = clusters_table %>% dplyr::filter(cluster %in% clusters)

  return(clusters_table)
}