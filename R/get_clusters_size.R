#' Extract size of clusters.
#'
#' @description Compute the size of CONGAS clusters, normalised or in absolute
#' cell numbers.
#'
#' @param x Input object with clusters.
#' @param normalised If \code{TRUE}, normalises the proportions.
#'
#' @return A named vector for cluster sizes.
#' @export
#'
#' @examples
#'
#' x = Rcongas::congas_example
#'
#' get_clusters_size(x, normalised = FALSE)
#' get_clusters_size(x, normalised = TRUE)
get_clusters_size = function(x, normalised = FALSE)
{
  if(!has_inference(x)) stop("Cannot extract clustering information if not avaiable, or not a CONGAS object.")

  best_model = get_best_model(x)

  # Counts from the assignments
  n = table(best_model$parameters$assignement)

  # Make it a named vector
  v = as.vector(n)
  names(v) = names(n)

  if (normalised)
    v = v / sum(v)

  return(v)
}