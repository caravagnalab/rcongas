#' Extract the number of mixture components.
#'
#' @description Returns \code{k>0} for the best mixture available in the
#' object.
#'
#' @param x Input object with clusters.
#'
#' @return A scalar.
#' @export
#'
#' @examples
#' x = Rcongas::congas_example
#' get_k(x)
get_k = function(x)
{
  if(!has_inference(x)) stop("Cannot extract clustering information if not avaiable, or not a CONGAS object.")

  best = get_best_model(x)
  best$parameters$assignement %>%
    unique() %>%
    length()
}