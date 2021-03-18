#' Get fit data.
#' 
#' @description General fit accessing getter function to return
#' segment_factors, CNAs (point estimates), posterior CNA (distribution),
#' mixing proportions, clustering assignments or the latent responsibilities.
#' 
#' @param x 
#' @param what Any of \code{"segment_factors"}, \code{"CNA"},  \code{"posterior_CNA"},
#'  \code{"mixing_proportions"},  \code{"cluster_assignments"},  or \code{"z_nk"}.
#'
#' @return A tibble.
#' @export
#'
#' @examples
get_fit = function(x, what = 'segment_factors')
{
  x %>% sanitize_obj()
  
  if(!('best_fit' %in% names(x)))
    stop("Missing fits for the input object, cannot extract ")
  
  if(what == 'segment_factors') return(x %>% get_segment_factors())
  if(what == 'CNA') return(x %>% get_CNA())
  if(what == 'posterior_CNA') return(x %>% get_posterior_CNA())
  if(what == 'mixing_proportions') return(x %>% get_mixing_proportions())
  if(what == 'cluster_assignments') return(x %>% get_cluster_assignments())
  if(what == 'z_nk') return(x %>% get_z_nk())

  stop("Unrecognised 'what': use any of 'segment_factors', 'CNA', 'posterior_CNA',
       'mixing_proportions', 'cluster_assignments, or 'z_nk'.")
}

get_segment_factors = function(x)
{
  x$best_fit$segment_factors
}

get_CNA = function(x)
{
  x$best_fit$CNA
}

get_posterior_CNA = function(x)
{
  x$best_fit$posterior_CNA
}

get_mixing_proportions = function(x)
{
  x$best_fit$mixing_proportions
}

get_cluster_assignments = function(x)
{
  x$best_fit$cluster_assignments
}

get_z_nk = function(x)
{
  x$best_fit$z_nk
}