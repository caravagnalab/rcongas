#' Extract fit data.
#'
#' @description This is like function \code{\link{get_input}}, but for fit information.
#' With this getter you can obtain:
#'
#' * segment parameters (scaling factors and component-specific
#' values for the used distributions);
#' * inferred Copy Number Alteration (CNA) values;
#' * the posterior distribution over CNAs;
#' * mixing proportions;
#' * clustering assignments;
#' * z_nk (the latent variables of the model)
#'
#' Like \code{get_input}, the function uses the \code{what} parameter to return the
#' appropriate type of information.
#'
#' @param x An object of class \code{rcongasplus}.
#' @param what Any of \code{"segment_parameters"}, \code{"CNA"},  \code{"posterior_CNA"},
#'  \code{"mixing_proportions"},  \code{"cluster_assignments"},  or \code{"z_nk"}.
#'
#' @return A tibble; its format depends on \code{what}. See the examples.
#' @export
#'
#' @examples
#' data(example_object)
#'
#' # Extract segment parameters
#' get_fit(example_object, what = 'segment_parameters')
#'
#' # Extract CNAs
#' get_fit(example_object, what = 'CNA')
#'
#' # Extract CNAs
#' get_fit(example_object, what = 'posterior_CNA')
#'
#' # Extract mixing proportions
#' get_fit(example_object, what = 'mixing_proportions')
#'
#' # Extract clustering assignments
#' get_fit(example_object, what = 'cluster_assignments')
#'
#' # Extract the clustering responsibilities
#' get_fit(example_object, what = 'z_nk')
get_fit = function(x, what = 'CNA')
{
  x %>% sanitize_obj()

  if(!('best_fit' %in% names(x)))
    stop("Missing fits for the input object, cannot extract ")

  if(what == 'segment_parameters') return(x %>% get_segment_parameters())
  if(what == 'CNA') return(x %>% get_CNA())
  if(what == 'CNA_real') return(x %>% get_CNA_real())
  if(what == 'posterior_CNA') return(x %>% get_posterior_CNA())
  if(what == 'mixing_proportions') return(x %>% get_mixing_proportions())
  if(what == 'cluster_assignments') return(x %>% get_cluster_assignments())
  if(what == 'z_nk') return(x %>% get_z_nk())

  stop("Unrecognised 'what': use any of 'segment_factors', 'CNA', 'posterior_CNA',
       'mixing_proportions', 'cluster_assignments, or 'z_nk'.")
}

get_segment_factors = function(x)
{
  x$best_fit$segment_parameters
}

get_CNA = function(x)
{
  x$best_fit$CNA
}

get_CNA_real = function(x)
{
  x$best_fit$CNA_real
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
