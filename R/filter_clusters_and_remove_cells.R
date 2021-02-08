#' Filter output clusters and remove cells beloging to those clusters 
#'
#'The function works in a similar way to \code{\link{filter_clusters}} but it also removes the cells belonging to those clusters 
#'
#' @param x Rcongas object
#' @param ncells minimum size of a valid cluster expressed as absolute number of cells
#' @param abundance minimum size of a cluster expressed as a percentage of the total cell number 
#'
#' @return
#' @export
#'
#' @examples
#' @examples
#'
#' x = Rcongas::congas_example
#'
#' print(x)
#'
#' x = x %>% Rcongas:::filter_clusters_and_remove_cells(ncells = 150)
#' print(x)
#' 

filter_clusters_and_remove_cells <- function(x, ncells = 10, abundance = 0.03){
  
  x$inference$models <-
    lapply(x$inference$models, function(x)
      filter_clusters_and_remove_cells_aux(x , ncells = ncells, abundance = abundance))
  #x$inference$model_selection$IC <- recalculate_information_criteria(x, x$inference$model_selection$IC_type)
  return(x)
  
}


filter_clusters_and_remove_cells_aux <- function(x, ncells, abundance){

  
  
  if(length(x$parameters$mixture_weights) == 1) return(x)
  
  ta <-  table(x$parameters$assignement)
  mask <-  (x$parameters$mixture_weights > abundance) & (ta > ncells)
  
  mask <-  (x$parameters$mixture_weights > abundance) & (ta > ncells)
  
  # This is important
  if(all(mask)) return(x)
  
  nremoved <-  sum(mask)
  cli::cli_alert_warning("Filtering {nremoved} cluster{?s} due to low cell counts or abudance")
  
  
  # This cannot work if you do not HAVE to remove at least one cluster 
  to_remove <- names(x$parameters$mixture_weights[!mask])
  
  # Removal
  return(x[-which(x$parameters$assignement %in% to_remove),])
  
}