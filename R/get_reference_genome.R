#' Reference genome used in a model.
#' 
#' @description Extracts the code for the reference genome used in a model.
#'
#' @param x Input CONGAS model.
#'
#' @return A code for any of the supported references
#' 
#' @export
#'
#' @examples
#' 
#' x = Rcongas::congas_example
#'
#' get_gene_annotations(x)
get_reference_genome = function(x)
{
  return(x$reference_genome)
}