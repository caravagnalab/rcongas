

#' Return gene mapping on segments
#' 
#' The function returns a mapping of all the genes associated to the currently used segments
#'
#' @param x an Rcongas object
#' @param chromosomes filter for specific chromosomes
#'
#' @return a table with gene names and segment coordinates
#' @export
#'
#' @examples
#' 
#' x = Rcongas::congas_example
#' get_mapped_genes(x)
#' 
get_mapped_genes = function(x,
                            chromosomes = paste0("chr", c(1:22, "X", "Y")))
{
  x$data$gene_locations %>%
    # dplyr::select(-segment_id) %>%
    dplyr::filter(chr %in% chromosomes)
}




