

#' Title
#'
#' @param x 
#' @param chromosomes 
#'
#' @return
#' @export
#'
#' @examples
get_mapped_genes = function(x,
                            chromosomes = paste0("chr", c(1:22, "X", "Y")))
{
  x$data$gene_locations %>%
    # dplyr::select(-segment_id) %>%
    dplyr::filter(chr %in% chromosomes)
}




