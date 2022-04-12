#' Retrieve the CN table
#'
#'A function to get the copy-number (CN) table actually stored in the Rcongas object.
#'
#' @param x Rcongas object
#' @param chromosomes filter for specific chromosomes
#'
#' @return The CN table stored by the Rcongas object
#' @export
#'
#' @examples
#' 
#' x = Rcongas::congas_example
#' 
#' get_input_segmentation(x)

get_input_segmentation = function(x,
                                  chromosomes = paste0("chr", c(1:22, "X", "Y")))
{
  df_segments = x$data$cnv %>% mutate(size = to - from)
  
  return(df_segments %>% dplyr::filter(chr %in% chromosomes))
}

