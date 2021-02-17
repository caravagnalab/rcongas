#' Title
#'
#' @param x
#' @param chromosomes
#'
#' @return
#' @export
#'
#' @examples
get_input_segmentation = function(x,
                                  chromosomes = paste0("chr", c(1:22, "X", "Y")))
{
  df_segments = x$data$cnv %>% mutate(size = to - from)
  
  return(df_segments %>% dplyr::filter(chr %in% chromosomes))
}

