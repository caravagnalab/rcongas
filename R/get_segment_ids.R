#' Extract the segment ids (for indexing).
#' 
#' @description For the functions that require a segment identifier as input,
#' for instance \code{plot_segment_density},
#' this function returns the list of segment ids used in the input object. If 
#' run with \code{highlight = TRUE} it will return only the segments that
#' are highlighted in the plots; this requires clustering to be computed
#'
#' @param x Input object (with or without clustering).
#' @param highlight If \code{TRUE}, returns return only the segments that are
#' highlighted in the plots. Use ellipsis to specify parameters such as \code{alpha}
#' to function \code{get_clusters_ploidy}.
#' @param ... Parameters forwarded to \code{get_clusters_ploidy} to get segment
#' ids.
#'
#' @return A vector of strings.
#' 
#' @export
#'
#' @examples
#' 
#' x = Rcongas::congas_example
#' 
#' # All segments
#' get_segment_ids(x)
#' 
#' # All segments that are highlighted with default value
#' get_segment_ids(x, highlight = TRUE)
#' 
#' # All segments that are highlighted with stricter value than default
#' get_segment_ids(x, highlight = TRUE, alpha = 0.01)
get_segment_ids = function(x, highlight = FALSE, ...)
{
  if (!highlight)
    return(get_input_segmentation(x) %>%
             idify() %>%
             pull(segment_id)) %>% 
    unique
  
  if(!has_inference(x)) stop("Input is missing  clustering information or not a CONGAS object.")
  
  get_clusters_ploidy(x, ...) %>%
    idify() %>% 
    filter(highlight) %>% 
    pull(segment_id) %>% 
    unique
}
