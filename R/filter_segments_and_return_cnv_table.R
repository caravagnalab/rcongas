#' Title
#'
#' @param x 
#' @param mu 
#' @param length 
#'
#' @return
#' @export
#'
#' @examples
filter_segments_and_return_cnv_table <- function(x, mu = 15, length = 2.e7){
  
  x$data$cnv$dist <- x$data$cnv$to - x$data$cnv$from
  
  to_remove <- (x$data$cnv$mu < mu) || (x$data$cnv$dist < length)
  
  x$data$cnv <- x$data$cnv[!to_remove,]
  
  x$data$cnv$tot <- x$data$cnv$ploidy_real
  x$data$cnv <-  Rcongas:::correct_bins(x$data$cnv, 0, FALSE)
  
  return(x$data$cnv)
  
}