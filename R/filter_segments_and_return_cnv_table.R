#' Filter a cnv table 
#' 
#'The function takes and rcongas object as input and returns a cnv table with the filtered segments. 
#'
#' @param x rcongas object
#' @param mu number of genes to filter
#' @param length length of the segment to filter
#' @param merge merge continous segments with the same CNV value 
#'
#' @return a cnv_table with only the filtered segments
#' @export
#'
#' @examples
filter_segments_and_return_cnv_table <- function(x, mu = 15, length = 2.e7, merge = TRUE){
  
  x$data$cnv$dist <- x$data$cnv$to - x$data$cnv$from
  
  to_remove <- (x$data$cnv$mu < mu) | (x$data$cnv$dist < length)
  
  x$data$cnv <- x$data$cnv[!to_remove,]
  
  x$data$cnv$tot <- x$data$cnv$ploidy_real
  if(merge)
    x$data$cnv <-  Rcongas:::correct_bins(x$data$cnv , 0, FALSE)

  return(x$data$cnv)
  
}