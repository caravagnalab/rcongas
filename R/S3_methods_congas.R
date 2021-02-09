#' @export
`[.congas` <- function(x,i,j) {
  
  
  x$parameters$cnv_probs <- x$parameters$cnv_probs[,j, drop = FALSE]
  x$parameters$norm_factor <- x$parameters$norm_factor[i]
  x$parameters$assignement <- x$parameters$assignement[i]
  if(is_MAP_Z(x, congas_obj = TRUE)){
    x$parameters$assignment_probs <-  x$parameters$assignment_probs[i,, drop = FALSE]
  }
  x <- clean_clusters(x)
  
  
}

 