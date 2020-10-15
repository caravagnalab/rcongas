filter_clusters  <- function(x, ncells = 10, abundance = 0.03) {

  x$inference$models <- lapply(x$inference$models, function(x) filter_clsuter_aux(x , ncells = ncells, abundance = abundance))

  return(x)
}


filter_clsuter_aux <- function(x, ncells, abundance) {

  if(length(x$parameters$mixture_weights) == 1) return(x)

  mask <-  (x$parameters$mixture_weights > abundance) & (table(x$parameters$assignement) > ncells)
  x$parameters$mixture_weights <- x$parameters$mixture_weights[mask]
  x$parameters$cnv_probs <- x$parameters$cnv_probs[mask,]
  if(is.null(x$parameters$assignment_probs)){
  } else {

    x$parameters$assignment_probs <- x$parameters$assignment_probs[,mask]
    x$parameters$assignment_probs <- x$parameters$assignment_probs / rowSums(x$parameters$assignment_probs)
    x$parameters$assignement <- apply(x$parameters$assignment_probs, 1, which.max)

  }

  return(x)

}
