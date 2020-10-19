filter_clusters  <- function(x, ncells = 10, abundance = 0.03) {

  x$inference$models <- lapply(x$inference$models, function(x) filter_clsuter_aux(x , ncells = ncells, abundance = abundance))

  return(x)
}


filter_clsuter_aux <- function(x, ncells, abundance) {

  if(length(x$parameters$mixture_weights) == 1) return(x)

  mask <-  (x$parameters$mixture_weights > abundance) & (table(x$parameters$assignement) > ncells)

  nremoved <-  sum(mask)

  cli::cli_alert_warning("Filtering {no(nremoved)} cluster{?s} due to low cell counts or abudance")
  if(nremoved == 0) return(x)

  x$parameters$mixture_weights <- x$parameters$mixture_weights[mask]
  cnv_probs_new <- x$parameters$cnv_probs[mask,]
  if(is.null(x$parameters$assignment_probs)){

    cli::cli_alert_info("No posterior probabilities, using euclidean distance to merge clusters")

    cnv_no_assign <- x$parameters$cnv_probs[!mask,, drop = FALSE]
    for(c in which(!mask)){

      ## !!! note, this works only because the labelling of clusters is in order of abudance (descending)
      distance <- distance(cnv_no_assign[c,], cnv_probs_new)
      x$parameters$assignement[x$parameters$assignement == c] <- which.min(distance)
    }


  } else {

    cli::cli_alert_info("Reculcating cluster assignement and renormalizing posterior probabilities")


    x$parameters$assignment_probs <- x$parameters$assignment_probs[,mask]
    x$parameters$assignment_probs <- x$parameters$assignment_probs / rowSums(x$parameters$assignment_probs)
    x$parameters$assignement <- apply(x$parameters$assignment_probs, 1, which.max)

  }

  x$parameters$cnv_probs <- cnv_probs_new

  return(x)

}
