filter_clusters  <- function(x, ncells = 10, abundance = 0.03) {

  x$inference$models <- lapply(x$inference$models, function(x) filter_clsuter_aux(x , ncells = ncells, abundance = abundance))

  return(x)
}


filter_clsuter_aux <- function(x, ncells, abundance) {

  if(length(x$parameters$mixture_weights) == 1) return(x)

  ta <-  table(x$parameters$assignement)

  diff_len <-  length(x$parameters$mixture_weights) - length(ta)

  if( diff_len != 0 ){ ta <- c(ta, rep(x = 0, diff_len))}


  mask <-  (x$parameters$mixture_weights > abundance) & (ta > ncells)

  nremoved <-  sum(mask)

  cli::cli_alert_warning("Filtering {nremoved} cluster{?s} due to low cell counts or abudance")
  if(nremoved == 0) return(x)

  x$parameters$mixture_weights <- x$parameters$mixture_weights[mask]
  cnv_probs_new <- x$parameters$cnv_probs[mask,]
  if(!x$run_information$posteriors){

    cli::cli_alert_info("No posterior probabilities, using euclidean distance to merge clusters")
    if(diff_len != 0){
      for(d in seq_len(diff_len)){
        x$parameters$assignment_probs <-  cbind(x$parameters$assignment_probs, rep(0, nrow(x$parameters$assignment_probs)))
      }
      colnames(x$parameters$assignment_probs) <- seq_along(mask)
    }
    x$parameters$assignment_probs <- x$parameters$assignment_probs[,mask]
    cnv_no_assign <- x$parameters$cnv_probs[!mask,, drop = FALSE]
    distance <- as.matrix(dist(as.matrix(x$parameters$cnv_probs), diag = T))

    for(c in which(!mask)){
      distance[c,c] <-  Inf
      ## !!! note, this works only because the labelling of clusters is in order of abudance (descending)
      x$parameters$assignement[x$parameters$assignement == c] <- which.min(distance[c,])

    }


  } else {

    cli::cli_alert_info("Reculcating cluster assignement and renormalizing posterior probabilities")


    x$parameters$assignment_probs <- x$parameters$assignment_probs[,mask]
    x$parameters$assignment_probs <- x$parameters$assignment_probs / rowSums(as.matrix(x$parameters$assignment_probs))
    x$parameters$assignement <- apply(as.matrix(x$parameters$assignment_probs), 1, which.max)

  }

  x$parameters$cnv_probs <- cnv_probs_new

  return(x)

}
