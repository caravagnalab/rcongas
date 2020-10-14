get_assignment_probs <- function(x) {
  best_model <-  get_best_model(x)
  if(is_MAP_Z(x)){

    return(best_model$parameters$assignment_probs)
  } else{
    cli::cli_alert_warning("Posteriors probabilities where not calculated in the previous step, please set posterior=TRUE")
    return(NULL)
  }
}

get_CNV_mean <- function(x) {
  best_model <-  get_best_model(x)
  return(best_model$parameters$cnv_probs)
}

get_CNV_sd <- function(x) {
  best_model <-  get_best_model(x)
  if(is_MAP_CNV){
    cli::cli_alert_warning("CNV values where learned with MAP, please re-run the analysis with MAP=FALSE if you want to do otherwise")
    return(NULL)
  }else{
    return(best_model$parameters$cnv_var)
  }
}


get_cluster_assignments <- function(x) {
  best_model <-  get_best_model(x)
  return(paste(best_model$parameters$assignement))

}

