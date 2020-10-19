has_DE = function(x)
{
  return(!all(is.null(x$DE)))
}

is_MAP_Z <-  function(x){
  best_model <- get_best_model(x)
  return(is.null(best_model$parameters$assignment_probs))
}

is_MAP_CN <-  function(x){
  best_model <- get_best_model(x)
  return(is.null(best_model$parameters$cnv_var))
}

has_inference <-  function(x){

  return(!is.null(x$inference))
}