has_DE = function(x)
{
  return(!all(is.null(x$DE)))
}

is_MAP_Z <-  function(x){
  best_model <- get_best_model(x)
  if(is.null(best_model$run_information_old))
    return(best_model$run_information$posteriors)
  else
    return(best_model$run_information_old$posteriors)
}

is_MAP_CN <-  function(x){
  best_model <- get_best_model(x)
  return(best_model$run_information$MAP)

}

has_inference <-  function(x){

  return(inherits(x, 'rcongas') && !is.null(x$inference))
}