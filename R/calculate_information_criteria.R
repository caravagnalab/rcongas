
#' Title
#'
#' @param x 
#' @param method 
#'
#' @return
#' @export
#'
#' @examples
#' 
recalculate_information_criteria <-  function(x, method){
  #if(!has_inference(x)) stop("Cannot calculate ICs prior to inference")
  x$inference$model_selection$IC <- calculate_information_criteria(x$inference$models, x, method)
  x$inference$model_selection$IC_type <- method
  return(x)  
}


calculate_information_criteria <- function(res, x, method){
  
  model <-  res[[1]]$run_information$model
  lik_fun <- choose_likelihood(model)
  
  if (method == "BIC"){
    IC <-
      sapply(res, function(y)
        calculate_BIC(y, x$data$counts, x$data$cnv$mu, llikelihood = lik_fun))
  } else if (method == "AIC") {
    IC <-
      sapply(res, function(y)
        calculate_AIC(y, x$data$counts, x$data$cnv$mu, llikelihood = lik_fun))
  } else if (method == "ICL") {
    IC <-
      sapply(res, function(y)
        calculate_ICL(y, x$data$counts, x$data$cnv$mu, llikelihood = lik_fun))
  } else if (method == "lk") {
    IC <-
      sapply(res, function(y)
        lik_fun(x$data$counts, x$data$cnv$mu, y$parameters) * -1)
  } else {
    stop("Information criterium not present in the package")
  }
  
  return(IC)
}


choose_likelihood <- function(model) {
  if (grepl(tolower(model), pattern = "norm")) {
    lik_fun <-  gauss_lik_norm
  } else if (grepl(tolower(model), pattern = "Categorical")) {
    lik_fun <- gauss_lik_with_means
  } else if (grepl(tolower(model), pattern = "Old")) {
    lik_fun <- gauss_lik_old
  } else {
    lik_fun <-  gauss_lik
  }
}