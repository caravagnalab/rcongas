#' Approximation of the library size prior
#'
#'This function fit a theta 
#'
#' @param x 
#'
#' @return
#' @export
#'
#' @examples
approx_theta_prior_params <-  function(x, plot = FALSE, mult_factor = 1){
  data_mle = tranform_data_for_theta_inf(x)
  
  LL <-  function(theta_shape, theta_rate){
    
    R = dgamma(data_mle,shape = theta_shape, rate = theta_rate, log = TRUE)
    return(-sum(R))
  }
  
  coeff <-  stats4::mle(minuslogl = LL, start = list(theta_shape = 1, theta_rate = 1), ,
                lower = c(0, 0), upper = c(Inf, Inf))@coef * mult_factor
  if(plot){
    
    hist(data_mle, prob = TRUE, col = "grey" ,main= "Inferred theta prior", xlab = "")
    curve(dgamma(x,shape = coeff[1], rate = coeff[2] ), from = min(data_mle), to = max(data_mle), 
          n = 1000, col = "blue", add = TRUE, main= "", xlab = "", lwd = 2)
  }
    
  return(coeff)
}



tranform_data_for_theta_inf <- function(x){
  
  data_t <- apply( x$data$counts, 1, function(y) y/x$data$cnv$mu)
  data_t <- apply(data_t, 2, function(y) y/ x$data$cnv$ploidy_real) %>% t
  return(data_t + 1e-6)
  
}