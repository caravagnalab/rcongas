estimate_segment_factors <- function(data, norm_factors, pld, plot = T){

  x <-  apply(data, 2, function(y) y/  norm_factors)
  x <-  apply(x, 1, function(y) y / pld) %>%  t


  res <-  lapply(1:ncol(x),function(y) {estimate_segment_factors_aux(x[,y], plot = plot)})

  names(res) <-  colnames(data$data$counts)
  return(res)
}

estimate_segment_factors_aux <- function(data_mle, plot){
  data_mle <-  ifelse(data_mle == 0, rnorm(1,1e2, 1e2), data_mle)

  LL <-  function(theta_shape, theta_rate){

    R = dgamma(data_mle,shape = theta_shape, rate = theta_rate, log = TRUE)
    return(-sum(R))
  }

  coeff <-  stats4::mle(minuslogl = LL, start = list(theta_shape = mean(data_mle), theta_rate = 1), ,
                        lower = c(1e-16, 1e-16), upper = c(Inf, Inf))@coef
  if(plot){

    par(mfrow = c(5,5))


    hist(data_mle, prob = TRUE, col = "grey" ,main= "Inferred segment prior", xlab = "", breaks = 50)
    curve(dgamma(x,shape = coeff[1], rate = coeff[2] ), from = min(data_mle), to = max(data_mle),
          n = 1000, col = "blue", add = TRUE, main= "", xlab = "", lwd = 2)
  }

  return(coeff)
}


