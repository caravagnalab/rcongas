estimate_segment_factors <-
  function(data, norm_factors, pld, plot = T) {
    
    cli::cli_h3("Estimating segment factors")
    
    x <-  apply(data, 2, function(y)
      y /  norm_factors) 
    
    if (ncol(x) > 1) {
      x <-  apply(x, 1, function(y)
        y / pld) %>% t
      
    } else{
      x <- x / pld
    }
    
    res <-
      lapply(1:ncol(x), function(y) {
        w = estimate_segment_factors_aux(data_mle = x[, y], plot = plot)
        cli::cli_alert("{y}: {.field {colnames(x)[y]}} theta_shape = {w['theta_shape']}, theta_rate = {w['theta_rate']}")
  
        w
      })
    
    names(res) <-  colnames(data)
    return(res)
  }

estimate_segment_factors_aux <- function(data_mle, plot) 
{
  # data_mle = x[, 2]
  # hist(data_mle, breaks = 100)
  
  # data_mle <-  ifelse(data_mle == 0, rnorm(1, 1e2, 1e2), data_mle)
  data_mle = data_mle[data_mle > 0]
  data_mle = data_mle[!is.na(data_mle)]
    
  quants = quantile(data_mle, probs = c(0.005, 0.995), na.rm = TRUE)
  
  data_mle <-  ifelse(data_mle < quants[1], quants[1], data_mle)
  data_mle <-  ifelse(data_mle > quants[2], quants[2], data_mle)
  
  LL <-  function(theta_shape, theta_rate) {
    R = dgamma(data_mle,
               shape = theta_shape,
               rate = theta_rate,
               log = TRUE)
    return(-sum(R))
  }
  
  # library(MASS)
  # fitdistr(data_mle, "gamma", start=list(shape=1, rate=1))$estimate
  
  
  MAXT = 10
  ct = 0
  while (ct < MAXT) {
    tryCatch({
      coeff <-
        stats4::mle(
          minuslogl = LL,
          start = list(
            theta_shape = mean(data_mle, na.rm = TRUE),
            theta_rate = sd(data_mle, na.rm = TRUE)
          ),
          lower = c(1e-16, 1e-16),
          upper = c(Inf, Inf)
        )@coef
      ct = MAXT
    }, ct = ct + 1)
  }
  
  if (plot) {
    hist(
      data_mle,
      prob = TRUE,
      col = "grey" ,
      main = "Inferred segment prior",
      xlab = "",
      breaks = 500
    )
    curve(
      dgamma(x, shape = coeff[1], rate = coeff[2]),
      from = min(data_mle),
      to = max(data_mle),
      n = 1000,
      col = "blue",
      add = TRUE,
      main = "",
      xlab = "",
      lwd = 2
    )
  }
  
  return(coeff)
}
