norm_factor <- function(data) {


  norm <- rowMeans((data$counts / data$cnv$mu) / data$cnv$ploidy_real)
  return(torch$tensor(norm))

}

log_sum_exp <- function(x) {

  c <-  max(x)
  x <-  x - c

  ret <- c + log(sum(exp(x)))

  return(ret)
}


calculate_entropy <- function(x) {
  tmp <- apply(x, 1, function(x) {
    x = x + 1e9
    x = x/sum(x)
    
    x * log(x) 
  }
              )
  -sum(tmp)
}


calculate_ICL <- function(inf, data, mu,llikelihood = gauss_lik) {

  BIC <- calculate_BIC(inf, data, mu,llikelihood)
  H <- calculate_entropy(inf$parameters$assignment_probs)
  return(BIC + H)
}


calculate_AIC <-  function(inf, data, mu,llikelihood = gauss_lik) {

  n_param <- param_total(inf$parameters)
  log_lik <- llikelihood(data,mu,inf$parameters)




  return(n_param * 2 - 2 * log_lik)

}


calculate_BIC <-  function(inf, data, mu,llikelihood = gauss_lik) {

  n_param <- param_total(inf$parameters)
  log_lik <- llikelihood(data,mu,inf$parameters)


  return(n_param * log(nrow(data)) - 2 * log_lik)


}

gauss_lik_norm <-  function(data,mu,par) {



  mixture_weights <- par$mixture_weights

  lambdas <- par$cnv_probs

  data <- as.matrix(data)


  sd <- par$norm_sd

  log_lk <- matrix(nrow = nrow(data), ncol = ncol(data))

  for(n in 1:nrow(data)){
    for(i in 1:ncol(data)){
      lk <-  vector(length = length(mixture_weights))
      for(k in 1:length(mixture_weights)){
        lk[k] <- (dnorm((data[n,i]),mean =  as.numeric(lambdas[k,i]), sd = as.numeric(sd[i]), log = TRUE)  + log(mixture_weights[k]))
      }
      log_lk[n,i] <- log_sum_exp(lk)
    }


  }

  return(sum(log_lk))
}


gauss_lik_with_means <-  function(data,mu,par) {

  mixture_weights <- par$mixture_weights


  lambdas <-  lapply(1:length(mixture_weights), function(x) {


    lambdas <- matrix(par$norm_factor, ncol = 1) %*% (as.numeric(par$cnv_probs[x,]) * mu * as.numeric(par$segment_factor[x]))

    return(lambdas)

  })

  log_lk <- matrix(nrow = nrow(data), ncol = ncol(data))

  for(n in 1:nrow(data)){
    for(i in 1:ncol(data)){
      lk <-  vector(length = length(mixture_weights))
      for(k in 1:length(mixture_weights))
        lk[k] <- dpois(as.numeric(data[n,i]), as.numeric(lambdas[[k]][n,i]), log = TRUE) + log(mixture_weights[k])
      log_lk[n,i] <- log_sum_exp(lk)
    }
  }


  log_lk <-  sum(log_lk)
}

param_total <-  function(param_list) {
  param_list <-  param_list[!(names(param_list) %in% c("assignment_probs", "assignement"))]
  res <- sapply(param_list,function(x) if(is.null(dim(x))) length(x) else prod(dim(x)))
  return(sum(res))
}

gauss_lik_old <-  function(data,mu,par) {

  mixture_weights <- par$mixture_weights


  lambdas <-  lapply(1:length(mixture_weights), function(x) {


    lambdas <- matrix(par$norm_factor, ncol = 1) %*% (as.numeric(par$cnv_probs[x,]) * mu)

    return(lambdas)

  })

  log_lk <- array(0, dim = c(length(mixture_weights), nrow(data), ncol(data)))

  for(n in 1:nrow(data)){
    for(i in 1:ncol(data)){
      lk <-  vector(length = length(mixture_weights))
      for(k in 1:length(mixture_weights))
        log_lk[k,n,i] <- dpois(as.numeric(data[n,i]), as.numeric(lambdas[[k]][n,i]), log = TRUE) 
    }
  }
  log_lk <- apply(log_lk,c(1,2), sum) + log(mixture_weights)
  log_lk <- log_sum_exp(log_lk)
  
  log_lk <-  sum(log_lk) / ncol(data)

  return(log_lk)
}


gauss_lik <-  function(data,mu,par) {

  mixture_weights <- par$mixture_weights


  lambdas <-  lapply(1:length(mixture_weights), function(x) {


    lambdas <- matrix(par$norm_factor, ncol = 1) %*% (as.numeric(par$cnv_probs[x,]) * mu) / (sum(as.numeric(par$cnv_probs[x,]) * mu ) / sum(mu))

    return(lambdas)

  })

  log_lk <- array(0, dim = c(length(mixture_weights), nrow(data), ncol(data)))
  
  for(n in 1:nrow(data)){
    for(i in 1:ncol(data)){
      lk <-  vector(length = length(mixture_weights))
      for(k in 1:length(mixture_weights))
        log_lk[k,n,i] <- dpois(as.numeric(data[n,i]), as.numeric(lambdas[[k]][n,i]), log = TRUE) 
    }
  }
  log_lk <- apply(log_lk,c(1,2), sum) 
  log_lk <- apply(log_lk, 2, function(x) x + log(mixture_weights))
  if(length(mixture_weights) > 1 )
   log_lk <- apply(log_lk, 2, function(x) log_sum_exp(x))
  
  log_lk <-  sum(log_lk) / ncol(data)

  return(log_lk)
}
