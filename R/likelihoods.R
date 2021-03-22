NB_likelihood <- function(x, segment_id, cluster, modality = "RNA"){

  bm <-  x$best_fit
  segment_factors <- bm$segment_parameters %>%
    filter(modality == !!modality, parameter=="segment_factor", segment_id == !!segment_id) %>%
      pull(value)
  size <- bm$segment_parameters %>%
    filter(modality == !!modality, parameter=="NB_size", segment_id == !!segment_id) %>%
      pull(value)
  CNA <- bm$CNA %>% filter(cluster == !!cluster, segment_id == !!segment_id) %>% pull()

  ret <- function(linspace){
    return(dnbinom(linspace, mu = CNA * segment_factors, size = size))
  }

  return(ret)

}


G_likelihood <- function(x, segment_id, cluster, modality = "RNA"){

  bm <-  x$best_fit

  sd <- bm$segment_parameters %>%
    filter(modality == !!modality, parameter=="normal_sd", segment_id == !!segment_id) %>%
    pull(value)
  CNA <- bm$CNA %>% filter(cluster == !!cluster, segment_id == !!segment_id) %>% pull()

  ret <- function(linspace){
    return(dnorm(linspace, mean = CNA, sd = sd))
  }

  return(ret)

}



P_likelihood <- function(x, segment_id, cluster, modality = "RNA"){

  bm <-  x$best_fit
  segment_factors <- bm$segment_parameters %>%
    filter(modality == !!modality, parameter=="segment_factor", segment_id == !!segment_id) %>%
    pull(value)
  CNA <- bm$CNA %>% filter(cluster == !!cluster, segment_id == !!segment_id) %>% pull()

  ret <- function(linspace){
    return(dpois(linspace, lambda = CNA * segment_factors))
  }

  return(ret)

}
