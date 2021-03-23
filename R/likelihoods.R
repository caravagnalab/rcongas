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



assemble_likelihood_tibble <- function(x, segments){

  norm_factors = get_input(x, what = 'normalisation')

  atac_segs <- c()
  if(x %>%  has_atac){
    what_ATAC = get_input(x, what = 'data') %>%
      filter(modality == "ATAC")

    if (which_likelihood(x, "ATAC") != "G"){
      what_ATAC = normalise_modality(what_ATAC, norm_factors %>% filter(modality == "ATAC")) %>%
        group_by(segment_id) %>%  summarize(min = min(value), max = max(value))
    } else {
      what_ATAC = what_ATAC %>%
        group_by(segment_id) %>%  summarize(min = min(value), max = max(value))
    }


    atac_segs <- lapply(segments, function(s) assemble_likelihood_tibble_aux(x,
                                                                             what_ATAC %>%  filter(segment_id == s) %>%  pull(min),
                                                                             what_ATAC %>%  filter(segment_id == s) %>%  pull(max),
                                                                             s,
                                                                             "ATAC"
    ) )
    atac_segs <-  do.call(rbind, atac_segs) %>%  unique()
  }

  rna_segs <-  c()
  if(x %>%  has_rna){
    what_RNA = get_input(x, what = 'data') %>%
      filter(modality == "RNA")

    if (which_likelihood(x, "RNA") != "G"){
      what_RNA = normalise_modality(what_RNA, norm_factors %>% filter(modality == "RNA")) %>%
        group_by(segment_id) %>%  summarize(min = min(value), max = max(value))
    } else {
      what_RNA = what_RNA %>%
        group_by(segment_id) %>%  summarize(min = min(value), max = max(value))
    }


    rna_segs <- lapply(segments, function(s) assemble_likelihood_tibble_aux(x,
                                                                            what_RNA %>%  filter(segment_id == s) %>%  pull(min),
                                                                            what_RNA %>%  filter(segment_id == s) %>%  pull(max),
                                                                            s,
                                                                            "RNA"
    ) )
    rna_segs <-  do.call(rbind, rna_segs) %>%  unique()
  }

  ret <-  rbind(rna_segs, atac_segs)
  return(ret)


}


assemble_likelihood_tibble_aux <-  function(x,linspace_min, linspace_max, segment, modality){

  clusters <- get_cluster_assignments(x) %>%  filter(modality == !!modality) %>%  pull(cluster) %>%  unique()


  ret <-  lapply(clusters,FUN =function(cc) evaluate_likelihood(x,linspace_min, linspace_max, segment, modality, cc))

  ret <- do.call(rbind, ret)

  ret$modality <-  modality
  ret$segment <-  segment

  return(ret)

}


evaluate_likelihood <-  function(x,linspace_min, linspace_max, segment, modality, cluster){

  linspace <-  seq(linspace_min, linspace_max, length.out = 1000)


  if(which_likelihood(x, modality) == "G"){
    lk <- G_likelihood(x, segment, cluster, modality)

  } else if(which_likelihood(x, modality) == "P"){
    lk <- P_likelihood(x, segment, cluster, modality)
    linspace <- round(linspace)
    mult = 1000
  } else {
    lk <- NB_likelihood(x, segment, cluster, modality)
    linspace <- round(linspace)
  }
  lk <-  as.data.frame(lk(linspace))
  colnames(lk) <- "value"
  lk$X <- linspace
  lk$cluster <-  cluster
  return(lk)


}
