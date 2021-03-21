tensorize <-  function(x){

  torch <- reticulate::import("torch")
  ret <-  torch$tensor(x)
  return(ret)
}

detensorize <-  function(x){

  if(typeof(x) == "environment"){
    return(x$detach()$numpy())
  }
  return(x)
}


input_data_from_rcongas <- function(x){
  ret <-  list()
  if(has_atac(x)){
    ret$data_atac <- get_data(x) %>%  filter(modality == "ATAC") %>%  reshape2::acast(segment_id  ~ cell, value.var = "value")
    ret$data_atac <- ret$data_atac[order(rownames(ret$data_atac)),order(colnames(ret$data_atac))]
    norm_tmp <- get_normalisation(x) %>%  filter(modality == "ATAC")
    norm_factor_atac <- norm_tmp$normalisation_factor
    names(norm_factor_atac) <- norm_tmp$cell
    ret$norm_factor_atac <- norm_factor_atac[colnames(ret$data_atac)]
  }

  if(has_rna(x)){
    ret$data_rna <- get_data(x) %>%  filter(modality == "RNA") %>%  reshape2::acast(segment_id  ~ cell, value.var = "value")
    ret$data_rna <- ret$data_rna[order(rownames(ret$data_rna)),order(colnames(ret$data_rna))]
    norm_tmp <- get_normalisation(x) %>%  filter(modality == "RNA")
    norm_factor_rna <- norm_tmp$normalisation_factor
    names(norm_factor_rna) <- norm_tmp$cell
    ret$norm_factor_rna <- norm_factor_rna[colnames(ret$data_rna)]
  }

  segs_tmp <- get_segmentation(x) %>% select(copies, segment_id)
  segs <- segs_tmp$copies
  names(segs) <- segs_tmp$segment_id

  ret$pld <- segs[order(names(segs))]

  ret <-  lapply(ret, tensorize)

  ret$segments <- as.integer(length(ret$pld))

  return(ret)

}



format_best_model <-  function(x, inf){
  ret <-  list()

  hyperpars <- lapply(inf$hyperparameters, detensorize)

  ret$parameters <-  list(ICs = inf$model_selection, hyperparameters = hyperpars)

  cluster_names <- paste0("C", 1:inf$hyperparameters$K)
  segment_names <-  get_segmentation(x) %>%  pull(segment_id) %>%  sort()


  mixing_proportions_atac <-  c()
  segment_parameters_atac <-  c()
  z_nk_atac <- c()
  cluster_assignments_atac <-  c()



  if(has_atac(x)){

    ### Cell names for ATAC

    cell_names_atac <- get_data(x) %>% filter(modality == "ATAC") %>%  pull(cell)  %>%  unique() %>%  sort()


    ### Mixing proportions for ATAC

    mixing_proportions_atac <- inf$inferred_params$mixture_weights_atac
    names(mixing_proportions_atac) <- cluster_names

    z_nk_atac <- detensorize(inf$inferred_params$assignment_probs_atac) %>%  t() %>%  as.data.frame()
    colnames(z_nk_atac) <- cluster_names

    z_nk_atac$cell <- cell_names_atac

    cluster_assignments_atac <-  detensorize(inf$inferred_params$assignment_atac)
    names(cluster_assignments_atac) <- cell_names_atac

    segment_factors_atac <- inf$inferred_params$segment_factor_atac
    names(segment_factors_atac) <- segment_names
    segment_factors_atac <-  segment_factors_atac  %>%  as.data.frame() %>%  tibble::rownames_to_column("segment_id")
    colnames(segment_factors_atac) <- c("segment_id", "value")
    segment_factors_atac$parameter <- "segment_factor"

    norm_sd_atac <-  c()
    if(which_likelihood(x, modality = "ATAC") == "G"){
      norm_sd_atac <- inf$inferred_params$norm_sd_atac
      names(norm_sd_atac) <- segment_names
      norm_sd_atac <-  norm_sd_atac %>% as.data.frame() %>%  tibble::rownames_to_column("segment_id")
      colnames(norm_sd_atac) <- c("segment_id", "value")
      norm_sd_atac$parameter <- "normal_sd"
    }
    NB_size_atac <- c()
    if(which_likelihood(x, modality = "ATAC") == "NB"){
      NB_size_atac <- inf$inferred_params$NB_size_atac
      names(NB_size_atac) <- segment_names
      NB_size_atac <-  NB_size_atac %>% as.data.frame() %>%  tibble::rownames_to_column("segment_id")
      colnames(NB_size_atac) <- c("segment_id", "value")
      NB_size_atac$parameter <- "NB_size"
    }

    segment_parameters_atac <- rbind(segment_factors_atac,norm_sd_atac,NB_size_atac) %>%  as_tibble()
    segment_parameters_atac$modality <-  "ATAC"

    mixing_proportions_atac <-  mixing_proportions_atac %>% as.data.frame()  %>%  tibble::rownames_to_column("cluster")
    colnames(mixing_proportions_atac)[2] <- "mixing"

    z_nk_atac <- z_nk_atac %>% as_tibble()  %>% mutate(modality = "ATAC")

    cluster_assignments_atac <-  cluster_assignments_atac %>% as.data.frame()  %>%  tibble::rownames_to_column("cell") %>%  as_tibble()
    colnames(cluster_assignments_atac)[2] <- "cluster"
    cluster_assignments_atac$cluster <-  paste0("C",cluster_assignments_atac$cluster + 1)
    cluster_assignments_atac$modality <-"ATAC"

  }

  if(has_rna(x)){

    ### Cell names for RNA

    cell_names_rna <- get_data(x) %>% filter(modality == "RNA") %>%  pull(cell)  %>%  unique() %>%  sort()


    ### Mixing proportions for RNA

    mixing_proportions_rna <- inf$inferred_params$mixture_weights_rna
    names(mixing_proportions_rna) <- cluster_names

    z_nk_rna <- detensorize(inf$inferred_params$assignment_probs_rna) %>%  t() %>%  as.data.frame()
    colnames(z_nk_rna) <- cluster_names

    z_nk_rna$cell <- cell_names_rna

    cluster_assignments_rna <-  detensorize(inf$inferred_params$assignment_rna)
    names(cluster_assignments_rna) <- cell_names_rna

    segment_factors_rna <- inf$inferred_params$segment_factor_rna
    names(segment_factors_rna) <- segment_names
    segment_factors_rna <-  segment_factors_rna  %>%  as.data.frame() %>%  tibble::rownames_to_column("segment_id")
    colnames(segment_factors_rna) <- c("segment_id", "value")
    segment_factors_rna$parameter <- "segment_factor"

    norm_sd_rna <-  c()
    if(which_likelihood(x, modality = "RNA") == "G"){
      norm_sd_rna <- inf$inferred_params$norm_sd_rna
      names(norm_sd_rna) <- segment_names
      norm_sd_rna <-  norm_sd_rna %>% as.data.frame() %>%  tibble::rownames_to_column("segment_id")
      colnames(norm_sd_rna) <- c("segment_id", "value")
      norm_sd_rna$parameter <- "normal_sd"
    }
    NB_size_rna <- c()
    if(which_likelihood(x, modality = "RNA") == "NB"){
      NB_size_rna <- inf$inferred_params$NB_size_rna
      names(NB_size_rna) <- segment_names
      NB_size_rna <-  NB_size_rna %>% as.data.frame() %>%  tibble::rownames_to_column("segment_id")
      colnames(NB_size_rna) <- c("segment_id", "value")
      NB_size_rna$parameter <- "NB_size"
    }

    segment_parameters_rna <- rbind(segment_factors_rna,norm_sd_rna,NB_size_rna) %>%  as_tibble()
    segment_parameters_rna$modality <-  "RNA"

    mixing_proportions_rna <-  mixing_proportions_rna %>% as.data.frame()  %>%  tibble::rownames_to_column("cluster")
    colnames(mixing_proportions_rna)[2] <- "mixing"

    z_nk_rna <- z_nk_rna %>% as_tibble()  %>% mutate(modality = "RNA")

    cluster_assignments_rna <-  cluster_assignments_rna %>% as.data.frame()  %>%  tibble::rownames_to_column("cell") %>%  as_tibble()
    colnames(cluster_assignments_rna)[2] <- "cluster"
    cluster_assignments_rna$cluster <-  paste0("C",cluster_assignments_rna$cluster + 1)
    cluster_assignments_rna$modality <-"RNA"

  }


  posterior_CNA <- inf$inferred_params$CNV_probabilities
  posterior_CNA <-  apply(posterior_CNA, 2, function(y) data.frame(y))
  names(posterior_CNA) <-  segment_names

  posterior_CNA <-  mapply(posterior_CNA,names(posterior_CNA), SIMPLIFY = F,FUN = function(y,z){
    colnames(y) <- paste(1:ncol(y))
    rownames(y) <-  cluster_names
    y %>% tibble::rownames_to_column("cluster") %>% tidyr::pivot_longer(!cluster, names_to = "value", values_to = "probability") %>%  mutate("segment_id"= z)
  })

  posterior_CNA <-  do.call(args = posterior_CNA, rbind)


  CNA <- inf$inferred_params$CNA
  if(any(dim(CNA) == NULL))
    CNA <-  as.data.frame(CNA) %>% t()
  colnames(CNA) <- segment_names
  rownames(CNA) <- cluster_names
  CNA <-  as.data.frame(CNA)  %>% tibble::rownames_to_column("cluster") %>% tidyr::pivot_longer(!cluster, names_to = "segment_id", values_to = "value")

  ret$CNA <-  CNA
  ret$posterior_CNA <-  posterior_CNA
  ret$mixing_proportions <- rbind(mixing_proportions_atac, mixing_proportions_rna)
  ret$cluster_assignments <- rbind(cluster_assignments_atac, cluster_assignments_rna)
  ret$z_nk <-  rbind(z_nk_atac, z_nk_rna)
  ret$segment_parameters <- rbind(segment_parameters_atac, segment_parameters_rna)


  return(ret)


}
