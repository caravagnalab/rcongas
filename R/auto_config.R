auto_config_run <- function(x, K, NB_size_atac = 150, NB_size_rna = 150, lambda = 0.3, a_sd = 1, b_sd = 100, prior_cn = c(0.2,0.6,0.2,0.05,0.05)){

    param_list <-  list()

    torch <-  reticulate::import("torch")

    param_list$probs <-  torch$tensor(prior_cn)

    if(has_atac(x)){
      params_atac <- list()

      if(startsWith(which_likelihood(x, "ATAC"), prefix = "G")){
        sds <- get_data(x) %>%  filter(modality == "ATAC") %>% group_by(segment_id) %>%  summarize(sd = sd(value)) %>%  pull(sd)
        params_atac$norm_init_sd_atac <- torch$tensor(sds / max(K))
        params_atac$likelihood_atac <-  "Gaussian"
      } else if(startsWith(which_likelihood(x, "ATAC"), prefix = "NB")){
        params_atac <- gamma_shape_rate(x, modality = "ATAC")
        names(params_atac) <-  c("theta_shape_atac", "theta_rate_atac")
        params_atac$nb_size_init_atac <- torch$tensor(rep(NB_size_atac, length(get_segmentation(x) %>%  pull(segment_id) %>%  unique())))
        params_atac$likelihood_atac <-  "NB"
      } else {
        params_atac <- gamma_shape_rate(x, modality = "ATAC", torch = torch)
        names(params_atac) <-  c("theta_shape_atac", "theta_rate_atac")
        params_atac$likelihood_atac <-  "Poisson"
      }

      param_list <-  c(param_list, params_atac)

    }

  if(has_rna(x)){

      params_rna <- list()

      if(startsWith(which_likelihood(x, "RNA"), prefix = "G")){
        sds <- get_data(x) %>%  filter(modality == "RNA") %>% group_by(segment_id) %>%  summarize(sd = sd(value)) %>%  pull(sd)
        params_rna$norm_init_sd_rna <- torch$tensor(sds / max(K))
        params_rna$likelihood_rna <-  "Gaussian"
      } else if (startsWith(which_likelihood(x, "RNA"), prefix = "NB")){
        params_rna <- gamma_shape_rate(x, modality = "RNA")
        names(params_rna) <-  c("theta_shape_rna", "theta_rate_rna")
        params_rna$nb_size_init_rna <- torch$tensor(rep(NB_size_rna, get_segmentation(x) %>%  pull(segment_id) %>%  unique()))
        params_rna$likelihood_rna <-  "NB"
      } else {
        params_rna <- gamma_shape_rate(x, modality = "RNA", torch = torch)
        names(params_rna) <-  c("theta_shape_rna", "theta_rate_rna")
        params_rna$likelihood_rna <-  "Poisson"
      }

      param_list <-  c(param_list, params_rna)

  }

  param_list$lambda <- lambda

  return(param_list)

}

gamma_shape_rate <-  function(x, modality = "RNA", torch = reticulate::import("torch")){


  inp = reshape2::acast(get_data(x) %>% filter(modality == !!modality), cell~segment_id, value.var="value")
  inp = inp[order(rownames(inp)),]
  norm_raw = get_normalisation(x) %>% filter(modality == !!modality) %>% select(normalisation_factor, cell)
  norm = norm_raw$normalisation_factor
  names(norm) = norm_raw$cell
  norm = norm[order(rownames(inp))]
  ploidy <- x$input$segmentation$copies

  theta_factors = estimate_segment_factors(inp,norm, ploidy,plot=F)

  theta_shape = torch$tensor(sapply(theta_factors, function(x) x[1]))
  theta_rate = torch$tensor(sapply(theta_factors, function(x) x[2]))

  ret = list(theta_shape, theta_rate)

  return(ret)
}


