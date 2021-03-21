auto_config_run <- function(x, K, NB_size_atac = 150, NB_size_rna = 150, lambda = 0.3, a_sd = 1, b_sd = 100, prior_cn = c(0.2,0.6,0.2,0.05,0.05)){

    param_list <-  list()

    torch <-  reticulate::import("torch")

    param_list$probs <-  torch$tensor(prior_cn)

    if(has_atac(x)){

      params_atac <- gamma_shape_rate(x, modality = "ATAC")
      names(params_atac) <-  c("theta_shape_atac", "theta_rate_atac")
      if(startsWith(which_likelihood(x, "ATAC"), prefix = "[G|N]")){
        sds <- get_data(x) %>%  filter(modality == "ATAC") %>% group_by(segment_id) %>%  summarize(sd = sd(value)) %>%  pull(sd)
        params_atac$norm_init_sd_atac <- sds / max(K)
        params_atac$likelihood_atac <-  "Gaussian"
      } else if(startsWith(which_likelihood(x, "ATAC"), prefix = "NB")){
        params_atac$nb_size_init_atac <- rep(NB_size_atac, length(K))
        params_atac$likelihood_atac <-  "NB"
      } else {
        params_atac$likelihood_atac <-  "Poisson"
      }

      param_list <-  c(param_list, params_atac)

    }

  if(has_rna(x)){

      params_rna <- gamma_shape_rate(x, modality = "RNA")
      names(params_rna) <-  c("theta_shape_rna", "theta_rate_rna")
      if(startsWith(which_likelihood(x, "RNA"), prefix = "[G|N]")){
        sds <- get_data(x) %>%  filter(modality == "RNA") %>% group_by(segment_id) %>%  summarize(sd = sd(value)) %>%  pull(sd)
        params_atac$norm_init_sd_rna <- sds / max(K)
        params_rna$likelihood_rna <-  "Gaussian"
      } else if (startsWith(which_likelihood(x, "RNA"), prefix = "NB")){
        params_rna$nb_size_init_rna <- rep(NB_size_rna, length(K))
        params_rna$likelihood_rna <-  "NB"
      } else {
        params_rna$likelihood_rna <-  "Poisson"
      }

      param_list <-  c(param_list, params_rna)

  }

  param_list$lambda <- lambda

}

gamma_shape_rate <-  function(x, modality = "RNA"){

  inp = reshape2::acast(get_data(x) %>% filter(modality == modality), cell~segment_id, value.var="value")
  inp = inp[order(rownames(inp)),]
  norm_raw = get_normalisation(x) %>% filter(modality == modality) %>% select(normalisation_factor, cell)
  norm = norm_raw$normalisation_factor
  names(norm) = norm_raw$cell
  norm = norm[order(rownames(inp))]
  ploidy <- x$input$segmentation$copies

  theta_factors = estimate_segment_factors(inp,norm, ploidy,plot=F)

  theta_shape = sapply(theta_factors, function(x) x[1])
  theta_rate = sapply(theta_factors, function(x) x[2])

  ret = list(theta_shape, theta_rate)

  return(ret)
}


