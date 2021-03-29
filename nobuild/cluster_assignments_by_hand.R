calculate_cluster_assignements_r <- function(x,run){

  data_atac <- get_data(x) %>%  filter(modality == "ATAC") %>%  reshape2::acast(segment_id  ~ cell, value.var = "value")
  data_atac[is.na(data_atac)] <-  0
  data_atac <- data_atac[order(rownames(data_atac)),order(colnames(data_atac))]
  norm_tmp <- get_normalisation(x) %>%  filter(modality == "ATAC")
  norm_factor_atac <- norm_tmp$normalisation_factor
  names(norm_factor_atac) <- norm_tmp$cell
  norm_factor_atac <- norm_factor_atac[colnames(data_atac)]

  data_rna <- get_data(x) %>%  filter(modality == "RNA") %>%  reshape2::acast(segment_id  ~ cell, value.var = "value")
  data_rna[is.na(data_rna)] <-  0
  data_rna <- data_rna[order(rownames(data_rna)),order(colnames(data_rna))]
  norm_tmp <- get_normalisation(x) %>%  filter(modality == "RNA")
  norm_factor_rna <- norm_tmp$normalisation_factor
  names(norm_factor_rna) <- norm_tmp$cell
  norm_factor_rna <- norm_factor_rna[colnames(data_rna)]

  CNV <-  run$inferred_params$CNA
  rna_segment_factors <- run$inferred_params$segment_factor_rna
  atac_segment_factors <- run$inferred_params$segment_factor_atac
  norm_factors_rna <- x$input$normalisation %>% filter(modality=="RNA") %>% pull(normalisation_factor)
  norm_factors_atac <- x$input$normalisation %>% filter(modality=="ATAC") %>% pull(normalisation_factor)
  mixture_weights_rna <- run$inferred_params$mixture_weights_rna
  mixture_weights_atac <- run$inferred_params$mixture_weights_atac

  NB_size_rna <- run$inferred_params$NB_size_rna
  NB_size_atac <- run$inferred_params$NB_size_atac

  ass_rna <- array(dim = c(length(mixture_weights_rna),length(rna_segment_factors) , length(norm_factors_rna)))

  ass_atac <- array(dim = c(length(mixture_weights_atac),length(atac_segment_factors) , length(norm_factors_atac)))

  for(k in 1:length(mixture_weights_atac)){
    for(i in 1:length(atac_segment_factors)){
      for(n in 1:length(norm_factors_atac)){

        mu <- atac_segment_factors[i] *  CNV[k,i] * norm_factors_atac[n]

        ass_atac[k,i,n] <- dnbinom(data_atac[i,n], mu = mu, size = 10000, log = T)

      }
    }
  }

  ass_atac <- ass_atac[,2,]
  ass_atac_summed <- apply(ass_atac,c(1,3), function(y) sum(y))

  ass_atac_summed <- ass_atac

  ass_atac_summed <-  apply(ass_atac_summed,2, function(y) y + log(mixture_weights_atac))

  ass_atac_summed_norm <- log_sum_exp(ass_atac_summed)

  ass_atac_norm <- exp(ass_atac_summed - ass_atac_summed_norm)

  ass_atac <- paste0("C", apply(ass_atac_norm,2,which.max)) %>%  as.data.frame()


  for(k in 1:length(mixture_weights_rna)){
    for(i in 1:length(rna_segment_factors)){
      for(n in 1:length(norm_factors_rna)){

        mu <- rna_segment_factors[i] *  CNV[k,i] * norm_factors_rna[n]

        ass_rna[k,i,n] <- dnbinom(data_rna[i,n], mu = mu, size = NB_size_rna[i], log = T)

      }
    }
  }

  ass_rna_summed <- apply(ass_rna,c(1,3), function(y) sum(y))

  ass_rna_summed <-  apply(ass_rna_summed,2, function(y) y + log(mixture_weights_rna))

  ass_rna_summed_norm <- log_sum_exp(ass_rna_summed)

  ass_rna_norm <- exp(ass_rna_summed - ass_rna_summed_norm)

  ass_rna <- paste0("C", apply(ass_rna_norm,2,which.max)) %>%  as.data.frame()



  ass_rna$modality <-  "RNA"
  ass_rna$cell <-  names(norm_factor_rna)
  colnames(ass_rna)[1] <-  "cluster"
  ass_atac$modality <-  "ATAC"
  ass_atac$cell <-  names(norm_factor_atac)
  colnames(ass_atac)[1] <-  "cluster"

  return(rbind(ass_rna, ass_atac))


}



log_sum_exp <-  function(val){

  res <- vector(length = length(val[1,]))

  for(n in 1:length(res)){
    c <-  max(val[,n])
    xx <- val[,n] - c
    res[n] <- c + log(sum(exp(xx)))
  }
  return(res)

}
