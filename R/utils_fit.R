tensorize <-  function(x){

  torch <- reticulate::import("torch")
  ret <-  torch$tensor(x)
  return(ret)
}


input_data_from_rcongas <- function(x){
  ret <-  list()
  if(has_atac(x)){
    ret$data_atac <- get_data(x) %>%  filter(modality == "ATAC") %>%  reshape2::acast(segment_id  ~ cell, value.var = "value")
    ret$data_atac <- ret$data_atac[order(rownames(ret$data_atac)),]
    norm_tmp <- get_normalisation(x) %>%  filter(modality == "ATAC")
    norm_factor_atac <- norm_tmp$normalisation_factor
    names(norm_factor_atac) <- norm_tmp$cell
    ret$norm_factor_atac <- norm_factor_atac[colnames(ret$data_atac)]
  }

  if(has_rna(x)){
    ret$data_rna <- get_data(x) %>%  filter(modality == "RNA") %>%  reshape2::acast(segment_id  ~ cell, value.var = "value")
    ret$data_rna <- ret$data_rna[order(rownames(ret$data_rna)),]
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



format_best_model <-  function(x){
  ret <-  list()

  ret$parameters <-  x$ICs

  #ret$cluster_assignements <-

  #ret$segment_factors <-


}
