has_atac <-  function(x){
  return("ATAC" %in% toupper(get_data(x) %>%  pull(modality) %>%  unique()))
}

has_rna <-  function(x){
  return("RNA" %in% toupper(get_data(x) %>%  pull(modality) %>%  unique()))
}

which_likelihood <-  function(x, modality = "RNA"){
  return(get_data(x) %>% filter(modality == !!modality) %>%  pull(value_type) %>%  unique())
}
