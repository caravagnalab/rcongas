# Key creation and decryption
idify = function(y) {
  y %>% dplyr::mutate(segment_id = paste(chr, as.integer(from), as.integer(to), sep = ":"))
}

deidify = function(y) {
  y %>% tidyr::separate(segment_id,
                        into = c('chr', 'from', 'to'),
                        sep = ":") %>%
    dplyr::mutate(from = as.integer(from), to = as.integer(to))
}

# Named vector of colours
modality_colors = function(what)
{
  # Colors come from this
  # v = c(ggsci::pal_jama()(2), 'indianred3')
  
  # RNA         ATAC         <NA> 
  v = c("#374E55FF", "#DF8F44FF", "indianred3") 
  names(v) = c("RNA", "ATAC", "segmentation")
  
  w_rna = grep("RNA", what)
  w_atac = grep("ATAC", what)
  w_segm = grep("segmentation", what)

  if(length(w_rna) > 0) names(v)[1] = what[w_rna]
  if(length(w_atac) > 0) names(v)[2] = what[w_atac]
  if(length(w_segm) > 0) names(v)[3] = what[w_segm]
  
  return(v)
}

# Normalise modality data by the corresponding normalisation factor
normalise_modality = function(modality_data, normalisation_factors)
{
  what_modality_data = modality_data$modality %>% unique()
  what_modality_norm = normalisation_factors$modality %>% unique()
  
  if(
    (what_modality_data %>% length() > 1) |
    (what_modality_norm %>% length() > 1) 
  )
    stop("Data contains multiple modalities: ", 
         paste(what_modality_data, collapse = ', '),
         '. Select only one at a time by subsetting the data.')
  
  if(what_modality_data != what_modality_norm)
    stop("Data and factors do contain the same modalities: ", 
         paste(what_modality_data, 'vs', what_modality_norm),
         '. Select only one at a time by subsetting the data.')
  
  cli::cli_alert("Normalising {.field {what_modality_norm}} counts using input normalisation factors.")
  
  factors = normalisation_factors %>% dplyr::select(-modality)
  
  # Divide counts by per-cell normalization_factor
  modality_data = modality_data %>% 
    left_join(factors, by = 'cell') %>% 
    mutate(value = value/normalisation_factor)
  
  return(modality_data %>% select(-normalisation_factor))
}

# Normalise modality data by the corresponding normalisation factor
# rescale_segment = function(modality_data, normalisation_factors)
# {
#   what_modality_data = modality_data$modality %>% unique()
#   what_modality_norm = normalisation_factors$modality %>% unique()
#   
#   if(
#     (what_modality_data %>% length() > 1) |
#     (what_modality_norm %>% length() > 1) 
#   )
#     stop("Data contains multiple modalities: ", 
#          paste(what_modality_data, collapse = ', '),
#          '. Select only one at a time by subsetting the data.')
#   
#   if(what_modality_data != what_modality_norm)
#     stop("Data and factors do contain the same modalities: ", 
#          paste(what_modality_data, 'vs', what_modality_norm),
#          '. Select only one at a time by subsetting the data.')
#   
#   cli::cli_alert("Normalising {.field {what_modality_norm}} counts using input normalisation factors.")
#   
#   factors = normalisation_factors %>% dplyr::select(-modality)
#   
#   # Divide counts by per-cell normalization_factor
#   modality_data = modality_data %>% 
#     left_join(factors, by = 'cell') %>% 
#     mutate(value = value/normalisation_factor)
#   
#   return(modality_data %>% select(-normalisation_factor))
# }
