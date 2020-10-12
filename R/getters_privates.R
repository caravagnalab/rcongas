
# PRIVATE GETTERS

get_karyotype <-  function(x) {
  if (x$reference_genome %in% c('hg19', 'GRCh37'))
    data('hg19_karyo')
  if (x$reference_genome %in% c('hg38', 'GRCh38'))
    data('hg38_karyo')
}

get_best_model <- function(X) {
  X$inference$models[[X$inference$model_selection$best_K]]
}





# Key creation and decrypt
idify = function(y) {
  y %>% dplyr::mutate(segment_id = paste(chr, as.integer(from), as.integer(to), sep = ":"))
}

deidify = function(y) {
  y %>% tidyr::separate(segment_id,
                        into = c('chr', 'from', 'to'),
                        sep = ":") %>%
    dplyr::mutate(from = as.integer(from), to = as.integer(to))
}
