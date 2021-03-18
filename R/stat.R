#' Title
#'
#' @param x 
#' @param what 
#'
#' @return
#' @export
#'
#' @examples
stat = function(x, what = 'data')
{
  x %>% sanitize_obj()
  
  if(what == 'data') return(x %>% stat_data())
  if(what == 'clusters') return(x %>% stat_clusters())

  stop("Unrecognised 'what': use any of 'data' or 'clusters'.")
}
  
stat_data = function(x)
{
  # number of cells, segments, genes, peaks and modalities
  ncells_RNA = x %>% get_data() %>% filter(modality == 'RNA') %>% distinct(cell) %>% nrow()
  ncells_ATAC = x %>% get_data() %>% filter(modality == 'ATAC') %>% distinct(cell) %>% nrow()
  
  nsegments = x %>% get_segmentation() %>% nrow()
  nmodalities = x %>% get_data() %>% pull(modality) %>% unique() %>% length()
  
  rna_events = x %>% get_segmentation() %>% pull(RNA_events) %>% sum()
  rna_genes = x %>% get_segmentation() %>% pull(RNA_genes) %>% sum()
  
  atac_events = x %>% get_segmentation() %>% pull(ATAC_events) %>% sum()
  atac_peaks = x %>% get_segmentation() %>% pull(ATAC_peaks) %>% sum()
  
  # data types and likelihood adopted
  rna_dtype = x %>% get_data() %>% filter(modality == 'RNA') %>% pull(value_type) %>% unique()
  atac_dtype = x %>% get_data() %>% filter(modality == 'ATAC') %>% pull(value_type) %>% unique()
  
  # Modalities
  modalities = x %>% get_data() %>% pull(modality) %>% unique()
  
  return(
    list(
      ncells_RNA = ncells_RNA,
      ncells_ATAC = ncells_ATAC,
      rna_events = rna_events,
      rna_genes = rna_genes,
      atac_events = atac_events,
      atac_peaks = atac_peaks,
      nsegments = nsegments,
      nmodalities = nmodalities,
      modalities = modalities,
      rna_dtype = rna_dtype,
      atac_dtype = atac_dtype
    )
  )
}

stat_clusters = function(x)
{
 return(NULL)
}
