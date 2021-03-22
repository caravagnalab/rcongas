#' Compute the statistics of a dataset.
#' 
#' @description The statistics of a dataset regard the input data, and the fits.
#' These are accessed by using different values for the \code{what} parameter
#' passed to this function (by default, data are used).
#' 
#' The result value is a named list with several information that are usefull
#' to summarise the model. The names should be self explicatory.
#' 
#' For data  the following information are reported:
#' 
#'  * it is reported the number of cells in each modality,
#'  * the number of non-zero ATAC/RNA input entries used to create the object
#'  * the number of RNA genes and ATAC peaks used
#'  * the number of segments
#'  * the number of modalities
#'  * the keyword for each  modality (\code{"RNA"} or \code{"ATAC"})
#'  * the likelihood for each  modality (\code{"G"} for Gaussian or \code{"NB"} for Negative Binomial)
#' 
#' For fits the following information are reported:
#' 
#'  * ...
#'  
#' If fits are not available, \code{NULL} is returned. 
#'
#' @md
#' 
#' @param x An object of class \code{rcongasplus}.
#' @param what The type of statistics that needs to be computed. It can be one of
#' \code{"data"} or \code{"fit"}. By deafult it is \code{"data"}.
#'
#' @return A named list with the available statistics. 
#' @export
#'
#' @examples
#' data(example_object)
#' 
#' # Compare the information here 
#' print(example_object)
#' 
#' # ... to the one reported here for data
#' stat(example_object, what = 'data')
#' 
#' # ... and here for fits
#' stat(example_object, what = 'fit')
stat = function(x, what = 'data')
{
  x %>% sanitize_obj()
  
  if(what == 'data') return(x %>% stat_data())
  if(what == 'fit') return(x %>% stat_fit())

  stop("Unrecognised 'what': use any of 'data' or 'clusters'.")
}
  
stat_data = function(x)
{
  # number of cells, segments, genes, peaks and modalities
  ncells_RNA = x %>% get_data() %>% filter(modality == 'RNA') %>% distinct(cell) %>% nrow()
  ncells_ATAC = x %>% get_data() %>% filter(modality == 'ATAC') %>% distinct(cell) %>% nrow()
  
  nsegments = x %>% get_segmentation() %>% nrow()
  nmodalities = x %>% get_data() %>% pull(modality) %>% unique() %>% length()
  
  rna_events = x %>% get_segmentation() %>% pull(RNA_nonzerovals) %>% sum()
  rna_genes = x %>% get_segmentation() %>% pull(RNA_genes) %>% sum()
  
  atac_events = x %>% get_segmentation() %>% pull(ATAC_nonzerovals) %>% sum()
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

stat_fit = function(x)
{
  if(!('best_fit' %in% names(x))) return(NULL)
  
  return("BINGOBONGO")
}
