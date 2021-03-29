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
  
  # Zero-counts cells - for sanitization
  zeroes_in_modality = function(modality)
  {
    zero_counts_cells = NULL
    
    M = x$input$dataset %>% 
      filter(modality == !!modality) %>% 
      select(segment_id, cell, value) %>% 
      pivot_wider(names_from = segment_id, values_from = value)
    
    # Cells stats - number of 0s
    cell_ids = which(is.na(M), arr.ind = TRUE)[, 1]
    
    if((cell_ids %>% length()) > 1)
    {
      cell_freq_z = M$cell[cell_ids] %>% table
      cell_freq_z = sort(cell_freq_z, decreasing = TRUE)
      
      # Stats
      ncells = length(cell_freq_z)
      nsegs = dim(M)[2] - 1
      
      # Proportions
      cell_freq_z_p = ((cell_freq_z/nsegs) * 100) %>% round
      cell_freq_z_p = paste0(cell_freq_z_p, '%')
      
      zero_counts_cells = cell_freq_z %>% 
        as_tibble() %>% 
        mutate(`%` = ((n/nsegs) * 100) %>% round, modality = !!modality)
      
      colnames(zero_counts_cells)[1] = 'cell'
    }
    
    return(zero_counts_cells)
  }
  
  zero_counts_cells_RNA = zeroes_in_modality("RNA")
  zero_counts_cells_ATAC = zeroes_in_modality("ATAC")
  
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
      atac_dtype = atac_dtype,
      zero_counts_cells_RNA = zero_counts_cells_RNA,
      zero_counts_cells_ATAC = zero_counts_cells_ATAC
    )
  )
}

stat_fit = function(x)
{
  if(!('best_fit' %in% names(x))) return(NULL)
  
  # Runs and tested k
  nruns = x$runs %>% length
  tested_k = x$runs %>% names
  
  # IC and score
  fit_IC = x$used_IC
  fit_score = x$best_fit$parameters$ICs[fit_IC] %>% as.numeric()
  
  # Hyperparameters 
  fit_model = x$best_fit$parameters$hyperparameters$model
  fit_loss = x$best_fit$parameters$hyperparameters$loss
  fit_optimizer = x$best_fit$parameters$hyperparameters$optimizer
  fit_lambda = x$best_fit$parameters$hyperparameters$lambda
  fit_k = x$best_fit$parameters$hyperparameters$K
  
  # Mixing
  fit_mixing = x$best_fit$mixing_proportions
  
  fit_mixing_RNA_n = x$best_fit$cluster_assignments %>% 
    filter(modality == 'RNA') %>% 
    pull(cluster) %>% 
    table %>% 
    as.numeric()
  names(fit_mixing_RNA_n) = x$best_fit$cluster_assignments %>% 
    filter(modality == 'RNA') %>% 
    pull(cluster) %>% 
    table %>% 
    names
  
  fit_mixing_ATAC_n = x$best_fit$cluster_assignments %>% 
    filter(modality == 'ATAC') %>% 
    pull(cluster) %>% 
    table %>% 
    as.numeric()
  
  names(fit_mixing_ATAC_n) = x$best_fit$cluster_assignments %>% 
    filter(modality == 'ATAC') %>% 
    pull(cluster) %>% 
    table %>% 
    names
  
  return(
    list(
      nruns = nruns,
      tested_k = tested_k,
      fit_IC = fit_IC,
      fit_score = fit_score,
      fit_model = fit_model,
      fit_loss = fit_loss,
      fit_optimizer = fit_optimizer,
      fit_lambda = fit_lambda,
      fit_k = fit_k,
      fit_mixing = fit_mixing,
      fit_mixing_RNA_n = fit_mixing_RNA_n,
      fit_mixing_ATAC_n = fit_mixing_ATAC_n
    )
  )
}
