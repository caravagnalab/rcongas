#' Return dataset statistics
#' 
#' @description From a CONGAS object returns the data statistics such ash
#' numnber of cells, semgents, clusters (if any), score used etc.
#'
#' @param x A CONGAS object
#'
#' @return A list of named values that contain the required statistics.
#' 
#' @export
#'
#' @examples
#' 
#' x = Rcongas::congas_example
#'
#' get_dataset_stats(x)
get_dataset_stats = function(x)
{
  # number of cells, segments
  ncells = x$data$counts %>% nrow()
  nsegments = x$data$counts %>% ncol()
  ngenes = x$data$gene_locations$gene %>% unique() %>% length()
  
  # clustering stats
  nclusters = sizes_cl = psizes_cl = NULL
  
  if(Rcongas:::has_inference(x))
  {
    nclusters =  Rcongas::get_k(x)
    sizes_cl = Rcongas::get_clusters_size(x)
    psizes_cl = Rcongas::get_clusters_size(x, normalised = TRUE)
  }
  
  # scores
  x$inference$model_selection
  
  
  return(
    list(
      ncells = ncells,
      ngenes = ngenes,
      nsegments = nsegments,
      clusters_k = nclusters,
      clusters_n = sizes_cl,
      clusters_pi = psizes_cl,
      score_type = x$inference$model_selection$IC_type,
      score = x$inference$model_selection$IC[x$inference$model_selection$best_K]
    )
  )
}
