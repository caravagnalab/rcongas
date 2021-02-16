#' Title
#'
#' @param x
#' @param transpose 
#' @param all_cells 
#' @param as_tibble 
#'
#' @return
#' @export
#'
#' @examples
get_input_raw_data = function(x, 
                              transpose = FALSE, 
                              only_clustered_cells = TRUE, 
                              only_mapped_genes = TRUE,
                              as_tibble = FALSE)
{
  if (all(is.null(x$data$gene_counts))) {
    cli::cli_alert_warning("Input data has not been stored in the object, re-run the analysis with XXX = TRUE ...")
    return(NULL)
  }
  
  y = x$data$gene_counts
  
  if(transpose) y = t(x$data$gene_counts)
  
  if(only_clustered_cells && has_inference(x))
     y = y[, colnames(y) %in% (get_clusters(x) %>% pull(cell)), drop = FALSE]
  
  # Retain only genes
  if(only_mapped_genes && has_inference(x))
    y = y[rownames(y) %in% (get_mapped_genes(x) %>% pull(gene)), , drop = FALSE]
  
  if(!as_tibble) return(y)
  
  y_tb = y %>% 
    t %>% 
    as_tibble 
  # %>% 
  #   mutate_all(~ replace(., . == 0, NA))
  
  y_tb$cell = colnames(y)
  y_tb_melt = reshape2::melt(y_tb, id = 'cell') %>% dplyr::as_tibble()
    
  y_tb_melt = y_tb_melt %>% 
    filter(value > 0) %>% 
    rename(gene = variable, count = value) 
  
  return(y_tb_melt)
}
