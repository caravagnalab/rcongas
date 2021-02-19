#' Title
#'
#' @param x
#' @param all_genes 
#' @param as_matrix 
#' @param all_cells 
#'
#' @return
#' @export
#'
#' @examples
get_input_raw_data = function(x, 
                              all_cells = FALSE, 
                              all_genes = FALSE,
                              add_locations = FALSE,
                              add_clusters = FALSE,
                              as_matrix = FALSE)
{
  if (all(is.null(x$data$raw))) {
    
    cli::cli_alert_warning("Input counts are missing, check the init function to see how you created the data")
    
    return(NULL)
  }
  
  if(as_matrix & add_locations){
    cli::cli_alert_warning("Incompatible TRUE values for 'add_locations' and 'as_matrix', using only 'as_matrix'.")
    add_locations = FALSE
  }

  if(as_matrix & add_clusters){
    cli::cli_alert_warning("Incompatible TRUE values for 'add_clusters' and 'as_matrix', using only 'as_matrix'.")
    add_clusters = FALSE
  }
  
  # What we return
  y = x$data$raw
  
  # Cells (retain only those used to compute the mappings)
  if(!all_cells) 
  {
    y = y %>% filter(cell %in% (x$data$counts %>% rownames))
  }
  
  # Cells (retain only those used to compute the mappings)
  if(!all_genes) 
  {
    y = y %>% filter(gene %in% (x$data$gene_locations$gene))
  }
  
  # Matrix transform only if required
  if(as_matrix) 
  {
    # Spread, change NAs to 0
    y_matrix = y %>% spread(gene, n, fill = 0)
    y_mmatrix = y_matrix %>% dplyr::select(-cell) %>% as.matrix()
    rownames(y_mmatrix) = y_matrix$cell
  
    return(y_mmatrix)
  }
  
  # add_locations
  if(add_locations)
    y = y %>% left_join(x %>% get_mapped_genes(), by = 'gene')
  
  # Cells (retain only those used to compute the mappings)
  if(add_clusters) 
  {
    y = y %>% 
      left_join(
        get_clusters(x) %>% 
          select(cell, cluster),
        by = 'cell'
      )
  }
  
  return(y)

}
