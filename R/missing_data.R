#' Filter cells with missing data.
#' 
#' @description After mapping counts data to segments, this function
#' can be used to determine cells with missing data, and remove
#' them The function requires and returns an (R)CONGAS+ object.
#' 
#' This filter works by a proportion, as reported by the \code{\link{stat}} function.
#' 
#' If these cells are not removed, during inference missing values
#' are imputed to be \code{0}. This can create an excess of mixture components
#' fitting 0-counts data.
#'
#' @param x An \code{rcongasplus} object.
#' @param proportion_RNA The RNA proportion cut for a cell to be removed, default 5%.
#' @param proportion_ATAC The ATAC proportion cut for a cell to be removed, default 5%.
#'
#' @return The object \code{x} where 0-counts cells have been removed.
#' @export
#'
#' @examples
#' data('example_object')
#' 
#' # Default
#' print(example_object)
#' 
#' example_object = filter_missing_data(example_object)
#' 
#' # After filtering
#' print(example_object)
filter_missing_data = function(x, proportion_RNA = 0.05, proportion_ATAC = 0.05)
{
  stat_x = stat(x, what = 'data')
  
  x %>% sanitize_zeroes()
  
  all_cells = x$input$dataset$cell %>% unique
  n_orig = all_cells %>% length
  
  # RNA to remove
  if(!is.null(stat_x$zero_counts_cells_RNA))
  {
    to_remove = stat_x$zero_counts_cells_RNA %>% 
      filter(`%` > 100 * proportion_RNA) %>% 
      pull(cell)
    
    n = stat_x$zero_counts_cells_RNA %>% nrow
    
    cli::cli_h3("{.field {'RNA'}} modality, proportions cut {.field {proportion_RNA}}")
    cli::cli_alert("Cells with missing data {.field {n}}, removing {.field {to_remove %>% length}}")
    
    if(length(to_remove) > 0) all_cells = setdiff(all_cells, to_remove)
  }
  
  # ATAC to remove
  if(!is.null(stat_x$zero_counts_cells_ATAC))
  {
    to_remove = stat_x$zero_counts_cells_ATAC %>% 
      filter(`%` > 100 * proportion_ATAC) %>% 
      pull(cell)
    
    n = stat_x$zero_counts_cells_ATAC %>% nrow
    
    cli::cli_h3("{.field {'ATAC'}} modality, proportions cut {.field {proportion_ATAC}}")
    cli::cli_alert("Cells with missing data {.field {n}}, removing {.field {to_remove %>% length}}")
    
    if(length(to_remove) > 0) all_cells = setdiff(all_cells, to_remove)
  }
  
  # Filter
  if(n_orig > (all_cells %>% length))
  {
    x$input$dataset = x$input$dataset %>% filter(cell %in% all_cells)
    x$input$normalisation = x$input$normalisation %>% filter(cell %in% all_cells)
    
    x$log = paste0(x$log, '\n- ', 
                   Sys.time(), " Filtered missing data: [", proportion_RNA, 
                   '|', proportion_ATAC, ']')
    
  }
  

  
  
  return(x)
}