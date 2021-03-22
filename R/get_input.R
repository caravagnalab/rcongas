#' Extract input data.
#'
#' @description General data accessing getter function to return
#' any of:
#' 
#' * data, 
#' * segmentation, 
#' * normalisation factors.
#' 
#' The function uses the \code{what} parameter to return the 
#' appropriate type of information
#' 
#' Since the input formats are different, the outputs are also
#' different based on the value of \code{what}, but they are 
#' always in tibble format. 
#' 
#' Besides obvious columns, these information will also be 
#' available in the returned tibbles.
#' 
#' * \code{what = "data"}: \code{modality} and \code{value_type}, the latter reporting the likelihood
#' code associated to the modality.
#' 
#' * \code{what = "segmentation"}: 
#'     * \code{ATAC_nonzerovals} and \code{ATAC_peaks}, reporting the number of ATAC entries mapped to a 
#'     segment, and the number of peaks these come from. Note that one non-zero entry is given by a cell that reported a 
#'     ATAC count in a peak mapping within the segment. 
#' 
#'     * \code{RNA_nonzerovals} and \code{RNA_genes}, reporting the number of RNA entries mapped to a segment, 
#'     and the number of genes these come from. For RNA, one non-zero entry is given by a cell that reported an RNA count
#'      in a gene mapping within the segment. 
#'      
#' * \code{what = "normalisation"}: just \code{modality}.
#'
#' @param x An object of class \code{rcongasplus}.
#' @param what Any of \code{"data"}, \code{"segmentation"} or
#' \code{"normalisation"}.
#'
#' @return A tibble; its format depends on \code{what}. See the examples.
#' @export
#'
#' @examples
#' data(example_object)
#' 
#' # Extract the input data, after mapping to segments
#' get_input(example_object, what = 'data') 
#' 
#' # e.g., in this way you can get input ATAC values
#' get_input(example_object, what = 'data') %>% dplyr::filter(modality == 'ATAC')
#'  
#' # Extract the input segmentation. 
#' get_input(example_object, what = 'segmentation')
#'  
#' # Extract the input normalisation factors
#' get_input(example_object, what = 'normalisation')
get_input = function(x, what = 'data')
{
  x %>% sanitize_obj()

  if(what == 'data') return(x %>% get_data())
  if(what == 'segmentation') return(x %>% get_segmentation())
  if(what == 'normalisation') return(x %>% get_normalisation())
  # if(what == 'modalities') return(x %>% get_modalities())

  stop("Unrecognised 'what': use any of 'data', 'segmentation', 'normalisation' or 'modalities'.")
}

get_data = function(x)
{
  x$input$dataset
}

get_segmentation = function(x)
{
  x$input$segmentation
}

get_normalisation = function(x)
{
  x$input$normalisation
}

# The stat(x) function returns this information. These functions
# should return only tibbles.
get_modalities <-  function(x){

  x$input$dataset$modality %>%  unique()
}
