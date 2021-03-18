#' Access data.
#' 
#' @description General data accessing getter function to return
#' 
#'
#' @param x 
#' @param what 
#'
#' @return
#' @export
#'
#' @examples
get_input = function(x, what = 'data')
{
  x %>% sanitize_obj()
  
  if(what == 'data') return(x %>% get_data())
  if(what == 'segmentation') return(x %>% get_segmentation())
  if(what == 'normalisation') return(x %>% get_normalisation())

  stop("Unrecognised 'what': use any of 'data', 'segmentation' or 'normalisation'.")
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