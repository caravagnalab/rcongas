#' Get input data.
#'
#' @description General data accessing getter function to return
#' input data, segmentation or normalisation.
#'
#' @param x
#' @param what Any of \code{"data"}, \code{"segmentation"} or
#' \code{"normalisation"}.
#'
#' @return A tibble.
#' @export
#'
#' @examples
get_input = function(x, what = 'data')
{
  x %>% sanitize_obj()

  if(what == 'data') return(x %>% get_data())
  if(what == 'segmentation') return(x %>% get_segmentation())
  if(what == 'normalisation') return(x %>% get_normalisation())
  if(what == 'modalities') return(x %>% get_modalities())

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

get_modalities <-  function(x){

  x$input$dataset$modality %>%  unique()
}
