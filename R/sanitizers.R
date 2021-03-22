sanitize_obj = function(x)
{
  # Compare colnames vs expected
  chk_cn = function(x, cn, where)
  {
    m = setdiff(cn, x %>% colnames)

    if(length(m) > 1)
    {

      message("Error with this tibble: ", where)
      print(x)

      stop(
        paste("Missing columns", paste(m, collapse = ', '))
      )
    }
  }

  # Check for NA things
  chk_na = function(x, cols)
  {
    what_incomplete = !complete.cases(x[, cols])
    if(what_incomplete %>% any)
    {
      message("Error with this tibble: ")
      print(x[what_incomplete, ])

      stop(
        paste("NA values should not be here.")
      )
    }
  }

  # General input
  stopifnot(inherits(x, 'rcongasplus'))

  # Input data
  if(!('input' %in% names(x))) stop("Missing input")
  if(!('dataset' %in% names(x$input))) stop("Missing input$dataset")
  if(!('normalisation' %in% names(x$input))) stop("Missing input$normalisation")
  if(!('segmentation' %in% names(x$input))) stop("Missing input$segmentation")

  # Check that all the tibbles have the required information
  chk_cn(x$input$dataset, c("cell", "segment_id", "value", "modality"), "$input$dataset")
  chk_cn(x$input$normalisation, c("cell", "normalisation_factor", "modality"), "$input$normalisation")
  chk_cn(x$input$segmentation, c("chr", "from", "to", "copies", "segment_id"), "$input$segmentation")

  # Some things CANNOT be NA
  chk_na(x$input$dataset, c("cell", "segment_id", "value", "modality"))
  chk_na(x$input$normalisation, c("cell", "normalisation_factor", "modality"))
  chk_na(x$input$segmentation, c("chr", "from", "to", "copies", "segment_id"))

  # If available, fit information
  if('runs' %in% names(x))
  {
    if(!('best_fit' %in% names(x))) stop("Missing best_fit")

    if(!('parameters' %in% names(x$best_fit))) stop("Missing best_fit$parameters")
    if(!('cluster_assignments' %in% names(x$best_fit))) stop("Missing best_fit$cluster_assigments")
    if(!('segment_parameters' %in% names(x$best_fit))) stop("Missing best_fit$segment_parameters")
    if(!('CNA' %in% names(x$best_fit))) stop("Missing best_fit$CNA")
    if(!('posterior_CNA' %in% names(x$best_fit))) stop("Missing best_fit$posterior_CNA")
    if(!('mixing_proportions' %in% names(x$best_fit))) stop("Missing best_fit$mixing_proportions")
    if(!('cluster_assignments' %in% names(x$best_fit))) stop("Missing best_fit$cluster_assignments")
    if(!('z_nk' %in% names(x$best_fit))) stop("Missing best_fit$z_nk")

    if(!('model_selection' %in% names(x))) stop("Missing model_selection")

    # Check that all the tibbles have the required information
    chk_cn(x$best_fit$segment_parameters, c("segment_id", "parameter", "value", "modality"))
    chk_cn(x$best_fit$CNA, c("cluster", "segment_id", "value"))
    chk_cn(x$best_fit$posterior_CNA, c("cluster", "segment_id", "probability", "value"))
    chk_cn(x$best_fit$mixing_proportions, c("cluster", "mixing", "modality"))
    chk_cn(x$best_fit$cluster_assignments, c("cluster", "cell", "modality"))
    chk_cn(x$best_fit$z_nk, c("cluster", "cell", "z_nk", "modality"))
  }

  return(x)
}

sanitize_input = function(x, required_input_columns, types_required)
{
  if(is.null(x)) return()

  # No NAs
  if(any(is.na(x)))
  {
      message("Error with this tibble")
      print(x)

      stop("NAs in the input")
  }

  # Missing columns
  m = setdiff(required_input_columns, x %>% colnames)

  if(length(m) > 1)
  {
    message("Error with this tibble")
    print(x)

    stop(paste("Missing columns", paste(m, collapse = ', ')))
  }

  # Types
  classes = sapply(x[required_input_columns], class)

  if(!all(classes == types_required))
  {
    message("Error with this tibble")
    print(x)

    which_mismatch = classes != types_required

    stop(paste(
      "Types mismatched expected",
      paste(types_required[which_mismatch], collapse = ', '),
      "for",
      paste(required_input_columns[which_mismatch], collapse = ', ')
    ))
  }

  invisible(1)
}





