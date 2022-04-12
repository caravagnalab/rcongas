

#' Return expression level per segment
#'
#' @param inf_obj an Rcongas object
#' @param chromosomes Chromosome to filter 
#' @param normalise should the counts be divided by the library size?
#' @param z_score returns the Z score for each segment instead of the actual expression value
#' @param should the the denominator of the likelihood be included? (usually TRUE unless you are using very old versions)
#'
#' @return A tibble with expression values for each segment and cluster
#' @export
#'
#' @examples
#' 
#' library(Rcongas)
#' 
#' data('congas_example', package = 'Rcongas')
#' 
#' get_counts(congas_example)
#' 
#' 
get_counts <-
  function(x,
           chromosomes = paste0("chr", c(1:22, "X", "Y")),
           normalise = TRUE,
           z_score = FALSE,
           sum_denominator = TRUE)
  {
    data_matrix = x$data$counts
    
    if (normalise) {
      
      if (!has_inference(x)) stop("To use normalise = TRUE you need to compute a model fit first.")
        
      best_model <- Rcongas:::get_best_model(x)
      
      if (is.null(best_model$parameters$norm_factor))
        best_model$parameters$norm_factor <-
          rep(1, length(best_model$parameters$assignement))
      
      normalisation_factors = best_model$parameters$norm_factor
      
      assignments = get_cluster_assignments(x) %>% gsub(pattern = "C",
                                                        replacement = "",
                                                        ignore.case = T)  %>% as.numeric
      mu = get_input_segmentation(x)$mu
      total_cn =  rowSums(best_model$parameters$cnv_probs * mu)
      
      
      # Handle this special case which happens for already normalised data
      if (is.null(normalisation_factors)) {
        warning(
          "normalisation_factors are NULL - replacing them with all 1s. Maybe you're using normalised data."
        )
        
        normalisation_factors = rep(1, ncol(data_matrix))
      }
      
      for (i in 1:ncol(data_matrix))
        if (sum_denominator)
          data_matrix[, i] = (data_matrix[, i] / normalisation_factors)  * (total_cn[assignments] / sum(mu))
      else
        data_matrix[, i] = (data_matrix[, i] / normalisation_factors)
    }
    
    if (z_score) {
      rnames <-  rownames(data_matrix)
      data_matrix <-
        apply(data_matrix, 2, scale, center = T, scale = T)
      rownames(data_matrix) <-  rnames
    }
    
    M <-
      long_counts(data_matrix) %>%
      select(chr, from, to, cell, n) %>%
      dplyr::arrange(chr, from, to)
    
    if (has_inference(x))
    {
      best_model <- Rcongas:::get_best_model(x)
      clts <- as.data.frame(best_model$parameters$assignement)
      colnames(clts) <- "cluster"
      clts$cell <- rownames(clts)
      
      M <-  dplyr::inner_join(M, clts, by = "cell") %>%
        dplyr::select(chr, from, to, cell, n, cluster) %>%
        dplyr::mutate(cluster = paste(cluster))
    }
    else
      M$cluster = 'Not Available'
    
    
    return(M %>% filter(chr %in% chromosomes))
    
  }
