#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
#' 
#' x = Rcongas::congas_example
#' 
#' plot_normalization_factors(x) 
plot_normalization_factors <- function(x)
{
  best_model <- Rcongas:::get_best_model(x)
  
  if(is_gaussian(x)){
    
    df = data.frame(
      mu = best_model$parameters$norm_sd,
      segs = colnames(best_model$parameters$cnv_probs),
      stringsAsFactors = FALSE
    )
    
    ggplot(
      data = df,
      aes(x = mu)
    ) + geom_histogram(bins = 60, color = "black", fill = "royalblue2") +
      ggtitle("Distribution of segments standard deviations") +
      xlab("SD") +
      CNAqc:::my_ggplot_theme()
    
  } else{
    
    df = data.frame(
      mu = best_model$parameters$norm_factor,
      cell = names(best_model$parameters$norm_factor),
      stringsAsFactors = FALSE
    ) %>%
      as_tibble() %>%
      left_join(
        Rcongas:::get_clusters(x),
        by = 'cell'
      )
    
    ggplot(
      data = df,
      aes(x = mu, fill = cluster)
    ) + geom_histogram(bins = 60) +
      ggtitle("Distribution of library size factors") +
      xlab("Library size factor") +
      CNAqc:::my_ggplot_theme() +
      scale_fill_manual(values = Rcongas:::get_clusters_colors(df$cluster))
    
  }
  
  
}