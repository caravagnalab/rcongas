#' Per segment cut of high and lower quantiles 
#'
#'The function simply calculates two quantiles for each segments and the caps all the values higher than `high_quant` and lower than 
#' `low_quant`
#'
#' @param x Rcongas object
#' @param low_quant lower quantile 
#' @param high_quant high quantile
#'
#' @return
#' @export
#'
#' @examples
#' 
#' x <-  Rcongas::congas_example
#' 
#' get_counts(x) %>% pull(n) %>% max
#' 
#' x <- filter_quantiles_per_segments(x)
#' 
#' get_counts(x) %>% pull(n) %>% max
#' 
filter_quantiles_per_segments <-  function(x, low_quant = 0.01, high_quant = 0.99){
  
  tmp <- x$data$counts 
  
  tmp <- apply(tmp,2, function(y) remove_high_low_quant(y, low_quant, high_quant))
  
  x$data$counts <-  tmp
  
  return(x)
  
}


remove_high_low_quant <-  function(x, low_quant, high_quant){
  
  quantiles <-  quantile(x, probs = c(low_quant, high_quant))
  x <-  ifelse(x > quantiles[1],x, round(quantiles[1]))
  x <-  ifelse(x < quantiles[2],x, round(quantiles[2]))
  return(x)
}