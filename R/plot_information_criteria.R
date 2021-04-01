plot_information_criteria <- function(x){
  
  ic <- c("lk", "AIC", "BIC", "ICL")
  res <- matrix(0,nrow = length(x$inference$models), ncol = length(ic))
  
  colnames(res) <- ic
  rownames(res) <- sapply(x$inference$models, function(y) y$run_information$input_hyper_params$K)
    
  for(i in ic){
    tmp <- recalculate_information_criteria(x, i)
    res[,i] <- tmp$inference$model_selection$IC
  }
  
  res_l <-  res %>% reshape2::melt()
  colnames(res_l) <- c("clusters", "IC", "value")
  p1 <- ggplot(data = res_l, aes(x = clusters, y = value, color = IC)) + geom_point(size = 1.5) + geom_line(size = 0.5) + CNAqc:::my_ggplot_theme() 
  return(p1)
    
}
