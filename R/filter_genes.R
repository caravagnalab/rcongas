

#' Removing highly variable and highly expressed genes
#'
#' @param cmat 
#' @param quantile_mean 
#' @param quantile_var 
#' @param plot 
#'
#' @return
#' @export
#'
#' @examples
#' 
filter_genes <- function(cmat, quantile_mean = 0.95, quantile_var = 0.95, plot = TRUE){
    
    joined <- cmat
    mean <-  log(joined[,1])
    variance <-  sqrt(joined[,2])
    quantile_mean_val <- quantile(mean, probs = quantile_mean)
    quantile_var_val <- quantile(variance, probs = quantile_var)
    
    joined$to_remove <- ifelse(mean > quantile_mean_val & variance > quantile_var_val, T, F)
    

    filtered_genes <-  rownames(joined)[which(joined$to_remove)]
    
    if(plot){
      
      p1 <-  ggplot(joined, aes(x = mean, y = variance, color = to_remove)) +
        geom_point(size = 1) +
        labs(title = paste0("Mean variance relationship")) +
        scale_color_manual("Removing?",values = c("grey", "red")) +
        CNAqc:::my_ggplot_theme() + 
        labs(
          x = "Log mean",
          y = "Standard deviation"
        ) +
        geom_vline(xintercept = quantile_mean_val, linetype = 'dashed', color = 'black') +
        geom_hline(yintercept = quantile_var_val, linetype = 'dashed', color = 'black') +
        guides(shape = guide_legend(
          paste0("Outside q={", paste(c(quantile_mean_val, quantile_var_val), collapse = ', '), '} n = ', sum(joined$to_remove))
        ))
      plot(p1)  
    }
    
    return(filtered_genes)
    
}