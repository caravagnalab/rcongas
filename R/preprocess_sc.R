preprocess_sc <- function(mat,filter_genes = T, filter_cells = T, filter_upper_quantile = F,perc_cells_gene_expr = 0.05, n_genes_in_cell = 3000, upper_quantile = 0.95) {



  if(filter_cells){
    cat("Filtering cells\n")
    filter_on_cell <- Matrix::rowSums(mat > 0) > n_genes_in_cell
    removed <- sum(!filter_on_cell)
    cli::cli_alert_success(paste0("removed ", removed, " cells from a total of ", nrow(mat), "\n"))
    mat <- mat[filter_on_cell,]
  }

  if(filter_genes){
    cat("Filtering genes\n")
    filter_on_genes <- Matrix::colSums(mat > 0) > (nrow(mat) * perc_cells_gene_expr)
    removed <- sum(!filter_on_genes)
    cli::cli_alert_success(paste0("removed ", removed, " genes from a total of ", ncol(mat), "\n"))
    mat <- mat[,filter_on_genes]

  }


  if(filter_upper_quantile){
    cat("Filtering genes on quantile\n")

    qq <- quantile(Matrix::colMeans(mat), c(upper_quantile))
    filter_on_quant <- Matrix::colMeans(mat) < qq
    removed <- sum(!filter_on_quant)
    cli::cli_alert_success(paste0("removed ", removed, " genes from a total of ", ncol(mat), "\n"))
    mat <- mat[,filter_on_quant]

  }

  return(mat)
}
