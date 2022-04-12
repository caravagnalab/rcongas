
#' Get a table of DE genes
#' 
#' The actual function for DE is implemented following 
#' the Seurat `FindMarkers` method and as such conserves the same output structure
#'
#' @param x an Rcongas object after DE analysis
#' @param chromosomes filter for specific chromosomes
#' @param cut_pvalue corrected (BH) p-value filter
#' @param cut_lfc log fold change filter
#'
#' @return a table with differentially expressed genes 
#' @export
#'
#' @examples
get_DE_table <- function(x,
                         chromosomes = paste0("chr", c(1:22, "X", "Y")),
                         cut_pvalue = 0.001,
                         cut_lfc = 0.25)
{
  if (!has_DE(x))
  {
    cli::cli_alert_danger("No DE table in this object, compute it first!")
    return(NULL)
  }

  table_data = x$DE$table %>%
    dplyr::filter(p_val_adj < cut_pvalue,
                  abs(avg_log2FC) > cut_lfc,
                  chr %in% chromosomes)

  return(table_data)
}

