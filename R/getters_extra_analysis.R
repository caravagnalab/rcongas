
#' Title
#'
#' @param x
#' @param chromosomes
#' @param cut_pvalue
#' @param cut_lfc
#'
#' @return
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

