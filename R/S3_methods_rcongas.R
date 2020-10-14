#' Summary for an object of class \code{'rcongas'} is a print.
#'
#' @param object An obj of class \code{'rcongas'}.
#' @param ...
#'
#' @return See \code{\link{print}}.
#' @exportS3Method summary rcongas
#'
#' @examples
#' data(fit_example)
#' summary(fit_example$best)
summary.rcongas = function(object, ...) {
  print.rcongas(object, ...)
}

#' Summaries for an object of class \code{'rcongas'} is like a print.
#'
#' @param x An obj of class \code{'rcongas'}.
#' @param ...
#'
#' @return nothing.
#' @exportS3Method print rcongas
#' @importFrom crayon white red green yellow black bgYellow blue bold
#' @importFrom cli cli_rule cli_text
#' @importFrom clisymbols symbol
#'
#' @examples
#' data(fit_example)
#' print(fit_example$best)
print.rcongas = function(x, ...)
{
  stopifnot(inherits(x, "rcongas"))

  stats_data = Rcongas::get_dataset_stats(x)

  cli::cli_rule(
    paste(
      crayon::bgYellow(
        crayon::black("[ Rcongas ] {.value {Rcongas:::get_model_description(x)}}")
      ),
      '{.field n = {stats_data$ncells}} cells with {.field k = {stats_data$nsegments}} segments, grouped in {.field k = {stats_data$clusters_k}} clusters.'
    )
  )

  myp = function (m, symbol = "clisymbols::symbol$pointer")
  {
    paste("{", symbol, "}", m)
  }


  pi = (stats_data$clusters_pi * 100) %>% round(2)

  for (i in names(stats_data$clusters_n))
    myp(
      paste0(
        "Cluster {.field {i}}, n = {.value {stats_data$clusters_n[i]}} [{.value { pi[i]}}% of total cells]."
      ),
      symbol = 'clisymbols::symbol$bullet'
    ) %>% cli::cli_text()

  # cat('\n')

  paste(
    "{crayon::white(clisymbols::symbol$info)} Model scored with {.field {stats_data$score_type}} = {.value {round(stats_data$score, 2)}}"
  ) %>%
    cli::cli_text()

  if (Rcongas:::has_DE(x))
  {
    DE_table = Rcongas::get_DE_table(x, cut_pvalue = 0.01, cut_lfc = 0.25)

    nde = DE_table %>% nrow()

    myp(
      paste0(
        "Differential Expression analysis for {.field {x$DE$params$clone1}} vs {.field {x$DE$params$clone2}} with the {.field {x$DE$params$method}} method:
            {.field {nde}} DE genes at level alpha 0.01, with |lfc| > 0.25."
      ),
      symbol = 'clisymbols::symbol$bullet'
    ) %>% cli::cli_text()
  }

}


#' Plot an Rcongas fit.
#'
#' @param ...
#'
#' @return A ggplot object for the plot.
#'
#' @import ggplot2
#'
#' @exportS3Method plot rcongas
#' @export
#'
#' @examples
#' x=1
plot.rcongas = function(x, ...)
{
  # default plot
  plot_gw_cna_profiles(x, whole_genome = TRUE)
}
