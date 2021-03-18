#' Print for an object of class \code{'rcongasplus'}.
#'
#' @param x Object of class \code{'rcongasplus'}.
#' @param ... Unused.
#'
#' @return Nothing.
#'
#' @exportS3Method print rcongasplus
#'
#' @importFrom crayon white red green yellow black bgYellow blue bold
#' @importFrom cli cli_rule cli_text
#' @importFrom clisymbols symbol
#'
#' @examples
print.rcongasplus = function(x, ...)
{
  stopifnot(inherits(x, "rcongasplus"))
  
  stats_data = stat(x, what = 'data')
  
  # Header
  cli::cli_rule(paste0(
    crayon::bgYellow(crayon::black("[ (R)CONGAS+ ]")),
    crayon::blue(" {.value {x$description}}")
  ))
  # cat('\n')
  
  cli::cli_h3("CNA segments (reference: {.field {x$reference_genome}})")


  # Segments
  meanp = get_segmentation(x) %>% pull(copies) %>% mean %>% round(2)
  cli::cli_alert('Input {.field {stats_data$nsegments}} CNA segments, mean ploidy {.field {meanp}}.')
    
  cli::cli_h3("Modalities")
  
  # RNA Data
  what_rna_lik = case_when(
    stats_data$rna_dtype == "NB" ~ "Negative Binomial",
    stats_data$rna_dtype == "P" ~ "Poisson",
    stats_data$rna_dtype == "G" ~ "Gaussian"
  )
  
  if("RNA" %in% stats_data$modalities)
    cli::cli_alert(
      'RNA: {.field {stats_data$ncells_RNA}} cells with \\
      {.field {stats_data$rna_genes}} mapped genes, \\
      {.field {stats_data$rna_events}} events annotated. \\
      Likelihood: {.field {what_rna_lik}}.'
    )
  else
    cli::cli_alert(
      ' RNA: {crayon::red("not available")}'
    )
  
  # ATAC data
  what_atac_lik = case_when(
    stats_data$atac_dtype == "NB" ~ "Negative Binomial",
    stats_data$atac_dtype == "P" ~ "Poisson",
    stats_data$atac_dtype == "G" ~ "Gaussian"
  )
  
  if("ATAC" %in% stats_data$modalities)
    cli::cli_alert(
      'ATAC: {.field {stats_data$ncells_ATAC}} cells with \\
      {.field {stats_data$atac_peak}} mapped peaks, \\
      {.field {stats_data$atac_events}} events annotated. \\
      Likelihood: {.field {what_atac_lik}}.'
    )
  else
    cli::cli_alert(
      'ATAC: {crayon::red("not available")}'
    )
  
  cat("\n")
  
  if (!is.null(stats_data$clusters_k))
    cli::cli_alert_info(
      'Clusters: {.field k = {stats_data$clusters_k}}, model with {.field {stats_data$score_type}} = {.value {round(stats_data$score, 2)}}.'
    )
  else
    cli::cli_alert_warning('Clusters: {crayon::red("not available")}.')
  
  # cat('\n')
  
  
  # myp = function (m, symbol = "clisymbols::symbol$pointer")
  # {
  #   paste("{", symbol, "}", m)
  # }
  # 
  # if (!has_inference(x))
  #   return()
  # pi = (stats_data$clusters_pi * 100) %>% round(2)
  # 
  # cat('\n')
  # 
  # for (i in names(stats_data$clusters_n))
  #   myp(
  #     paste0(
  #       "Cluster {.field {i}}, n = {.value {stats_data$clusters_n[i]}} [{.value { pi[i]}}% of total cells]."
  #     ),
  #     symbol = 'clisymbols::symbol$bullet'
  #   ) %>% cli::cli_text()
}
