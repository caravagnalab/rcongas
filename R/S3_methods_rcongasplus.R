#' Print for an object of class \code{'rcongasplus'}.
#'
#' @param x An object of class \code{'rcongasplus'}.
#' @param ... Unused ellipsis.
#'
#' @return Nothing.
#'
#' @exportS3Method print rcongasplus
#'
#' @examples
#' data("example_object")
#' print(example_object)
print.rcongasplus = function(x, ...)
{
  log = function()
  {
    # cat('\n')
    # cli::cli_rule(crayon::red('LOG'))
    
    cli::cli_h2(crayon::bgMagenta(crayon::white(' LOG ')))
    cat(x$log)
  }
  
  inline_segments_printer = function(x, what = 'segmentation')
  {
    # Color coding
    colors = c(
      `0` = crayon::bgCyan(' '),
      `1` = crayon::bgBlue(' '),
      `2` = crayon::bgGreen(' '),
      `3` = crayon::bgYellow(' '),
      `4` = crayon::bgRed(' '),
      `5` = crayon::bgMagenta(' '),
      `*` = crayon::bgWhite(' ')
    )
    
    # Segments
    if(what == 'segmentation')
      segments_vals = x %>% 
        get_input(what = 'segmentation') %>% 
        group_by(chr) %>% 
        mutate(copies = ifelse(copies > 4, "*", paste(copies))) %>% 
        summarise(
          copies = paste(copies, collapse = '')
        ) %>% 
        summarise(
          copies = paste(copies, collapse = '|')
        ) %>% 
        pull(copies)
    
    # It is a cluster, eg. what = 'C1'
    if(what != 'segmentation')
      segments_vals = get_fit(x, what = 'CNA') %>% 
        filter(cluster == what) %>% 
        deidify() %>% 
        mutate(value = ifelse(value > 4, "*", paste(value))) %>% 
        group_by(chr) %>% 
        summarise(
          value = paste(value, collapse = '')
        ) %>% 
        summarise(
          value = paste(value, collapse = '|')
        ) %>% 
        pull(value)
  
    # Different print strategies
    if(what == 'segmentation')
    {
      cat("\n\t")
      for(char in strsplit(segments_vals,'') %>% unlist)
      {
        if(char %in% names(colors)) cat(colors[char])
        else cat(char)
      }
      cat("\n\n\t", crayon::underline("Ploidy:"), ' ')
      for(char in colors %>% names)
        cat(colors[char], char, '  ')
      cat('\n')
    }
    else{
      cat("\n  ", what, '  ')
      for(char in strsplit(segments_vals,'') %>% unlist)
      {
        if(char %in% names(colors)) cat(colors[char])
        else cat(char)
      }
    }
  }
  
  ############################################################################# 
  stopifnot(inherits(x, "rcongasplus"))
  
  # Sanitize zeroes - first thing to show
  x %>% sanitize_zeroes()
  
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
  
  x %>% inline_segments_printer(what = 'segmentation')
    
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
      {.field {stats_data$rna_events}} non-zero values. \\
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
      {.field {stats_data$atac_events}} non-zero values. \\
      Likelihood: {.field {what_atac_lik}}.'
    )
  else
    cli::cli_alert(
      'ATAC: {crayon::red("not available")}'
    )
  
  # cat("\n")
  
  # Handle clusters
  stats_fit = stat(x, what = 'fit')
  
  if (!is.null(stats_fit))
    cli::cli_h3(
      'Clusters: {.field k = {stats_fit$fit_k}}, model with {.field {stats_fit$fit_IC}} = {.value {round(stats_fit$fit_score, 2)}}.'
    )
  else
    {
      cli::cli_alert_warning('Clusters: {crayon::red("not available")}.')
      
      log()
      invisible(return(0))
    }
  
  cluster_labels = get_fit(x, what = 'CNA') %>% pull(cluster) %>% unique
  
  for(cluster in cluster_labels)
    x %>% inline_segments_printer(what = cluster)
  
  cat('\n\n')
  
  if(x %>% has_rna)
  {
    L_n = stats_fit$fit_mixing_RNA_n 
    L = (25 * (L_n/sum(L_n))) %>% round()
    
    cat('   ', crayon::underline("RNA"), '\n')
    for(cl in names(L))
      cat('\t', cl, ': ', rep("\u25A0", L[cl]) %>% paste(collapse = ''), 'n =', L_n[cl], '\n')  
  }
  
  if(x %>% has_atac)
  {
    L_n = stats_fit$fit_mixing_ATAC_n 
    L = (25 * (L_n/sum(L_n))) %>% round()
    
    cat('   ', crayon::underline("ATAC"), '\n')
    for(cl in names(L))
      cat('\t', cl, ': ', rep("\u25A0", L[cl]) %>% paste(collapse = ''), 'n =', L_n[cl], '\n')  
  }
  
  
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
  
  
  log()
 
  invisible(return(0))
}

#' Plot for an object of class \code{'rcongasplus'}.
#'
#' @param x An object of class \code{'rcongasplus'}.
#' @param ... Unused ellipsis.
#'
#' @return One of the plots available in the package, computed with
#' \code{\link{plot_data}} or \code{\link{plot_data}} functions.
#'
#' @exportS3Method plot rcongasplus
#'
#' @examples
#' data("example_object")
#' plot(example_object)
plot.rcongasplus = function(x, ...)
{
  stopifnot(inherits(x, "rcongasplus"))
  
  if(!('best_fit' %in% (x %>% names))) 
    return(x %>% plot_fit(what = 'CNA'))
  else
    return(x %>% plot_data(what = 'histogram'))
}



