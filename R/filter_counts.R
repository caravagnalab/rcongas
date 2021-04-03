#' Filter genes from an RNA tibble.
#'
#' @param x The tibble for RNA.
#' @param what A code string. If this contains "r" ribosomal genes are removed.
#' If this contains "m" mitocondrial genes are removed. Both are matched by
#' a regular expression.
#' @param specials A list of genes to remove a priori.
#'
#' @return The RNA tibble without the required gene entries.
#' @export
#'
#' @examples
#' data('example_input')
#' filter_known_genes(example_input$rna)
filter_known_genes = function(x, what = 'rm', specials = 'MALAT1')
{
  if(!('gene' %in% names(x)))
  {
    cli::cli_alert_danger("The {.field {'gene'}} coulumn is missing from input, returning input data.")
    return(x)
  }
    
  all_genes = x$gene %>% unique
  codes = strsplit(what, split = '')[[1]]
  
  if ('r' %in% codes)
  {
    # Ribosomal genes
    cli::cli_h3("Ribosomal")
    
    ribo = all_genes %>% grepl(pattern = 'RP[SL]')
    nribo = ribo %>% sum
    
    cli::cli_alert_info("{.field n = {nribo}} genes found.")
    
    if (nribo > 0) {
      cat(all_genes[ribo] %>% head %>% paste(collapse = ', '), ', ...\n')
      all_genes = setdiff(all_genes, all_genes[ribo])
    }
  }
  
  if ('m' %in% codes)
  {
    # Mitochondrial genes
    cli::cli_h3("Mitochondrial")
    
    mito = all_genes %>% grepl(pattern = 'MT')
    nmito = mito %>% sum
    
    cli::cli_alert_info("{.field {nmito}} genes found.")
    
    if (nmito > 0)
    {
      cat(all_genes[mito] %>% head %>% paste(collapse = ', '), ', ...\n')
      all_genes = setdiff(all_genes, all_genes[mito])
    }
  }
  
  if (length(specials) >= 1)
  {
    # Ribosomal genes
    cli::cli_h3("Specials")
    
    spec = all_genes %>% intersect(specials)
    nspec = spec %>% length
    
    cli::cli_alert_info("{.field n = {nspec}} genes found.\n")
    
    if(nspec > 0) all_genes = setdiff(all_genes, spec)
  }
  
  return(x %>% filter(gene %in% all_genes))
}

#' Filter observed counts by quantile.
#'
#' @param x An input RNA/ATAC dataset where entries are indexeable by genomic 
#' coordinate: "chr", "from" and "to".
#' @param upper_quantile The maximum quantile to determine cuts. If a value
#' is above the quantile it is removed
#'
#' @return The input data with removed entries.
#' @export
#'
#' @examples
#' data('example_input')
#' filter_counts_by_quantile(example_input$rna, upper_quantile = .98)
#' filter_counts_by_quantile(example_input$atac, upper_quantile = .98)
filter_counts_by_quantile = function(x, upper_quantile = .98)
{
  map_quantiles = x %>% 
    group_by(chr, from, to) %>% 
    summarise(
      q_max = quantile(value, upper_quantile),
      # sd = sd(value),
      # median = median(value),
      .groups = 'drop'
    )
  
  x = x %>% 
    left_join(map_quantiles, by = c('chr', 'from', 'to')) %>% 
    mutate(del = value > q_max)
  
  # x %>% arrange(desc(sd))
  # x %>% arrange(desc(value))
  # hist(x %>% filter(chr == 'chr9', to == '140135734', from== '140135234') %>% pull(value), breaks = 100) 
  # 
  # m_atac = quantile(x$value, .98)
  # x %>% filter
  
  # hist(norm_atac$normalisation_factor, breaks = 100)
  
  
  # hist(log(x$sd), breaks = 100) 
  
  # ggplot(x %>% filter(log(x$sd) > 2), aes(x = median, y = sd, color = del)) +
  #   geom_point()
  
  r_cap = x %>% filter(del)
  
  if(nrow(r_cap) > 0)
  {
    cli::cli_h3("Upper quantile {.field {upper_quantile}}")
    
    cli::cli_alert_info("{.field n = {nrow(r_cap)}} entries to remove")
    
    cat("\n")
    r_cap  %>% print
    
    x = x %>% 
      filter(!del)
  }
  
  return(x %>% select(-q_max, -del))
}

#' #' Cap observed values by quantile.
#' #'
#' #' @param x An input RNA/ATAC dataset where entries are indexeable by genomic 
#' #' coordinate: "chr", "from" and "to".
#' #' @param upper_quantile The maximum quantile to determine cuts. If a value
#' #' is above the quantile it is capped.
#' #'
#' #' @return The input data with capped entries.
#' #' @export
#' #'
#' #' @examples
#' #' data('example_input')
#' #' cap_values_by_quantile(example_input$x_rna, upper_quantile = .98)
#' #' cap_values_by_quantile(example_input$x_atac, upper_quantile = .98)
#' cap_values_by_quantile = function(x, upper_quantile = .98)
#' {
#'   map_quantiles = x %>% 
#'     group_by(chr, from, to) %>% 
#'     summarise(
#'       q_max = quantile(value, upper_quantile),
#'       .groups = 'drop'
#'     )
#'   
#'   x = x %>% 
#'     left_join(map_quantiles, by = c('chr', 'from', 'to')) %>% 
#'     mutate(cap = value > q_max)
#'   
#'   r_cap = x %>% filter(cap)
#'   
#'   if(nrow(r_cap) > 0)
#'   {
#'     cli::cli_h3("Upper quantile {.field {upper_quantile}}")
#'     
#'     cli::cli_alert_info("{.field n = {nrow(r_cap)}} entries to cap")
#'     
#'     cat("\n")
#'     r_cap  %>% print
#'     
#'     x = x %>% 
#'       mutate(value = ifelse(cap, q_max, value))
#'   }
#'   
#'   return(x %>% select(-q_max, -cap))
#' }