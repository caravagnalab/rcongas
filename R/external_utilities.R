smoothie = function(segments, reference = 'GRCh38')
{
  reference = CNAqc:::get_reference(ref = reference)
  
  # Use the CNAqc implementation, requires some tweaks
  # segments$Major = segments$minor = segments$CCF = 1
  
  # CNAqc:::split_cna_to_arms(x = list(reference_genome = 'GRCh38'), segments)
  
  smoothed_segments = NULL
  
  for (chr in unique(segments$chr))
  {
    chr_segments = segments %>% filter(chr == !!chr)
    
    if (nrow(chr_segments) == 1) {
      smoothed_segments = bind_rows(smoothed_segments,
                                    chr_segments)
      next
    }
    
    index = 1
    
    repeat {
      template = chr_segments[index, ]
      
      j = index
      repeat {
        if (j == nrow(chr_segments))
          break
        
        copies_match = chr_segments$copies[j + 1] == chr_segments$copies[j]
        
        if (copies_match)
          j = j + 1
        else
          break
      }
      
      template$to = chr_segments$to[j]
      template$length = template$to - template$from
      smoothed_segments = bind_rows(smoothed_segments,
                                    template)
      if (j == nrow(chr_segments))
        break
      index = j + 1
    }
    
  }
  
  
}

#' Filter genes from an RNA tibble.
#'
#' @param x The tibble for RNA.
#' @param what A code string. If this contains "r" ribosomal genes are removed.
#' If this contains "m" mitocondrial genes are removed. Both are matched by
#' a regex.
#' @param specials A list of genes to remove a priori.
#'
#' @return The RNA tibble withouth the required gene entries.
#' @export
#'
#' @examples
#' data('example_input')
#' filter_known_genes(example_input$x_rna)
filter_known_genes = function(x, what = 'rm', specials = 'MALAT1')
{
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

#' Cap observed values by quantile.
#'
#' @param x An input RNA/ATAC dataset where entries are indexeable by genomic 
#' coordinate: "chr", "from" and "to".
#' @param upper_quantile The maximum quantile to determine cuts. If a value
#' is above the quantile it is capped.
#'
#' @return The input data with capped entries.
#' @export
#'
#' @examples
#' data('example_input')
#' cap_values_by_quantile(example_input$x_rna, upper_quantile = .98)
#' cap_values_by_quantile(example_input$x_atac, upper_quantile = .98)
cap_values_by_quantile = function(x, upper_quantile = .98)
{
  map_quantiles = x %>% 
    group_by(chr, from, to) %>% 
    summarise(
      q_max = quantile(value, upper_quantile),
      .groups = 'drop'
    )
  
  x = x %>% 
    left_join(map_quantiles, by = c('chr', 'from', 'to')) %>% 
    mutate(cap = value > q_max)
  
  r_cap = x %>% filter(cap)
  
  if(nrow(r_cap) > 0)
  {
    cli::cli_h3("Upper quantile {.field {upper_quantile}}")
    
    cli::cli_alert_info("{.field n = {nrow(r_cap)}} entries to cap")
    
    cat("\n")
    r_cap  %>% print
    
    x = x %>% 
      mutate(value = ifelse(cap, q_max, value))
  }
  
  return(x %>% select(-q_max, -cap))
}

#' Filter observed values by quantile.
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
#' filter_values_by_quantile(example_input$x_rna, upper_quantile = .98)
#' filter_values_by_quantile(example_input$x_atac, upper_quantile = .98)
filter_values_by_quantile = function(x, upper_quantile = .98)
{
  map_quantiles = x %>% 
    group_by(chr, from, to) %>% 
    summarise(
      q_max = quantile(value, upper_quantile),
      .groups = 'drop'
    )
  
  x = x %>% 
    left_join(map_quantiles, by = c('chr', 'from', 'to')) %>% 
    mutate(del = value > q_max)
  
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
  
  return(x %>% select(-q_max, -cap))
}


# rna %>%
#   group_by(gene) %>%
#   summarise(sd = sd(value)) %>%
#   arrange(desc(sd))
#
# hist(x %>% filter(gene == "FTL") %>% pull(value), breaks = 100)
#
# IQR(x %>% filter(gene == "FTL") %>% pull(value))
