smoothie = function(segments, reference = 'GRCh38')
{
  stop("TODO")
  
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

#' Filter cells by upper and lower quantile of their normalisation factors.
#'
#' @param x A tibble with a cell and normalisation_factor column.
#' @param lower_quantile The minimum quantile to determine cuts. 
#' @param upper_quantile The maximum quantile to determine cuts. 
#'
#' @return The input data with normalisation factors in between the quantiles range.
#' @export
#'
#' @examples
#' data('example_input')
filter_cells_by_quantile = function(x, lower_quantile = 0.05, upper_quantile = .95)
{
  m_quantile = quantile(x$normalisation_factor, upper_quantile)
  l_quantile = quantile(x$normalisation_factor, lower_quantile)
  
  x = x %>% 
    mutate(del = normalisation_factor > m_quantile | normalisation_factor < l_quantile)
  
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
  
  return(x %>% select(-del))
}






outliers_inspector = function(x, x_original)
{
  all_plots = lapply(x$input$segmentation$segment_id, function(s_id)
  {
    # x
    what_x_rna = x$input$dataset %>% filter(segment_id == s_id, modality == 'RNA')
    what_x_atac = x$input$dataset %>% filter(segment_id == s_id, modality == 'ATAC')
    
    if((x %>% has_rna) && what_x_rna$value_type[1] == 'NB')
      what_x_rna = normalise_modality(what_x_rna, x$input$normalisation %>% filter(modality == 'RNA'))
    
    if((x %>% has_atac) && what_x_atac$value_type[1] == 'NB')
      what_x_atac = normalise_modality(what_x_atac, x$input$normalisation %>% filter(modality == 'ATAC'))
    
    what_x = bind_rows(what_x_rna, what_x_atac)
    what_x_cellids = what_x %>% filter(outlier_cap) %>% pull(cell)
    
    # x original
    what_x_rna2 = x_original$input$dataset %>% filter(segment_id == s_id, modality == 'RNA')
    what_x_atac2 = x_original$input$dataset %>% filter(segment_id == s_id, modality == 'ATAC')
    
    if((x_original %>% has_rna) && what_x_rna2$value_type[1] == 'NB')
      what_x_rna2 = normalise_modality(what_x_rna2, x_original$input$normalisation %>% filter(modality == 'RNA'))
    
    if((x_original %>% has_atac) && what_x_atac2$value_type[1] == 'NB')
      what_x_atac2 = normalise_modality(what_x_atac2, x_original$input$normalisation %>% filter(modality == 'ATAC'))
    
    what_x_original = bind_rows(what_x_rna2, what_x_atac2)
    what_x_original$outlier = what_x_original$cell %in% what_x_cellids
    
    
    cowplot::plot_grid(
      ggplot(what_x) +
        geom_histogram(aes(value, fill = outlier_cap), bins = 70) +
        facet_wrap(segment_id~modality, scales = 'free') +
        theme_linedraw(base_size = 9) +
        labs(title = "After filtering (red: capped)") +
        scale_fill_manual(values = c(`TRUE` = 'indianred3', `FALSE` = 'gray'))  +
        guides(fill = FALSE),
      ggplot(what_x_original) +
        geom_histogram(aes(value, fill = outlier), bins = 70) +
        facet_wrap(segment_id~modality, scales = 'free') +
        theme_linedraw(base_size = 9) +
        labs(title = "Before filtering (red: outliers)") +
        scale_fill_manual(values = c(`TRUE` = 'indianred3', `FALSE` = 'gray'))  +
        guides(fill = FALSE),
      ncol = 1, 
      nrow = 2
    )
  })
  
  return(all_plots)
}





outliers_inspector = function(x, x_original)
{
  all_plots = lapply(x$input$segmentation$segment_id, function(s_id)
  {
    # x
    what_x_rna = x$input$dataset %>% filter(segment_id == s_id, modality == 'RNA')
    what_x_atac = x$input$dataset %>% filter(segment_id == s_id, modality == 'ATAC')
    
    if((x %>% has_rna) && what_x_rna$value_type[1] == 'NB')
      what_x_rna = normalise_modality(what_x_rna, x$input$normalisation %>% filter(modality == 'RNA'))
    
    if((x %>% has_atac) && what_x_atac$value_type[1] == 'NB')
      what_x_atac = normalise_modality(what_x_atac, x$input$normalisation %>% filter(modality == 'ATAC'))
    
    what_x = bind_rows(what_x_rna, what_x_atac)
    what_x_cellids = what_x %>% filter(outlier_cap) %>% pull(cell)
    
    # x original
    what_x_rna2 = x_original$input$dataset %>% filter(segment_id == s_id, modality == 'RNA')
    what_x_atac2 = x_original$input$dataset %>% filter(segment_id == s_id, modality == 'ATAC')
    
    if((x_original %>% has_rna) && what_x_rna2$value_type[1] == 'NB')
      what_x_rna2 = normalise_modality(what_x_rna2, x_original$input$normalisation %>% filter(modality == 'RNA'))
    
    if((x_original %>% has_atac) && what_x_atac2$value_type[1] == 'NB')
      what_x_atac2 = normalise_modality(what_x_atac2, x_original$input$normalisation %>% filter(modality == 'ATAC'))
    
    what_x_original = bind_rows(what_x_rna2, what_x_atac2)
    what_x_original$outlier = what_x_original$cell %in% what_x_cellids
    
    
    cowplot::plot_grid(
      ggplot(what_x) +
        geom_histogram(aes(value, fill = outlier_cap), bins = 70) +
        facet_wrap(segment_id~modality, scales = 'free') +
        theme_linedraw(base_size = 9) +
        labs(title = "After filtering (red: capped)") +
        scale_fill_manual(values = c(`TRUE` = 'indianred3', `FALSE` = 'gray'))  +
        guides(fill = FALSE),
      ggplot(what_x_original) +
        geom_histogram(aes(value, fill = outlier), bins = 70) +
        facet_wrap(segment_id~modality, scales = 'free') +
        theme_linedraw(base_size = 9) +
        labs(title = "Before filtering (red: outliers)") +
        scale_fill_manual(values = c(`TRUE` = 'indianred3', `FALSE` = 'gray'))  +
        guides(fill = FALSE),
      ncol = 1, 
      nrow = 2
    )
  })
  
  return(all_plots)
}
