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


#' Filter per-segment outliers by quantiles.
#' 
#' @description After mapping counts data to segments, this function
#' can be used to determine quantiles of mapped data, and identify
#' outliers in each segment and modality. 
#' 
#' An outlier can then be removed or capped to the median cell value.
#' The former option introduced 0-counts in the data, which we suggest
#' to check with the \code{stat} function, and possibly remove by using
#' the \code{filter_missing_data} function. Removal can be important as
#' an excess of 0-counts cells (_missing data_) will drive the fit
#' to use 0-mean components.
#' 
#' Capping does not introduce any 0-count cell, and is the suggested choice.
#' The capped values is either a count value or a z-score, depending on
#' the modality type of likelihood.
#' 
#' In both cases pre-filtering normalisation factors are no longer adequate
#' after filtering, and have to be recomputed. If the modality adopts a 
#' Gaussian likelihood this is not a problem, since those are set to 1
#' when the object is created, and remain 1 afterwards. In the case of counts
#' based likelihood like Negative Binomials these are re-computed for all
#' input cells by using the \code{auto_normalisation_factor} function.
#' 
#' Therefore, if custom factors have been computing this function might
#' affect the general signal in the data, and factors should be handled
#' explicitly by the user.
#' 
#' The function requires and returns an (R)CONGAS+ object.
#'
#' @param x An \code{rcongasplus} object.
#' @param lower_quantile The lower quantile, default 1%.
#' @param upper_quantile The upper quantile, default 99%.
#' @param action If \code{"remove"}, outliers will be set to 0. If \code{"cap"}, 
#' outliers will be capped at the median per-cell counts.
#'
#' @return The object \code{x} where outliers have been identified and removec
#' or capped according to the parameters.
#' 
#' @export
#'
#' @examples
#' data('example_object')
#' 
#' # Default
#' print(example_object)
#' 
#' example_object %>% 
#'   filter_outliers(, action = 'cap') %>% 
#'   print()
#' 
#' example_object %>% 
#'   filter_outliers(, action = 'remove') %>% 
#'   print()
# filter_outliers = function(x, 
#                            lower_quantile = 0.01, 
#                            upper_quantile = .99,
#                            action = 'cap')
# {
#   # Manipulate normalisation factors
#   # aux_fun_nf = function(data, modality)
#   # {
#   #   cli::cli_h3("{.field {modality}} recomputing normalisation factors.")
#   #   
#   #   input = data %>% 
#   #     filter(modality == !!modality) 
#   # 
#   #   # Post-mapping we need to update all the factors, not just those that
#   #   # changed because some entries might not have been mapped succesfully
#   #   # to the input segments
#   #   # cells_to_update = input %>% filter(capped_outlier) %>% pull(cell) %>% unique
#   #   cells_to_update = input %>% pull(cell) %>% unique
#   #   
#   #   norm_factors = get_input(x, what = 'normalisation') %>% 
#   #     filter(modality == !!modality) 
#   #   
#   #   norm_factors_unaffected = norm_factors %>% filter(!(cell %in% cells_to_update)) 
#   #   
#   #   if(input$value_type[1] == "G")
#   #     cli::cli_alert("Gaussian likelihood, factors already coherced to 1 will not be changed.")
#   #   else 
#   #   {
#   #     norm_factors_affected = auto_normalisation_factor(
#   #       input %>% filter(cell %in% cells_to_update) # Update those required
#   #       ) %>% mutate(modality = !!modality)
#   #     
#   #     norm_factors = rbind(norm_factors_unaffected, norm_factors_affected)
#   #   }
#   #     
#   #   norm_factors %>% return
#   # }
# 
#   retained_rna = retained_atac = norm_factors_rna = norm_factors_atac = NULL
#   
#   if(x %>% has_rna)
#   {
#     retained_rna = compute_outliers_per_segment(data = get_input(x, what = 'data'), modality = "RNA", lower_quantile, upper_quantile)
#     norm_factors_rna = aux_fun_nf(data = retained_rna, modality = "RNA")
#   }
#   
#   if(x %>% has_atac)
#   {
#     retained_atac = compute_outliers_per_segment(data = get_input(x, what = 'data'), modality = "ATAC", lower_quantile, upper_quantile)
#     norm_factors_atac = aux_fun_nf(data = retained_atac, modality = "ATAC")
#   }
#   
#   retained = bind_rows(retained_rna, retained_atac)
#   norm_factors = bind_rows(norm_factors_rna, norm_factors_atac)
#   
#   # cli::cli_h1("Rebuilding (R)CONGAS+ object")
#   # 
#   # # Rebuild object
#   # 
#   # retained_cells = retained$cell %>% unique()
#   # retained_rna = retained %>% filter(modality == "RNA") %>% deidify()
#   # retained_atac = retained %>% filter(modality == "ATAC") %>% deidify()
#   # retained_factors = x %>% 
#   #   get_input(what = 'normalisation') %>% 
#   #   filter(cell %in% retained_cells)
#   # 
#   # new_x = init(
#   #   rna = retained_rna,
#   #   ata = retained_atac,
#   #   segmentation = x %>% get_input(what = 'segmentation'), # Same
#   #   normalisation_factors = retained_factors,
#   #   rna_likelihood = stat(x)$rna_dtype,
#   #   atac_likelihood = stat(x)$atac_dtype,
#   #   reference_genome = x$reference_genome,
#   #   description =  paste0(x$description, "; post-map q[", lower_quantile, ', ', upper_quantile, ']')
#   # )
#   
#   # Copy input x object, before making a new one
#   x_original = x
#   
#   # Rebuild loses some info, not super important
#   x$input$dataset = retained
#   x$input$normalisation = norm_factors
# 
#   x$input$normalisation = x %>% 
#     get_input(what = 'normalisation') %>% 
#     filter(cell %in% retained$cell %>% unique())
#   
#   x$description = paste0(
#     x$description, 
#     "; post-map q[", lower_quantile, ', ', upper_quantile, '] ', action)
#   
#   if(action == 'remove')
#   {
#     cli::cli_alert_warning("With {.field {action}} this filter will \\
#                            introduce 0-counts for outliers, consider removing them \\
#                            with {.field {'filter_missing_data'}}.")
#   }
#   
#   # if(action == 'cap')
#   # {
#   #   cli::cli_alert_warning("With {.field {action}} we are required to \\
#   #                          recompute normalisation factors per-cell. This will \\
#   #                          be done with {.field {'auto_normalisation_factor'}}.")
#   #   
#   #   rna_norm = ifelse(
#   #     x %>% has_rna,
#   #     auto_normalisation_factor(retained %>% filter(modality == 'RNA'))
#   #   )
#   #   auto_normalisation_factor(
#   #     x$input$dataset %>% m
#   #   )
#   # }
#   
#   # auto_normalisation_factor
#   
#   x$log = paste0(x$log, '\n- ', 
#                       Sys.time(), " Filtered outliers: [", action, 
#                       '|', lower_quantile, '|', upper_quantile, ']')
#   
#   
#   return(x)
# }

#' Filter cells with missing data.
#' 
#' @description After mapping counts data to segments, this function
#' can be used to determine cells with missing data, and remove
#' them The function requires and returns an (R)CONGAS+ object.
#' 
#' This filter works by a proportion, as reported by the \code{\link{stat}} function.
#' 
#' If these cells are not removed, during inference missing values
#' are imputed to be \code{0}. This can create an excess of mixture components
#' fitting 0-counts data.
#'
#' @param x An \code{rcongasplus} object.
#' @param proportion_RNA The RNA proportion cut for a cell to be removed, default 5%.
#' @param proportion_ATAC The ATAC proportion cut for a cell to be removed, default 5%.
#'
#' @return The object \code{x} where 0-counts cells have been removed.
#' @export
#'
#' @examples
#' data('example_object')
#' 
#' # Default
#' print(example_object)
#' 
#' example_object = filter_missing_data(example_object)
#' 
#' # After filtering
#' print(example_object)
filter_missing_data = function(x, proportion_RNA = 0.05, proportion_ATAC = 0.05)
{
  stat_x = stat(x, what = 'data')
  
  x %>% sanitize_zeroes()
  
  all_cells = x$input$dataset$cell %>% unique
  n_orig = all_cells %>% length
  
  # RNA to remove
  if(!is.null(stat_x$zero_counts_cells_RNA))
  {
    to_remove = stat_x$zero_counts_cells_RNA %>% 
      filter(`%` > 100 * proportion_RNA) %>% 
      pull(cell)
    
    n = stat_x$zero_counts_cells_RNA %>% nrow
    
    cli::cli_h3("{.field {'RNA'}} modality, proportions cut {.field {proportion_RNA}}")
    cli::cli_alert("Cells with missing data {.field {n}}, removing {.field {to_remove %>% length}}")

    if(length(to_remove) > 0) all_cells = setdiff(all_cells, to_remove)
  }
  
  # ATAC to remove
  if(!is.null(stat_x$zero_counts_cells_ATAC))
  {
    to_remove = stat_x$zero_counts_cells_ATAC %>% 
      filter(`%` > 100 * proportion_ATAC) %>% 
      pull(cell)
    
    n = stat_x$zero_counts_cells_ATAC %>% nrow
    
    cli::cli_h3("{.field {'ATAC'}} modality, proportions cut {.field {proportion_ATAC}}")
    cli::cli_alert("Cells with missing data {.field {n}}, removing {.field {to_remove %>% length}}")
  
    if(length(to_remove) > 0) all_cells = setdiff(all_cells, to_remove)
  }
  
  # Filter
  if(n_orig < (all_cells %>% length))
  {
    x$input$dataset = x$input$dataset %>% filter(cell %in% all_cells)
    x$input$normalisation = x$input$normalisation %>% filter(cell %in% all_cells)
    
    x$description = paste0(x$description, "; post-map 0-cells >", proportion_RNA, ', >', proportion_ATAC)
  }
  
  return(x)
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


# Compute segment outliers for one modality
compute_outliers_per_segment = function(x, modality, lower_quantile, upper_quantile)
{
  data = get_input(x, what = 'data') %>% filter(modality == !!modality)
  
  cli::cli_h3("{.field {modality}} outliers detection via quantiles: lower {.field {lower_quantile}}, upper {.field {upper_quantile}}.")
  
  input =  data %>% 
    filter(modality == !!modality) %>% 
    mutate(original_value = value)
  
  if(input$value_type[1] == "NB")
    scaled_data = normalise_modality(
      input,
      get_input(x, what = 'normalisation') %>% filter(modality == !!modality)
    )
  else scaled_data = input
  
  map_quantiles = scaled_data %>% 
    group_by(segment_id, modality) %>% 
    summarise(
      q_min = quantile(value, lower_quantile),
      q_max = quantile(value, upper_quantile),
      .groups = 'drop'
    )
  
  x_mapped = scaled_data  %>% 
    left_join(map_quantiles, by = c('segment_id', 'modality')) %>% 
    mutate(outlier = value > q_max | value < q_min)
  
  x_mapped %>% return
}

filter_outliers = function(x, 
                           frequency_cutoff = 0.2 * stat(x)$nsegments,
                           lower_quantile = 0.01, 
                           upper_quantile = .99)
{
  # Manipulate normalisation factors
  aux_fun_nf = function(data, modality)
  {
    cli::cli_h3("{.field {modality}} recomputing normalisation factors.")

    input = data %>%
      filter(modality == !!modality)

    # Post-mapping we need to update all the factors, not just those that
    # changed because some entries might not have been mapped succesfully
    # to the input segments
    # cells_to_update = input %>% filter(capped_outlier) %>% pull(cell) %>% unique
    cells_to_update = input %>% pull(cell) %>% unique

    norm_factors = get_input(x, what = 'normalisation') %>%
      filter(modality == !!modality)

    norm_factors_unaffected = norm_factors %>% filter(!(cell %in% cells_to_update))

    if(input$value_type[1] == "G")
      cli::cli_alert("Gaussian likelihood, factors already coherced to 1 will not be changed.")
    else
    {
      norm_factors_affected = auto_normalisation_factor(
        input %>% filter(cell %in% cells_to_update) # Update those required
        ) %>% mutate(modality = !!modality)

      norm_factors = rbind(norm_factors_unaffected, norm_factors_affected)
    }

    norm_factors %>% return
  }

  retained_rna = retained_atac = norm_factors_rna = norm_factors_atac = NULL
  
  if(x %>% has_rna)
  {
    retained_rna = compute_outliers_per_segment(x, modality = "RNA", lower_quantile, upper_quantile) %>% 
      group_by(cell) %>% 
      summarise(n_outlier = sum(outlier)) %>% 
      arrange(desc(n_outlier)) %>% 
      mutate(remove = n_outlier > frequency_cutoff)
    
    ncell = sum(retained_rna$remove)
    nprop = ((ncell/stat(x)$ncells_RNA) * 100) %>% round
    cli::cli_alert("{.field {ncell}} out of {.field {nrow(retained_rna)}} will be removed ({.field {nprop}%})")
  }
  
  if(x %>% has_atac)
  {
    retained_atac = compute_outliers_per_segment(x, modality = "ATAC", lower_quantile, upper_quantile) %>% 
      group_by(cell) %>% 
      summarise(n_outlier = sum(outlier)) %>% 
      arrange(desc(n_outlier)) %>% 
      mutate(remove = n_outlier > frequency_cutoff)
    
    ncell = sum(retained_atac$remove)
    nprop = ((ncell/stat(x)$ncells_ATAC) * 100) %>% round
    cli::cli_alert("{.field {ncell}} out of {.field {nrow(retained_atac)}} will be removed ({.field {nprop}%})")
  }
  
  retained = bind_rows(retained_rna, retained_atac) %>% 
    filter(!remove) %>% 
    pull(cell)
  
  # Rebuild loses some info, not super important
  x$input$dataset = x$input$dataset %>% filter(cell %in% retained)
  x$input$normalisation = x$input$normalisation %>% filter(cell %in% retained)
  
  x$log = paste0(x$log, '\n- ', 
                 Sys.time(), " Filtered outliers: [", frequency_cutoff, 
                 '|', lower_quantile, '|', upper_quantile, ']')
  
  return(x)
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
