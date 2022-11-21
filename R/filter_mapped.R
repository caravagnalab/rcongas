#' Filter segments.
#' 
#' @description After mapping data to segments, this function
#' can be used to remove segments that are too short, or that have too
#' few mapped RNA genes or ATAC peaks. A parameters allows to input a list
#' of segments that will be retained no matter what.
#' 
#' The function requires and returns an (R)CONGAS+ object.
#'
#' @param x An \code{rcongasplus} object.
#' @param length Mimium size in number of bases for each segment. This is \code{0}
#' by default, and therefore non effective.
#' @param RNA_genes Required number of RNA genes, if RNA is available. This is
#' \code{50} by default.
#' @param ATAC_peaks Required number of ATAC peaks, if ATAC is available. This is
#' \code{50} by default.
#' @param segments A list of segment ids to retain regardless of the filters.
#' 
#' @return The object \code{x} where segments have been identified and 
#' removed.
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
#'   filter_segments() %>% 
#'   print()
#' 
#' # Minimum size
#' example_object %>% 
#'   filter_segments(length = 1e7) %>% 
#'   print()
filter_segments = function(x, 
                           length = 0,
                           RNA_genes = 50,
                           ATAC_peaks = 50,
                           segments = list()
                           )
{
  x %>% sanitize_obj()
  
  retained_segments = get_input(x, 'segmentation') %>%
    mutate(
      force = segment_id %in% segments,
      RNA = ifelse(has_rna((x)),
                   RNA_genes > !!RNA_genes,
                   TRUE),
      ATAC = ifelse(has_atac((x)),
                    ATAC_peaks > !!ATAC_peaks,
                    TRUE),
      L = (as.double(to) - as.double(from)) > !!length
    ) %>%
    filter(
      (force) | (RNA & L & ATAC))
  
  cli::cli_h3("Segments filter")
  cli::cli_alert("{.field {nrow(retained_segments)}} retained segments out of {.field {stat(x)$nsegments}}.")
  
  if(nrow(retained_segments) == 0) stop("All segments removed with these parameters.")
  
  x$input$segmentation = retained_segments %>% dplyr::select(-RNA, -ATAC)
  
  x$input$dataset = x$input$dataset %>% 
    filter(segment_id %in% retained_segments$segment_id)

  x$log = paste0(x$log, '\n- ', 
                 Sys.time(), " Filtered segments: [", length, 
                 '|', RNA_genes, '|', ATAC_peaks, ']')
  
  return(x)
}

#' Filter per-segment outliers by quantiles.
#' 
#' @description After mapping counts data to segments, this function
#' can be used to determine quantiles of mapped data, and identify
#' outliers in each segment and modality. 
#' 
#' An outliers is an entry for a cell/segment pair; with this function
#' we compute how often a certain cell is marked as an outlier. Then
#' the function removes cells that are flagged as containing outliers
#' more then a certain input cutoff. This helps picking up which cells
#' are often showing counts that seem deviating from the main signal
#' in the data.
#' 
#' The function requires and returns an (R)CONGAS+ object.
#'
#' @param x An \code{rcongasplus} object.
#' @param lower_quantile The lower quantile, default 3%.
#' @param upper_quantile The upper quantile, default 97%.
#' @param frequency_cutoff The cutoff to determine if a cell should be removed
#' or not from the data because it has too many outliers. By default, this cut
#' is 20% of the input number of segments.
#' 
#' @return The object \code{x} where outlier cells have been identified and 
#' removed.
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
#'   filter_outliers() %>% 
#'   print()
#' 
#' example_object %>% 
#'   filter_outliers(, action = 'remove') %>% 
#'   print()
filter_outliers = function(x, 
                           frequency_cutoff = 0.2 * stat(x)$nsegments,
                           lower_quantile = 0.03, 
                           upper_quantile = .97)
{
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

#   retained_rna = retained_atac = norm_factors_rna = norm_factors_atac = NULL
  
#   if(x %>% has_rna)
#   {
#     retained_rna = compute_outliers_per_segment(data = get_input(x, what = 'data'), modality = "RNA", lower_quantile, upper_quantile)
#     norm_factors_rna = aux_fun_nf(data = retained_rna, modality = "RNA")
#   }
  
#   if(x %>% has_atac)
#   {
#     retained_atac = compute_outliers_per_segment(data = get_input(x, what = 'data'), modality = "ATAC", lower_quantile, upper_quantile)
#     norm_factors_atac = aux_fun_nf(data = retained_atac, modality = "ATAC")
#   }
  
#   retained = bind_rows(retained_rna, retained_atac)
#   norm_factors = bind_rows(norm_factors_rna, norm_factors_atac)
  
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
  
#   # Copy input x object, before making a new one
#   x_original = x
  
#   # Rebuild loses some info, not super important
#   x$input$dataset = retained
#   x$input$normalisation = norm_factors

#   x$input$normalisation = x %>% 
#     get_input(what = 'normalisation') %>% 
#     filter(cell %in% retained$cell %>% unique())
  
#   x$description = paste0(
#     x$description, 
#     "; post-map q[", lower_quantile, ', ', upper_quantile, '] ', action)
  
#   if(action == 'remove')
#   {
#     cli::cli_alert_warning("With {.field {action}} this filter will \\
#                            introduce 0-counts for outliers, consider removing them \\
#                            with {.field {'filter_missing_data'}}.")
#   }
  
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
  
#   # auto_normalisation_factor
  
#   x$log = paste0(x$log, '\n- ', 
#                       Sys.time(), " Filtered outliers: [", action, 
#                       '|', lower_quantile, '|', upper_quantile, ']')
  
  
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
  if(n_orig > (all_cells %>% length))
  {
    x$input$dataset = x$input$dataset %>% filter(cell %in% all_cells)
    x$input$normalisation = x$input$normalisation %>% filter(cell %in% all_cells)
    
    x$log = paste0(x$log, '\n- ', 
                   Sys.time(), " Filtered missing data: [", proportion_RNA, 
                   '|', proportion_ATAC, ']')
    
  }
  
  
  
  
  return(x)
}