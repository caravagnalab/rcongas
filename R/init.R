init = function(
  rna,
  atac,
  segmentation,
  normalisation_factors,
  rna_likelihood = "G",
  atac_likelihood = "NB",
  reference_genome = 'GRCh38',
  description = "CONGAS+ model"
)
{
  if(is.null(rna) & is.null(atac))
    stop("Cannot have both assays null.")
  
  # Output object
  ret_obj = list()
  class(ret_obj) = 'rcongasplus'
  
  ret_obj$description = description
  ret_obj$reference_genome = reference_genome
  
  cli::boxx(paste("(R)CONGAS+:", description), 
            background_col = 'orange',
            padding = 1,
            float = 'center') %>% cat
  cat("\n")

  # Sanitise by data type and required columns
  sanitize_input(rna, 
                 required_input_columns = c("chr", "from", "to", "count", "cell"),
                 types_required = c("character", "integer", "integer", "integer", "character")
  )
  
  sanitize_input(atac, 
                 required_input_columns = c("chr", "from", "to", "count", "cell"),
                 types_required = c("character", "integer", "integer", "integer", "character")
  )
  
  sanitize_input(segmentation, 
                 required_input_columns = c("chr", "from", "to", "copies"),
                 types_required = c("character", "integer", "integer", "integer")
  )
  
  sanitize_input(normalisation_factors, 
                 required_input_columns = c("cell", "normalisation_factor", "modality"),
                 types_required = c("character", "numeric", "character")
  )
  
  # Reference
  if(!(reference_genome %in% c("hg19", 'GRCh37', 'hg38', "GRCh38")))
    stop("Unsupported reference, use any of 'hg19'/'GRCh37' or 'hg38'/'GRCh38'")
  
  # Non unique cell ids
  nu_ids = intersect(rna$cell, atac$cell)
  
  if(!is.null(nu_ids) & (length(nu_ids) > 0))
  {
    stop("ATAC and RNA cells have shared ids.")
  }
  
  # Check that normalisation factors are available for all cells
  all_sample_cells = c(rna$cell, atac$cell) %>% unique
  
  if(!all(all_sample_cells %in% normalisation_factors$cell))
  {
    message("Error with this tibble")
    normalisation_factors %>% print()
    
    stop("Missing normalisation factors for some input cells.")
  }
  
  # Check likelihood
  if(!(rna_likelihood %in% c("NB", "P", "G")))
  {
    stop("Unsupported RNA likelihood, use:
      - NB (Negative-Binomial), 
      - P (Poisson),
      - G (Gaussian, with z-score).")
  }
  
  if(!(atac_likelihood %in% c("NB", "P", "G")))
  {
    stop("Unsupported RNA likelihood, use:
      - NB (Negative-Binomial), 
      - P (Poisson),
      - G (Gaussian, with z-score).")
  }
  
  # Prepare segment
  segmentation = segmentation %>% idify()
  
  segmentation$RNA_genes = 
    segmentation$RNA_events =   
    segmentation$ATAC_peaks = 
    segmentation$ATAC_events = 0
    
  # Create RNA modality data
  rna_modality_data = create_modality(
    modality = "RNA",
    data = rna,
    segmentation = segmentation,
    normalisation_factors = normalisation_factors %>% filter(modality == "RNA"),
    likelihood = rna_likelihood)
  
  if(!is.null(rna))
  {
    rna = rna_modality_data$data %>% select(segment_id, cell, value, modality, value_type)
    segmentation = rna_modality_data$segmentation
  }
  
  # Create ATAC modality data
  atac_modality_data = create_modality(
    modality = "ATAC",
    data = atac,
    segmentation = segmentation,
    normalisation_factors = normalisation_factors %>% filter(modality == "ATAC"),
    likelihood = atac_likelihood)
  
  if(!is.null(atac))
  {
    atac = atac_modality_data$data %>% select(segment_id, cell, value, modality, value_type)
    segmentation = atac_modality_data$segmentation
  }

  # Add to the return object
  ret_obj$input$dataset = bind_rows(rna, atac)
  ret_obj$input$normalisation = normalisation_factors
  ret_obj$input$segmentation = segmentation
  
  return(ret_obj %>% sanitize_obj)
}

create_modality = function(modality, data, segmentation, normalisation_factors, likelihood)
{
  # Special case, data are missing
  if(is.null(data)) {
    
    return(list(data = NULL, segmentation = segmentation))
  }
  
  # # Report input information for RNA
  cli::cli_h3("{.value {modality}} modality")
  cat("\n")
  
  cli::cli_alert("Input events: {.field {data %>% nrow}}")
  cli::cli_alert("Input cells: {.field {data %>% distinct(cell) %>% nrow}}")
  cli::cli_alert("Input locations: {.field {data %>% distinct(chr, from, to) %>% nrow}}")
  
  # Computing mapping for RNA 
  segmentation = segmentation %>% idify()
  
  evt_lbs = paste0(modality, '_events')
  loc_lbs = ifelse(
    modality == 'RNA',
    paste0(modality, '_genes'),
    paste0(modality, '_peaks')
  )
  
  segmentation[[evt_lbs]] = 0
  segmentation[[loc_lbs]] = 0
  
  data$segment_id = NA

  for(i in 1:nrow(segmentation))
  {
    what_maps = which(
      data$chr == segmentation$chr[i] &
        data$from >= segmentation$from[i] &
        data$to <= segmentation$to[i]
    )
    
    if(length(what_maps) == 0) next;
    
    data$segment_id[what_maps] = segmentation$segment_id[i]
    
    segmentation[[evt_lbs]][i] = what_maps %>% length
    segmentation[[loc_lbs]][i] = data[what_maps, ] %>% 
      distinct(chr, from, to) %>% 
      nrow()
  }
  
  n_na = is.na(data$segment_id) %>% sum()
  nn_na = (data %>% nrow) - n_na
  
  cli::cli_alert("Entries mapped: {.field {nn_na}}, with {.field {n_na}} outside input segments.")
  
  # Mapped RNA counts 
  mapped = data %>% 
    group_by(segment_id, cell) %>% 
    summarise(count = sum(count), .groups = 'keep') %>% 
    ungroup() %>% 
    mutate(modality = !!modality) %>% 
    rename(value = count)
  
  # Handle data type request: convert to z-score if required
  mapped$value_type = likelihood
  
  what_lik = case_when(
    likelihood == "NB" ~ "Negative Binomial",
    likelihood == "P" ~ "Poisson",
    likelihood == "G" ~ "Gaussian"
  )
  
  cli::cli_alert("Using RNA likelihood: {.field {what_lik}}.")
  
  if(likelihood %in% c("G"))
  {
    cli::cli_alert("Required z-score representation, computing after normalising data with input factors.")
    
    # Normalise by factor - divide counts by per-cell normalization_factor
    mapped = normalise_modality(mapped, normalisation_factors)
    
    zscore_params = mapped %>% 
      group_by(segment_id) %>% 
      summarise(counts_mean = mean(value), counts_sd = sd(value), .groups = 'keep')
    
    mapped = mapped %>% 
      left_join(zscore_params, by = 'segment_id') %>% 
      mutate(
        value = (value - counts_mean)/counts_sd # z-score
      ) 
    
    n_na = mapped$value %>% is.na %>% sum
    
    if(n_na > 0)
    {
      cli::cli_alert_warning("There are {.field {n_na}} z-scores that are NA, will be removed.")
      mapped %>% filter(is.na(value)) %>% print
      
      mapped = mapped %>% filter(!is.na(value))
    }
  }
  
  return(
    list(
      data = mapped,
      segmentation = segmentation
    )
  )
}

