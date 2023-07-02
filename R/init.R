#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
#' @param lhs A value or the magrittr placeholder.
#' @param rhs A function call using the magrittr semantics.
#' @return The result of calling `rhs(lhs)`.
NULL


#' Create a dataset.
#'
#' @description
#'
#' This function creates a dataset (an object of class \code{rcongasplus}) by assembling multiple single-cell input measurements
#' (ATAC and/or RNA data modalities), the input segmentation (from bulk DNA sequencing),
#' and the per-cell normalisation factors for the data.
#'
#' All input data are passed as tibbles; the input formats are as follows:
#'
#' * for single-cell ATAC/RNA data, the \code{cell} identifier, the genomic coordinates
#' (\code{chr}, \code{from}, \code{to}) which refer either to an ATAC peak, or an RNA gene
#' identifier, and a \code{value} reporting the reads mapped.
#'
#' * for the input segmentation, the genomic coordinates
#' (\code{chr}, \code{from}, \code{to}) which refer to the segment, and the number of
#' \code{copies} (i.e., DNA ploidy) of the segment.
#'
#' * for normalization factors the \code{cell} identifier, the actual \code{normalisation_factor}
#' and the \code{modality} to wihch the factor refers to
#'
#' This function receives also other parameters - e.g., the models likelihoods - which
#' will determine the overall behaviour of the underlying model, and how data are preared for inference.
#'
#' * A Negative Binomial likelihood (\code{"NB"}), which works directly from raw counts data
#'
#' * A Gaussian likelihood (\code{"G"}), which requires a z-score transformation of the data. This consists
#' in :
#'     * scaling raw counts by the input normalization factors;
#'     * computing z-scores per cell;
#'     * summing up z-scores per segment;
#'     * computing z-scores per segment;
#'     * center the z-scores mean to the input ploidy.
#'
#' @param rna A tibble with single-cell RNA data.
#' @param atac A tibble with single-cell ATAC data.
#' @param segmentation A tibble with the input segmentation.
#' @param rna_normalisation_factors The RNA tibble with the input per-cell normalisation factors.
#' By default these are computed by function \code{auto_normalisation_factor}.
#' @param atac_normalisation_factors The ATAC tibble with the input per-cell normalisation factors.
#' By default these are computed by function \code{auto_normalisation_factor}.
#' @param rna_likelihood Type of likelihood used for RNA data (\code{"G"} for Gaussian and
#' \code{""NB} for Negative Binomial). The RNA default is \code{"G"}.
#' @param atac_likelihood Type of likelihood used for ATAC data, with default \code{"NB"}.
#' @param reference_genome Either \code{"GRCh38"} or \code{"hg19"}.
#' @param description A model in-words description.
#' @param smooth If yes, input segments are smootheed by joining per chromosome segments that
#' have the same ploidy.
#'
#' @return An object of class \code{rcongasplus}
#'
#' @importFrom tidyr separate pivot_longer
#' @importFrom tibble rownames_to_column
#' @importFrom crayon bgCyan bgBlue bgGreen bgRed bgMagenta bgWhite underline bgYellow blue red
#' @importFrom cli cli_rule cli_h3 cli_alert cli_alert_info cli_alert_warning cli_alert_danger
#' @importFrom reshape2 acast melt
#' @importFrom stats4 mle
#' @importFrom cowplot plot_grid
#' @importFrom progress progress_bar
#' @importFrom tidyr unite
#' @importFrom graphics curve hist
#' @importFrom stats complete.cases dgamma quantile rnorm sd
#' @importFrom utils head object.size
#' @importFrom gtools mixedsort
#' @import dplyr
#' @import ggplot2
#' @import CNAqc
#'
#' @export
#'
#' @examples
#' data("example_input")
#'
#' # For instance, RNA data
#' example_input$x_rna %>% print
#'
#' # .. or ATAC data
#' example_input$x_atac %>% print
#'
#' # .. and segmentation
#' example_input$x_segmentation %>% print
#'
#' # .. and normalisation factors can be computed (default)
#' example_input$x_rna %>% auto_normalisation_factor()
#'
#' x = init(
#'   rna = example_input$x_rna,
#'   atac = example_input$x_atac,
#'   segmentation = example_input$x_segmentation,
#'   rna_likelihood = "G",
#'   atac_likelihood = 'NB',
#'   description = 'My model')
#'
#' print(x)
init = function(
  rna,
  atac,
  segmentation,
  rna_normalisation_factors = rna %>% auto_normalisation_factor(),
  atac_normalisation_factors = atac %>% auto_normalisation_factor(),
  rna_likelihood = "NB",
  atac_likelihood = "NB",
  reference_genome = 'GRCh38',
  description = "(R)CONGAS+ model",
  smooth = FALSE,
  out.rm = T
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
  rna = rna %>% filter(!is.na(chr))
  sanitize_input(rna,
                 required_input_columns = c("chr", "from", "to", "value", "cell"),
                 types_required = c("character", "integer", "integer", "integer", "character")
  )

  sanitize_input(atac,
                 required_input_columns = c("chr", "from", "to", "value", "cell"),
                 types_required = c("character", "integer", "integer", "integer", "character")
  )

  sanitize_input(segmentation,
                 required_input_columns = c("chr", "from", "to", "copies"),
                 types_required = c("character", "integer", "integer", "integer")
  )

  sanitize_input(rna_normalisation_factors,
                 required_input_columns = c("cell", "normalisation_factor"),
                 types_required = c("character", "numeric")
  )
  
  sanitize_input(atac_normalisation_factors,
                 required_input_columns = c("cell", "normalisation_factor"),
                 types_required = c("character", "numeric")
  )

  # Reference
  if(!(reference_genome %in% c("hg19", 'GRCh37', 'hg38', "GRCh38")))
    stop("Unsupported reference, use any of 'hg19'/'GRCh37' or 'hg38'/'GRCh38'")

  # Non unique cell ids
  nu_ids = intersect(rna$cell, atac$cell)

  if(!is.null(nu_ids) & (length(nu_ids) > 0))
  {
    stop("ATAC and RNA cells have shared ids, this is not possibile.")
  }

  # Check that normalization factors are available for all cells
  if(!all(rna$cell %in% rna_normalisation_factors$cell))
  {
    message("Error with this RNA tibble")
    rna_normalisation_factors %>% print()
    
    stop("Missing normalisation factors for some input cells.")
  }
  
  if(!all(atac$cell %in% atac_normalisation_factors$cell))
  {
    message("Error with this ATAC tibble")
    atac_normalisation_factors %>% print()

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
  if(smooth)
    segmentation = smoothie(segments = segmentation, reference = reference_genome)
  
  segmentation = segmentation %>% idify()

  segmentation$RNA_genes =
    segmentation$RNA_nonzerovals =
    segmentation$ATAC_peaks =
    segmentation$ATAC_nonzerovals = 
    segmentation$ATAC_nonzerocells =
    segmentation$RNA_nonzerocells = 0
  

  # Create RNA modality data
  rna_modality_data = create_modality(
    modality = "RNA",
    data = rna,
    segmentation = segmentation,
    normalisation_factors = rna_normalisation_factors,
    likelihood = rna_likelihood,
    out.rm = out.rm)

  if(!is.null(rna))
  {
    rna = rna_modality_data$data %>% select(segment_id, cell, value, modality, value_type)
    segmentation = rna_modality_data$segmentation
    rna_normalisation_factors = rna_modality_data$normalisation
  }

  # Create ATAC modality data
  atac_modality_data = create_modality(
    modality = "ATAC",
    data = atac,
    segmentation = segmentation,
    normalisation_factors = atac_normalisation_factors,
    likelihood = atac_likelihood,
    out.rm = out.rm)

  if(!is.null(atac))
  {
    atac = atac_modality_data$data %>% select(segment_id, cell, value, modality, value_type)
    segmentation = atac_modality_data$segmentation
    
    atac_normalisation_factors = atac_modality_data$normalisation
  }
  
  # The last thing we do is to sanitise the segmentation, this removes useless segments
  cli::cli_h3("Checking segmentation")
  cat("\n")

  # First, test if some segments have NO events associated
  s_subset_zero = segmentation %>% filter(RNA_nonzerovals == 0 & ATAC_nonzerovals == 0)
  segmentation = segmentation %>% filter(RNA_nonzerovals > 0 | ATAC_nonzerovals > 0)

  if (nrow(s_subset_zero) > 0)
  {
    cli::cli_alert_warning("These segments have no events associated and will be removed - we suggest you to check if these can be further smoothed.")
    cli::cli_alert_info("{crayon::bgCyan(crayon::white(' Hint '))} {crayon::italic('Check if the reduced segments can be further smoothed!')}")
    s_subset_zero %>% print
  }
  else cli::cli_alert_success("Nothing to process.")

  # Add to the return object
  ret_obj$input$dataset = bind_rows(rna, atac)
  ret_obj$input$normalisation = bind_rows(rna_normalisation_factors,
                                          atac_normalisation_factors)
  ret_obj$input$segmentation = segmentation

  ret_obj$log = paste('-', Sys.time(), "Created input object.")
  
  return(ret_obj %>% sanitize_obj %>% sanitize_zeroes
  )
}

create_modality = function(modality, data, segmentation, normalisation_factors, likelihood, out.rm=T)
{
  # Special case, data are missing
  if(is.null(data)) {
    return(list(data = NULL, segmentation = segmentation))
  }
  
  normalisation_factors = normalisation_factors %>% mutate(modality = !!modality)

  os = object.size(data)/1e6

  # Report input information for RNA
  cli::cli_h3("{.field {modality}} modality ({.value {os} Mb})")
  cat("\n")

  cli::cli_alert("Input events: {.field {data %>% nrow}}")
  cli::cli_alert("Input cells: {.field {data %>% distinct(cell) %>% nrow}}")
  cli::cli_alert("Input locations: {.field {data %>% distinct(chr, from, to) %>% nrow}}")

  # Computing mapping for RNA
  segmentation = segmentation %>% idify()

  evt_lbs = paste0(modality, '_nonzerovals')
  loc_lbs = ifelse(
    modality == 'RNA',
    paste0(modality, '_genes'),
    paste0(modality, '_peaks')
  )
  cells_lbs = paste0(modality, '_nonzerocells')

  segmentation[[evt_lbs]] = 0
  segmentation[[loc_lbs]] = 0
  segmentation[[cells_lbs]] = 0

  data$segment_id = NA

  pb <- progress::progress_bar$new(total = nrow(segmentation))

  pb$tick(0)

  for(i in 1:nrow(segmentation))
  {
    pb$tick()

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

    segmentation[[cells_lbs]][i] = data[what_maps, ] %>%
      pull(cell) %>% unique %>% length
  }

  n_na = is.na(data$segment_id) %>% sum()
  nn_na = (data %>% nrow) - n_na
  
  if (out.rm)
    data = clean_outliers_persegment(modality, data, normalisation_factors)$data_cleaned
  
  cli::cli_alert("Entries mapped: {.field {nn_na}}, with {.field {n_na}} outside input segments that will be discarded.")
  if(n_na > 0) data = data %>% filter(!is.na(segment_id))

  cli::cli_alert("Using likelihood: {.field {likelihood}}.")

  # Compute z-score per mapped gene/peak according to the required likelihood
  if(likelihood %in% c("G"))
  {
    cli::cli_alert_warning("Gaussian likelihood requires z-score representation of input values.")

    # Normalise by factor - divide counts by per-cell normalization_factor
    data = normalise_modality(data %>% mutate(modality = !!modality), normalisation_factors)

    if ('gene' %in% colnames(data)){
      # data  = data %>% mutate(idFeature = paste0(gene,chr,from,to))
      features = data %>% 
        select(gene, chr, from, to)  %>% 
        distinct() %>%
        unite(idFeature, c(gene, chr, from, to), remove = F) 
      
      data = data %>% left_join(features)
    }
    else{
     # data = data %>% mutate(idFeature = paste0(chr,from,to))
      features = data %>% 
        select(chr, from, to)  %>% 
        distinct() %>%
        unite(idFeature, c(chr, from, to), remove = F) 
      
      data = data %>% left_join(features)
    }
      
    # Compute z-score
    cli::cli_alert("Computing z-score.")

    gene_segments = data %>% group_by(idFeature, segment_id) %>% summarise(n_cells = n()) %>%
    select(-n_cells) %>% column_to_rownames('idFeature')

    inp = reshape2::acast(data,
                          cell ~ idFeature,
                          value.var = "value")
    inp[is.na(inp)] <- 0

    data_scaled = scale(inp)
    data_scaled[is.na(data_scaled)] <- 0
    mapped = sapply(segmentation$segment_id, function (x)  {
      curr_genes = rownames(gene_segments %>% filter(segment_id == x))
      tmp = data_scaled[,curr_genes,drop=F] 
      tmp = rowSums(tmp)
      return(tmp)
    }, USE.NAMES = T )

    mapped = as.data.frame(mapped) %>% rownames_to_column(var = 'cell') %>%
      pivot_longer(cols = setdiff(colnames(mapped), 'cell'), names_to = 'segment_id')

    mapped$modality = modality

    # # zscore_params = data %>%
    # #   group_by(chr, from, to) %>%
    # #   summarise(value_mean = mean(value), value_sd = sd(value), .groups = 'keep')
    # # NEW CODE FOR Z-SCORES

    # data_ = data %>% select(cell,value, idFeature) %>% 
    #   pivot_wider(values_from = value, names_from = idFeature)  %>%
    #   column_to_rownames('cell') %>% 
    #   replace(is.na(.), 0)

    # means = colMeans(data_)
    # sds  = apply(data_, 2, sd)

    # zscore_params = tibble(value_mean = means, value_sd = sds, idFeature = names(means))

    # # Set to 1 normalisation factors
    cli::cli_alert_warning("With Gaussian likelihood normalisation factors are changed to 1.")
    normalisation_factors$normalisation_factor = 1
    
    # # data = data %>%
    # #   left_join(zscore_params, by = c('chr', 'from', 'to')) %>%
    # #   mutate(
    # #     value = (value - value_mean)/value_sd # z-score
    # #   )
    # data = data %>%
    #   left_join(zscore_params) %>%
    #   mutate(
    #     value = (value - value_mean)/value_sd # z-score
    #   )

    # n_na = data$value %>% is.na %>% sum

    # if(n_na > 0)
    # {
    #   cli::cli_alert_warning("There are {.field {n_na}} z-scores that are NA, will be removed.")
    #   data %>% filter(is.na(value)) %>% print

    #   data = data %>% filter(!is.na(value))
    # }
    # min_value = min(data$value)
    # data = data %>% mutate(value = value + (min_value * -1))
    # # data = data %>% select(-value_mean, -value_sd, -idFeature)
    
  } else {
  # Mapped counts
  mapped = data %>%
    group_by(segment_id, cell) %>%
    summarise(value = sum(value), .groups = 'keep') %>%
    ungroup() %>%
    mutate(modality = !!modality)
  }

  # Center the new scores to the ploidy value
  if(likelihood %in% c("G")){
    cli::cli_alert("Centering the new scores around input ploidy values.")

    # data_ = mapped %>% select(cell,value, segment_id) %>% 
    #   pivot_wider(values_from = value, names_from = segment_id)  %>%
    #   column_to_rownames('cell') %>% 
    #   replace(is.na(.), 0)

    # means = colMeans(data_)
    # sds  = apply(data_, 2, sd)

    inp = reshape2::acast(mapped,
                          cell ~ segment_id,
                          value.var = "value")
    inp[is.na(inp)] <- 0

    data_scaled = scale(inp)
    data_scaled[is.na(data_scaled)] <- 0

    #zscore_params = tibble(value_mean = means, value_sd = sds, segment_id = names(means))

    # zscore_params = mapped %>%
    #   group_by(segment_id) %>%
    #   summarise(value_mean = mean(value), value_sd = sd(value), .groups = 'keep')

    # mapped = mapped %>%
    #   left_join(zscore_params, by = c("segment_id")) %>%
    #   mutate(
    #     value = (value - value_mean)/value_sd # z-score
    #   )
    mapped = as.data.frame(data_scaled) %>% rownames_to_column(var = 'cell') %>%
      pivot_longer(cols = setdiff(colnames(data_scaled), 'cell'), names_to = 'segment_id') %>%
      mutate(modality = !!modality)

    print(colnames(mapped))

    mapped =  mapped %>%
      left_join(segmentation, by = c("segment_id")) %>%
      mutate(
        value = value + copies
      ) %>%  select(segment_id, cell, value, modality)
  }

  # Handle data type request: convert to z-score if required
  mapped$value_type = likelihood

  what_lik = case_when(
    likelihood == "NB" ~ "Negative Binomial",
    likelihood == "P" ~ "Poisson",
    likelihood == "G" ~ "Gaussian"
  )

  return(
    list(
      data = mapped,
      normalisation = normalisation_factors,
      segmentation = segmentation
    )
  )
}

# Segments smoothing function
smoothie = function(segments, reference = 'GRCh38')
{
  cli::cli_alert("Smoothing {.field {nrow(segments)}} input segments.")
  reference = CNAqc:::get_reference(ref = reference)
  
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
  
  cli::cli_alert("After smoothing there are {.field {nrow(smoothed_segments)}} input segments.")
  
  return(smoothed_segments)
}
