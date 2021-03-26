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
#' @param normalisation_factors A tibble with the input per-cell normalisation factors.
#' @param rna_likelihood Type of likelihood used for RNA data (\code{"G"} for Gaussian and
#' \code{""NB} for Negative Binomial). The RNA default is \code{"G"}.
#' @param atac_likelihood Type of likelihood used for ATAC data, with default \code{"NB"}.
#' @param reference_genome Either \code{"GRCh38"} or \code{"hg19"}.
#' @param description A model in-words description.
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
#' # .. and x_normalisation factors
#' example_input$x_normalisation_factors %>% print
#'
#' x = init(
#'   rna = example_input$x_rna,
#'   atac = example_input$x_atac,
#'   segmentation = example_input$x_segmentation,
#'   normalisation_factors = example_input$x_normalisation_factors,
#'   rna_likelihood = "G",
#'   atac_likelihood = 'NB',
#'   description = 'My crazy model')
#'
#' print(x)
init = function(
  rna,
  atac,
  segmentation,
  normalisation_factors,
  rna_likelihood = "G",
  atac_likelihood = "NB",
  reference_genome = 'GRCh38',
  description = "A (R)CONGAS+ model"
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
    stop("ATAC and RNA cells have shared ids, this is not possibile.")
  }

  # Check that normalization factors are available for all cells
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
    segmentation$RNA_nonzerovals =
    segmentation$ATAC_peaks =
    segmentation$ATAC_nonzerovals = 0

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

  segmentation[[evt_lbs]] = 0
  segmentation[[loc_lbs]] = 0

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
  }



  n_na = is.na(data$segment_id) %>% sum()
  nn_na = (data %>% nrow) - n_na

  cli::cli_alert("Entries mapped: {.field {nn_na}}, with {.field {n_na}} outside input segments that will be discarded.")
  if(n_na > 0) data = data %>% filter(!is.na(segment_id))

  cli::cli_alert("Using likelihood: {.field {likelihood}}.")

  # Compute z-score per mapped gene/peak according to the required likelihood
  if(likelihood %in% c("G"))
  {
    cli::cli_alert_warning("Gaussian likelihood requires z-score representation of input values.")

    # Normalise by factor - divide counts by per-cell normalization_factor
    data = normalise_modality(data %>% mutate(modality = !!modality), normalisation_factors)

    # Compute z-score
    cli::cli_alert("Computing z-score.")
    zscore_params = data %>%
      group_by(chr, from, to) %>%
      summarise(value_mean = mean(value), value_sd = sd(value), .groups = 'keep')

    data = data %>%
      left_join(zscore_params, by = c('chr', 'from', 'to')) %>%
      mutate(
        value = (value - value_mean)/value_sd # z-score
      )

    n_na = data$value %>% is.na %>% sum

    if(n_na > 0)
    {
      cli::cli_alert_warning("There are {.field {n_na}} z-scores that are NA, will be removed.")
      data %>% filter(is.na(value)) %>% print

      data = data %>% filter(!is.na(value))
    }

    data = data %>% select(-value_mean, -value_sd)
  }

  # Mapped counts
  mapped = data %>%
    group_by(segment_id, cell) %>%
    summarise(value = sum(value), .groups = 'keep') %>%
    ungroup() %>%
    mutate(modality = !!modality)

  # Center the new scores to the ploidy value
  if(likelihood %in% c("G")){
    cli::cli_alert("Centering the new scores around input ploidy values.")

    zscore_params = mapped %>%
      group_by(segment_id) %>%
      summarise(value_mean = mean(value), value_sd = sd(value), .groups = 'keep')

    mapped = mapped %>%
      left_join(zscore_params, by = c("segment_id")) %>%
      mutate(
        value = (value - value_mean)/value_sd # z-score
      )

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
      segmentation = segmentation
    )
  )
}

