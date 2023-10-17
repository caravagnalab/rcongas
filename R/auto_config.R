#' Determine data-based hyperparameters for the Bayesian priors.
#'
#' @description This function determines the hyperparameters for the Bayesian priors based on the input data.
#'
#' @param x (required) CONGAS+ object.
#' @param K (required) Number of clusters that will be tested during the inference.
#' @param NB_size_atac Float (optional). Default 150. Value used to inizialize the size hyperparametr of the Negative Binomial for ATAC in case the likelihood for ATAC is set to NB
#' @param NB_size_rna Float (optional). Default 150. Value used to inizialize the size hyperparametr of the Negative Binomial for RNA in case the likelihood for RNA is set to NB
#' @param a_sd Float (optional). Default 0.1. Lower bound of the Uniform prior for the Gussian standard deviation. Used when one of RNA or ATAC likelihoods are gaussian.  
#' @param b_sd Float (optional). Default 1. Upper bound of the Uniform prior for the Gussian standard deviation. Used when one of RNA or ATAC likelihoods are gaussian.  
#' @param NB_size_priors (optional) Default c(15, 1000). Lower and upper bound of the Uniform prior on the Negative binomial size hyperparameter.
#' @param hidden_dim (Optional) defualt 5. Number of discrete copy number states to model. By default this is from 1 to 5.
#' @param prior_cn (Optional) Default c(0.1, 0.6, 0.1, 0.1, 0.1). Prior for the copy number state of every segment.
#' @param init_importance (Optional) Default 0.6. Value used to initialize the distribution over possible copy number states for every cluster and every segment. \code{init_importance}
#' is used to initialize the value corresponding to the copy number state inferred from bulk and \code{1-init_importance / (hidden_dim - 1)} is used to initialize the other ploidy states.
#' This distribution is initialized based on the value of init importance, and then its prior is defined in \code{prior_cn}.
#' @param purity Optional (default NULL). This hyperparameter can be used to inject prior knowledge about the purity of the sample. It can be set as the purity inferred from bulk Whole Genome Sequencing. Leave this to NULL
#' in case the purity of the sample is expected to be 100%.
#' @param multiome Default to FALSE. Flag indicating whether the RNA and ATAC observations are the result of a matched RNA-ATAC sequencing assay such as 10x multiome assay.
#' @param CUDA Defualt FALSE. Flag indicating whether to use GPU computation.
#' @param normal_cells Default to FALSE. Flag that can be used to inject prior knowledge about the presence of normal cells in the sample. In case this is set to TRUE, the copy number
#' distribution for one of the clusters will be initialized with values skewed the diploid state in every segment.
#'
#' @return CONAGS+ object
#' @export
#'
auto_config_run <-
  function(x,
           K = c(1:3),
           NB_size_atac = 150,
           NB_size_rna = 150,
           # lambda = 0.3,
           a_sd = 0.1,
           b_sd = 1,
           prior_cn = c(0.1, 0.6, 0.1, 0.1, 0.1), #, 0.025),
           hidden_dim = 5,#length(prior_cn),
           init_importance = 0.6,
           NB_size_priors = c(15, 1000),
           purity = NULL,
	         multiome = FALSE, 
           CUDA = FALSE,
           normal_cells = FALSE
           )

    {
    cli::cli_h1("(R)CONGAS+ hyperparameters auto-config")
    param_list <-  list()

    torch <-  reticulate::import("torch")
    
    if(CUDA){
      torch$set_default_tensor_type('torch.cuda.FloatTensor')
    } else{
      torch$set_default_tensor_type('torch.FloatTensor')
    }

    param_list$init_probs <- init_importance
    hidden_dim = if (! is.character(prior_cn)) length(prior_cn) else hidden_dim

    mode = if (is.character(prior_cn)) 'segment_specific' else 'same'
    param_list$probs = torch$tensor(prior_cn)
    #init_dirichlet_cna(x, hidden_dim = hidden_dim, mode = mode)
    

    param_list$purity <- purity
    param_list$multiome <- multiome
    param_list$normal_cells = normal_cells

    if (has_atac(x)) {

      cli::cli_h2("ATAC modality")

      params_atac <- list()

      if (startsWith(which_likelihood(x, "ATAC"), prefix = "G"))
      {
        cli::cli_alert("Gaussian likelihood")

        sds <-
          get_data(x) %>%
          filter(modality == "ATAC") %>%
          group_by(segment_id) %>%
          summarize(sd = sd(value)) %>%
          pull(sd)

        params_atac$norm_init_sd_atac <- torch$tensor(sds / max(K))
        params_atac$likelihood_atac <-  "G"
        I <-  length(get_segmentation(x) %>%  pull(segment_id) %>%  unique())
        params_atac$theta_shape_atac <-  torch$tensor(rep(1, I))
        params_atac$theta_rate_atac <-  torch$tensor(rep(1,I))

      } else if (startsWith(which_likelihood(x, "ATAC"), prefix = "NB"))
        {

        cli::cli_alert("Negative Binomial likelihood, estimating Gamma shape and rate")

        # cat('\n')
        # cli::cli_process_start("Estimating Gamma shape and rate")
        # cat('\n')

        params_atac <- gamma_shape_rate(x, modality = "ATAC")
        names(params_atac) <-
          c("theta_shape_atac", "theta_rate_atac")

        # cat('\n')
        # cli::cli_process_done()
        # cat('\n')

        params_atac$nb_size_init_atac <-
          torch$tensor(rep(
            NB_size_atac,
            length(get_segmentation(x) %>%  pull(segment_id) %>%  unique())
          ))

        params_atac$likelihood_atac <-  "NB"

      } else {
        params_atac <- gamma_shape_rate(x, modality = "ATAC", torch = torch)
        names(params_atac) <-
          c("theta_shape_atac", "theta_rate_atac")
        params_atac$likelihood_atac <-  "P"
      }

      param_list <-  c(param_list, params_atac)

    }

    if (has_rna(x)) {
      params_rna <- list()

      cli::cli_h2("RNA modality")


      if (startsWith(which_likelihood(x, "RNA"), prefix = "G")) {

        cli::cli_alert("Gaussian likelihood")

        sds <-
          get_data(x) %>%  filter(modality == "RNA") %>% group_by(segment_id) %>%  summarize(sd = sd(value)) %>%  pull(sd)
        params_rna$norm_init_sd_rna <- torch$tensor(sds / max(K))
        params_rna$likelihood_rna <-  "G"
        I <-  length(get_segmentation(x) %>%  pull(segment_id) %>%  unique())
        params_rna$theta_shape_rna <-  torch$tensor(rep(1, I))
        params_rna$theta_rate_rna <-  torch$tensor(rep(1,I))

      } else if (startsWith(which_likelihood(x, "RNA"), prefix = "NB"))
      {
        cli::cli_alert("Negative Binomial likelihood, estimating Gamma shape and rate")

        params_rna <- gamma_shape_rate(x, modality = "RNA")
        names(params_rna) <-  c("theta_shape_rna", "theta_rate_rna")

        params_rna$nb_size_init_rna <-
          torch$tensor(rep(
            NB_size_rna,
            length(get_segmentation(x) %>%  pull(segment_id) %>%  unique())
          ))

        params_rna$likelihood_rna <-  "NB"
      } else {
        params_rna <- gamma_shape_rate(x, modality = "RNA", torch = torch)
        names(params_rna) <-  c("theta_shape_rna", "theta_rate_rna")
        params_rna$likelihood_rna <-  "P"
      }

      param_list <-  c(param_list, params_rna)

    }

    param_list$a <-  a_sd
    param_list$b <- b_sd

    # param_list$lambda <- lambda
    param_list$hidden_dim <- as.integer(hidden_dim)
    param_list$binom_prior_limits <- NB_size_priors


    return(param_list)

  }

gamma_shape_rate <-
  function(x,
           modality = "RNA",
           torch = reticulate::import("torch"))
    {

    inp = reshape2::acast(Rcongas:::get_data(x) %>% filter(modality == !!modality),
                          cell ~ segment_id,
                          value.var = "value")
    inp[is.na(inp)] <- 0

    inp = inp[order(rownames(inp)), order(colnames(inp)), drop = FALSE]

    norm_raw = Rcongas:::get_normalisation(x) %>%
      filter(modality == !!modality) %>%
      select(normalisation_factor, cell)

    norm = norm_raw$normalisation_factor
    names(norm) = norm_raw$cell
    norm = norm[rownames(inp)]

    ploidy <- x$input$segmentation$copies

    names(ploidy) <- x$input$segmentation$segment_id

    ploidy <-  ploidy[colnames(inp)]

    theta_factors = Rcongas:::estimate_segment_factors(data = inp, norm_factors = norm, pld = ploidy, plot = F)

    theta_shape = torch$tensor(sapply(theta_factors, function(x)
      x[1]))
    theta_rate = torch$tensor(sapply(theta_factors, function(x)
      x[2]))

    ret = list(theta_shape, theta_rate)

    return(ret)
  }


#' Compute empirical normalisation factors (library size)
#'
#' @description For a tibble dataset with a \code{cell} and \code{value}
#' columns, it sums up all values per cell (total library size), and scales
#' the value for \code{10^x} where \code{x} is the median nuimber of digits
#' in each total value size.
#'
#' For instance for numbers of the order of ~1000 it will rescale the library
#' size by a factor 1000, taking all values around ~1.
#'
#' @param x A tibble with a \code{cell} and \code{value} columns.
#'
#' @return A tibble, aggregated by cell with total values (sum), rescaled
#' by the estimated constant.
#'
#' @export
#'
#' @examples
#' data('example_input')
#' example_input$x_rna %>% auto_normalisation_factor
auto_normalisation_factor = function(x)
{
  if(is.null(x)) return(NULL)
  if(!("cell" %in% colnames(x))) stop("Missing cell column.")
  if(!("value" %in% colnames(x))) stop("Missing value column.")

  cli::cli_alert("Computing library size factors as total counts.")

  norm_x = x %>%
    group_by(cell) %>%
    summarise(normalisation_factor = sum(value, na.rm = TRUE))

  ndigits = median(nchar(round(norm_x$normalisation_factor)))
  scaling = `^`(10,ndigits)

  cli::cli_alert("Median digits per factor {.field {ndigits}}, scaling by {.field {scaling}}.")

  # norm_x$normalisation_factor %>% summary %>% print()
  norm_x$normalisation_factor = norm_x$normalisation_factor / scaling
  # norm_x$normalisation_factor %>% summary
  norm_x$normalisation_factor %>% summary %>% print

  norm_x %>% return
}


init_dirichlet_cna = function(x, hidden_dim, mode = 'same') {
  if (mode == 'same') {
    segs = get_input(x, what = 'segmentation')  %>% arrange(segment_id)
    dirichlet_prior = lapply(segs$copies, function(x) {
      dir_conc = rep((1-init_importance) / (hidden_dim-1), hidden_dim)
      dir_conc[x] = init_importance
      return(dir_conc)})

    names(dirichlet_prior) = segs$segment_id
    dirichlet_prior = torch$tensor(matrix(unlist(dirichlet_prior), nrow = length(dirichlet_prior), byrow = T))
    param_list$probs = dirichlet_prior
  } else {
    # Use one vector for all segments.
    param_list$probs <-  torch$tensor(prior_cn)

  }
}


