#' Determine data-based hyperparameters for the Bayesian priors.
#'
#' @description
#'
#' @param x
#' @param K
#' @param NB_size_atac
#' @param NB_size_rna
#' @param lambda
#' @param a_sd
#' @param b_sd
#' @param prior_cn
#' @param purity
#'
#' @return
#' @export
#'
#' @examples
auto_config_run <-
  function(x,
           K,
           NB_size_atac = 150,
           NB_size_rna = 150,
           # lambda = 0.3,
           a_sd = 0.1,
           b_sd = 100,
           prior_cn = c(0.2, 0.6, 0.1, 0.05, 0.05), #, 0.025),
           hidden_dim = 5,#length(prior_cn),
           init_importance = 0.6,
           NB_size_priors = c(15, 1000),
           purity = 1
           )

    {
    cli::cli_h1("(R)CONGAS+ hyperparameters auto-config")
    param_list <-  list()

    torch <-  reticulate::import("torch")

    #param_list$probs <-  torch$tensor(prior_cn)
    if (is.character(prior_cn)) {
      param_list$init_probs <- init_importance
      segs = get_input(x, what = 'segmentation')  %>% arrange(segment_id)
      dirichlet_prior = lapply(segs$copies, function(x) {
        dir_conc = rep((1-init_importance) / hidden_dim, hidden_dim)
        dir_conc[x] = init_importance
        return(dir_conc)})

      names(dirichlet_prior) = segs$segment_id
      dirichlet_prior = torch$tensor(matrix(unlist(dirichlet_prior), nrow = length(dirichlet_prior), byrow = T))
      param_list$probs = dirichlet_prior
    } else {
      hidden_dim = length(prior_cn)
      # Use one vector for all segments.
      param_list$probs <-  torch$tensor(prior_cn)
      param_list$init_probs <- init_importance
    }

    param_list$purity <- purity

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
