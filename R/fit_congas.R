#' Fit a list of CONGAS model
#'
#' This function is general interface for fitting a \href{https://github.com/Militeee/congas}{congas} model in R. The model briefly consist in a joint
#' mixture model over two modalities, currently scATAC and scRNA-seq. For more information about the theoretical fundations of the approach refer
#' to the vignette. This function performs modele selection over a specified number of clusters, using a specific information criterium (IC). ICs
#' and results for all the runs are, however, reported in the object.
#'
#' The functions assume a list of model hyperparameters. As the the model formulation isquite complex, and those hyperparameters are extremely difficult to
#' set by hand we suggest the usage of the function \code{\link[Rcongas::auto_config_run]{Rcongas::auto_config_run()}}
#'
#' @param x An \code{rcongasplus} object with the input dataset, constructed with \code{\link[Rcongas::init]{Rcongas::init}}.
#' @param K a vector of integers with the number of clusters we want to test
#' @param learning_rate a learning rate for the Adam optimizer
#' @param model_parameters a list with model hyperparameters. As errors coming from wrong hyperparameters initialization
#' are quite hard to troubleshoot is higly suggested to use \code{\link[Rcongas::auto_config_run]{Rcongas::auto_config_run()}} to generate
#' a template and eventually modify it.
#' @param latent_variables specify the nature of the latent variable modelling the copy number profile. Currently only "D" (discrete) is available
#' @param compile use JIT compiler for the Pyro ba
#' @param steps number of steps of optimization
#' @param model_selection information criteria to which perform the model selection (one of ICL, NLL, BIC, AIC)
#'
#' @return An object ot class \code{rcongasplus} with a slot \code{bset_fit} with the learned parameters for the selected model in tiblle format. A slot \code{runs}
#' with all the runs performed ordered by the selectde IC and a slot \code{model_selection} with all the information to perform model selection.
#' @export
#'
#' @examples
#' library(Rcongas)
#'\dontrun{
#' K <-  1:4
#' hyperparams <- auto_config_run(example_object, 1:4)
#'
#' fit <- fit_congas(example_object, K = 1:4,learning_rate = 0.05, model_parameters = hyperparams)
#' }
fit_congas <-
  function(x,
           K,
           model_parameters,
           learning_rate = 0.05,
           latent_variables = "D",
           compile = FALSE,
           steps = 500,
           model_selection = "ICL",
           temperature = 10) {
    if (!inherits(x, "rcongasplus")) {
      stop("Input object needs to be an rcongas instance!")
    }

    runs <-
      lapply(K, function(k)
        fit_congas_single_run(
          x,
          model_parameters,
          k,
          learning_rate,
          latent_variables,
          compile,
          steps,
          temperature
        ))
    names(runs) <-  paste(K)

    cli::cli_h3("Inference completed, choosing best model.")

    model_selection_df <-  lapply(runs, function(y)
      y$ICs)
    model_selection_df <-
      do.call(rbind, model_selection_df) %>%  as_tibble()

    if(model_selection_df %>% is.na %>% any)
    {
      cli::cli_alert_danger("NA found in model scores, aborting")
      model_selection_df[!complete.cases(model_selection_df), ] %>% print
      stop("Cannot select best model!")
    }

    model_selection_df$K <- K

    bms_idx <- order(model_selection_df %>%  pull(!!model_selection))

    runs <-  runs[bms_idx]

    best_fit  <-  format_best_model(x, runs[[1]])

    x$runs <-  runs
    x$best_fit <-  best_fit
    x$model_selection <-  model_selection_df
    x$used_IC <- model_selection

    return(x)

  }


fit_congas_single_run <-
  function(x,
           parameters,
           K,
           learning_rate,
           latent_variables,
           compile,
           steps,
           temperature) {

    cli::cli_h3("Fit with k = {.field {K}}.")

    cg <- reticulate::import("congas")
    cg_mod <- reticulate::import("congas.models")
    pyro_optim <- reticulate::import("pyro.optim")

    data <- input_data_from_rcongas(x)
    param_optimizer <-  list()
    parameters$K <-  as.integer(K)
    parameters$Temperature <-  temperature
    param_optimizer$lr <- learning_rate

    ### LOAD MODEL AND GUIDE ###
    if (latent_variables == "D") {
      model <- cg_mod$LatentCategorical
      model_string <-  "LatentCategorical"
      parameters$latent_type <- "D"
    } else if (latent_variables == "C"){

      model <- cg_mod$LatentCategorical
      model_string <-  "LatentCategorical"
      parameters$latent_type <- "C"

    }else {
      stop("Continous latent variable model not yet implemented.")
    }

    ### LOAD THE ELBO-LOSS ###

    elbo_p <- reticulate::import("pyro.infer")
    if (compile) {
      elbo <- elbo_p$TraceGraph_ELBO
      elbo_string <- "TraceGraph_ELBO"
    } else {
      elbo <-  elbo_p$JitTraceGraph_ELBO
      elbo_string <- "JitTraceGraph_ELBO"
    }
    int <- cg$Interface()

    int$set_model(model)
    int$set_optimizer(pyro_optim$ClippedAdam)
    int$set_loss(elbo)
    int$initialize_model(data)
    int$set_model_params(parameters)

    loss = int$run(steps = as.integer(steps), param_optimizer = param_optimizer)

    fit_params = int$learned_parameters()
    ICs = int$calculate_ICs()


    hyperparams = c(
      list(
        "model" = model_string,
        "loss" = elbo_string,
        "optimizer" = "ClippedAdam"
      ),
      parameters
    )

    ret <-  list(
      "inferred_params" = fit_params,
      "hyperparameters" = lapply(hyperparams, detensorize),
      "ICs" = ICs %>% as.data.frame(),
      "loss" = loss
    )
    class(ret) <- "congas"
    return(ret)


  }
