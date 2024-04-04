#' Fit an (R)CONGAS+ model
#'
#' @description
#'
#' This function is a general interface for fitting a \href{https://github.com/Militeee/congas}{congas} Python model in R. The model briefly consist in a joint
#' mixture model over two modalities, currently scATAC and scRNA-seq. For more information about the theoretical fundations of the approach refer
#' to the vignette. This function performs modele selection over a specified number of clusters, using a specific information criterium (IC). ICs
#' and results for all the runs are, however, reported in the object.
#'
#' The functions assume a list of model hyperparameters. As the the model formulation isquite complex, and those hyperparameters are extremely difficult to
#' set by hand we suggest the usage of the function \code{\link[Rcongas::auto_config_run]{Rcongas::auto_config_run()}}
#'
#' @param x An \code{rcongasplus} object with the input dataset, constructed with \code{\link[Rcongas::init]{Rcongas::init}}.
#' @param K a vector of integers with the number of clusters we want to test
#' @param lambdas Float (Optional). Default 0.5. Value of the hyperparameter that controls the weight given to RNA and ATAC modalities during the inference. Values closer to 0 give more weight to
#' the ATAC likelihood, while values closer to 1 result in higher weight given to the RNA likelihood.

#' @param learning_rate a learning rate for the Adam optimizer
#' @param model_parameters a list with model hyperparameters. As errors coming from wrong hyperparameters initialization
#' are quite hard to troubleshoot is higly suggested to use \code{\link[Rcongas::auto_config_run]{Rcongas::auto_config_run()}} to generate
#' a template and eventually modify it.
#' @param latent_variables specify the nature of the latent variable modelling the copy number profile. Currently only "G" is available, 
#' @param CUDA use GPU if avilable for training
#' @param steps number of steps of optimization
#' @param samples Number of times a model is fit for each value of \code{K}.
#' @param model_selection information criteria to which perform the model selection (one of ICL, NLL, BIC, AIC)
#' @param same_mixing boolean that indicates whether to use the same mixing proportions for both RNA and ATAC or use different vectors for the two
#' modalities. Default is FALSE. When the data is multi-omic, the mixing proportions are forced to be the same. 
#' @param threshold. Float, default is \code{learning_rate * 0.1}. It corresponds to the threshold that determines the early stopping of the training procedure. When the difference between parameters in step t and step t+1 is
#' lower than this threshold for a number of steps equal to the parameter \code{patience} the inference is stopped.
#' @param patience Integer. Number of steps to wait before stopping the inference. See \code{threshold} for more details.
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
           lambdas,
           model_parameters,
           learning_rate = 0.01,
           latent_variables = "G",
           CUDA = FALSE,
           steps = 500,
           samples = 1,
           parallel = FALSE,
           model_selection = "ICL",
           temperature = 10,
           equal_variance = TRUE,
           threshold = learning_rate * 0.1,
           patience = 5,
           same_mixing = FALSE
           ) {

    if (!inherits(x, "rcongasplus")) {
      stop("Input object needs to be an rcongas instance!")
    }

    # Sanitizers obj and zeroes
    x %>% sanitize_obj()
    x %>% sanitize_zeroes()
    x = sort_multiome(x)
    model_parameters$equal_sizes_sd <- equal_variance

    if (x$input$multiome) {
      same_mixing = TRUE
    }
    cli::cli_h1("{crayon::bgYellow(' (R)CONGAS+ ')} Variational Inference")

    # Multi-run fit
    one_k = function(k_lambda)
    {
      k = k_lambda['k']
      l = k_lambda['l']
      cli::cli_h3("Fit with k = {.field {k}} and lambda = {.field {l}}.")

      lapply(1:samples, function(w)
        fit_congas_single_run(
          x,
          model_parameters,
          l,
          k,
          learning_rate,
          latent_variables,
          CUDA,
          steps,
          temperature,
          threshold,
          patience,
          same_mixing
        ))
    }

    TIME = as.POSIXct(Sys.time(), format = "%H:%M:%S")

    k_lambda_pairs = expand.grid(k = K, l = lambdas)

    runs <- easypar::run(
      FUN = one_k,
      PARAMS = apply(k_lambda_pairs, 1, list), #lapply(K, list),
      parallel = parallel,
      progress_bar = FALSE
    )

    # Unroll everything
    runs = Reduce(append, runs)

    # Report timing to screen
    TIME = difftime(as.POSIXct(Sys.time(), format = "%H:%M:%S"), TIME, units = "mins")

    cat('\n\n')
    cli::cli_h2("{crayon::bold('(R)CONGAS+ fits')} completed in {.field {prettyunits::pretty_dt(TIME)}}.")
    cat('\n')

    # names(runs) <-  paste(K)

    # runs <-
    #   lapply(K, function(k)
    #     fit_congas_single_run(
    #       x,
    #       model_parameters,
    #       k,
    #       learning_rate,
    #       latent_variables,
    #       compile,
    #       steps,
    #       temperature
    #     ))
    # names(runs) <-  paste(K)

    cat("\n")
    # cli::cli_h3("Inference completed, choosing best model.")

    model_selection_df <-  lapply(runs, function(y)
      y$ICs)
    model_selection_df <-
      do.call(rbind, model_selection_df) %>%  as_tibble()

    # TODO: sistemare il codice in modo che i cluster che scompaiono non vengano proprio ritornati.
    model_selection_df$hyperparameter_K = sapply(runs, function(w) w$hyperparameters$K)
    model_selection_df$K = sapply(runs, function(w) {
      ass_atac = Rcongas:::detensorize(w$inferred_params$assignment_atac, CUDA)
      ass_rna = Rcongas:::detensorize(w$inferred_params$assignment_rna, CUDA)
      return(length(unique(c(ass_atac, ass_rna))))
    })
    model_selection_df$K_rna = sapply(runs, function(w) {
      ass_rna = Rcongas:::detensorize(w$inferred_params$assignment_rna, CUDA)
      return(length(unique(ass_rna)))
    })
    model_selection_df$K_atac = sapply(runs, function(w) {
      ass_atac = Rcongas:::detensorize(w$inferred_params$assignment_atac, CUDA)
      return(length(unique(ass_atac)))
    })
    model_selection_df$lambda = sapply(runs, function(w) Rcongas:::detensorize(w$hyperparameters$lambda, CUDA))

    if(model_selection_df %>% is.na %>% any)
    {
      cli::cli_alert_danger("NA found in model scores, aborting")
      model_selection_df[!complete.cases(model_selection_df), ] %>% print
      stop("Cannot select best model!")
    }

    bms_idx <- order(model_selection_df %>%  pull(!!model_selection))

    runs <-  runs[bms_idx]

    best_fit  <-  format_best_model(x, runs[[1]], same_mixing)

    x$runs <-  runs
    x$best_fit <-  best_fit
    x$model_selection <-  model_selection_df
    x$used_IC <- model_selection

    x %>%  print

    return(x)

  }


fit_congas_single_run <-
  function(x,
           parameters,
           l,
           K,
           learning_rate,
           latent_variables,
           CUDA,
           steps,
           temperature,
           threshold,
           patience,
           same_mixing) {

    # cli::cli_h3("Fit with k = {.field {K}}.")

    cg <- reticulate::import("congas")
    cg_mod <- reticulate::import("congas.models")
    pyro_optim <- reticulate::import("pyro.optim")

    data <- input_data_from_rcongas(x, CUDA = CUDA)
    param_optimizer <-  list()
    parameters$K <-  as.integer(K)
    parameters$lambda <- as.double(l)
    parameters$Temperature <-  temperature
    parameters$equal_mixture_weights = same_mixing
    parameters$CUDA = CUDA
    param_optimizer$lr <- learning_rate

    ### LOAD MODEL AND GUIDE ###
    if (latent_variables == "M") {
      model <- cg_mod$LatentCategorical
      model_string <-  "LatentCategoricalMarginalized"
      parameters$latent_type <- "M"
    } else if (latent_variables == "G"){

      model <- cg_mod$LatentCategorical
      model_string <-  "LatentCategoricalGumble"
      parameters$latent_type <- "G"

    } else if (latent_variables == "B"){

      model <- cg_mod$LatentCategorical
      model_string <-  "LatentCategoricalBaseline"
      parameters$latent_type <- "B"

   } else {
      stop("Continous latent variable model of this type not yet implemented. Use one of c('M', 'G', 'B')")
    }

    ### LOAD THE ELBO-LOSS ###

    elbo_p <- reticulate::import("pyro.infer")
    elbo <- elbo_p$Trace_ELBO
    elbo_string <- "Trace_ELBO"
 
    int <- cg$Interface(CUDA = CUDA)

    int$set_model(model)
    int$set_optimizer(pyro_optim$ClippedAdam)
    int$set_loss(elbo)
    int$initialize_model(data)
    int$set_model_params(parameters)

    warnings <-  reticulate::import("warnings")

    #ignore tracer warnings
    wcont = warnings$catch_warnings()

    wcont$`__enter__`()
    warnings$filterwarnings("ignore")
    loss = int$run(steps = as.integer(steps), param_optimizer = param_optimizer, e = threshold, patience = patience)
    wcont$`__exit__`()


    fit_params = int$learned_parameters()
    ICs = int$calculate_ICs()
    params_history = int$params_history


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
      "hyperparameters" = lapply(hyperparams, detensorize, CUDA = CUDA),
      "ICs" = ICs %>% as.data.frame(),
      "loss" = loss,
      "params_history" = params_history
    )

    class(ret) <- "congas"
    return(ret)


  }
