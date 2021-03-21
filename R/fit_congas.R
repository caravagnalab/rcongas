fit_congas <-  function(x, K, learning_rate, model_parameters, latent_variables = "D", compile = FALSE, steps = 500, model_selection = "ICL"){

  if(!inherits(x, "rconags")) {
    stop("Input object needs to be an rcongas instance!")
  }

  runs <- lapply(K, function(k) fit_congas_single_run(x,model_parameters, k, learning_rate, latent_variables, compile, steps))
  names(run) <-  paste(K)

  model_selection <-  sapply(runs, function(y) y$ICs)
  model_selection$K <- K

  bms_idx <- order(model_selection[[model_selection]])

  runs <-  runs[bms_idx]

  best_model  <-  format_best_model(runs[1])

  x$runs <-  runs
  x$best_model <-  best_model
  x$model_selection <-  model_selection

  return(x)

}


fit_congas_single_run <-  function(x, parameters,K, learning_rate,latent_variables, compile, steps){

  cg <- reticulate::import("congas")
  cg_mod <- reticulate::import("congas.models")
  pyro_optim <- reticulate::import("pyro.optim")

  data <- input_data_from_rcongas(x)
  param_optimizer <-  list()
  parameters$K <-  as.integer(K)
  param_optimizer$lr <- learning_rate

  ### LOAD MODEL AND GUIDE ###
  if(latent_variables == "D"){
    model <- cg_mod$LatentCategorical
    model_string <-  "LatentCategorical"
  } else {
    stop("Continous latent variable model not yet implemented.")
  }

  ### LOAD THE ELBO-LOSS ###

  elbo_p <- reticulate::import("pyro.infer")
  if(compile){
    elbo <- elbo_p$TraceGraph_ELBO
    elbo_string <- "TraceGraph_ELBO"
  } else {
    elbo <-  elbo$JitTraceGraph_ELBO
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

  hyperparams = c(list("model" = model_string, "loss" = elbo_string, "optimizer" = "ClippedAdam"),parameters )

  ret <-  list("inferred_params" = fit_params,
               "hyperparameters" = hyperparams,
               "ICs" = ICs %>% as.data.frame(),
               "loss" = loss)
  class(ret) <- "congas"
  return(ret)


}
