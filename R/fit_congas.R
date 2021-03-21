fit_congas <-  function(x, K, learning_rate, model_parameters, latent_variables = "D", compile = FALSE, steps = 500){

  if(!inherits(x, "rconags")) {
    stop("Input object needs to be an rcongas instance!")
  }




}


fit_congas_single_run <-  function(x, parameters,K, learning_rate,latent_variables, compile, steps){

  cg <- reticulate::import("congas")
  cg_mod <- reticulate::import("congas.models")

  ### LOAD MODEL AND GUIDE ###
  if(latent_variables == "D"){
    model <- cg_mod$La
  } else {
    stop("Continous latent variable model not yet implemented.")
  }

  ### LOAD THE ELBO-LOSS ###

  elbo <- reticulate::import("pyro.infer")
  if(compile){
    elbo <- elbo$TraceGraph_ELBO
  } else {
    elbo <-  elbo$JitTraceGraph_ELBO
  }

  int$run_analysis(steps = as.integer(steps), lr = learning_rate)





}
