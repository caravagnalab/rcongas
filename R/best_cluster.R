
#' Model selection with CONGAS.
#' 
#' @description Selects the best number of clusters \code{k} to fit a dataset.
#' Options allow to select the type of model available in the package, parameters
#' for the Pyro optimizer and a number of different configurations to learn
#' Bayesian posteriors or MAP estimates.
#' 
#' @param X 
#' @param model Any of the model keywords available via \code{list_models()}.
#' @param clusters A vector of values for \code{k}.
#' @param optim 
#' @param elbo 
#' @param inf_type 
#' @param steps 
#' @param lr 
#' @param param_list 
#' @param MAP 
#' @param posteriors 
#' @param seed 
#' @param step_post 
#' @param mixture 
#' @param method 
#'
#' @return
#' @export
#'
#' @examples
best_cluster <-
  function(X ,
           model,
           clusters ,
           optim = "ClippedAdam",
           elbo = "TraceEnum_ELBO",
           inf_type = "SVI",
           steps = 300,
           lr = 0.05,
           param_list = list(),
           MAP = TRUE,
           seed = 3,
           mixture = NULL,
           method = "AIC") {
    if (is.null(mixture)) {
      mixture <- rep(NULL, length(clusters))
    }
    
    res <-
      lapply(clusters, function(x)
        run_inference(
          X =  X,
          model = model,
          optim = optim,
          elbo = elbo,
          inf_type = inf_type,
          steps = steps,
          lr = lr,
          param_list = c(param_list, list('K' = x, 'mixture' = mixture[[x]])),
          MAP = MAP,
          seed = seed,
        ))

    IC <- calculate_information_criteria(res,X,method)
    ind <-  which.min(IC)
    print(paste0("Best number of cluster is " , ind))
    ret <-
      list(
        models = res,
        model_selection = list(
          best_K = ind,
          IC = IC,
          IC_type = method,
          clusters = clusters
        )
      )
    X$inference <- ret
    return(X)
  }
