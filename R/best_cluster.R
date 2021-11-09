
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
#' @param optim a pyro.optim object, in string format
#' @param elbo  a pyro.inference object, in string format
#' @param inf_type currently just "SVI" is supported
#' @param steps number of step of the gradient-based VI algorithm
#' @param lr learning rate of the gradient-based VI algorithm
#' @param param_list named list of parameters for the specific model
#' @param MAP perform MAP or full bayesian inference (MAP is suggested as it converges better)
#' @param seed seed for torch random number generator
#' @param mixture initial mixture weights, if NULL then they are equal to 1/K
#' @param method informtion criteria to score the best solution (currently supports "AIC", "BIC" and "ICL)
#'
#' @return A congas (described in [\fun{init}]) object with following elements:
#' - *inference* :
#'   - *models* with all the runs performed
#'   - *model_selection* with information about the IC and run statistics 
#' 
#' @details  While this function has a lot of personalizable arguments, the default
#' ones are usually good for most of the applications. It is better to be confident with Pyro
#' before changing them
#' 
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
