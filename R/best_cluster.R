
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
           posteriors = FALSE,
           seed = 3,
           step_post = 300,
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
          posteriors = posteriors,
          seed = seed,
          step_post = step_post
        ))
    
    if (grepl(tolower(model), pattern = "norm")) {
      lik_fun <-  gauss_lik_norm
    } else if (grepl(tolower(model), pattern = "EXP")) {
      lik_fun <- gauss_lik_with_means
    } else if (grepl(tolower(model), pattern = "Old")) {
      lik_fun <- gauss_lik_old
    } else {
      lik_fun <-  gauss_lik
    }
    
    if (method == "BIC")
    {
      IC <-
        sapply(res, function(x)
          calculate_BIC(x, X$data$counts, X$data$cnv$mu, llikelihood = lik_fun))
    } else if (method == "AIC") {
      IC <-
        sapply(res, function(x)
          calculate_AIC(x, X$data$counts, X$data$cnv$mu, llikelihood = lik_fun))
    } else if (method == "ICL") {
      IC <-
        sapply(res, function(x)
          calculate_ICL(x, X$data$counts, X$data$cnv$mu, llikelihood = lik_fun))
    } else if (method == "lk") {
      IC <-
        sapply(res, function(x)
          lik_fun(X$data$counts, X$data$cnv$mu, x$parameters) * -1)
    } else {
      stop("Information criterium not present in the package")
    }
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
