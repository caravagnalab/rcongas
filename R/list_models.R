#' List models available in the package.
#' 
#' @description Returns a simple description of the models that can be 
#' choosen to fit the data.
#'
#' @return Nothing
#' @export
#'
#' @examples
#' list_models()
list_models <-  function() {
  
  cli::cli_h3(paste0("Models available in version: ", packageVersion("Rcongas"))) 
  
  data.frame(
      model = c(
        "HmmSimple",
        "MixtureGaussian", # Shoule it have Poisson in the name?
        "MixtureGaussianDMP",
        "MixtureDirichlet",
        "HmmMixtureRNA",
        "HmmSegmenter"
      ),
      synopsis = c(
        NA,
        "Poisson-based mixture with Lognormal copy-numbers",
        NA,
        "Dirichlet-process version of XXXXXXXXXXXXX",
        NA,
        NA
        ),
      vignette = 
        c(
          "link to vignette for the model",
          
          NA
        )
      ) %>% 
    # as_tibble() %>% 
    print(row.names = F)
}