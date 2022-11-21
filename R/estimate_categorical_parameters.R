estimate_segment_factors <-
  function(data, norm_factors, pld, plot = T) {

    cli::cli_h3("Estimating segment factors")

    x <-  apply(data, 2, function(y)
      y /  norm_factors)

    if (ncol(x) > 1) {
      x <-  apply(x, 1, function(y)
        y / pld) %>% t

    } else{
      x <- x / pld
    }

    res <-
      lapply(1:ncol(x), function(y) {
        w = estimate_segment_factors_aux(data_mle = x[, y], plot = plot)
        cli::cli_alert("{y}: {.field {colnames(x)[y]}} theta_shape = {w['theta_shape']}, theta_rate = {w['theta_rate']}")

        w
      })

    names(res) <-  colnames(data)
    return(res)
  }

estimate_segment_factors_aux <- function(data_mle, plot)
{
  # data_mle = x[, 2]
  # hist(data_mle, breaks = 100)

  # data_mle <-  ifelse(data_mle == 0, rnorm(1, 1e2, 1e2), data_mle)
  data_mle = data_mle[data_mle > 0]
  data_mle = data_mle[!is.na(data_mle)]

  quants = quantile(data_mle, probs = c(0.005, 0.995), na.rm = TRUE)

  data_mle <-  ifelse(data_mle < quants[1], quants[1], data_mle)
  data_mle <-  ifelse(data_mle > quants[2], quants[2], data_mle)

  MAXT = 10
  success = FALSE
  while (!success) {
    tryCatch(expr = {

      MAXT <-  MAXT -1

      if(MAXT == 0){
        success = TRUE
      }

      if (MAXT < 9) {
        coeff <-
        fitdistrplus::fitdist(data_mle
          ,distr = "gamma", method = "mle", lower = c(0,0)
        )
      } else {
        coeff <-
          fitdistrplus::fitdist(data_mle
          ,distr = "gamma", method = "mle"
        )
      }


      summary(coeff)


      if(is.null(coeff))
        stop()

      coeff <- coeff$estimate

      names(coeff) <-  c("theta_shape", "theta_rate")

      success = TRUE



    }, error = function(e) {cli::cli_alert_warning("Iteration not successfull, retrying by setting a lower bound to the parameters")})
  }

  if(MAXT == 0){
    cli::cli_alert_danger("Max number of iteration reached, returning NULL!")
    return(c("thetha_shape" = NA, "theta_rate" = NA))
  }

  if (plot) {
    hist(
      data_mle,
      prob = TRUE,
      col = "grey" ,
      main = "Inferred segment prior",
      xlab = "",
      breaks = 500
    )
    curve(
      dgamma(x, shape = coeff[1], rate = coeff[2]),
      from = min(data_mle),
      to = max(data_mle),
      n = 1000,
      col = "blue",
      add = TRUE,
      main = "",
      xlab = "",
      lwd = 2
    )
  }

  return(coeff)
}
