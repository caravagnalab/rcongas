


segment_selector <- function(x, what="congas", K=1:3, score="ICL", mod="ATAC", CUDA = FALSE, ...){


  if(what=="congas"){

       obj=segments_selector_congas(x, K = K, score = score, CUDA = CUDA, ...)

  }else if(what=="nbmix"){

       obj=segments_selector_nbmix(x,K=K,score=score,mod=mod)

  }else{

      stop("Unknown method")
  }

     return(obj)
}



segments_selector_nbmix <- function(obj, K = 1:3, score = "ICL", modality="ATAC") {

  fit_mixture = fit_nbmix(obj, K = K, score = score, mod=modality) %>% filter(rank == 1)

  plot = fit_mixture$plot

  x_red = fit_mixture %>% filter(status == "Polyclonal")

  polyclonal_segments = x_red$segment_id %>% unique()

  filt_cyto <- obj %>% get_input(what = "segmentation") %>%
    filter(segment_id %in% polyclonal_segments)

  filt_dataset <- obj %>% get_input(what = "data") %>%
    filter(segment_id %in% polyclonal_segments)

  if (!(length(rownames(filt_cyto)) == 0)) {
    obj$input$dataset = filt_dataset

    obj$input$segmentation = filt_cyto

  }

  obj$polyclonal_segments_plot <- plot

  obj$polyclonal_segments <- polyclonal_segments

  return(obj)

}


fit_nbmix = function(x, K = 1:3, score="ICL", mod="ATAC")
{
  library(purrr)

  # EM
  fit_em = function(x, k, segment_id)
  {
    # dnbinom(x, size, prob, mu, log = FALSE) is the likelihood function for a NB using size, mu as parameters
    fit_nb_mle <- function(x) {
      ll <- function(size, mu) {
        -sum(dnbinom(x, size = size, mu = mu, log = TRUE))
      }

      m <-
        stats4::mle(
          ll,
          start = list(size = 1, mu = mean(x) %>% round),
          method = "L-BFGS-B",
          lower = c(.001, .001),
          upper = c(Inf, max(x))
        )

      ab <- stats4::coef(m)

      tibble(size = ab[1], mu = ab[2], number = length(x))
    }

    # EM algorithm NB likelihood
    iterate_em <- function(state, ...)
    {
      fits <- state$assignments %>%
        group_by(cluster) %>%
        do(mutate(fit_nb_mle(.$value), number = nrow(.))) %>%
        ungroup() %>%
        mutate(prior = number / sum(number))

      all_model <- state$assignments %>%
        select(cell:value) %>%
        crossing(fits) %>%
        mutate(likelihood = prior * dnbinom(x = value, mu = mu, size = size, log = FALSE)) # Lik

      assignments = all_model %>%
        group_by(cell) %>%
        top_n(1, likelihood) %>%
        ungroup()

      list(assignments = assignments, all_model = all_model, fits = fits)
    }

    # NLL
    loglik = function(x) { -sum(x$assignments$likelihood) }

    # initialize assignments
    km = kmeans(x$value, centers=k, nstart = 100,  trace=FALSE)

    index = 1:k
    labels = sample(index, round(length(km$cluster) * 0.3), replace = T)
    km$cluster[1:round(length(km$cluster) * 0.3)] = labels
    x$cluster = km$cluster


    #Starting point random
    x = x %>%
      rowwise() %>%
      mutate(value = value %>% round) %>%
      ungroup()

    y = NULL
    repeat {
      # cat('.')

      iterations <-
        accumulate(1:5, iterate_em, .init = list(assignments = x))

      y = append(y, iterations)
      delta_ll = abs((iterations[[5]] %>% loglik) - (iterations[[4]] %>% loglik))

      # cat(iterations[[5]] %>% loglik, '\n')
      # cat(iterations[[4]] %>% loglik, '\n')
      # cat(delta_ll < epsilon, '\n')
      # cat(delta_ll, '\n \n')

      if (delta_ll < 1e-4)
        break

       x = iterations[[5]]$assignments
    }

    # cat("FINISHED")

    final_step = y[[length(y)]]

    H = final_step$all_model %>%
      mutate(posterior = prior*likelihood) %>%
      group_by(cell) %>%
      mutate(posterior = posterior/sum(posterior)) %>%
      mutate(entropy= -sum(posterior*log(posterior))) %>%
      distinct(cell, entropy) %>%
      ungroup()


    # Values returned
    NLL = final_step$all_model %>%
      group_by(cell) %>%
      summarise(lik = log(sum(likelihood))) %>%
      summarise(lik = sum(lik))

    K_effective = final_step$fits %>% filter(number > 0) %>% nrow

    N_params = K_effective + 2 * K_effective
    BIC = N_params * log(x %>% nrow) + 2 * NLL
    ICL = BIC + sum(H$entropy)


    # Plots assignments
    P = ggplot(final_step$assignments,
               aes(value, fill = cluster %>% paste)) +
      geom_histogram(bins = 50) +
      guides(fill = FALSE) +
      labs(title = paste0("K= ", K_effective, ", segment_id = ", segment_id )) +
      theme_linedraw(base_size = 8)

    # Density
    PD = NULL
    for(i in final_step$fits$cluster %>% seq){

      PD_df = data.frame(
        x = seq(min(x$value), max(x$value), 10),
        cluster = final_step$fits$cluster[i]
      ) %>%
        mutate(
          y = dnbinom(x, mu = final_step$fits$mu[i], size = final_step$fits$size[i]) * final_step$fits$prior[i]
        )
      PD = bind_rows(PD, PD_df)
    }

    PD = ggplot(PD,
                aes(x = x, y = y, cluster = cluster %>% paste)) +
      geom_line() +
      guides(fill = FALSE)  +
      theme_linedraw(base_size = 8)

    P = cowplot::plot_grid(P, PD, ncol = 2, align = 'h', axis = 'tb')

    tibble(
      k = K_effective,
      NLL = NLL%>% as.numeric,
      BIC = BIC %>% as.numeric,
      ICL = ICL %>% as.numeric,
      plot = P %>% list,
      nsteps = length(y)
    ) %>%
      return
  }

  # Fit for k
  runner_k = function(x,segment_id)
  {
    trials = expand.grid(k = K, r = 1:10)

    # Model selection
    # trials_fit = lapply(1:nrow(trials), function(i){
    #   fit_em(x = x, k = trials$k[i], epsilon)
    # })

    # Model selection
    trials_fit = easypar::run(
      FUN = function(i){
        fit_em(x = x, k = trials$k[i], segment_id=segment_id)
      },
      PARAMS = lapply(1:nrow(trials), list),
      parallel = FALSE
    )

    if(score=="BIC"){

    trials_fit = Reduce(bind_rows, trials_fit) %>%
      arrange(BIC)

    }else if(score=="ICL"){

      trials_fit = Reduce(bind_rows, trials_fit) %>%
        arrange(ICL)

    }else{

         stop("unkonwn score")
    }
}

  # Modality data, normalised by factors
  modality_data = x %>%
    get_data() %>%
    filter(modality == mod)

  atac_norm_factors = x %>% get_input(what = "normalisation") %>% filter(modality == mod)

  modality_data = Rcongas:::normalise_modality(modality_data, atac_norm_factors) %>% group_split(segment_id)
  names(modality_data) = x$input$segmentation$segment_id %>% unique()

  all_fits = lapply(modality_data %>% names, function(x)
  {
    values = modality_data[[x]]

    cli::cli_h3("NB mixture: {.field {x}}")

    results = runner_k(values %>% select(cell, value),x) %>%
      mutate(segment_id = x)

    status = ifelse(results$k[1] == 1, "Monoclonal", "Polyclonal")

    cat("\n")
    cli::cli_alert("Fit K = {.field {results$k[1]}}, {.field {status}}")

    return(results)
  })

  all_fits = Reduce(bind_rows, all_fits)

  clonal_status = all_fits %>%
    group_by(segment_id) %>%
    mutate(rank = row_number()) %>%
    mutate(status = ifelse(k > 1, "Polyclonal", "Monoclonal"))

  cli::cli_h1("Status")

  cli::cli_alert("Polyclonal")
  print(clonal_status %>% filter(rank == 1, status == "Polyclonal"))

  cli::cli_alert("Monoclonal")
  print(clonal_status %>% filter(rank == 1, status == "Monoclonal"))

  return(clonal_status)
}

#' Filter segments where the optimal number of clusters is 1.
#' 
#' @description This function can be used to exclude from the inference step those segments that are unimodal.
#' 
#' The function requires and returns an (R)CONGAS+ object.
#'
#' @param obj An \code{rcongasplus} object.
#' @param K_max Max Number of clusters to test for the congas runs on single segments
#' @param score Model selection score to be used to select the number of clusters
#' @param lambda lambda to use for each congas run
#' @param cores_ratio fraction of cores that we be used to run the single CONGAS+ inferences in parallel.
#' 
#' @return The object \code{obj} where segments have been identified and 
#' removed.
#' 
#' @export
#'
segments_selector_congas <- function(obj, K_max = 3, score = "BIC", lambda = 0.5, cores_ratio = 0.5, CUDA = F, binom_limits = c(40,1000)){
  print('segment_selector')
  congas_single_segment <- function(obj, binom_limits, CUDA, lambda){
    # seg_id = segment_ids[i]
    K_max = 3#params$K_max
    lr = 0.01#params$lr
    steps = 2000#params$steps
    score = 'BIC'#params$score
    temperature = 20#params$temperature

    # obj = Rcongas:::select_segments(obj, segment_ids = c(seg_id))

    model_params = Rcongas:::auto_config_run(obj,
                                            K = 1:K_max,
                                            prior_cn = c(0.2, 0.6, 0.05, 0.025, 0.025),
                                            CUDA = CUDA)

    model_params$lambda=lambda

    model_params$binom_prior_limits = binom_limits

    fit_obj = Rcongas:::fit_congas(
      obj,
      K = 1:K_max,
      model_parameters = model_params,
      lambdas = lambda,
      learning_rate = lr,
      steps = steps,
      temperature = temperature,
      model_selection = score,
      threshold = 0.001,
      CUDA = CUDA,
      same_mixing = FALSE
    )

    # p=Rcongas:::plot_fit_density(fit_obj, highlight=F)
    bms_idx <- order(fit_obj$model_selection %>% pull(!!score))[1]

    model = tibble::tibble(
      k = fit_obj$model_selection$K[bms_idx],
      segment_id = unique(obj$input$segmentation$segment_id)
      # plot=p
    )

    return(model)

  }

  segment_ids = unique(Rcongas:::get_data(obj)$segment_id)

  # params = list(K_max = 3,
  #              lambda=0.5,
  #              lr=0.01,
  #              steps=2000, 
  #              score="BIC", 
  #              purity = NULL, 
  #              temperature=20)

  report = easypar::run(
    FUN = congas_single_segment,
    PARAMS = lapply(segment_ids, function (x) {
      list(obj = Rcongas:::select_segments(obj, segment_ids = c(x)) %>%
        adjust_multiome_cells(),
          binom_limits = binom_limits,
          CUDA = CUDA,
          lambda = lambda)
      }),
    parallel = FALSE,
    cores.ratio = cores_ratio,
    filter_errors = FALSE
  )
  print(report)
  report = Reduce(bind_rows, report)

  polyclonal_segments <-
    filter(report, k > 1) %>% dplyr::select(segment_id) %>% unique()

  if (!(length(rownames(polyclonal_segments)) == 0)) {

    obj = Rcongas:::select_segments(obj, polyclonal_segments$segment_id)

  }
  cat(nrow(polyclonal_segments), " segments are found polyclonal")
  # obj$segments_plot <- segments_plot

  obj$polyclonal_segments <- polyclonal_segments

  return(obj)

}


select_segments <- function(x, segment_ids) {
  x_filt = x %>% get_input(what = "data") %>%
      filter(segment_id %in% segment_ids)
  cyto = x %>% 
      get_input(what = "segmentation") %>%
      filter(segment_id %in% segment_ids)

  x$input$dataset = x_filt
  x$input$segmentation = cyto

  return(x)

}


adjust_multiome_cells = function(x) {
  if(!all(x$input$multiome)) return(x)
  cells_keep = x$input$dataset %>% select(modality, multiome_barcode) %>% distinct() %>%
	group_by(multiome_barcode) %>% summarise(nmodalities = n()) %>% 	
	filter(nmodalities == 2) %>% pull(multiome_barcode)

  x$input$dataset = x$input$dataset %>% filter(multiome_barcode %in% cells_keep)
  x$input$normalisation = x$input$normalisation %>% filter(multiome_barcode %in% cells_keep)
  return(x)
}




