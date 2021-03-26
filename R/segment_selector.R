# x <- readRDS("3.rcongas_tumour.rds")

model_selection = function(x, K = 1:3, samples = 10, epsilon = 1e-4)
{
  library(purrr)
  
  # EM
  fit_em = function(x, k, epsilon = 1e-4)
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
          lower = c(0.001, .001),
          upper = c(100, max(x))
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
    
    # Starting point random
    x = x %>% 
      rowwise() %>% 
      mutate(value = value %>% round, cluster = sample(paste(1:k), 1)) %>% 
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
      
      if (delta_ll < epsilon)
        break
      
      x = iterations[[5]]$assignments
    }
    
    # cat("FINISHED")
    
    final_step = y[[length(y)]]
    
    # Values returned
    NLL = final_step$all_model %>% 
      group_by(cell) %>% 
      summarise(lik = log(sum(likelihood))) %>% 
      summarise(lik = sum(lik))
    
    K_effective = final_step$fits %>% filter(number > 0) %>% nrow
    
    N_params = K_effective + 2 * K_effective
    BIC = N_params * log(x %>% nrow) + 2 * NLL
    
    # Plots assignmenta
    P = ggplot(final_step$assignments,
               aes(value, fill = cluster %>% paste)) +
      geom_histogram(bins = 50) +
      guides(fill = FALSE) +
      labs(title = paste0("k = ", K_effective, ', NLL = ', NLL %>% round, " BIC = ", BIC %>% round)) +
      theme_linedraw(base_size = 8)
    
    # Density
    PD = NULL
    for(i in final_step$fits$cluster %>% seq)
    {
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
      plot = P %>% list,
      nsteps = length(y)
    ) %>% 
      return
  }
  
  # Fit for k
  runner_k = function(x)
  {
    trials = expand.grid(k = K, r = 1:samples)
    
    # Model selection 
    # trials_fit = lapply(1:nrow(trials), function(i){
    #   fit_em(x = x, k = trials$k[i], epsilon)
    # })
    
    # Model selection 
    trials_fit = easypar::run(
      FUN = function(i){
        fit_em(x = x, k = trials$k[i], epsilon)
      },
      PARAMS = lapply(1:nrow(trials), list),
      parallel = FALSE
    )
    
    trials_fit = Reduce(bind_rows, trials_fit) %>% 
      arrange(BIC)
  }
  
  # Modality data, normalised by factors
  modality_data = x %>% 
    get_data() %>% 
    filter(modality == "ATAC") 
  
  atac_norm_factors = x %>% get_input(what = "normalisation") %>% filter(modality == "ATAC")
  
  modality_data = Rcongas:::normalise_modality(modality_data, atac_norm_factors) %>% 
    group_split(segment_id)
  names(modality_data) = x$input$segmentation$segment_id %>% unique()
  
  all_fits = lapply(modality_data %>% names, function(x)
  {
    values = modality_data[[x]]
    
    cli::cli_h3("NB mixture: {.field {x}}")
    
    results = runner_k(values %>% select(cell, value)) %>% 
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







Message Riccardo Bergamin















