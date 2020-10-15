
# torch <- reticulate::import("torch")
# congas <- reticulate::import("congas")


choose_model <-  function(model_string){


  tryCatch({model <-  reticulate::import(paste0("congas.models.", model_string))

  return(reticulate::py_get_attr(model,model_string))

  }
  , error = function(e) print("Model not present, please run list_models() to see what is available"))
}


list_models <-  function(){
  cnv_mod <- c("HmmSimple", "MixtureGaussian", "MixtureGaussianDMP", "MixtureDirichlet", "HmmMixtureRNA", "HmmSegmenter")
  cat(cnv_mod, sep = "\n")
}

choose_optim <-  function(optim_string){
  tryCatch({
  optim <-  reticulate::import("pyro.optim")

  return(reticulate::py_get_attr(optim,optim_string))

  }
  , error = function(e) print("Optimizer not present, please look at congas documentation to see what is available"))
}

choose_loss <-  function(loss_string){
  tryCatch({
    infer <- reticulate::import("pyro.infer")

    return(reticulate::py_get_attr(infer,loss_string))

  }
  , error = function(e) print("Loss not present, please look at congas documentation to see what is available"))
}

choose_type <-  function(type_string){
  tryCatch({
    infer <- reticulate::import("pyro.infer")

    return(reticulate::py_get_attr(infer,type_string))

  }
  , error = function(e) print("Inferrence type not present, please look at congas documentation to see what is available"))
}


from_simulation_to_data_list <- function(x){



  return(list(data = t(x$data$counts), mu= x$data$cnv$mu, pld = as.vector(x$data$cnv$ploidy_real), segments = nrow(x$data$cnv)))
}


tensorize <- function(list){

  torch <- reticulate::import("torch")
  np <- reticulate::import("numpy")
  res <- lapply(list, function(x) torch$tensor(x, dtype = torch$float32))
  return(res)
}




set_names <-  function(an){

  names(an$parameters) <- gsub(names(an$parameters), pattern = "param_", replacement = "")

  if(length(as.vector(an$parameters$assignement)) == 1) {

    an$parameters$assignement <- rep(x = an$parameters$assignement, length(an$dim_names$cell_names) )

  }

  uniq_clus_ids <- names(sort(table(an$parameters$assignement),TRUE))


  tmp <- vector(length = length(an$parameters$assignement))
  new_clusters <-  1:length(uniq_clus_ids)

  for(j in new_clusters){
    tmp[which(an$parameters$assignement == as.integer(uniq_clus_ids[j]))] <- j
  }

  an$parameters$assignement <-  tmp



  mix_order <- order(an$parameters$mixture_weights, decreasing = TRUE)

  an$parameters$mixture_weights <-  an$parameters$mixture_weights[mix_order]
  if(!is.null(an$parameters$assignment_probs)) an$parameters$assignment_probs <- an$parameters$assignment_probs[,mix_order]
  names(an$parameters$mixture_weights) <-  new_clusters
  an$parameters$cnv_probs <- data.frame(an$parameters$cnv_probs)
  an$parameters$cnv_probs <-  an$parameters$cnv_probs[mix_order,, drop = FALSE]
  #rownames(an$parameters$cnv_probs) <- new_clusters

  colnames(an$parameters$cnv_probs) <- an$dim_names$seg_names
  if(!is.null(an$parameters$norm_factor))
    names(an$parameters$norm_factor) <-  an$dim_names$cell_names



  names(an$parameters$assignement) <-  an$dim_names$cell_names

  return(an)
}

#' Title
#'
#' @param data_list
#' @param model
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
#'@param verbose
#'
#' @return
#'
#' @importFrom dplyr %>% filter mutate select arrange desc pull row_number group_by
#' @importFrom dplyr summarise bind_cols rename bind_rows left_join distinct
#' @importFrom dplyr ungroup full_join right_join
#' @importFrom tidyr spread gather tibble tribble as_tibble
#' @importFrom magrittr %>%
#' @importFrom gtools mixedorder mixedsort
#' @importFrom purrr rbernoulli
#' @importFrom extraDistr rmvhyper
#' @export
#'
#' @examples
run_inference <-  function(X , model, optim = "ClippedAdam", elbo = "TraceEnum_ELBO", inf_type = "SVI", steps = 300, lr = 0.05,
                            param_list = list(), MAP = TRUE, posteriors = FALSE, seed = 3, step_post=300, verbose= FALSE, filt_merge = 0.5){

  an <- reticulate::import("congas")


  if(inherits(X, "CNVSimulation") | inherits(X, "rcongas")) {

    data_list <- from_simulation_to_data_list(X)
  }

  model_name <- model
  cell_names <- colnames(data_list$data)
  seg_names <- rownames(data_list$data)
  data_list <- tensorize(data_list)
  seed <- as.integer(seed)
  steps <-  as.integer(steps)
  step_post <-  as.integer(step_post)
  model <-  choose_model(model)
  optim <-  choose_optim(optim)
  elbo <-  choose_loss(elbo)
  inf_type <- choose_type(inf_type)
  lr <- as.double(lr)

  int <- an$Interface(model,optim,elbo,inf_type)

  int$initialize_model(data_list)
  int$set_model_params(param_list)



  loss <- int$run(steps=steps, seed = seed, param_optimizer=list('lr'= lr), verbose = verbose, MAP = MAP)
  parameters <- int$learned_parameters(posterior=posteriors, verbose=verbose, steps=step_post)

  if(posteriors){
    parameters$assignment_probs <- parameters$assignment_probs$numpy()
    parameters$assignement <- apply(parameters$assignment_probs, 1, which.max)
    rownames(parameters$assignment_probs) <- cell_names
  }

  dim(parameters$assignment_probs)

  dim_names <- list(cell_names = cell_names, seg_names = seg_names)

  an <-  set_names(list(loss = loss, parameters = parameters, dim_names = dim_names))

  # if(model_name == "MixtureGaussianDMP") {
  #
  #   an$parameters <- merge_clusters(an$parameters, "DMP", posterior=posteriors)
  #
  # } else {
  #
  #   an$parameters <- merge_clusters(an$parameters, type = "NONE", filt = filt_merge, posterior=posteriors)
  # }

  an$run_information <-  list(model = model,optim = optim, elbo = elbo, inf_type = inf_type,
                              steps = steps, lr = lr, input_hyper_params = param_list, MAP = MAP,
                              posteriors = posteriors, seed = seed, step_post=step_post)
  X$inference$models <- list(an)

  return(structure(an, class = "congas"))

}

best_cluster <- function(X , model, clusters ,optim = "ClippedAdam", elbo = "TraceEnum_ELBO", inf_type = "SVI", steps = 300, lr = 0.05,
                         param_list = list(), MAP = TRUE, posteriors = FALSE, seed = 3, step_post=300, verbose= FALSE, mixture = NULL, method = "AIC", filt_merge = 0.5){


  if(is.null(mixture)) {
    mixture <- rep(NULL, length(clusters))
  }

  res <- lapply(clusters, function(x) run_inference(X =  X, model = model,optim = optim, elbo = elbo, inf_type = inf_type,
                                        steps = steps, lr = lr, param_list = c(param_list, list('K' = x, 'mixture' = mixture[[x]])), MAP = MAP, posteriors = posteriors,
                                        seed = seed, step_post=step_post, verbose= verbose, filt_merge = filt_merge))

  lik_fun <-  ifelse(grepl(tolower(model), pattern = "norm"), gauss_lik_norm, gauss_lik)

  if(method == "BIC")
  {
    IC <- sapply(res, function(x) calculate_BIC(x, X$data$counts, X$data$cnv$mu, llikelihood = lik_fun))
  } else if (method == "AIC"){
    IC <- sapply(res, function(x) calculate_AIC(x, X$data$counts, X$data$cnv$mu, llikelihood = lik_fun))
  } else if (method == "ICL") {
    IC <- sapply(res, function(x) calculate_ICL(x, X$data$counts, X$data$cnv$mu, llikelihood = lik_fun))
  }else {
    stop("Information criterium not present in the package")
  }
  ind <-  which.min(IC)
  print(paste0("Best number of cluster is " , ind))
  ret <-  list(models = res, model_selection = list(best_K = ind, IC = IC, IC_type = method, clusters = clusters))
  X$inference <- ret
  return(X)
}


merge_clusters <-  function(parameters, type = "DMP", filt = 0.5, posterior=F) {

  curr <-  1
  counter <-  1
  res <- data.frame()
  weights <- vector()
  if(type == "DMP")
    betas <- vector()
  if(posterior){
    assignment_probs <- data.frame()
  assignments <-  rep(1,length(parameters$assignement) )

  unique_idx <- which(parameters$mixture_weights > 0.05)
  parameters$cnv_probs <-  parameters$cnv_probs[unique_idx,, drop = FALSE]
  parameters$mixture_weights <- parameters$mixture_weights[unique_idx]
  if(type == "DMP")
    parameters$betas <-  parameters$betas[unique_idx]
  max_rows <- nrow(parameters$cnv_probs)
  iter <-  seq(curr, max_rows)

  while(TRUE){


    mask <- sapply(iter, function(x) all(abs(parameters$cnv_probs[curr,] - parameters$cnv_probs[x,]) < filt))
    idx <- iter[which(mask)]
    res <- rbind(res, colMeans(parameters$cnv_probs[idx,, drop = F]))
    assignments[which(parameters$assignement %in% idx)] <- counter
    weights[counter] <- sum(parameters$mixture_weights[idx])
    if(type == "DMP")
      betas[counter] <- sum(parameters$betas[idx])
    if(posterior)
      assignment_probs <- rbind(assignment_probs, rowSums(parameters$assignment_probs[,idx, drop = F]))
    if(length(which(!mask)) == 0)
      break

    max_rows <- iter[max(which(!mask))]

    curr <- iter[min(which(!mask))]

    iter <-  iter[which(!mask)]
    counter <- counter + 1
  }

  colnames(res) <-  colnames(parameters$cnv_probs)
  parameters$cnv_probs <- res
  parameters$assignment_probs <- t(assignment_probs)
  names(weights) <- seq(1, length(weights))
  if(type == "DMP")
    names(betas) <- seq(1, length(betas))
  parameters$mixture_weights <- weights
  if(type == "DMP")
    parameters$betas <- betas
  names(assignments) <-  names(parameters$assignement)
  parameters$assignement <-  paste0("c",assignments)

  }

  return(parameters)

}

write_result.congas <-  function(an, new_dir = FALSE,dir_name,  out_prefix, sep = ","){
  if(new_dir){
    dir.create(file.path(".", new_dir), showWarnings = FALSE)
    out_prefix <- file.path(".", new_dir, out_prefix)
  } else {
    out_prefix <-  file.path(".", out_prefix)
  }
  write.table(x = an$paramters$mixture_weights , file = paste0(out_prefix,"_mixture_weights.csv"), sep =sep, row.names = FALSE)
  write.table(x = an$paramters$cnv_probs , file = paste0(out_prefix,"_cnv_probs.csv"), sep = sep)
  write.table(x = an$paramters$norm_factor , file = paste0(out_prefix,"_norm_factor.csv"), sep = sep)
  write.table(x = an$paramters$assignement , file = paste0(out_prefix,"_assignement.csv"), sep = sep)


}

param_total <-  function(param_list) {
  param_list <-  param_list[!(names(param_list) %in% c("assignment_probs", "assignement"))]
  res <- sapply(param_list,function(x) if(is.null(dim(x))) length(x) else prod(dim(x)))
  return(sum(res))
}

gauss_lik <-  function(data,mu,par) {

  lambdas <-  lapply(1:length(par$mixture_weights), function(x) {


    lambdas <- (matrix(par$norm_factor, ncol = 1) %*% (as.numeric(par$cnv_probs[x,]) * mu)) + 1e-8

    return(lambdas)

  })

  lk <- logLikePoisMix(data, lambdas, par$mixture_weights)

  return(lk$ll)
}


gauss_lik_with_means <-  function(data,mu,par) {

  lk <-  sapply(1:nrow(data), function(x) {

    lambdas <- (par$norm_factor[x] * round(par$cnv_probs[par$assignement[x],]) * mu) + 1e-6

    lk <- dpois(as.numeric(data[x,]), as.numeric(lambdas), log = TRUE)
    return(sum(lk))

  })

  log_lk <-  sum(lk)
  return(log_lk)
}



gauss_lik_norm <-  function(data,mu,par) {

  lk <-  sapply(1:nrow(data), function(x) {

    lambdas <- par$cnv_probs[par$assignement[x],]

    prec <- par$norm_prec[par$assignement[x],]

    lk <- dnorm(as.numeric(data[x,]),mean =  as.numeric(lambdas), sd = 1/prec, log = TRUE)
    return(sum(lk))

  })

  log_lk <-  sum(lk)

  return(log_lk)
}



calculate_BIC <-  function(inf, data, mu,llikelihood = gauss_lik) {

  n_param <- param_total(inf$parameters)
  log_lik <- llikelihood(data,mu,inf$parameters)



  return(n_param * log(nrow(data)) - 2 * log_lik)


}

calculate_AIC <-  function(inf, data, mu,llikelihood = gauss_lik) {

  n_param <- param_total(inf$parameters)
  log_lik <- llikelihood(data,mu,inf$parameters)

  return(n_param * 2 - 2 * log_lik)

}

calculate_entropy <- function(x) {
  -sum(x * log(x))
}


calculate_ICL <- function(inf, data, mu,llikelihood = gauss_lik) {

  BIC <- calculate_BIC(inf, data, mu,llikelihood)
  H <- calculate_entropy(inf$parameters$assignment_probs)
  return(BIC + H)
}


segment_genome <-  function(MAF , optim = "Adam", elbo = "TraceEnum_ELBO", inf_type = "SVI", steps = 300, lr = 0.05,
                           param_list = list(), median_filtering = TRUE, filter_size = 51, seed = 3) {

  model <-  "HmmSegmenter"



  res <- list()

  MAF_L <- MAF  %>% group_by(chr)

  MAF_splitted <-  split(MAF_L, f = MAF_l$chr)

  if(median_filtering) {

    MAF_splitted <- filter_MAF_seg(MAF_splitted, k = filter_size)

  }

  an <- reticulate::import("congas")

  seed <- as.integer(seed)
  steps <-  as.integer(steps)
  model <-  choose_model(model)
  optim <-  choose_optim(optim)
  elbo <-  choose_loss(elbo)
  inf_type <- choose_type(inf_type)

  for(chr in MAF_splitted) {
    val_mat <- matrix(chr$value,ncol = 1)
    val_mat <- val_mat + 0.001
    data_list <- list(data = torch$tensor(val_mat)$float())

    print(unique(chr$chr))
    int <- an$Interface(model,optim,elbo,inf_type)

    int$initialize_model(data_list)
    int$set_model_params(param_list)



    loss <- int$run(steps=steps, seed = seed, param_optimizer=list('lr'= lr), MAP = TRUE, verbose = FALSE)
    parameters <- int$learned_parameters(posterior=FALSE, verbose=FALSE)

    res <- c(res, list(parameters))

  }

  names(res) <- names(MAF_splitted)

  segments <-  res

  return(segments)

}


from_inference_to_segments <-  function(inf,maf) {

    maf_s <- maf %>% split(., .$chr)
    maf_s <-  maf_s[names(inf)]
    res <- list()
    for(j in seq_along(inf)){
      s_e <-  find_start_and_end(inf[[j]]$z)
      res[[j]] <- cbind(chr = as.numeric(maf_s[[j]]$chr[s_e$start]), start = as.numeric(maf_s[[j]]$start[s_e$start]),
                        end = as.numeric(maf_s[[j]]$start[s_e$end]), tot = as.numeric(s_e$tot) + 1)
    }
    return(as.data.frame(Reduce(f = rbind, res)))
}


find_start_and_end <- function(z){

  N <- length(z)
  start <- 1
  end <-  1
  start_vec <- c()
  end_vec <-  c()
  value_vec <- c()

  while(end <= N){
    if(end != N & z[end] == z[end + 1] ){
      end <- end + 1
    } else {
      start_vec <- c(start_vec, start)
      end_vec <- c(end_vec, end)
      value_vec <-  c(value_vec, z[start])
      start <-  end + 1
      end <- start
    }
  }

  return(list(start = start_vec, end = end_vec, tot = value_vec))


}


norm_factor <- function(data) {


  norm <- rowMeans((data$counts / data$cnv$mu) / data$cnv$ploidy_real)
  return(torch$tensor(norm))

}



