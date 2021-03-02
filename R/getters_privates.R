
# PRIVATE GETTERS

get_karyotype <-  function(x) {
  if (x$reference_genome %in% c('hg19', 'GRCh37'))
    return(Rcongas::hg19_karyo)
  if (x$reference_genome %in% c('hg38', 'GRCh38'))
    return(Rcongas::hg38_karyo)
  if (x$reference_genome %in% c('mm10', 'GRCm38'))
    return(Rcongas::mm10_karyo)
}

get_best_model <- function(X) {
  X$inference$models[[X$inference$model_selection$best_K]]
}

get_congas_model_used = function(x)
{
  if(any(is.double(x$data$counts))) return("MixtureGaussianNorm")
  if(all(is.integer(x$data$counts))) return("MixtureGaussian")

  return("HmmSegmenter")
  # MixtureGaussian
  # MixtureGaussianNorm
  # HmmSegmenter
}


# Key creation and decrypt
idify = function(y) {
  y %>% dplyr::mutate(segment_id = paste(chr, as.integer(from), as.integer(to), sep = ":"))
}

deidify = function(y) {
  y %>% tidyr::separate(segment_id,
                        into = c('chr', 'from', 'to'),
                        sep = ":") %>%
    dplyr::mutate(from = as.integer(from), to = as.integer(to))
}


get_poisson_parameters <-  function(x) {

  bm <- Rcongas:::get_best_model(x)


  lambdas <- apply(bm$parameters$cnv_probs, 1,function(y) (y * x$data$cnv$mu))

  ret  <-  reshape2::melt(lambdas %>%  as.matrix)

  colnames(ret) <- c("segment_id", "cluster", "lambda")

  ret <-  ret %>% Rcongas:::deidify() %>% Rcongas:::idify()

  return(ret)

}


get_gaussian_parameters <-  function(x) {

  bm <- get_best_model(x)

  mean <- bm$parameters$cnv_probs

  rownames(mean) <- 1:nrow(mean)

  sd <-  bm$parameters$norm_sd

  ret  <-  reshape2::melt(mean %>%  as.matrix)

  ret$sd <- rep(sd, nrow(mean))

  colnames(ret) <- c("cluster", "segment_id", "mean", "sd")


  ret <-  ret %>% deidify() %>% idify()

  return(ret)

}

is_gaussian <-  function(x) {

  bm <-  get_best_model(x)
  return(grepl("Norm", bm$run_information$model, ignore.case = T))
}

is_categorical <-  function(x){
  bm <-  get_best_model(x)
  
  return(grepl("Categorical", bm$run_information$model, ignore.case=T))
  
}

