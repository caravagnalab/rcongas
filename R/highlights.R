#' Title
#'
#' @param x 
#' @param alpha 
#'
#' @return
#' @export
#'
#' @examples
highlights <- function(x, alpha = 0.05)
{
  stopifnot(has_inference(x))
  
  # Best model
  best_model <- get_best_model(x)
  
  # CN table
  cn_table =
    data.frame(
      t(best_model$parameters$cnv_probs),
      stringsAsFactors = FALSE,
      check.names = FALSE
    ) %>%
    dplyr::mutate(segment_id = rownames(.)) %>%
    deidify() %>%
    reshape2::melt(
      id.vars = c("chr", "from", "to"),
      variable.name = "cluster",
      value.name = "CN"
    ) %>%
    dplyr::mutate(cluster = as.character(cluster)) %>% mutate(cluster = ifelse(
      str_starts(cluster, pattern = "c|C"),
      cluster,
      paste0("c", cluster)
    )) %>%
    dplyr::as_tibble() %>% idify()
  
  if(length(unique(cn_table$cluster)) == 1) 
  {
    cn_table$versus =
      cn_table$diff = 
      cn_table$shape =
      cn_table$rate = NA
  
    cn_table$highlight = FALSE
    
    return(cn_table)
  }
    
  
  # Pairwise testing
  clusters = combn(unique(cn_table$cluster), 2)
  
  tests_table = NULL
  
  for (i in 1:ncol(clusters))
  {
    tests_table = tests_table %>%
      bind_rows(get_density_gamma_highlights(cn_table, clusters[1, i], clusters[2, i], alpha))
  }
  
  return(tests_table)
  
  # tests_table %>% filter(highlight)
  
  # tests_table
  
  
  #
  #
  # highlights = NULL
  #
  # for (i in clusters)
  #   for (j in clusters)
  #     highlights = c(highlights, compare_clones(x, i, j, alpha))
  #
  # highlights = highlights %>% unique
  #
  
  # ggplot(x, aes(lognorm_CN, fill = cluster)) +
  #   geom_histogram() +
  #   facet_wrap(~segment_id)
  #
  # nclones <-
  #
  # ret <-  list(vector(length = nrow(x)))
  #
  # for (i in nclones)
  #   ret[[i]] <- compare_clones(x, i, nclones, alpha)
  
  
  # x$segment_id %in% highlights
}

# compare_clones <-  function(x, i, nclones, alpha = 0.05)
# {
#   x_i <-  x %>% filter(cluster == i)
#
#   nclones <-  nclones[nclones != i]
#   ret <-  logical(length = nrow(x_i))
#
#   for (j in nclones)
#   {
#     x_j <- x %>% filter(cluster == j)
#     diff <- abs(as.numeric(x_i$lognorm_CN - x_j$lognorm_CN))
#     fit <-  fitdistrplus::fitdist(diff, "gamma")
#
#     # p-value
#     cut_p = qgamma(1-alpha, fit$estimate[[1]], fit$estimate[[2]])
#
#     # p_val <- dgamma(diff, fit$estimate[[1]], fit$estimate[[2]])
#
#     ret[which(diff > cut_p)] <-  TRUE
#
#
#     # hist(diff, breaks = 30)
#     # points(diff, dgamma(diff, fit$estimate[[1]], fit$estimate[[2]]))
#     # title(sub = pgamma(1-alpha, fit$estimate[[1]], fit$estimate[[2]]))
#     # ret[which(p_val < alpha)] <-  TRUE
#   }
#
#   return(ret)
# }

get_density_gamma_highlights = function(x, i, j, alpha)
{
  x_i <-  x %>% filter(cluster == i)
  x_j <-  x %>% filter(cluster == j)
  
  diff <- abs(as.numeric(x_i$CN - x_j$CN))
  
  fit_ij = fitdistrplus::fitdist(diff, "gamma")
  
  cut_p = qgamma(1 - alpha, fit_ij$estimate[[1]], fit_ij$estimate[[2]])
  
  x_i$versus = j
  x_i$diff = diff
  
  return(
    x_i %>%
      mutate(
        highlight = diff >= cut_p,
        shape = fit_ij$estimate[[1]],
        rate = fit_ij$estimate[[2]]
      )
  )
}


