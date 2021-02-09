create_highlights <- function(x, alpha = 0.05)
{
  nclones <- 1:length(unique(x$cluster))

  ret <-  list(vector(length = nrow(x)))

  for (i in nclones)
    ret[[i]] <- compare_clones(x, i, nclones, alpha)

  return(do.call(c, ret))
}

compare_clones <-  function(x, i, nclones, alpha = 0.05)
{
  x_i <-  x %>% filter(cluster == i)

  nclones <-  nclones[nclones != i]
  ret <-  logical(length = nrow(x_i))

  for (j in nclones) {
    x_j <- x %>% filter(cluster == j)
    diff <- abs(as.numeric(x_i$lognorm_CN - x_j$lognorm_CN))
    fit <-  fitdistrplus::fitdist(diff, "gamma")

    # p-value
    cut_p = qgamma(1-alpha, fit$estimate[[1]], fit$estimate[[2]])

    # p_val <- dgamma(diff, fit$estimate[[1]], fit$estimate[[2]])

    ret[which(diff > cut_p)] <-  TRUE


    # hist(diff, breaks = 30)
    # points(diff, dgamma(diff, fit$estimate[[1]], fit$estimate[[2]]))
    # title(sub = pgamma(1-alpha, fit$estimate[[1]], fit$estimate[[2]]))
    # ret[which(p_val < alpha)] <-  TRUE
  }

  return(ret)
}
