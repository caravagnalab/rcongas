# create_highlights <- function(x, alpha = 0.05)
# {
#   x = x %>% idify()
#   
#   clusters = unique(x$cluster)
#   
#   highlights = NULL
#   
#   for(i in clusters)
#     for(j in clusters)
#       highlights = c(highlights, compare_clones(x, i,j,alpha))
#   
#   highlights = highlights %>% unique
#   
# 
#   # ggplot(x, aes(lognorm_CN, fill = cluster)) +
#   #   geom_histogram() +
#   #   facet_wrap(~segment_id)
#   # 
#   # nclones <- 
#   # 
#   # ret <-  list(vector(length = nrow(x)))
#   # 
#   # for (i in nclones)
#   #   ret[[i]] <- compare_clones(x, i, nclones, alpha)
# 
#   
#   x$segment_id %in% highlights
# }
# 
# # compare_clones <-  function(x, i, nclones, alpha = 0.05)
# # {
# #   x_i <-  x %>% filter(cluster == i)
# # 
# #   nclones <-  nclones[nclones != i]
# #   ret <-  logical(length = nrow(x_i))
# # 
# #   for (j in nclones) 
# #   {
# #     x_j <- x %>% filter(cluster == j)
# #     diff <- abs(as.numeric(x_i$lognorm_CN - x_j$lognorm_CN))
# #     fit <-  fitdistrplus::fitdist(diff, "gamma")
# # 
# #     # p-value
# #     cut_p = qgamma(1-alpha, fit$estimate[[1]], fit$estimate[[2]])
# # 
# #     # p_val <- dgamma(diff, fit$estimate[[1]], fit$estimate[[2]])
# # 
# #     ret[which(diff > cut_p)] <-  TRUE
# # 
# # 
# #     # hist(diff, breaks = 30)
# #     # points(diff, dgamma(diff, fit$estimate[[1]], fit$estimate[[2]]))
# #     # title(sub = pgamma(1-alpha, fit$estimate[[1]], fit$estimate[[2]]))
# #     # ret[which(p_val < alpha)] <-  TRUE
# #   }
# # 
# #   return(ret)
# # }
# 
# get_density_gamma_highlights = function(x, i, j)
# {
#   x_i <-  x %>% filter(cluster == i)
#   x_j <-  x %>% filter(cluster == j)
#   
#   diff <- abs(as.numeric(x_i$lognorm_CN - x_j$lognorm_CN))
#   
#   fitdistrplus::fitdist(diff, "gamma")
# }
# 
# compare_clones <-  function(x, i, j, alpha = 0.05)
# {
#   if(i == j) return(NULL)
#   
#   # x_i <-  x %>% filter(cluster == i)
#   # x_j <-  x %>% filter(cluster == j)
#   # 
#   # diff <- abs(as.numeric(x_i$lognorm_CN - x_j$lognorm_CN))
#   fit <-  get_density_gamma_highlights(x, i, j)
#     
#   # p-value
#   cut_p = qgamma(1-alpha, fit$estimate[[1]], fit$estimate[[2]])
#     
#   x$segment_id[which(diff > cut_p)]  
# }
# 
