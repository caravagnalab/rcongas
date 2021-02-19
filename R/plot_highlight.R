
#' plot_highlights = function(x, alpha = 0.05)
#' {
#'   # Best model
#'   best_model <- Rcongas:::get_best_model(x)
#'   
#'   res <-
#'     data.frame(
#'       t(best_model$parameters$cnv_probs),
#'       stringsAsFactors = FALSE,
#'       check.names = FALSE
#'     ) %>%
#'     dplyr::mutate(segment_id = rownames(.)) %>%
#'     deidify() %>%
#'     reshape2::melt(
#'       id.vars = c("chr", "from", "to"),
#'       variable.name = "cluster",
#'       value.name = "CN"
#'     ) %>%
#'     dplyr::mutate(cluster = as.character(cluster)) %>% mutate(cluster = ifelse(
#'       str_starts(cluster, pattern = "c|C"),
#'       cluster,
#'       paste0("c", cluster)
#'     )) %>%
#'     dplyr::as_tibble() %>%
#'     dplyr::mutate(lognorm_CN = CN)
#'   
#'   res$highlights = create_highlights(res, alpha)
#'   
#'   # Normalise CNA values for comparisons
#'   # means = res %>%
#'   #   dplyr::group_by(chr, from, to) %>%
#'   #   dplyr::summarise(segment_mean = mean(CN), .groups = 'keep') %>%
#'   #   dplyr::ungroup()
#'   #
#'   # # Bind everything
#'   # res = res %>%
#'   #   dplyr::full_join(means, by = c('chr', 'from', 'to')) %>%
#'   #   dplyr::mutate(lognorm_CN = CN,
#'   #                 CN = CN - segment_mean)
#'   
#'   res$highlight = create_highlights(res, alpha)
#'   
#'   # Clones
#'   clones = Rcongas::get_clusters(x)$cluster %>% unique %>% sort
#'   
#'   
#'   res_split = res %>% group_split(cluster)
#'   names(res_split) = sapply(res_split, function(x)
#'     x$cluster[1])
#'   
#'   df_density = NULL
#'   
#'   for (i in clones)
#'     for (j in clones)
#'     {
#'       if (i != j)
#'       {
#'         f  = get_density_gamma_highlights(res, i, j)
#'         x_ij = x %>%
#'           filter(cluster %in% c(i, j)) %>%
#'           group_by(chr, from, to) %>%
#'           summarise(diff_ij = abs(diff(lognorm_CN)))
#'         
#'         # Fit density function
#'         den_ij = data.frame(
#'           x = seq(0, max(x_ij$diff_ij) * 1.25, length.out = 30),
#'           y = dgamma(
#'             seq(0, max(x_ij$diff_ij) * 1.25, length.out = 30),
#'             fit_ij$estimate[[1]],
#'             fit_ij$estimate[[2]]
#'           ),
#'           stringsAsFactors = FALSE
#'         )
#'         
#'         df_density = df_density %>%
#'           bind_rows(den_ij %>% mutate(pair = paste0(i, ' vs ', j)))
#'       }
#'     }
#'   
#'   ggplot(df_density, aes(x=x,y=y, color=pair))+geom_point()
#'   
#' }
#' 
#' # 
#' # # Plot
#' # plots_list = apply(combn(clones, 2, simplify = TRUE),
#' #                    2,
#' #                    function(w)
#' #                    {
#' #                      i = w[1]
#' #                      j = w[2]
#' #                      
#' #                      x_i = res_split[[i]]
#' #                      
#' #                      x_j = res_split[[j]]
#' #                      
#' #                      x_ij = full_join(x_i,
#' #                                       x_j,
#' #                                       by = c("chr", "from", "to"),
#' #                                       suffix = c(".i", ".j")) %>%
#' #                        mutate(diff_ij = abs(lognorm_CN.i - lognorm_CN.j))
#' #                      
#' #                      fit_ij <-
#' #                        fitdistrplus::fitdist(x_ij$diff_ij, "gamma")
#' #                      
#' #                      # Fit density function
#' #                      den_ij = data.frame(
#' #                        x = seq(0, max(x_ij$diff_ij) * 1.25, length.out = 30),
#' #                        y = dgamma(
#' #                          seq(0, max(x_ij$diff_ij) * 1.25, length.out = 30),
#' #                          fit_ij$estimate[[1]],
#' #                          fit_ij$estimate[[2]]
#' #                        ),
#' #                        stringsAsFactors = FALSE
#' #                      )
#' #                      
#' #                      # p-value
#' #                      cut_p = qgamma(1 - alpha, fit_ij$estimate[[1]], fit_ij$estimate[[2]])
#' #                      
#' #                      # cut_p = 1 - pgamma(alpha, fit_ij$estimate[[1]], fit_ij$estimate[[2]])
#' #                      x_ij$highlight = x_ij$diff_ij > cut_p
#' #                      
#' #                      pl = ggplot(x_ij,
#' #                                  aes(diff_ij)) +
#' #                        geom_histogram(aes(fill = highlight), bins = 30) +
#' #                        scale_fill_manual("Highlight",
#' #                                          values = c(`TRUE` = 'forestgreen', `FALSE` = 'darkred')) +
#' #                        geom_line(data = den_ij,
#' #                                  aes(x = x, y = y),
#' #                                  inherit.aes = FALSE) +
#' #                        geom_point(
#' #                          data = den_ij,
#' #                          aes(x = x, y = y),
#' #                          inherit.aes = FALSE,
#' #                          color = 'black',
#' #                          fill = 'orange',
#' #                          pch = 21
#' #                        ) +
#' #                        CNAqc:::my_ggplot_theme() +
#' #                        geom_vline(xintercept = cut_p,
#' #                                   color = 'darkred',
#' #                                   linetype = 'dashed') +
#' #                        labs(
#' #                          x = paste("Log-normal difference (alpha = ", alpha, ')'),
#' #                          title = paste0(i, " vs ", j)
#' #                        ) +
#' #                        theme(legend.position = c(.8, .8))
#' #                      
#' #                      # if(any(x_ij$highlight)) pl = pl +
#' #                      #   ggrepel::geom_label_repel(
#' #                      #     data = x_ij %>% filter(highlight),
#' #                      #     aes(x = diff_ij,  label = chr),
#' #                      #     inherit.aes = FALSE,
#' #                      #     y = max(den_ij$y),
#' #                      #     # y = -0.1,
#' #                      #     ylim = c(max(den_ij$y) - 1,max(den_ij$y)-1),
#' #                      #     size = 2
#' #                      #   )
#' #                      
#' #                      pl
#' #                    })
#' # 
#' # ggarrange(plotlist = plots_list, nrow = 1)
#' # }


#' Title
#'
#' @param x
#' @param alpha
#'
#' @return
#' @export
#'
#' @examples
#'
#' x = Rcongas::congas_example
#'
#' print(x)
#'
#' # Default view
#' plot_highlights(x)
#'
#' # More stringent
#' plot_highlights(x, alpha = 0.01)
#'
#' # Even more stringent
#' plot_highlights(x, alpha = 0.001)
plot_highlights = function(x, alpha = 0.05)
{
  stopifnot(has_inference(x))
  
  tests_table = highlights(x, alpha) %>%
    mutate(label = paste(cluster , 'vs', versus))
  
  nh = tests_table %>% filter(highlight) %>% nrow
  
  # Density points
  params = tests_table %>% distinct(cluster, versus, shape, rate)
  
  df_points = lapply(1:nrow(params), function(i) {
    x_ij = tests_table %>%
      filter(cluster == params$cluster[i], versus == params$versus[i])
    
    domain_x = seq(0, max(x_ij$diff) * 1.25, length.out = 30)
    image_y = dgamma(domain_x,
                     params$shape[i],
                     params$rate[i])
    
    data.frame(
      x = domain_x,
      y = image_y,
      cluster = params$cluster[i],
      versus = params$versus[i]
    )
  }) %>%
    Reduce(f = bind_rows) %>%
    as_tibble() %>%
    mutate(label = paste(cluster , 'vs', versus))
  
  # Plot
  ggplot(tests_table,
         aes(diff, fill = highlight)) +
    geom_histogram(bins = 30) +
    facet_wrap(~ label) +
    CNAqc:::my_ggplot_theme() +
    scale_fill_manual(bquote("Above " * q[alpha] * ""),
                      values = c(`TRUE` = 'indianred3', `FALSE` = 'gray')) +
    labs(
      x = "Absolute CN difference",
      y = "Observations",
      title = bquote(.(nh)~"highlight(s) (" * q[alpha] * " = " * .(1 - alpha) *
                       ")")
    ) +
    geom_line(
      data = df_points,
      aes(x = x, y = y),
      color = 'black',
      inherit.aes = F,
      size = .5
    )
}