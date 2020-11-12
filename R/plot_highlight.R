#' Title
#'
#' @param x
#' @param alpha
#'
#' @return
#' @export
#'
#' @examples
plot_highlights = function(x, alpha = 0.05)
{
  # Best model
  best_model <- Rcongas:::get_best_model(x)

  res <-
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
    dplyr::mutate(cluster = as.character(cluster)) %>%
    dplyr::as_tibble()

  # Normalise CNA values for comparisons -- z_score_alike via offset_amplitude
  means = res %>%
    dplyr::group_by(chr, from, to) %>%
    dplyr::summarise(segment_mean = mean(CN), .groups = 'keep') %>%
    dplyr::ungroup()

  # Bind everything
  res = res %>%
    dplyr::full_join(means, by = c('chr', 'from', 'to')) %>%
    dplyr::mutate(lognorm_CN = CN,
                  CN = CN - segment_mean)

  # Clones
  clones = Rcongas::get_clusters(x)$cluster %>% unique %>% sort
  res_split = res %>% group_split(cluster)
  names(res_split) = sapply(res_split, function(x)
    x$cluster[1])


  # Plot
  plots_list = apply(combn(clones, 2, simplify = TRUE),
                     2,
                     function(w)
                     {
                       i = w[1]
                       j = w[2]

                       x_i = res_split[[i]]

                       x_j = res_split[[j]]

                       x_ij = full_join(x_i,
                                        x_j,
                                        by = c("chr", "from", "to"),
                                        suffix = c(".i", ".j")) %>%
                         mutate(diff_ij = abs(lognorm_CN.i - lognorm_CN.j))

                       fit_ij <-
                         fitdistrplus::fitdist(x_ij$diff_ij, "gamma")

                       # Fit density function
                       den_ij = data.frame(
                         x = seq(0, max(x_ij$diff_ij) * 1.25, length.out = 30),
                         y = dgamma(
                           seq(0, max(x_ij$diff_ij) * 1.25, length.out = 30),
                           fit_ij$estimate[[1]],
                           fit_ij$estimate[[2]]
                         ),
                         stringsAsFactors = FALSE
                       )

                       # p-value
                       cut_p = qgamma(1 - alpha, fit_ij$estimate[[1]], fit_ij$estimate[[2]])

                       # cut_p = 1 - pgamma(alpha, fit_ij$estimate[[1]], fit_ij$estimate[[2]])
                       x_ij$highlight = x_ij$diff_ij > cut_p

                       pl = ggplot(x_ij,
                              aes(diff_ij)) +
                         geom_histogram(aes(fill = highlight), bins = 30) +
                         scale_fill_manual("Highlight",
                                           values = c(`TRUE` = 'forestgreen', `FALSE` = 'darkred')) +
                         geom_line(data = den_ij,
                                   aes(x = x, y = y),
                                   inherit.aes = FALSE) +
                         geom_point(
                           data = den_ij,
                           aes(x = x, y = y),
                           inherit.aes = FALSE,
                           color = 'black',
                           fill = 'orange',
                           pch = 21
                         ) +
                         CNAqc:::my_ggplot_theme() +
                         geom_vline(xintercept = cut_p,
                                    color = 'darkred',
                                    linetype = 'dashed') +
                         labs(
                           x = "Log-normal difference",
                           title = paste0(i, " versus ", j),
                           subtitle = paste0("Alpha-level ", alpha)
                         )

                       if(any(x_ij$highlight)) pl = pl +
                         ggrepel::geom_label_repel(
                           data = x_ij %>% filter(highlight),
                           aes(x = diff_ij, y = Inf, label = chr),
                           inherit.aes = FALSE,
                           ylim = c(1, Inf),
                           size = 2
                         )

                       pl
                     })

  ggarrange(plotlist = plots_list, nrow = 1)
}
