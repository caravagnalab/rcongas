plot_inference_report <-  function(x) {
  Z_post <- Rcongas:::is_MAP_Z(x)
  CN_post <- Rcongas:::is_MAP_CN(x)

  if (Z_post & CN_post) {
    return(plot_inference_report_1(x))
  }
  else if (!Z_post & CN_post) {
    return(plot_inference_report_2(x))
  }
  else if (Z_post & !CN_post) {
    return(plot_inference_report_3(x))
  }
  else if (!Z_post & !CN_post) {
    return(plot_inference_report_4(x))
  }
}


plot_inference_report_1 <-  function(x) {
  return(
    cowplot::plot_grid(
        cowplot::plot_grid(
          Rcongas:::plot_loss(x),
          Rcongas:::plot_normalization_factors(x),
          Rcongas::plot_latent_variables(x),
          nrow = 1,
          axis = 'b',
          align = 'h',
          labels = c("a", "b", "c")
        ),
      Rcongas:::plot_CNV_distribution(x),
      nrow = 2,
      labels = c("", "d"),
      rel_heights = c(1, 3)
    )
  )
}

plot_inference_report_2 <-  function(x) {
  return(
    cowplot::plot_grid(
      cowplot::plot_grid(
        Rcongas:::plot_loss(x),
        Rcongas:::plot_normalization_factors(x),
        nrow = 1,
        axis = 'b',
        align = 'h',
        labels = c("a", "b")
      ),
      Rcongas:::plot_CNV_distribution(x),
      nrow = 2,
      labels = c("", "c"),
      rel_heights = c(1, 3)
    )
  )

}

plot_inference_report_3 <-  function(x) {
  cowplot::plot_grid(
    cowplot::plot_grid(
      Rcongas:::plot_loss(x),
      Rcongas:::plot_normalization_factors(x),
      Rcongas::plot_latent_variables(x),
      nrow = 1,
      axis = 'b',
      align = 'h',
      labels = c("a", "b", "c")
    ),
    Rcongas:::plot_CNV_MAP(x),
    nrow = 2,
    labels = c("", "d"),
    rel_heights = c(1, 3)
  )
}

plot_inference_report_4 <-  function(x) {
    cowplot::plot_grid(
      Rcongas:::plot_loss(x),
      Rcongas:::plot_normalization_factors(x),
      Rcongas::plot_latent_variables(x),
      nrow = 1,
      axis = 'b',
      align = 'h',
      labels = c("a", "b", "c")
    )
}


plot_loss <- function(x) {
  best_model <- Rcongas:::get_best_model(x)
  ggplot(data =  best_model$loss  %>% as_tibble, aes(y = value, x = seq_along(value))) +
    geom_line(color = "steelblue") +
    geom_point(size = .4) +
    ggtitle("ELBO") +
    ylab("ELBO") + xlab("VI step") +
    CNAqc:::my_ggplot_theme()
}


#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
plot_normalization_factors <- function(x)
{
  best_model <- Rcongas:::get_best_model(x)

  df = data.frame(
    mu = best_model$parameters$norm_factor,
    cell = names(best_model$parameters$norm_factor),
    stringsAsFactors = FALSE
  ) %>%
    as_tibble() %>%
    left_join(
      Rcongas:::get_clusters(x),
      by = 'cell'
    )

  ggplot(
    data = df,
    aes(x = mu, fill = cluster)
  ) + geom_histogram(bins = 60) +
    ggtitle("Distribution of library size factors") +
    xlab("Library size factor") +
    CNAqc:::my_ggplot_theme() +
    scale_fill_manual(values = Rcongas:::get_clusters_colors(df$cluster))
}




plot_CNV_distribution <-
  function(x,
           density_points = 500,
           xlims = c(0.1, 5)) {
    best_model <- Rcongas:::get_best_model(x)

    cnv_mean <-
      best_model$parameters$cnv_probs %>% as.matrix %>% reshape2::melt()

    if(!is.null(best_model$parameters$cnv_var))
    {
      cnv_sd <- best_model$parameters$cnv_var %>% reshape2::melt()
      cnv_mean$sd <- cnv_sd$value
    }
    else
    {
      message("Did you call the right function?")
    }



    density_matrix <-
      matrix(nrow = nrow(cnv_mean) , ncol = density_points)

    for (i in seq_len(nrow(cnv_mean))) {
      density_matrix[i, ] <-
        dlnorm(
          seq(xlims[1], xlims[2], length.out = density_points),
          log(cnv_mean$value[i]),
          cnv_mean$sd[i]
        )
    }

    colnames(density_matrix) <-
      seq(xlims[1], xlims[2], length.out = density_points)
    density_matrix <-  as.data.frame(density_matrix)
    density_matrix$cluster <- cnv_mean$Var1
    density_matrix$segment <- cnv_mean$Var2

    density_long <-
      reshape2::melt(density_matrix, id.vars = c("cluster", "segment"))

    quant <-
      cnv_mean %>%  dplyr::mutate(
        low = qlnorm(0.1, log(value), sd),
        up = qlnorm(0.90, log(value), sd),
        xlow  = qlnorm(0.0001, log(value), sd),
        xup = qlnorm(0.9999, log(value), sd)
      )
    colnames(quant)[1:3] <- c("cluster", "segment", "vl")
    density_long <- density_long %>% dplyr::inner_join(., quant)

    density_long <-
      density_long %>% dplyr::mutate(quantile = dplyr::case_when(
        as.numeric(as.character(variable)) >= up ~ "up",
        as.numeric(as.character(variable)) <= low ~ "low",
        TRUE ~ "mid"
      ))

    density_long <-
      density_long %>% dplyr::filter(as.numeric(as.character(variable)) >= xlow,
                                     as.numeric(as.character(variable)) <= xup)

    ggplot(density_long,
           aes(
             x = as.numeric(as.character(variable)),
             y = paste(cluster),
             height = value,
             fill = quantile
           )) +
      ggridges::geom_density_ridges_gradient(stat = "identity", size = 0.01) + facet_wrap( ~
                                                                                             segment, scales = "free_y") + xlab("CNV mean") + ylab("cluster") +
      CNAqc:::my_ggplot_theme() + scale_fill_manual(
        name = "Probability",
        values = c("#FF0000A0", "#A0A0A0A0", "#FF0000A0"),
        labels = c("(0, 0.1]", "(0.1, 0.90]", "(0.90, 1]")
      ) +
      labs(title = "Posterior for CNA real values")


  }

plot_CNV_MAP <-
  function(x,
           density_points = 500,
           xlims = c(0.1, 5)) {
    best_model <- Rcongas:::get_best_model(x)

    cnv_mean <-
      best_model$parameters$cnv_probs %>% as.matrix %>% reshape2::melt()
    colnames(cnv_mean) = c('cluster', 'segment', 'MAP')
    cnv_mean$cluster = paste(cnv_mean$cluster)

    ggplot(cnv_mean, aes(x = cluster, y = MAP, ymax = MAP, ymin = 0, color = cluster, fill = MAP)) +
      geom_linerange(size = .4) +
      geom_point(pch = 21, size = 3) +
      scale_color_manual(values = Rcongas:::get_clusters_colors(cnv_mean$cluster)) +
      # scale_color_distiller(palette = 'Spectral', direction = -1) +
      scale_fill_distiller(palette = 'Spectral', direction = -1) +
      CNAqc:::my_ggplot_theme() +
      facet_wrap(~segment) +
      labs(title = 'CNA MAP estimates')
  }

get_baf <- function(vcf) {
  tb1 <- vcf@gt %>% as_tibble()
  colnames(tb1) <-  c("FORMAT", "tmp")
  tb1 <-
    tb1 %>%  tidyr::separate(tmp, into = c("GT", "PL", "DP", "AD"), sep = ":")
  tb2 <-
    tb1 %>% select(AD) %>% tidyr::separate(AD, into = "ref", "alt") %>%  mutate(cov = ref + alt,
                                                                                chr = vcf@fix[, 1],
                                                                                start = vcf@fix[, 2])
}
