plot_inference_report <-  function(x){

  best_model <- get_best_model(x)


}


plot_inference_report_1 <-  function() {

}

plot_inference_report_2 <-  function() {

}

plot_inference_report_3 <-  function() {

}

plot_inference_report_4 <-  function() {

}


plot_loss <- function(x){

  best_model <- get_best_model(x)
  ggplot(data =  best_model$loss  %>% as_tibble,aes(y = value, x = seq_along(value))) + geom_line(color = "blue") +
    ggtitle("ELBO convergence") + ylab("ELBO") + xlab("Step") + CNAqc:::my_ggplot_theme()
}


plot_normalization_factors <- function(x){
  best_model <- get_best_model(x)
  ggplot(data = best_model$parameters$norm_factor  %>% as_tibble %>% dplyr::mutate(clust = get_cluster_assignments(x)), aes(x = value, fill = clust)) + geom_histogram(bins = 60) +
    ggtitle("Distribution of library size factors") + xlab("Library size factor") + CNAqc:::my_ggplot_theme()

}


plot_posterior_probabilities <- function(x) {
  best_model <- get_best_model(x)
  posteriors <- get_assignment_probs(x) %>% reshape2::melt() %>% arrange(Var2, value)

  ggplot(data = posteriors, aes(x = factor(Var1), y = paste(Var2), fill = value)) +
    geom_tile() + ggtitle("Posterior assignment probability") + xlab("cells")
    CNAqc:::my_ggplot_theme() + theme(axis.title.x=element_blank(),
                                              axis.text.x=element_blank(),
                                              axis.ticks.x=element_blank())

}


plot_CNV_distribution <- function(x, density_points = 500, xlims = c(0.1,5)) {

  best_model <- get_best_model(x)

  cnv_mean <-  best_model$parameters$cnv_probs %>% as.matrix %>% reshape2::melt()

  cnv_sd <- best_model$parameters$cnv_var %>% reshape2::melt()


  cnv_mean$sd <- cnv_sd$value



  density_matrix <- matrix(nrow = nrow(cnv_mean) , ncol = density_points)

  for(i in seq_len(nrow(cnv_mean))) {

    density_matrix[i,] <- dlnorm(seq(xlims[1], xlims[2], length.out = density_points), log(cnv_mean$value[i]), cnv_mean$sd[i])
  }

  colnames(density_matrix) <- seq(xlims[1], xlims[2], length.out = density_points)
  density_matrix <-  as.data.frame(density_matrix)
  density_matrix$cluster <- cnv_mean$Var1
  density_matrix$segment <- cnv_mean$Var2

  density_long <- reshape2::melt(density_matrix, id.vars = c("cluster", "segment"))

  quant <- cnv_mean %>%  dplyr::mutate(low = qlnorm(0.1, log(value), sd), up = qlnorm(0.90, log(value), sd),
                                xlow  = qlnorm(0.0001, log(value), sd), xup = qlnorm(0.9999, log(value), sd) )
  colnames(quant)[1:3] <- c("cluster", "segment", "vl")
  density_long <- density_long %>% dplyr::inner_join(., quant)

  density_long <-  density_long %>% dplyr::mutate(quantile = dplyr::case_when(as.numeric(as.character(variable)) >= up ~ "up", as.numeric(as.character(variable)) <= low ~ "low", TRUE ~ "mid") )

  density_long <-  density_long %>% dplyr::filter(as.numeric(as.character(variable)) >= xlow, as.numeric(as.character(variable)) <= xup)

  ggplot(density_long, aes(x = as.numeric(as.character(variable)), y = paste(cluster), height = value, fill = quantile)) +
    ggridges::geom_density_ridges_gradient(stat = "identity", size = 0.01) + facet_wrap(~segment, scales = "free_x") + xlab("CNV mean") + ylab("cluster") +
     CNAqc:::my_ggplot_theme() + scale_fill_manual(
      name = "Probability", values = c("#FF0000A0", "#A0A0A0A0", "#FF0000A0"),
      labels = c("(0, 0.1]", "(0.1, 0.90]", "(0.90, 1]")
    )


}

get_baf <- function(vcf) {
  tb1 <- vcf@gt %>% as_tibble()
  colnames(tb1) <-  c("FORMAT", "tmp")
  tb1 <- tb1 %>%  tidyr::separate(tmp, into = c("GT", "PL", "DP", "AD"), sep = ":")
  tb2 <- tb1 %>% select(AD) %>% tidyr::separate(AD, into = "ref", "alt") %>%  mutate(cov = ref + alt, chr = vcf@fix[,1], start = vcf@fix[,2])
}
