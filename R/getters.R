get_clones_ploidy <-  function(x, chromosomes = paste0("chr", c(1:22, "X", "Y"))) {

  inf_obj = x


  best_model <- get_best_model(inf_obj)
  res <-
    data.frame(t(best_model$parameters$cnv_probs), stringsAsFactors = FALSE)
  res <-
    res %>% mutate(tmp = rownames(res)) %>% tidyr::separate(col = tmp,
                                                            into = c("chr", "from", "to"),
                                                            sep = ":")
  colnames(res)[seq_along(best_model$parameters$mixture_weights)] <-
    paste(seq_along(best_model$parameters$mixture_weights))
  res <-
    reshape2::melt(
      res,
      id.vars = c("chr", "from", "to"),
      variable.name = "cluster",
      value.name = "CN"
    ) %>%
    dplyr::mutate(
      cluster = as.character(cluster),
      from = as.numeric(from),
      to = as.numeric(to)
    ) %>%
    dplyr::as_tibble()

  # Normalise CNA values for comparisons
  means = res %>%
    dplyr::group_by(chr, from, to) %>%
    dplyr::summarise(segment_mean = mean(CN))

  joined_res = res %>%
    dplyr::full_join(means, by = c('chr', 'from', 'to')) %>%
    dplyr::mutate(
      CN = CN - segment_mean
    )

  if(!grepl('chr', joined_res$chr[1]))
  {
    cli::cli_alert_warning("Missing `chr` prefix in chromosomes labels, added now.")

    joined_res = joined_res %>% dplyr::mutate(chr = paste0("chr", chr))
  }

  return(joined_res %>% dplyr::filter(chr %in% chromosomes))
}



get_counts <-  function(inf_obj, input,  chromosomes = paste0("chr", c(1:22, "X", "Y")), norm = T) {
  best_model <- get_best_model(inf_obj)
  if(norm) input$counts <- input$counts / best_model$parameters$norm_factor

  M <-
    long_counts(input$counts) %>%  select(chr, from, to, cell, n) %>%
    dplyr::arrange(chr, from, to)

  clts <- as.data.frame(best_model$parameters$assignement)
  colnames(clts) <- "cluster"
  clts$cell <- rownames(clts)

  ret <-  dplyr::inner_join(M, clts, by = "cell") %>%
    dplyr::select(chr, from, to, cell, n, cluster) %>%
    dplyr::mutate(cluster = paste(cluster))

  if (!grepl('chr', ret$chr[1]))
  {
    cli::cli_alert_warning("Missing `chr` prefix in chromosomes labels, added now.")

    ret <-  ret %>% dplyr::mutate(chr = paste0("chr", chr))
  }


  return(ret %>% filter(chr %in% chromosomes))

}



