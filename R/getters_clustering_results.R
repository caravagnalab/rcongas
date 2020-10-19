#' Title
#'
#' @param x
#' @param chromosomes
#'
#' @return
#' @export
#'
#' @examples
get_clones_ploidy <- function(x,
                              chromosomes = paste0("chr", c(1:22, "X", "Y")),
                              clusters = NULL)
{
  best_model <- get_best_model(x)
  res <-
    data.frame(t(best_model$parameters$cnv_probs), stringsAsFactors = FALSE)
  res <-
    res %>% dplyr::mutate(tmp = rownames(res)) %>% tidyr::separate(col = tmp,
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
    dplyr::summarise(segment_mean = mean(CN), .groups = 'keep') %>%
    dplyr::ungroup()

  joined_res = res %>%
    dplyr::full_join(means, by = c('chr', 'from', 'to')) %>%
    dplyr::mutate(lognorm_CN = CN,
                  CN = CN - segment_mean)

  if (!grepl('chr', joined_res$chr[1]))
  {
    cli::cli_alert_warning("Missing `chr` prefix in chromosomes labels, added now.")

    joined_res = joined_res %>% dplyr::mutate(chr = paste0("chr", chr))
  }

  # Apply filters
  joined_res = joined_res %>% dplyr::filter(chr %in% chromosomes)

  if (!is.null(clusters))
    joined_res = joined_res %>% dplyr::filter(cluster %in% clusters)


  return(joined_res)
}



#' Title
#'
#' @param x
#' @param normalised
#'
#' @return
#' @export
#'
#' @examples
get_clusters_size = function(x, normalised = FALSE)
{
  best_model = get_best_model(x)

  # Counts from the assignments
  n = table(best_model$parameters$assignement)

  # Make it a named vector
  v = as.vector(n)
  names(v) = names(n)

  if (normalised)
    v = v / sum(v)

  return(v)
}


#' Title
#'
#' @param x
#' @param cluster_label
#'
#' @return
#' @export
#'
#' @examples
get_clusters = function(x,
                        clusters = NULL,
                        cut_znk = 0)
{
  # TODO - take all assingments with z_nk > c, c >= 0
  best_model = Rcongas:::get_best_model(x)

  clusters_table = data.frame(
    cell = names(best_model$parameters$assignement),
    cluster = paste(best_model$parameters$assignement),
    stringsAsFactors = FALSE
  ) %>%
    as_tibble()

  # Latent variables
  z_nk = best_model$parameters$assignment_probs %>%  as.data.frame()
  colnames(z_nk) = Rcongas::get_clusters_size(x) %>% names

  z_nk$p_assignment = unlist(apply(z_nk, 1, function(x) {
    max(x, na.rm = TRUE)
  }))

  z_nk$cell = rownames(z)

  clusters_table = clusters_table %>%
    full_join(z_nk, by = "cell") %>%
    dplyr::mutate(
      cluster = ifelse(p_assignment > cut_znk, cluster, NA)
      ) %>%
    as_tibble()

  if (!is.null(clusters))
    clusters_table = clusters_table %>% dplyr::filter(cluster %in% clusters)

  return(clusters_table)
}



#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
get_k = function(x)
{
  best = get_best_model(x)
  best$parameters$assignement %>%
    unique() %>%
    length()
}


#' Title
#'
#' @param x
#' @param group1
#' @param group2
#' @param cutoff_p
#'
#' @return
#' @export
#'
#' @examples
get_segment_test_counts = function(x, group1, group2, cutoff_p = 0.01)
{
  counts = get_counts(x) %>% Rcongas:::idify() %>% dplyr::group_split(segment_id)
  names_counts = sapply(counts, function(x)
    x$segment_id[1])

  ntests = length(counts)

  # Test every segment
  tests = sapply(counts, function(count_segment) {
    gr1 = count_segment %>% dplyr::filter(cluster %in% !!group1)
    gr2 = count_segment %>% dplyr::filter(cluster %in% !!group2)

    p = wilcox.test(gr1$n, gr2$n)$p.value * ntests
    p = ifelse(p > 1, 1, p)

    p
  })

  data.frame(segment_id = names_counts,
             p = tests,
             stringsAsFactors = FALSE) %>%
    deidify() %>%
    as_tibble() %>%
    arrange(p) %>%
    dplyr::mutate(sign = p < cutoff_p)
}
