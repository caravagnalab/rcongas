#' Extract the clusters ploidy profiles.
#'
#' @description Extract the clone ploidy profiles from a fit object.
#' It can subset by chromosomes, cluster and select what to highlight
#' based on a parameter (alpha)
#'
#' @param x Input object with clusters.
#' @param chromosomes Chromosome id to subset.
#' @param clusters Cluster id to subset.
#' @param offset_amplitude If TRUE, normalise CNA values for comparisons (z-score alike)
#' @param alpha The parameter to select what to highlight.
#'
#' @return A tibble with CN values for each segment, cluster membership and if a segment is highlighted
#' @export
#'
#' @examples
#'
#' x = Rcongas::congas_example
#'
#' # Default view
#' x %>% get_clusters_ploidy()
#'
#' # Subset by chromosome
#' x %>% get_clusters_ploidy(chromosomes = 'chr1')
#'
#' # Subset by cluster id and chromosome
#' x %>% get_clusters_ploidy(chromosomes = 'chr1', clusters = "c1")
#'
#' # Change parameter to find what is most relevant
#' x %>% get_clusters_ploidy(alpha = 0.1)
get_clusters_ploidy <- function(x,
                                chromosomes = paste0("chr", c(1:22, "X", "Y")),
                                clusters = NULL,
                                offset_amplitude = TRUE,
                                alpha = 0.05)
{
  if (!has_inference(x))
    stop("Cannot extract clustering information if not avaiable, or not a CONGAS object.")

  best_model <- Rcongas:::get_best_model(x)
  res <-
    data.frame(t(best_model$parameters$cnv_probs), stringsAsFactors = FALSE)
  res <-
    res %>% dplyr::mutate(tmp = rownames(res)) %>% tidyr::separate(col = tmp,
                                                                   into = c("chr", "from", "to"),
                                                                   sep = ":")
  colnames(res)[seq_along(best_model$parameters$mixture_weights)] <-
    names(best_model$parameters$mixture_weights)

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

  joined_res = res

  # Normalise CNA values for comparisons -- z_score alike 
  means = res %>%
    dplyr::group_by(chr, from, to) %>%
    dplyr::summarise(segment_mean = mean(CN), .groups = 'keep') %>%
    dplyr::ungroup()
  
  joined_res = joined_res %>%
    dplyr::full_join(means, by = c('chr', 'from', 'to')) %>%
    dplyr::mutate(
      # lognorm_CN = CN,
      offset_CN = CN - segment_mean)
  
  joined_res = joined_res %>% idify
  
  # Make default CN to be the offset
  if (offset_amplitude)
    joined_res = joined_res %>% dplyr::mutate(CN = offset_CN)
  
  # Determine highlights
  table_highlight <- highlights(x, alpha) %>% 
    filter(highlight) %>% 
    pull(segment_id) %>% 
    unique
  
  joined_res$highlight = 
    joined_res$segment_id %in% table_highlight
  


  # if (!grepl('chr', joined_res$chr[1]))
  # {
  #   cli::cli_alert_warning("Missing `chr` prefix in chromosomes labels, added now.")
  # 
  #   joined_res = joined_res %>% dplyr::mutate(chr = paste0("chr", chr))
  # }

  # Apply filters
  joined_res = joined_res %>% dplyr::filter(chr %in% chromosomes)

  if (!is.null(clusters))
    joined_res = joined_res %>% dplyr::filter(cluster %in% clusters)


  return(joined_res)
  }