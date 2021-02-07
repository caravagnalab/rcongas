#' Wilcoxon test for difference in RNA counts among groups.
#'
#' @description Performs a \code{wilcox.test} for the counts in each segments,
#' split by pair of groups. A group can be a list of cluster ids; the test
#' p-value can be passed as a parameter, and is adjusted for Family-Wise Error
#' Rate (FWER) via Bonferroni correction.
#'
#' @param x Input object with clusters.
#' @param group1 Group of cluster ids (vs \code{group2}).
#' @param group2 Group of cluster ids (vs \code{group1}).
#' @param cutoff_p P-value cutoff
#'
#' @return A tibble with segment coordinates, p-values and a \code{sign} status
#' for the significance of the test.
#' @export
#'
#' @examples
#'
#' x = Rcongas::congas_example
#'
#' # Compare groups
#' get_segment_test_counts(x, "c1", "c2")
#'
#' # More stringent test
#' get_segment_test_counts(x, "c1", "c2", cutoff_p = 1e-4)
get_segment_test_counts = function(x, group1, group2, cutoff_p = 0.01)
{
  if(!has_inference(x)) stop("Cannot extract clustering information if not avaiable, or not a CONGAS object.")
  if(cutoff_p < 0) stop("Negative p-value")
  if(cutoff_p > 1) stop(">1 p-value")

  # Type of model used
  # - MixtureGaussian
  # - MixtureGaussianNorm
  # - HmmSegmenter
  type_of_model = Rcongas:::get_congas_model_used(x)

  if (type_of_model == "HmmSegmenter")
    stop("Cannot do any test for the HMM segementer..")

  # Check input
  cl_labels = get_clusters_ploidy(x)$cluster %>% unique()

  if(!all(group1 %in% cl_labels)) stop("Cluster labels in group1 are not correct.")
  if(!all(group2 %in% cl_labels)) stop("Cluster labels in group2 are not correct.")
  if(length(intersect(group1, group2)) > 0) stop("Some cluster labels appear in both groups.")

  # ..
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