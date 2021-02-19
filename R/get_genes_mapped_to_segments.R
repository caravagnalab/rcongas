#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
get_genes_mapped_to_segments = function(x, ...)
{
  stopifnot(inherits(x, 'rcongas'))

  # Get input raw data in the object
  input_rna = get_input_raw_data(x, as_matrix = FALSE, ...) %>%
    left_join(x %>% get_mapped_genes(), by = 'gene')

  if(
    input_rna[is.na(input_rna$chr), ] %>% nrow > 0
  )
    warning("Do not find a map for some of the genes - this should not happen.")

  input_rna
}