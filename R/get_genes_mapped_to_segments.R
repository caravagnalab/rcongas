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
  input_rna = get_input_raw_data(x, ...)
  
  # Unify gene names if required
  unique_gene_names = rownames(input_rna) %>% make.unique()
  nrow(input_rna) == length(unique_gene_names)
  rownames(input_rna) = unique_gene_names
  
  # Tibble it
  input_rna = input_rna %>% as_tibble()
  rn = unique_gene_names %>% as_tibble() %>% rename(gene = value)
  
  input_rna = input_rna %>%
    bind_cols(rn) %>%
    left_join(get_gene_annotations(x), by = 'gene') %>%
    select(gene, chr, from, to)
  
  if(
    input_rna[complete.cases(input_rna), ] %>% nrow <
    input_rna %>% nrow
  )
    warning("get_gene_annotations() is missing annotations for some input genes ...")
  
  # Per segment
  df_genes = easypar::run(function(i)
  {
    # Get genes mapped to the input segments
    segment = get_input_segmentation(x)[i,] %>% Rcongas:::idify()

    mapped_genes = input_rna %>%
      filter(chr == segment$chr,
             from >= segment$from,
             to <= segment$to)
    
    mapped_genes %>% 
      bind_cols(
        segment %>% select(-start, -end, dist, -chr, -from, -to)
      )
  },
  lapply(1:nrow(get_input_segmentation(x)), list),
  parallel = FALSE)
  
  Reduce(bind_rows, df_genes)
}