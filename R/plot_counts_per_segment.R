#' Title
#'
#' @param x
#' @param chromosomes
#'
#' @return
#' @export
#'
#' @examples
plot_counts_per_segment = function(x,
                                     chromosomes = paste0("chr", c(1:22, "X", "Y")),
                                     annotate_from = 3,
                                     ...)
{
  stopifnot(inherits(x, 'rcongas'))
  
  # Get input raw data in the object
  input_rna = get_input_raw_data(x, ...)
  
  if (nrow(input_rna) == 0)
  {
    cli::cli_alert_danger("RNA counts not found in the input data -- see ?get_input_raw_data.")
    return(ggplot())
  }
  
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
    select(gene, chr, from, to, everything())
  
  # Retain only required "chromosomes"
  input_rna = input_rna %>%
    filter(chr %in% chromosomes)
  
  cli::cli_alert_info(
    "Working with n = {.field {nrow(input_rna)}} genes, and {.field {ncol(input_rna) - 3}} cells; computing cutoffs per segment"
  )
  
  # Creation of the required tibble - specific behaviour if there are clusters
  with_clustering_results = Rcongas:::has_inference(x)
  clusters = get_clusters(x) %>% select(cell, cluster)
  
  df_genes = df_levels = NULL
  if (!with_clustering_results)
  {
    cli::cli_h3("Mapping data for all cells at once (no clusters)\n")
    
    df_genes = aux_plot_histocount_per_segment(x, input_rna) %>%
      mutate(cluster = "No clusters available")
  }
  else
  {
    df_genes = easypar::run(
      function(cl) {
        cli::cli_h3("Mapping data for cluster {.field {cl}}\n")
        
        cell_ids = clusters %>% filter(cluster == cl) %>% pull(cell)
        
        # Retain cluster-specific cells
        aux_plot_histocount_per_segment(x,
                                      input_rna %>%
                                        select(gene, chr, from, to,!!cell_ids)
                                      )  %>%
          mutate(cluster = cl)
      },
      PARAMS = lapply(clusters$cluster %>% unique %>% sort, list),
      parallel = FALSE,
      progress_bar = FALSE
    )
    
    df_genes = Reduce(bind_rows, df_genes)
  }
  
  
  ggplot(
    df_genes,
    aes(counts, fill = cluster)
  ) +
    geom_histogram(bins = 100) +
    facet_wrap(~segment_id, ncol = 6) +
    scale_y_log10() +
    CNAqc:::my_ggplot_theme() +
    scale_fill_brewer(palette = 'Set1')
}



aux_plot_histocount_per_segment = function(x, input_rna)
{
  # Data shaping
  df_genes = easypar::run(function(i)
  {
    # Get genes mapped to the input segments
    segment = get_input_segmentation(x)[i,] %>% Rcongas:::idify()
    
    mapped_genes = input_rna %>%
      filter(chr == segment$chr,
             from >= segment$from,
             to <= segment$to) %>%
      select(-gene,-chr,-from,-to) %>%
      as_vector()
    
    mapped_genes = mapped_genes[mapped_genes > 0]
    
    mapped_genes %>%
      enframe %>%
      select(-name) %>% 
      rename(counts = value) %>% 
      mutate(segment_id = segment$segment_id) 
  },
  lapply(1:nrow(get_input_segmentation(x)), list),
  parallel = FALSE)
  
  # Assembly of the df
  Reduce(bind_rows, df_genes)
}
