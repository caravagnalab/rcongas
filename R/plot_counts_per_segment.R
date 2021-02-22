#' Title
#'
#' @param x
#' @param chromosomes
#' @param genes
#' @param quantiles
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plot_counts_per_segment = function(x,
                                   chromosomes = paste0("chr", c(1:22, "X", "Y")),
                                   genes = NULL,
                                   quantiles = c(0.01, 0.99),
                                   normalize_library_size = FALSE,
                                   ...)
{
  stopifnot(inherits(x, 'rcongas'))
  
  # Get input raw data in the object
  if (has_inference(x))
  {
    input_RNA = get_input_raw_data(x, add_locations = TRUE, add_clusters = TRUE,  ...) 
  }
  else
  {
    input_RNA = get_input_raw_data(x, add_locations = TRUE, add_clusters = FALSE,  ...)
    
    input_RNA$cluster = "Not available"
  }
  
  if (length(quantiles) != 2) {
    cli::cli_alert_info("Quantiles {.field {quantiles}}; should have 2 values, using default.")
    quantiles = c(0.01, .99)
  }
  
  # Normalisation by library size
  if (normalize_library_size)
  {
    cli::cli_alert_info("Normalising for library size..")
    
    input_RNA = input_RNA %>%
      left_join(input_RNA %>%
                  group_by(cell) %>%
                  summarise(lib_size_factor = sum(n)),
                by = 'cell') %>% 
      mutate(n = n/lib_size_factor)
    
    # cell_ids = input_rna %>% colnames
    # 
    # for (cell in cell_ids)
    # {
    #   input_rna[, cell] = input_rna[, cell] / sum(input_rna[, cell])
    # }
  }
  
  if (!is.null(genes))
  {
    ng = input_rna %>% pull(gene) %>% unique %>% length()
    
    input_rna = input_rna %>% filter(gene %in% !!genes)
    
    n2g = input_rna %>% pull(gene) %>% unique %>% length()
    
    cli::cli_alert_info(
      "Subsetting {.field {ng}} available genes leaves {.field {n2g}} genes (out of {.field {length(genes)}} in the list)."
    )
  }
  
  # Used genes
  used_genes = input_RNA %>% pull(gene) %>% unique %>% length
  
  if (nrow(input_RNA) == 0)
  {
    cli::cli_alert_danger("No RNA counts in the data, returning empty plot.")
    return(ggplot())
  }
  
  # Retain only required "chromosomes"
  input_RNA = input_RNA %>%
    filter(chr %in% chromosomes)
  
  
  # # Unify gene names if required
  # unique_gene_names = rownames(input_RNA) %>% make.unique()
  # nrow(input_rna) == length(unique_gene_names)
  # rownames(input_rna) = unique_gene_names
  # 
  # # Tibble it
  # input_rna = input_rna %>% as_tibble()
  # rn = unique_gene_names %>% as_tibble() %>% rename(gene = value)
  # 
  # input_rna = input_rna %>%
  #   bind_cols(rn) %>%
  #   left_join(get_gene_annotations(x), by = 'gene') %>%
  #   select(gene, chr, from, to, everything())
  # 
  
  

  
  # cli::cli_alert_info(
  #   "Working with n = {.field {nrow(input_rna)}} genes, and {.field {ncol(input_rna) - 3}} cells; computing cutoffs per segment"
  # )
  
  # Creation of the required tibble - specific behaviour if there are clusters
  # with_clustering_results = Rcongas:::has_inference(x)
  
  
  #   
  # df_genes = df_levels = NULL
  # if (!with_clustering_results)
  # {
  #   cli::cli_h3("Mapping data for all cells at once (no clusters)\n")
  #   
  #   df_genes = aux_plot_histocount_per_segment(x, input_rna) %>%
  #     mutate(cluster = "No clusters available")
  # }
  # else
  # {
  #   df_genes = easypar::run(
  #     function(cl) {
  #       cli::cli_h3("Mapping data for cluster {.field {cl}}\n")
  #       
  #       cell_ids = clusters %>% filter(cluster == cl) %>% pull(cell)
  #       
  #       # Retain cluster-specific cells
  #       aux_plot_histocount_per_segment(x,
  #                                       input_rna %>%
  #                                         select(gene, chr, from, to, !!cell_ids))  %>%
  #         mutate(cluster = cl)
  #     },
  #     PARAMS = lapply(clusters$cluster %>% unique %>% sort, list),
  #     parallel = FALSE,
  #     progress_bar = FALSE
  #   )
  #   
  #   df_genes = Reduce(bind_rows, df_genes)
  # }
  
  # Segment size
  sz_lb = round(get_input_segmentation(x, chromosomes = chromosomes)$size /
                  1e6)
  pl_lb = get_input_segmentation(x, chromosomes = chromosomes)$ploidy_real
  
  names(sz_lb) = names(pl_lb) = get_input_segmentation(x, chromosomes = chromosomes) %>%
    Rcongas:::idify() %>%
    pull(segment_id)
  
  # Custom breaks
  base_breaks <- function(n = 1) {
    function(x) {
      # print(range(1 + x, na.rm = TRUE))
      # print(      axisTicks(log10(range(1 + x, na.rm = TRUE)), log = TRUE, n = n))
      v = axisTicks(log10(range(1 + x, na.rm = TRUE)), log = TRUE, n = n)
      if (length(v) > 3)
        v = v[seq(1, length(v), 2)]
      if (length(v) > 3)
        v = v[2:length(v)]
      
      v
    }
  }
  
  # Add ploidy information
  input_RNA = input_RNA  %>%
    left_join(
      get_input_segmentation(x, chromosomes = chromosomes) %>%
        Rcongas:::idify() %>%
        select(segment_id, size, ploidy_real) %>%
        mutate(ploidy_real = paste0('Ploidy ', ploidy_real,
                                    " (", round(size/1e6), " Mb)")),
      by = "segment_id"
    ) 
  
  # Computing quantiles per segment
  quantiles_df = input_RNA %>%
    group_by(segment_id) %>%
    mutate(q_1 = quantile(n, quantiles[1]),
           q_2 = quantile(n, quantiles[2])) %>%
    distinct(segment_id, ploidy_real, q_1, q_2)
  
  # Plot
  ggplot(input_RNA,
         aes(n, fill = cluster)) +
    geom_histogram(bins = 100) +
    facet_wrap(ploidy_real ~ segment_id, ncol = 6, scales = 'free_x') +
    scale_y_continuous(trans = scales::log1p_trans(),
                       breaks = base_breaks(4)) +
    CNAqc:::my_ggplot_theme() +
    scale_fill_brewer(palette = 'Set1')  +
    labs(
      x = "Counts",
      title = paste0("Log counts per segment"),
      subtitle = paste0("Quantiles ", paste(quantiles, collapse = ','), " annotated."),
      caption = ifelse(
        !is.null(genes),
        paste0("Genes: all available input (n = ", used_genes, ')'),
        paste0("Genes: custom list (n =", used_genes, ')')
      )
    ) +
    geom_vline(
      data = quantiles_df,
      aes(xintercept = q_1),
      color = 'black',
      linetype = 'dashed',
      size = .3
    ) +
    geom_vline(
      data = quantiles_df,
      aes(xintercept = q_2),
      color = 'black',
      linetype = 'dashed',
      size = .3
    )
}


aux_plot_histocount_per_segment = function(x, input_rna)
{
  # Data shaping
  df_genes = easypar::run(function(i)
  {
    # Get genes mapped to the input segments
    segment = get_input_segmentation(x)[i, ] %>% Rcongas:::idify()
    
    mapped_genes = input_rna %>%
      filter(chr == segment$chr,
             from >= segment$from,
             to <= segment$to) %>%
      select(-gene, -chr, -from, -to) %>%
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
