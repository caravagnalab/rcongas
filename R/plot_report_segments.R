#' Title
#'
#' @param x
#' @param chromosomes
#'
#' @return
#' @export
#'
#' @examples
plot_report_segments = function(x,
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
  
  # Computation of cutoffs
  cutoffs_table = aux_compute_cutoffs_per_segment(x, input_rna)
  
  # Creation of the required tibble - specific behaviour if there are clusters
  with_clustering_results = Rcongas:::has_inference(x)
  clusters = get_clusters(x) %>% select(cell, cluster)
  
  df_genes = df_levels = NULL
  if (!with_clustering_results)
  {
    cli::cli_h3("Mapping data for all cells at once (no clusters)\n")
    
    df_genes = aux_plot_coverage_per_segment(x, input_rna) %>% 
      mutate(cluster = "No clusters available")
  }
  else
  {
    df_genes = easypar::run(
      function(cl) {
        cli::cli_h3("Mapping data for cluster {.field {cl}}\n")
        
        cell_ids = clusters %>% filter(cluster == cl) %>% pull(cell)
        
        # Retain cluster-specific cells
        aux_plot_coverage_per_segment(x,
                                      input_rna %>%
                                        select(gene, chr, from, to, !!cell_ids),
                                      cutoffs_table = cutoffs_table)  %>%
          mutate(cluster = cl)
      },
      PARAMS = lapply(clusters$cluster %>% unique %>% sort, list),
      parallel = FALSE,
      progress_bar = FALSE
    )
    
    df_genes = Reduce(bind_rows, df_genes)
  }
  
  # First panel - segments length
  pl_size = aux_plot_panel_segments(df_genes)
  
  # Barplot of counts per segment
  pl_counts = aux_plot_panel_counts(df_genes)
  
  # Inform about the maximas
  pl_maxima = aux_plot_panel_maxima(df_genes)
  
  cowplot::plot_grid(
    pl_size + labs(title = 'Coverage per gene per segment '),
    pl_maxima,
    pl_counts,
    rel_heights = c(1, 1, 3),
    ncol = 1,
    align = 'v',
    axis = 'lr'
  )
  
  # cowplot::plot_grid(
  #   pl_size + labs(title = 'Coverage per gene per segment '),
  #   pl_maxima,
  #   # pl_counts,
  #   # rel_heights = c(1, 1, 3),
  #   ncol = 1,
  #   align = 'v',
  #   axis = 'lr'
  # )
  
}


aux_compute_cutoffs_per_segment = function(x, input_rna)
{
  # Data shaping
  df_genes = easypar::run(function(i)
  {
    # Get genes mapped to the input segments
    segment = get_input_segmentation(x)[i, ] %>% Rcongas:::idify()
    
    mapped_genes = input_rna %>%
      filter(chr == segment$chr,
             from >= segment$from,
             to <= segment$to)
    
    df_counts = mapped_genes %>%
      select(-gene, -chr, -from, -to) %>%
      as_vector()
    
    df_counts = df_counts[df_counts > 0]
    
    min_df_counts = min(df_counts)
    max_df_counts = max(df_counts)
    
    seq(min_df_counts,
        max_df_counts,
        (max_df_counts - min_df_counts) / 10) %>%
      as_tibble() %>%
      rename(min_range = value) %>%
      mutate(level = row_number(), segment_id = segment$segment_id)
  },
  lapply(1:nrow(get_input_segmentation(x)), list),
  parallel = FALSE)
  
  Reduce(bind_rows, df_genes)
}


aux_plot_coverage_per_segment = function(x, input_rna, cutoffs_table)
{
  sid = get_input_segmentation(x) %>% Rcongas:::idify() %>% pull(segment_id)
  
  # Data shaping
  df_genes = easypar::run(function(i)
  {
    # Get genes mapped to the input segments
    segment = get_input_segmentation(x)[i, ] %>% Rcongas:::idify()
    cutoffs_table_segment = cutoffs_table %>% filter(segment_id == sid[i]) %>%
      rename(value = level)
    
    mapped_genes = input_rna %>%
      filter(chr == segment$chr,
             from >= segment$from,
             to <= segment$to)
    
    df_counts = mapped_genes %>%
      select(-gene, -chr, -from, -to) %>%
      as_vector()
    
    df_counts = df_counts[df_counts > 0]
    
    cut_df_counts = cut(df_counts,
                        breaks = cutoffs_table_segment$min_range,
                        include.lowest = TRUE)
    ncut_df_counts = cut_df_counts %>% as.numeric()
    
    ncut_df_counts %>%
      as_tibble() %>%
      group_by(value) %>%
      summarise(n = n(), ln = log(n() + 1)) %>%
      mutate(chr = segment$chr,
             from = segment$from,
             to = segment$to) %>%
      Rcongas:::idify() %>%
      ungroup() %>%
      left_join(cutoffs_table_segment, by = c("segment_id", "value"))
    
  },
  lapply(1:nrow(get_input_segmentation(x)), list),
  parallel = FALSE)
  
  # Assembly of the df
  df_genes = Reduce(bind_rows, df_genes)
  
  # Segment metadata
  sm = df_genes %>%
    left_join(
      get_input_segmentation(x) %>%
        Rcongas:::idify() %>%
        select(segment_id, size, ploidy_real, mu),
      by = 'segment_id'
    )
  
  return(sm)
}

aux_plot_panel_segments = function(df_genes)
{
  df_genes_x = df_genes %>% distinct(segment_id, .keep_all = TRUE)
  
  ggplot(df_genes_x) +
    geom_bar(aes(
      fill = mu,
      x = factor(segment_id,
                 levels = gtools::mixedsort(segment_id) %>% unique()),
      size / 1e6
    ),
    stat = 'identity') +
    scale_fill_distiller(palette = "Reds", direction = 1) +
    CNAqc:::my_ggplot_theme() +
    theme(# axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.x = element_blank(),
      legend.position = 'right') +
    labs(x = NULL,
         y = 'Segment size (Mb)') +
    guides(fill = guide_colorbar("Mapped genes", barheight = unit(4, 'cm'))) +
    facet_grid(. ~ ploidy_real, scales = "free", space = 'free')
}

aux_plot_panel_counts = function(df_genes)
{
  df_genes_x = df_genes[complete.cases(df_genes), ]
  m_value = df_genes$value %>% unique() %>% rev
  
  # Barplot of counts per segment
  ggplot(df_genes_x) +
    geom_bar(aes(
      fill = factor(value, levels = m_value),
      x = factor(segment_id,
                 levels = gtools::mixedsort(segment_id) %>% unique()),
      ln
    ),
    # color = ploidy_real %>% paste()),),),
    stat = 'identity') +
    scale_fill_viridis_d(direction = -1) +
    CNAqc:::my_ggplot_theme() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = NULL) +
    guides(# color = guide_legend(
      #   "Input ploidy",
      #   override.aes = aes(fill = NA),
      #   nrow = 1
      # ),
      fill = guide_legend(
        "Log of counts prevalence",
        nrow = 1,
        reverse = TRUE
      )) +
    labs(y = "Log of counts") +
    facet_grid(cluster ~ ploidy_real, scales = "free", space = 'free')
}

aux_plot_panel_maxima = function(df_genes)
{
  df_genes_x = df_genes %>%
    group_by(segment_id, cluster) %>%
    arrange(desc(value)) %>%
    filter(row_number() == 1)
  
  df_genes_x %>%
    ggplot() +
    CNAqc:::my_ggplot_theme() +
    geom_bar(aes(
      x = interaction(
        cluster,
        factor(segment_id,
               levels = gtools::mixedsort(segment_id) %>% unique())
      ),
      y = min_range,
      fill = cluster
    ), stat = 'identity')  +
    # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(axis.text.x = element_blank(), legend.position = 'right') +
    labs(x = NULL,
         y = 'Maximum counts per cluster') +
    scale_fill_brewer(palette = 'Set1') +
    guides(fill = guide_legend("Cluster", ncol = 1)) +
    facet_grid(. ~ ploidy_real, scales = "free", space = 'free')
}
