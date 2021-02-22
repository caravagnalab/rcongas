#' Title
#'
#' @param x
#' @param chromosomes
#'
#' @return
#' @export
#'
#' @examples
plot_coverage_per_segment = function(x,
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
    "Working with n = {.field {nrow(input_rna)}} genes, and {.field {ncol(input_rna) - 3}} cells."
  )
  
  # Creation of the required tibble - specific behaviour if there are clusters
  with_clustering_results = Rcongas:::has_inference(x)
  clusters = get_clusters(x) %>% select(cell, cluster)
  
  df_genes = df_levels = NULL
  if (!with_clustering_results)
  {
    cli::cli_h3("Mapping data for all cells at once (no clusters)\n")
    
    df_genes = aux_plot_coverage_per_segment(x, input_rna)
    df_genes = df_genes[[2]] %>%
      rename(ranges = value, value = level) %>%
      left_join(df_genes[[1]], by = c("segment_id", 'value')) %>%
      mutate(cluster = "No clusters available")
  }
  else
  {
    df_genes = easypar::run(
      function(cl) {
        cli::cli_h3("Mapping data for cluster {.field {cl}}\n")
        
        cell_ids = clusters %>% filter(cluster == cl) %>% pull(cell)
        
        # Retain cluster-specific cells
        xx = aux_plot_coverage_per_segment(x,
                                           input_rna %>%
                                             select(gene, chr, from, to,!!cell_ids))
        
        xx[[2]] %>%
          rename(ranges = value, value = level) %>%
          left_join(xx[[1]], by = c("segment_id", 'value')) %>%
          mutate(cluster = cl)
        
        # df_levels = df_genes[[2]] %>%
        #   mutate(cluster = "No clusters available")
        # df_genes = df_genes[[1]] %>%
        #   mutate(cluster = "No clusters available")
        #
        #
        # %>%
        #   mutate(cluster = cl)
        
        
      },
      PARAMS = lapply(clusters$cluster %>% unique %>% sort, list),
      parallel = FALSE,
      progress_bar = FALSE
    )
    
    df_genes = Reduce(bind_rows, df_genes)
  }
  
  
  
  #
  #
  # # Data shaping
  # df_genes = easypar::run(function(i)
  # {
  #   # Get genes mapped to the input segments
  #   segment = get_input_segmentation(x)[i, ] %>% Rcongas:::idify()
  #
  #   mapped_genes = input_rna %>%
  #     filter(chr == segment$chr,
  #            from >= segment$from,
  #            to <= segment$to)
  #   nrow(mapped_genes)
  #
  #   # cli::cli_alert(
  #   #   "n = {.field {nrow(mapped_genes)}} genes mapping to segment {.field {segment$segment_id}}"
  #   # )
  #
  #   df_counts = mapped_genes %>%
  #     select(-gene, -chr, -from, -to) %>%
  #     as_vector()
  #
  #
  #   mapped_genes %>%
  #     select(-gene, -chr, -from, -to) %>%
  #     gather(key = var_name, value = value, 2:ncol(df)) %>%
  #     spread_(key = names(df)[1], value = 'value')
  #
  #
  #
  #   mapped_genes %>%
  #
  #
  #     df_counts = df_counts[df_counts > 0]
  #   cell_ids = df_counts %>% names
  #
  #   df_counts = df_counts %>% as_tibble()
  #   df_counts$cell = cell_ids
  #
  #   # Add a cluster label
  #   if (with_clustering_results)
  #   {
  #     df_counts %>%
  #       left_join(clusters)
  #
  #   }
  #
  #   %>%
  #     group_by(value) %>%
  #     summarise(n = n(), ln = log(n() + 1)) %>%
  #     mutate(chr = segment$chr,
  #            from = segment$from,
  #            to = segment$to) %>%
  #     ungroup()
  # },
  # lapply(1:nrow(get_input_segmentation(x)), list),
  # parallel = FALSE)
  #
  # # Assembly of the df
  # df_genes = Reduce(bind_rows, df_genes)
  #
  # ndf_genes = format(nrow(df_genes), scientific = TRUE)
  # # cli::cli_alert_info("Creating barplot from {.field {ndf_genes}} points.")
  #
  #
  # df_genes = df_genes %>% Rcongas:::idify()
  #
  # # df_genes = df_genes %>%
  # #   group_by(segment_id, value) %>%
  # #   summarise(ln = log(1 + n()), n = n(), .groups = 'keep')
  #
  # df_genes = df_genes %>%
  #   left_join(segment_len %>%
  #               select(segment_id, size, ploidy_real, mu),
  #             by = 'segment_id')
  
  # First panel - segments length
  pl_size = ggplot(df_genes) +
    geom_bar(aes(fill = mu, x = segment_id, size / 1e6), stat = 'identity') +
    scale_fill_distiller(palette = "Reds", direction = 1) +
    CNAqc:::my_ggplot_theme() +
    theme(# axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.x = element_blank(),
      legend.position = 'right') +
    labs(x = NULL,
         y = 'Segment size (Mb)') +
    guides(fill = guide_colorbar("Mapped genes", barheight = unit(4, 'cm')))
  # scale_x_discrete(limits = df_genes$segment_id %>% unique)
  
  # Barplot of counts per segment
  m_value = df_genes$value %>% unique() %>% rev
  pl_counts = ggplot(df_genes[complete.cases(df_genes),]) +
    geom_bar(aes(fill = factor(value, levels = m_value),
                 x = segment_id,
                 ln),
             # color = ploidy_real %>% paste()),),),
             stat = 'identity') +
    scale_fill_viridis_d(direction = -1) +
    # scale_fill_brewer(palette = "Reds", direction = 1) +
    CNAqc:::my_ggplot_theme() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = NULL) +
    # scale_color_brewer(palette = 'Set2')  +
    ggsci::scale_color_igv() +
    guides(
      color = guide_legend(
        "Input ploidy",
        override.aes = aes(fill = NA),
        nrow = 1
      ),
      fill = guide_legend(
        "Log of counts prevalence",
        nrow = 1,
        reverse = TRUE
      )
    ) +
    labs(y = "Log of counts") +
    facet_grid(cluster ~ ploidy_real, scales = "free", space = 'free')
  
  # pl_counts
  
  # pl_counts + ggforce::facet_row(vars(ploidy_real), scales = 'free', space = 'free')
  
  top_ranges = df_genes %>%
    filter(value > annotate_from) %>%
    arrange(desc(value), desc(n)) %>%
    pull(segment_id)
  
  df_cumsum <- plyr::ddply(
    df_genes %>%
      arrange(cluster, segment_id, value),
    c("segment_id", "cluster"),
    transform,
    label_ypos = cumsum(ln)
  ) %>%
    filter(value > annotate_from)
  
  # Filter what is not visible
  y_labels = ggplot_build(pl_counts)$layout$panel_params[[1]]$y$get_labels()
  y_labels = y_labels[!is.na(y_labels)] %>% as.numeric()
  delta = (y_labels[2] - y_labels[1]) / 10
  
  df_cumsum = df_cumsum %>%
    group_by(segment_id, cluster) %>%
    mutate(offset = label_ypos - lag(label_ypos, default = 0)) %>%
    filter(offset >= delta)
  
  # pl_counts = pl_counts +
  #   geom_text(
  #     data = df_cumsum,
  #     aes(y = label_ypos, label = n, x = segment_id),
  #     vjust = 1.6,
  #     color = "black",
  #     size = 2
  #   )
  
  # Inform about the maximas
  pl_maxima = df_genes %>%
    group_by(segment_id, cluster) %>%
    arrange(desc(value)) %>%
    filter(row_number() == 1) %>%
    separate(ranges, into  = c('m_1', 'm_2'), sep = ',') %>%
    mutate(m_2 = as.numeric(gsub(']', '', m_2))) %>%
    group_by(segment_id) %>% 
    filter(segment_id == "chr18:55950001:78150000")
    
    mutate(m_2 = m_2/sum(m_2)) %>% 
    
    
    ggplot() +
    CNAqc:::my_ggplot_theme() +
    geom_bar(aes(x = segment_id, y = m_2, fill = cluster), stat = 'identity')  +
    theme(axis.text.x = element_blank(), legend.position = 'right') +
    labs(x = NULL,
         y = '% counts per clone') +
    scale_fill_brewer(palette = 'Set1') +
    guides(fill = guide_legend("Cluster", ncol = 1)) +
    facet_grid(. ~ ploidy_real, scales = "free", space = 'free')
  
  
  # y_labels = ggplot_build(pl_counts)$layout$panel_params[[1]]$y$get_labels()
  # y_labels = y_labels[!is.na(y_labels)] %>% as.numeric()
  
  # pl_counts +
  #   scale_y_continuous(
  #     "Log of counts",
  #     breaks = y_labels,
  #     sec.axis = sec_axis(~ ., name = "mpg (UK)", labels = exp(y_labels))
  #   )
  
  cowplot::plot_grid(
    pl_size + labs(title = 'Coverage per gene per segment '),
    pl_counts,
    pl_maxima,
    rel_heights = c(1, 3, 1),
    ncol = 1,
    align = 'v',
    axis = 'lr'
  )
}


aux_plot_coverage_per_segment = function(x, input_rna)
{
  sid = get_input_segmentation(x) %>% Rcongas:::idify() %>% pull(segment_id)
  
  # Data shaping
  df_genes = easypar::run(function(i)
  {
    # Get genes mapped to the input segments
    segment = get_input_segmentation(x)[i,]
    
    mapped_genes = input_rna %>%
      filter(chr == segment$chr,
             from >= segment$from,
             to <= segment$to)
    nrow(mapped_genes)
    
    df_counts = mapped_genes %>%
      select(-gene,-chr,-from,-to) %>%
      as_vector()
    
    df_counts = df_counts[df_counts > 0]
    
    # hist(df_counts)
    # quantile(df_counts, seq(0, 1, .1))
    cut_df_counts = cut(df_counts, breaks = 10)
    ncut_df_counts = cut_df_counts %>% as.numeric()
    
    r = ncut_df_counts %>%
      as_tibble() %>%
      group_by(value) %>%
      summarise(n = n(), ln = log(n() + 1)) %>%
      mutate(chr = segment$chr,
             from = segment$from,
             to = segment$to) %>%
      ungroup()
    
    return(list(r, levels(cut_df_counts)))
  },
  lapply(1:nrow(get_input_segmentation(x)), list),
  parallel = FALSE)
  
  # Assembly of the df for levels
  df_levels = Reduce(bind_rows,
                     lapply(seq(df_genes), function(x)
                       df_genes[[x]][[2]] %>% as_tibble() %>%
                         mutate(segment_id = sid[x],
                                level = row_number())))
  
  # Assembly of the df for counts
  df_counts = Reduce(bind_rows,
                     lapply(seq(df_genes), function(x)
                       df_genes[[x]][[1]]))
  
  df_genes = df_counts
  
  # df_genes = Reduce(bind_rows, df_genes)
  df_genes = df_genes %>% Rcongas:::idify()
  
  # Segment metadata
  sm = df_genes %>%
    left_join(
      get_input_segmentation(x) %>%
        Rcongas:::idify() %>%
        select(segment_id, size, ploidy_real, mu),
      by = 'segment_id'
    )
  
  return(list(sm, df_levels))
}
