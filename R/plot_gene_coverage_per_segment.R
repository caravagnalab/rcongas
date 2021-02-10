#' Title
#'
#' @param x 
#' @param chromosomes 
#'
#' @return
#' @export
#'
#' @examples
plot_coverage_per_segment = function(x, chromosomes = paste0("chr", c(1:22, "X", "Y")), annotate_from = 3)
{
  stopifnot(inherits(x, 'rcongas'))
  
  # Get input raw data in the object
  input_rna = get_input_raw_data(x) %>% as_tibble()
  rn = get_input_raw_data(x) %>% rownames(input_rna) %>% as_tibble() %>% rename(gene = value)
  
  if (nrow(input_rna) == 0)
  {
    cli::cli_alert_danger("RNA counts not found in the input data -- see ?get_input_raw_data.")
    return(ggplot())
  }
  
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
  
  # Data shaping
  df_genes = easypar::run(function(i)
  {
    # Get genes mapped to the input segments
    segment = get_input_segmentation(x)[i,] %>% Rcongas:::idify()
    
    mapped_genes = input_rna %>%
      filter(chr == segment$chr,
             from >= segment$from,
             to <= segment$to)
    nrow(mapped_genes)
    
    # cli::cli_alert(
    #   "n = {.field {nrow(mapped_genes)}} genes mapping to segment {.field {segment$segment_id}}"
    # )
    
    df_counts = mapped_genes %>%
      select(-gene,-chr,-from,-to) %>%
      as_vector()
    
    df_counts = df_counts[df_counts > 0]
    df_counts %>% as_tibble() %>% 
      group_by(value) %>% 
      summarise(n = n(), ln = log(n() + 1)) %>%
      mutate(chr = segment$chr,
             from = segment$from,
             to = segment$to) %>% 
      ungroup()
  },
  lapply(1:nrow(
    get_input_segmentation(x)
  ), list),
  parallel = FALSE)
  
  # Assembly of the df
  df_genes = Reduce(bind_rows, df_genes)
  
  ndf_genes = format(nrow(df_genes), scientific = TRUE)
  # cli::cli_alert_info("Creating barplot from {.field {ndf_genes}} points.")
  
  
  df_genes = df_genes %>% Rcongas:::idify()
  
  # df_genes = df_genes %>%
  #   group_by(segment_id, value) %>%
  #   summarise(ln = log(1 + n()), n = n(), .groups = 'keep')
  
  df_genes = df_genes %>%
    left_join(segment_len %>%
                select(segment_id, size, ploidy_real, mu),
              by = 'segment_id')
  
  # First panel - segments length
  pl_size = ggplot(df_genes) +
    geom_bar(aes(fill = mu, x = segment_id, size / 1e6), stat = 'identity') +
    scale_fill_distiller(palette = "Reds", direction = 1) +
    CNAqc:::my_ggplot_theme() +
    theme(# axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.x = element_blank(),
      legend.position = 'top') +
    labs(x = NULL,
         y = 'Segment size (Mb)') +
    guides(fill = guide_colorbar("Mapped genes", barwidth = unit(4, 'cm')))
  # scale_x_discrete(limits = df_genes$segment_id %>% unique)
  
  # Barplot of counts per segment
  pl_counts = ggplot(df_genes) +
    geom_bar(aes(
      fill = factor(value, levels = rev(1:6)),
      x = segment_id,
      ln,
      color = ploidy_real %>% paste()
    ),
    stat = 'identity') +
    scale_fill_viridis_d(direction = -1) +
    # scale_fill_brewer(palette = "Reds", direction = 1) +
    CNAqc:::my_ggplot_theme() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = NULL) +
    scale_color_brewer(palette = 'Set1')  +
    guides(
      color = guide_legend("Input ploidy"),
      fill = guide_legend(
        "Log of counts prevalence",
        nrow = 1,
        reverse = TRUE
      )
    ) +
    labs(y = "Log of counts")
  
  top_ranges = df_genes %>% 
    filter(value > annotate_from) %>% 
    arrange(desc(value), desc(n)) %>% 
    pull(segment_id)
    
  df_cumsum <- plyr::ddply(df_genes %>%
                             arrange(segment_id, value),
                           "segment_id",
                           transform,
                           label_ypos = cumsum(ln)) %>% 
    filter(value > annotate_from)
  
  # Filter what is not visible
  y_labels = ggplot_build(pl_counts)$layout$panel_params[[1]]$y$get_labels()
  y_labels = y_labels[!is.na(y_labels)] %>% as.numeric()
  delta = (y_labels[2] - y_labels[1])/10
  
  df_cumsum = df_cumsum %>% 
    group_by(segment_id) %>% 
    mutate(offset = label_ypos - lag(label_ypos, default = 0)) %>% 
    filter(offset >= delta)
  
  
  pl_counts = pl_counts +
    geom_text(data = df_cumsum, 
              aes(y=label_ypos, label = n, x = segment_id), 
              vjust = 1.6, 
              color="black",
              size = 2)
    
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
    rel_heights = c(1, 3),
    ncol = 1,
    align = 'v',
    axis = 'lr'
  )
}
