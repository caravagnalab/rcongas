#' Title
#'
#' @param x
#' @param chromosomes
#' @param top
#' @param ...
#' @param by_segment
#'
#' @return
#' @export
#'
#' @examples
plot_gene_rank = function(x,
                          fun = sd,
                          chromosomes = paste0("chr", c(1:22, "X", "Y")),
                          top = 300,
                          by_segment = FALSE,
                          assembly = TRUE,
                          ...)
{
  if (has_inference(x))
  {
    input_RNA = get_input_raw_data(x, add_locations = TRUE, add_clusters = TRUE,  ...) %>%
      filter(chr %in% chromosomes)
  }
  else
  {
    input_RNA = get_input_raw_data(x, add_locations = TRUE, add_clusters = FALSE,  ...) %>%
      filter(chr %in% chromosomes)
    
    input_RNA$cluster = "Not available"
  }
  
  #
  #
  # # Gene maps
  # genes_map = get_genes_mapped_to_segments(x,
  #                                          only_clustered_cells = TRUE,
  #                                          only_mapped_genes = TRUE,
  #                                          ...) %>%
  #   filter(chr %in% chromosomes)
  #
  # # Map to CNAs
  # mapped_coordinates_cnas = get_genes_mapped_to_segments(x, ...) %>%
  #   select(gene, segment_id)
  #
  # # Input counts
  # input_RNA = get_input_raw_data(
  #   x,
  #   as_tibble = TRUE,
  #   only_clustered_cells = TRUE,
  #   only_mapped_genes = TRUE,
  #   ...
  # ) %>%
  #   left_join(get_clusters(x) %>%
  #               select(cell, cluster),
  #             by = 'cell') %>%
  #   left_join(genes_map  %>%
  #               select(gene, chr, from),
  #             by = 'gene') %>%
  #   left_join(mapped_coordinates_cnas, by = 'gene') %>%
  #   filter(chr %in% chromosomes)
  
  # Auxisliary function
  aux_fun = function(chr)
  {
    input_RNA = input_RNA %>% filter(chr %in% !!chr)
    
    if (nrow(input_RNA) == 0)
      return(CNAqc:::eplot())
    
    # Largest cluster
    if (has_inference(x))
      largest_cluster = which.max(get_clusters_size(x)) %>% names
    else
      largest_cluster = "Not available"
    
    # Stdev for each cluster, then sorted by largest_cluster
    cluster_stdev = input_RNA %>%
      group_by(gene, from, chr, cluster, .drop = 'none') %>%
      summarise(
        aggregated = fun(n),
        segment_id = segment_id[1],
        .groups = 'keep'
      ) %>%
      ungroup()
    
    ngenes = cluster_stdev %>% pull(gene) %>% unique %>% length
    
    order_genes = cluster_stdev %>%
      filter(cluster == largest_cluster) %>%
      arrange(desc(aggregated)) %>%
      filter(row_number() < top) %>%
      pull(gene)
    
    cluster_stdev = cluster_stdev %>% filter(gene %in% order_genes)
    
    
    # cluster_stdev = cluster_stdev %>% filter(stdev > cutoff)
    # order_genes = order_genes[order_genes %in% cluster_stdev$gene]
    
    clr = 'black'
    if (has_inference(x))
      clr =
      pio:::nmfy(get_clusters_size(x) %>% names,
                 RColorBrewer::brewer.pal(n = get_k(x), "Set1"))
    
    if (has_inference(x))
      cluster_stdev$cluster = factor(cluster_stdev$cluster, levels = get_clusters_size(x) %>% sort %>% names)
    
    
    pl = ggplot(cluster_stdev, aes(x = gene,
                                   aggregated,
                                   fill = cluster)) +
      geom_bar(stat = 'identity') +
      CNAqc:::my_ggplot_theme() +
      scale_x_discrete(limits = order_genes) +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      # scale_fill_brewer(palette = 'Set1') +
      scale_fill_manual(values = clr) +
      guides(fill = guide_legend("")) +
      # theme(legend.position = c(.9, .9)) +
      labs(x = paste0("Genes (top ", top, "/", ngenes, ")"),
           y = 'Statistics') +
      theme(legend.background = element_rect(colour = NA, fill = NA))
    
    if (by_segment)
      pl = pl + facet_grid(segment_id ~ cluster_stdev$chr[1])
    else
      pl = pl + facet_wrap(~ cluster_stdev$chr[1])
    
    pl
  }
  
  plots = NULL
  for (chr in chromosomes)
    plots = append(plots, list(aux_fun(chr)))
  
  if(!assembly) return(plots)
  
  if (length(chromosomes) > 1)
  {
    if (length(chromosomes) == 4)
      b = c(2, 2)
    else
      b = c(ceiling(length(plots) / 3), 3)
    
    plots = ggarrange(
      plotlist = plots,
      ncol = b[2],
      nrow = b[1],
      common.legend = TRUE,
      legend = 'bottom'
    )
  }
  else
    plots = plots[[1]]
  
  return(plots)
}
