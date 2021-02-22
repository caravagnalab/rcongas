plot_gene_importance_per_segments <- function(x, segments_id = get_highlights_names(x), ntop = 10){
  
  plts <- lapply(segments_id, function(y) plot_gene_importance_per_segments_aux(x,y, ntop))
  names(plts) <- segments_id
  return(plts)
}

plot_gene_importance_per_segments_aux <-function(x, segment, top){
  
  genes <- get_genes_mapped_to_segments(x) %>% filter(segment_id == segment) %>%  pull(gene)
  counts <- x$data$gene_counts[genes,]
  norm <-  colSums(counts)
  #CPM normalization
  counts_norm <- counts / norm * 10^6
  counts_norm <- log(counts_norm + 1)
  clusters <-  get_clusters(x) %>%  select(cell, cluster)
  counts_norm_l <-  counts_norm %>%  reshape2::melt()
  colnames(counts_norm_l) <-  c("gene", "cell", "value")
  counts_c <- dplyr::inner_join(counts_norm_l, clusters)
  plot_df <-  counts_c %>%  dplyr::group_by(cluster, gene) %>%  dplyr::summarize(mean = mean(value), sd = sd(value))
  sd_by_cluster <- plot_df %>%  dplyr::group_by(gene) %>%  dplyr::summarise(sd_gr = sd(mean), mean_gr = mean(mean)) %>% 
    dplyr::arrange(-sd_gr) %>% dplyr::top_n(top)
  plot_df_filt <- plot_df %>%  filter(gene %in% sd_by_cluster$gene)
  colors <- Rcongas:::get_clusters_colors(plot_df_filt$cluster %>%  unique())
  ggplot(aes(y = mean, x= gene, fill = cluster), data = plot_df_filt) + geom_bar(stat = "identity", position = "dodge") +
    CNAqc:::my_ggplot_theme() + coord_flip() + ggtitle(paste(segment)) + scale_fill_manual("Cluster", values = colors)
  
}


get_highlights_names <-  function(x) {
  hn <- highlights(x) %>% filter(highlight) %>% pull(segment_id) %>%  unique()
  return(hn)
}
