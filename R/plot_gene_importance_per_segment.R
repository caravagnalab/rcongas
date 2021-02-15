plot_gene_importance_per_segments <- function(x, segments_id = get_highlights_names(x), ntop = 10){
  
  plts <- lapply(segment_id, function(y) plot_gene_importance_per_segments_aux(x,y, ntop))
  return(plts)
}

plot_gene_importance_per_segments_aux <-function(x, segmentn, top){
  
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
  plot_df <-  counts_c %>%  group_by(cluster, gene) %>%  summarize(mean = mean(value), sd = sd(value))
  
}


get_highlights_names <-  function(x) {
  hn <- highlights(x) %>% filter(highlight) %>% pull(segment_id) %>%  unique()
  return(hn)
}
