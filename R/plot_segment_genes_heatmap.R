plot_segment_genes_heatmap <- function(x, segment_id = NULL,type = c("scaled", "relative", "lognorm"),
                                       filename = "heatmap_counts_segments.pdf", cluster_genes = FALSE,
                                       min = NULL, max = NULL, median_filt = FALSE, k = 5, filt_counts = 0){

  if(is.null(segment_id)){
    high <- get_highlights_names(x)
  }else {
    high <- segment_id
  }

  g_names <- get_genes_mapped_to_segments(x) %>% filter(segment_id %in% high) %>% select(gene, segment_id) %>%  arrange(segment_id) %>% unique()
  #To change with the new version
  counts <-  get_input_raw_data(x)

  counts <- counts[rowSums(counts) > filt_counts,]

  if(type == "relative"){
    counts <- counts / colSums(counts) * 1e6
  } else if(type == "lognorm") {
    counts <- log1p(counts / colSums(counts) * 1e6)
  }else if(type == "scaled"){
    counts <- log1p(counts / colSums(counts) * 1e6) %>% t %>%  scale %>%  t
  } else {
    stop("Normalization type not imllemented, use one among scaled, relative or lognorm")
  }

  counts <- counts[intersect(g_names$gene, rownames(counts)),]

  annot_row <- get_clusters(x) %>% select(cell,cluster) %>%  arrange(cluster) %>% tibble::column_to_rownames("cell")
  annot_col <- g_names %>% tibble::column_to_rownames("gene")
  annot_col <-  annot_col[gtools::mixedorder(annot_col)]
  colors <- rev(RColorBrewer::brewer.pal("RdBu", n = 11))
  if(is.null(min)) min = min(counts)
  if(is.null(max)) max = max(counts)
  if(median_filt) {
    rnames <- rownames(counts)
    counts <-  apply(counts, 2, function(x) runmed(x, k = k))
    rownames(counts) <-  rnames
  }
  pdf(filename, width = 20, height = 12)
  pheatmap::pheatmap(counts[, rownames(annot_row)] %>%  t, cluster_rows = F, color = colors,
                     cluster_cols = cluster_genes,breaks = seq(min, max, length.out = 11),
                     annotation_row = annot_row, annotation_col = annot_col,
                     show_rownames = F, show_colnames = F, gaps_row = cumsum(table(annot_row))[-length(table(annot_row))] ,
                     gaps_col = cumsum(table(annot_col)) [-length(table(annot_col))])
  dev.off()

}
