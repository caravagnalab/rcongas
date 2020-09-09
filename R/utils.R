
relative_to_absolute_coordinates <- function(df, genome = "hg38"){

  karyo <- eval(parse(text = paste0(genome, "_karyo")))
  karyo_sum <- cumsum(c(0,karyo))
  karyo_sum <- karyo_sum[-length(karyo_sum)]
  names(karyo_sum) <- names(karyo)


  df <-  df %>% mutate(to = to + karyo_sum[chr], from = from + karyo_sum[chr])

  return(df)
}


absolute_t_relative_coordinates <- function(df, genome = "hg38"){

  karyo <- load(paste0(genome, "_karyotype.rda"))
  karyo_sum <- cumsum(c(0,karyo))
  karyo_sum <- karyo_sum[-length(karyo_sum)]
  names(karyo_sum) <- names(karyo)

  df <-  df %>% mutate(to = to - karyo_sum[chr], from = from - karyo_sum[chr])

  return(df)
}


rntocl <- function(df, name = "rownames") {

  df[name] <- rownames(df)
  return(df)
}

long_counts <- function(M) {

  res <- as.matrix(M) %>%
    reshape2::melt() %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(seg = Var2) %>%
    tidyr::separate(Var2, sep = ':', into = c('chr', 'from', 'to')) %>%
    dplyr::mutate(
      Var1 = paste(Var1),
      from = as.integer(from),
      to = as.integer(to),
      value = as.numeric(value)
    ) %>%
    dplyr::rename(cell = Var1, n = value)

  return(res)

}

csidx <-  function(cumsm, idx){
  return((cumsm[idx]+1) : cumsm[idx+1])
}




gene_zscore = function(M)

{

  M <- M %>% reshape2::melt()
  colnames(M) <- c("cell", "gene", "counts")

  msd = M %>%
    group_by(gene) %>%
    summarise(m = mean(counts), sd = sd(counts))

  cm <- msd$m
  names(cm) <- msd$gene
  sd <- msd$sd
  sd <- ifelse(sd == 0, 1, sd)
  names(sd) <- msd$gene

  M <-  M %>%
            mutate(z = (counts - cm[gene])/sd[gene])

  return(M)

}

plot_singlecell_gwide = function(M, cells = unique(M$cell), ordering = unique(M$cell), legend = "value", by_seg = F, space = "free", title = "scRNA-seq over CNV segments")
{
  N = M %>% filter(cell %in% !!cells)
  N$chr <-  factor(N$chr,levels= gtools::mixedsort(unique(N$chr)))


  # Cells by n
  res <- ggplot(data = N) +
    geom_segment(
      aes(
        x = from,
        xend = to,
        y = cell,
        yend = cell,
        color = n
      )
    )  +
    scale_color_viridis_c(direction = -1, legend) + ggtitle(title) +
    geom_vline(data = N, aes(xintercept = from), color = 'white', size = .2) +
    scale_y_discrete(limits = ordering) + scale_x_continuous(expand = c(0,0)) +
    theme(panel.spacing = unit(0.01, "lines"), axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank(),
          panel.background = element_rect(fill = NA),  axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(),
          strip.text.x = element_text(size = 5))

  if(!by_seg){
    res <-  res + facet_grid(.~chr,scales = "free", space = space)
  } else{
    res <-  res + facet_grid(.~seg,scales = "free", space = space)
  }

  return(res)
}

cluster_cells = function(N, stat = 'n', cells = unique(M$cell))
{
  N = M %>% filter(cell %in% !!cells)

  # .. or cluster the groups
  datamat = N %>%
    arrange(cell) %>%
    mutate(L = paste0(chr, from, to)) %>%
    dplyr::select(cell, L, !!stat) %>%
    pivot_wider(names_from = L,  values_from = stat)

  cell_names = datamat$cell

  datamat = datamat %>% dplyr::select(-cell) %>% as.matrix() %>% apply(2, as.numeric)

  d <- dist(datamat)   # find distance matrix
  h = hclust(d)
  ordering = cell_names[h$order]


  ordering

}

plot_bulk = function(cnv, chromosomes = paste0(c(1:22, 'X', 'Y')), title = "")

  {

  cnv$chr <-  factor(cnv$chr,levels= gtools::mixedsort(unique(cnv$chr)))
  cnv <- cnv %>% select(matches("ploidy"), chr, from, to) %>% select(-ploidy_real)
  cnv <-  reshape2::melt(cnv, id.vars = c("chr", "from", "to"))
  colnames(cnv) <- c("chr", "from", "to", "cluster", "ploidy_real")
  cnv$cluster <- substring(cnv$cluster, first = 7)
  pl_segm_bulk <- ggplot() +
    geom_segment(
      data = cnv,
      aes(
        x = from,
        xend = to,
        y = ploidy_real + 0.02 * (as.numeric(cluster) - 1 ),
        yend = ploidy_real + 0.02 * (as.numeric(cluster) - 1),
        color = cluster,
      ),
      size = 1.5,
    ) + facet_grid(.~chr,scales = "free_x", space = "free_x")+
    geom_vline(data = cnv, aes(xintercept = from), linetype = 'dotted', color = 'grey50', size = .3, alpha = 0.8) +
    ylab("Ploidy") + theme(panel.spacing = unit(0.01, "lines"), axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank(),
                          panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank())+
    scale_y_continuous(breaks = seq(floor(min(cnv$ploidy_real)), ceiling(max(cnv$ploidy_real)) , by = 1)) + scale_x_continuous(expand = c(0,0))+ ggtitle(title)


}


filter_MAF_seg <-  function(MAF_splitted, k){

  lapply(MAF_splitted, function(x) {
    x$value <-  runmed(x$value, k = k, endrule = "constant")
    x$value <-  x$value
    return(x)
  })
}

plot_MAF <- function(MAF_l, k = 31) {

  g <-  as.factor(MAF_l$chr)
  MAF_splitted <-  split(MAF_l, f = g)

  MAF_splitted <- filter_MAF_seg(MAF_splitted, k = k)
  MAF_l <- do.call(rbind, MAF_splitted)

  MAF_l$chr <- factor(MAF_l$chr,levels = gtools::mixedsort(unique(MAF_l$chr)))
  ggplot(data = MAF_l, aes(x = as.numeric(start), y = runmed(value, k = k, endrule = "constant"), color = chr)) +
    geom_point() + xlab("") + facet_wrap(~chr, scales = "free_x") + ylab("MAF") +
    theme(axis.title.x=element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

}


simplify_segs <-  function(df) {

  df %>% group_by(chr,tot) %>%  dplyr::summarize(start = min(start), end = max(end))

}


bind_confettti <-  function(a,b){

  res <- a
  a$counts <- cbind(a$counts, b$counts)
  a$cnv <-   rbind(a$cnv, b$cnv)
  a$bindim <-   cbind(a$bindim, b$bindim)
  return(a)

}

MAF_clust_plot <-  function(ref, cov, idx1, idx2, chrs, filt = 30, labels = c("c1", "c2"), both = F){

  c1_cells <- idx1
  c2_cells <- idx2

  alt <-  cov - ref

  ref_c1 <- rowSums(ref[,c1_cells])
  alt_c1 <- rowSums(alt[,c1_cells])
  min_c1 <- pmin(ref_c1, alt_c1)
  maj_c1 <- pmax(ref_c1, alt_c1)
  cov_c1 <- ref_c1 + alt_c1

  ref_c2 <- rowSums(ref[,c2_cells])
  alt_c2 <- rowSums(alt[,c2_cells])
  min_c2 <- pmin(ref_c2, alt_c2)
  maj_c2 <- pmax(ref_c2, alt_c2)
  cov_c2 <- ref_c2 + alt_c2

  if(both)
    filter <- which(pmin(cov_c1,cov_c2) > filt)
  else
    filter <- which(cov_c1 > filt)

  MAF_c1 <- min_c1[filter] / cov_c1[filter]
  MAF_c2 <- min_c2[filter] / cov_c2[filter]

  df <- rbind("c1" = MAF_c1, "c2" = MAF_c2 )
  df_l <-  reshape2::melt(df)

  colnames(df_l) <- c("cluster", "SNP", "MAF")

  df_l <- df_l[grep(chrs,df_l$SNP),]

  print(df_l)

  ggplot(data = df_l, aes(x = SNP, y = MAF, color = cluster)) + geom_point() + scale_color_discrete(labels = labels) + theme_bw()  + theme(axis.text.x = element_text(angle = 90))


}


plot_confusion_matrix <-function(True_vec, Pred_vec){

  comm_cells <- intersect(names(True_vec), names(Pred_vec))

  True_vec <- True_vec[comm_cells]
  Pred_vec <- Pred_vec[comm_cells]

  cm <-  table(True_vec,Pred_vec)

  print(cm)

  colnames(cm) <- paste0("clone_", seq(ncol(cm)))
  #rownames(cm) <- paste0("real_", seq(nrow(cm)))

  cm_df <- reshape2::melt(cm)

  colnames(cm_df) <- c("True", "Pred", "Tot")

  ggplot(data =  cm_df, mapping = aes(x = True, y =Pred)) +
    geom_tile(aes(fill = Tot), colour = "white") +
    geom_text(aes(label = sprintf("%1.0f", Tot)), vjust = 1, size = 12.9) + xlab("Clonealign") + ylab("CONGAS") +
    scale_fill_gradient(low = "blue", high = "red") +
    theme_minimal() + theme(legend.position = "none") + coord_flip()

}



calculate_DGE <-  function(inference,gene_counts,clone1, clone2, method = "wilcox", normalize = T){

    so <- CreateSeuratObject(gene_counts)
    if(normalize)
      so <- NormalizeData(so)
    so@meta.data$membership <- factor(paste0(inference$parameters$assignement))
    so <- SetIdent(object = so, value = "membership")
    df_genes <- Seurat::FindMarkers(so, ident.1 = clone1, ident.2 = clone2, test.use = method, logfc.threshold
 = 0.25)

    return(df_genes)

}



calculate_GSEA <- function(inference,fc_df,clone1, clone2) {

  if(require("org.Hs.eg.db")){
    organism <- org.Hs.eg.db::org.Hs.eg.db
  } else {
    BiocManager::install("org.Hs.eg.db", character.only = TRUE)
    organism <- org.Hs.eg.db::org.Hs.eg.db
  }

  change <- log2(unlist(gtools::foldchange(fc_df[clone1,], fc_df[clone2,])))
  names(change) <-  colnames(fc_df)
  gene_list <- sort(change, decreasing = TRUE)
  gse <- gseGO(geneList=gene_list,
               ont ="ALL",
               keyType = "SYMBOL",
               minGSSize = 3,
               maxGSSize = 800,
               pvalueCutoff = 0.05,
               verbose = TRUE,
               OrgDb = organism,
               pAdjustMethod = "none", eps = 0)


  ids<-bitr(names(change), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)
  dedup_ids <-  ids[!duplicated(ids[c("SYMBOL")]),]
  change2 <-  change[names(change) %in% dedup_ids$SYMBOL]
  names(change2) <-  dedup_ids$ENTREZID
  kegg_gene_list <- change2
  kegg_gene_list<-na.omit(kegg_gene_list)
  kegg_gene_list <-  sort(kegg_gene_list, decreasing = TRUE)

  kegg_organism <-  "hsa"
  kk2 <- gseKEGG(geneList     = kegg_gene_list,
                 organism     = kegg_organism,
                 minGSSize    = 3,
                 maxGSSize    = 800,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "none",
                 keyType       = "ncbi-geneid", eps = 0)

  return(list(kegg = kk2, go = gse))

}



