






load_gene_annotations <- function(genome) {
  return(data(text = paste0(genome, "_genes")))
}

get_gene_annotations = function(x)
{
  if (x$reference_genome %in% c('hg19', 'GRCh37')) {
    data('hg19_gene_coordinates')
    return(hg19_gene_coordinates)
  }

  if(x$reference_genome %in% c('hg38', 'GRCh38')) {
    data('hg38_gene_coordinates')
    return(hg38_gene_coordinates)
  }

  if(x$reference_genome %in% c('mm10', 'GRCm38')) {
    data('mm10_gene_coordinates')
    return(mm10_gene_coordinates)
  }

  stop("reference unknown?")
}


relative_to_absolute_coordinates <- function(df, genome = "hg38"){

  karyo <- load_genome(genome)
  karyo_sum <- cumsum(c(0,karyo))
  karyo_sum <- karyo_sum[-length(karyo_sum)]
  names(karyo_sum) <- names(karyo)


  df <-  df %>% mutate(to = to + karyo_sum[chr], from = from + karyo_sum[chr])

  return(df)
}


absolute_t_relative_coordinates <- function(df, genome = "hg38"){

  karyo <- load_genome(genome)
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

plot_singlecell_gwide = function(M, cells = unique(M$cell), ordering = unique(M$cell), legend = "value", by_seg = F, space = "free", title = "scRNA-seq over CNV segments", font_size = 12)
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
          panel.background = element_rect(fill = NA),  axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank()) +
    theme(text = element_text(size=font_size))

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
    scale_y_continuous(seq(1,5))+ ylim(c(1,5)) + scale_x_continuous(expand = c(0,0))+ ggtitle(title)


}


filter_MAF_seg <-  function(MAF_splitted, k){

  lapply(MAF_splitted, function(x) {
    x$value <-  runmed(x$value, k = k, endrule = "constant")
    x$value <-  x$value
    return(x)
  })
}

plot_MAF <- function(MAF_l, k = 31, fsize = 20) {

  g <-  as.factor(MAF_l$chr)
  MAF_splitted <-  split(MAF_l, f = g)

  MAF_splitted <- filter_MAF_seg(MAF_splitted, k = k)
  MAF_l <- do.call(rbind, MAF_splitted)

  MAF_l$chr <- factor(MAF_l$chr,levels = gtools::mixedsort(unique(MAF_l$chr)))
  p = ggplot(data = MAF_l, aes(x = as.numeric(start), y = runmed(value, k = k, endrule = "constant"), color = chr)) +
    geom_point() + xlab("") + facet_wrap(~chr, scales = "free_x", nrow = 3) + ylab("MAF") +
    theme(axis.title.x=element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), text = element_text(size=fsize))
  return(p)
}

plot_HMM <- function(segs, MAF_l= NULL, k = 31, fsize = 20, padding = 0.05) {

  g <-  as.factor(segs$chr)
  segs_splitted <-  split(segs, f = g)
  segs_l <- do.call(rbind, segs_splitted)

  segs_l$chr <- factor(segs_l$chr,levels = gtools::mixedsort(unique(segs_l$chr)))


  if(is.null(MAF_l)) {
    p = ggplot(data = segs_l, aes(xmin = as.numeric(start), xmax = as.numeric(end), ymin = tot - padding, ymax = tot + padding, fill = chr)) + geom_rect()+
      xlab("") + facet_wrap(~chr, scales = "free_x", nrow = 3) + ylab("HMM state") +
      theme(axis.title.x=element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), text = element_text(size=fsize))


  }else {
    g <-  as.factor(MAF_l$chr)
    MAF_splitted <-  split(MAF_l, f = g)

    MAF_splitted <- filter_MAF_seg(MAF_splitted, k = k)
    MAF_l <- do.call(rbind, MAF_splitted)

    MAF_l$chr <- factor(MAF_l$chr,levels = gtools::mixedsort(unique(MAF_l$chr)))

    p = p + geom_point(x = as.numeric(MAF_l$start), y = runmed(MAF_l$value, k = k, endrule = "constant"), color = "grey65", alpha = 0.6) + facet_wrap(~MAF_l$chr, scales = "free_x", nrow = 3)
  }

  return(p)
}

## to finish

simplify_segs <-  function(df) {

  res <- data.frame()

  old_chr = df$chr[1]
  old_start = df$start[1]
  old_tot = df$tot[1]

  for(n in 2:nrow(df)) {
      curr_chr = df$chr[n]
      curr_tot = df$tot[n]
      curr_start = df$start[n]

      if(curr_chr != old_chr) {

        res <- rbind(res, data.frame(t1 = old_chr, t2 = old_start, t3 = df$end[n-1], t4 = old_tot))
        old_chr = curr_chr
        old_start = curr_start
        old_tot = curr_tot

      } else if(curr_tot != old_tot & ! (n == nrow(df))) {
        res <- rbind(res, data.frame(t1 = old_chr, t2 = old_start, t3 = df$end[n-1], t4 = old_tot))
        old_chr = curr_chr
        old_start = curr_start
        old_tot = curr_tot
      }

      if(n == nrow(df)){
        if(curr_tot != old_tot){
          res <- rbind(res, data.frame(t1 = old_chr, t2 = old_start, t3 = df$end[n-1],t4 =  old_tot))
        }
        res <- rbind(res, data.frame(t1 = curr_chr, t2 = curr_start, t3 = df$end[n],t4 =  curr_tot))
      }
  }
  colnames(res) <- colnames(df)
  return(res)
}


bind_rcongas_data <-  function(a,b){

  res <- a
  a$data$counts <- cbind(a$data$counts, b$data$counts)
  a$data$cnv <-   rbind(a$data$cnv, b$data$cnv)
  a$data$bindims <-   cbind(a$data$bindims, b$data$bindims)
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




gene_hist <-  function(counts, clusters,gname = "SLC14A1", fsize = 20){
  data <- as.data.frame(counts[,gname])
  data$clust <-  clusters
  colnames(data) <-  c("expr", "cluster")
  ggplot(aes(x = expr, y = ..count..), data = data) + geom_histogram(aes(fill = paste(cluster))) + ggtitle(gname) +
    xlab("") + ylab("") + scale_fill_discrete("Clone") + theme_bw(base_size = fsize)

}

# Obtain a set of colours for a liste of cluster names
get_clusters_colors = function(labels, palette = 'Set1')
{
  cols = RColorBrewer::brewer.pal(n = 9, palette)
  labels = sort(labels %>% unique)

  cols = cols[1:length(labels)]
  names(cols) = labels

  return(cols)
}



filter_segments.rcongas <- function(X, filter_mu = 30, filter_fixed = FALSE, filter_ploidy = 1:12, filter_cells = rownames(X$data$counts)) {

  if(filter_fixed)
    mask_seg <- X$data$cnv$fixed_mu > filter_mu & !is.na(X$data$cnv$fixed_mu)
  else
    mask_seg <- X$data$cnv$mu > filter_mu & !is.na(X$data$cnv$mu)

  mask_ploidy <- X$data$cnv$ploidy_real %in% filter_ploidy & !is.na(X$data$cnv$ploidy_real)

  print(mask_seg)
  X$data$cnv <- X$data$cnv[mask_seg & mask_ploidy,]

  X$data$counts <- X$data$counts[filter_cells,mask_seg & mask_ploidy]

  X$data$bindims <- X$data$bindims[filter_cells,mask_seg & mask_ploidy]

  X$data$gene_locations <- X$data$gene_locations %>% dplyr::filter(segment_id %in% X$data$cnv$segment_id)

  X$data$gene_counts <- X$data$gene_counts[which(rownames(X$data$gene_counts) %in% X$data$gene_locations$gene),]

  return(X)
}

lighten <- function(color, factor=0.5){
  if ((factor > 1) | (factor < 0)) stop("factor needs to be within [0,1]")
  col <- col2rgb(color)
  col <- col + (255 - col)*factor
  col <- rgb(t(as.matrix(apply(col, 1, function(x) if (x > 255) 255 else x))), maxColorValue=255)
  col
}


p_value_format = function(p)
{
  lb = paste0("p = ", format.pval(p, eps = 0.001))

  if(p < 0.001) return(paste0("p < 0.001 (***)"))
  if(p < 0.01) return(paste0(lb, " (**)"))
  if(p < 0.05) return(paste0(lb, " (*)"))

  lb
}

# p_value_format(0.005)

from_MAP_to_post <- function(df) {

  ncols <- unique(df$value)
  res <- matrix(ncol = length(ncols), nrow = nrow(df))
  colnames(res) <-  ncols
  for(cols in ncols)
    res[,cols] <-  ifelse(df$value == cols, 1, 0)
  rownames(res) <- rownames(df)
  return(res)

}
