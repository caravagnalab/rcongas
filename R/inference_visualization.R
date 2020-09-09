plot_loss <- function(anneal_inf,...) UseMethod("plot_loss", anneal_inf)

hist_plot <- function(anneal_inf,...) UseMethod("hist_plot", anneal_inf)

plot_clusters <- function(anneal_inf,...) UseMethod("plot_clusters", anneal_inf)



plot_loss.anneal <- function(anneal_inf) {

  loss <- anneal_inf$loss


  ggplot() + geom_line(aes(x= seq_along(loss),y=loss),size = 1.2, color = 'blue') + xlab("step") + ggtitle("Loss decay")

}


plot_counts.anneal <- function(anneal_inf, counts, diff = FALSE,chrs = c(1:22, "X") ,
                               order = names(sort(anneal_inf$parameters$assignement)), norm = F, by_seg = F, space = "free", legend = "Counts"){

  if(inherits(counts, "CNVSimulation")){
    real_ids <- counts$clust_ids
    counts_mat <-  counts$counts
    plot_dim <- 12
  } else {
    counts_mat <-  counts
    plot_dim <- 18
  }

  if(norm){
    counts_mat <- apply(counts_mat, MARGIN = 2, function(x) x / anneal_inf$parameters$norm_factor)
  }




  M <- long_counts(counts_mat)
  M <- M %>% dplyr::as_tibble() %>% filter(chr %in% chrs)

  if(diff){
    segments <- colnames(anneal_inf$parameters$cnv_probs[,which(abs(anneal_inf$parameters$cnv_probs[1,] - anneal_inf$parameters$cnv_probs[2,]) > 0.6)])
    M <- M %>% dplyr::as_tibble %>% filter(seg %in% segments)
  }

  counts_plot <-  plot_singlecell_gwide(M, legend = legend,ordering = order, by_seg = by_seg, space = space)

  if(inherits(counts, "CNVSimulation")){

    clus_inf <-  cbind(anneal_inf$parameters$assignement[order], real_ids[order])
    print(clus_inf)
    colnames(clus_inf) <- c("inferred", "real")

  } else {

    clus_inf <-  as.data.frame(anneal_inf$parameters$assignement[order])
    print(clus_inf)
    colnames(clus_inf) <- "inferred"

  }



  clus_inf <- clus_inf %>% mutate(cell_id = order, inferred = paste0("c",inferred)) %>% melt(id.vars = "cell_id")

  clus_inf$cell_id <-  factor(clus_inf$cell_id, levels = order)

  clust_plot <- ggplot(data = clus_inf) +
    geom_tile(aes(x = variable, y = cell_id, fill=paste(value) ))  + scale_fill_discrete() +
    scale_y_discrete(limits = gtools::mixedsort(unique(clus_inf$cell_id))) +
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1),
          axis.ticks.y=element_blank(),axis.ticks.x=element_blank(), axis.title.x = element_blank(), legend.position="none") +
    ggtitle("")

  cell_counts_and_clust_plot <- cowplot::plot_grid(clust_plot, counts_plot,ncol = 2, rel_widths = c(1,plot_dim), axis = "tb", align = "hv")


  return(cell_counts_and_clust_plot)

}

plot_clusters.anneal <- function(anneal_inf, diff = F, chrs = c(1:22, "X"), round = F, low_diff_rem = F){

  CNV <- anneal_inf$parameters$cnv_probs %>% as.data.frame
  if(diff){
    chrs <- CNV[,which(abs(CNV[1,] - CNV[2,]) > 0.6)]
  }

  CNV <- CNV %>% select(starts_with(chrs)) %>% rntocl("cid")  %>% reshape2::melt()

  if(round){
    CNV$value <-  round(CNV$value)
  }

  pl_segm <- ggplot() +
    geom_segment(
      data = CNV  %>% mutate(value = value + 0.01 * as.numeric(cid)),
      aes(
        x = 0,
        xend = 50,
        y = value,
        yend = value,
        color = paste(cid)
      ),
      size = 1,
    ) + facet_wrap(.~variable,scales = "free_x")+ scale_color_discrete("Cluster_id") +
    ylab("Ploidy") + theme(panel.spacing = unit(0.03, "lines"), axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank(),
                           panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), legend.position = "bottom") + scale_x_continuous(expand = c(0,0))

  return(pl_segm)
}


hist_plot.default <- function(counts, CNV = NULL,norm_f = NULL, cid = rep(1, nrow(counts)), norm = F, diff = F,
                              segments = c(1:22, "X"), sim = F, real_cid = NULL, norm_genes = F, norm_genes_df = NULL)
{
  if(norm) counts <- apply(counts,2, function(x) x / norm_f)

  if(norm_genes) counts <- t(apply(counts,1, function(x) x / norm_genes_df))


  if(diff){
    segments <- colnames(counts) [which(abs(CNV[1,] - CNV[2,]) > 0.6)]
  }

  counts <-  counts %>% dplyr::as_tibble() %>% select(dplyr::matches(segments)) %>%  mutate(cid = cid)


  if(sim) {
    counts <-  counts %>%  mutate(real = real_cid[[1]])
    counts <- counts %>% melt(id.vars = c("cid", "real")) %>%  mutate(diff = (cid != real))
  }
  else{
    counts <- counts %>% melt(id.vars = "cid")
  }

  counts$variable <-  factor(counts$variable, levels = gtools::mixedsort(unique(counts$variable)))

  diff_distr <- ggplot(data = counts) + geom_histogram(aes(y = ..density..,  x = value, fill = paste(cid)), alpha = 0.92) +
    facet_wrap(.~variable, scales="free") +
    theme_bw()+theme(legend.position = "bottom") + ggtitle("Counts distribution in segments") + scale_fill_discrete("Cluster Id") + theme(legend.position = "none")

  if(sim) {
    diff_distr <-  diff_distr +
      geom_point(data = counts,aes(x =value, y = 0, color = diff), alpha = 0.6, size = 0.6) + scale_color_manual("missclassified?", values = c("grey69", "red"))

  }


  return(diff_distr)


}

hist_plot.anneal <-  function(anneal_inf, counts, norm_genes_df = NULL, segments = c(1:22, "X"), diff = F, norm = F, norm_genes = F){

  sim <-  FALSE

  if(inherits(counts, "CNVSimulation")){
    real_cid <-  counts$clust_ids
    norm_f <- counts$norm_factors
    counts <- counts$counts
    sim <-  TRUE
    norm_genes_df <- counts$cnv$mu


  } else {
    norm_f <- anneal_inf$parameters$norm_factor

  }


  CNV <- anneal_inf$parameters$cnv_probs

  cid <- anneal_inf$parameters$assignement


  return(hist_plot.default(counts, CNV =CNV, norm_f = norm_f, cid = cid, norm=norm, diff=diff, segments = segments, sim = sim,
                           real_cid = real_cid,  norm_genes_df = norm_genes_df, norm_genes= norm_genes ))





}




plot_umap.anneal <-  function(anneal_inf, counts, norm = T, chrs =  c(1:22, "X"), custom_col = F, color = NULL, legend_title = "Clusters"){


  counts <- counts %>%  as_tibble %>% select(dplyr::starts_with(chrs))
  if(norm) counts <- apply(counts,2, function(x) x / anneal_inf$parameters$norm_factor)
  ump <- umap::umap(counts)$layout
  if(custom_col)
    ump <-  ump %>% as_tibble %>% mutate(clust = paste(color))
  else
    ump <-  ump %>% as_tibble %>% mutate(clust = paste(anneal_inf$parameters$assignement))
  colnames(ump) <-  c("dim1", "dim2", "clust")


  res <- ggplot(aes(x = dim1, y = dim2, col = clust), data = ump) + geom_point(size = 0.5)+ scale_color_discrete(legend_title) +theme(
                                                                      axis.text.x=element_blank(), axis.text.y = element_blank()
                                                                      ) + ggtitle("UMAP") + theme_bw()
  return(list(plot = res,dims = ump))
}