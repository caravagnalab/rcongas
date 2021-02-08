

plot_report <-
  function(cnv_sim, ...)
    UseMethod("plot_report", cnv_sim)


plot_counts <-
  function(cnv_sim, ...)
    UseMethod("plot_counts", cnv_sim)

plot_genome <-
  function(cnv_sim, ...)
    UseMethod("plot_genome", cnv_sim)


plot.CNVSimulation <- function(cnv_sim, plot = T) {
  cnv_mat <- cnv_sim$cnv_mat
  
  N <- long_counts(cnv_mat)
  
  
  cell_cm_pl = plot_singlecell_gwide(N %>% mutate(n = n), legend = "ploidy")
  
  if (plot) {
    plot(cell_cm_pl)
  }
  return(cell_cm_pl)
  
}


plot_counts.CNVSimulation <-
  function(cnv_sim,
           plot = F,
           order = rownames(cnv_sim$clust_ids),
           by_seg = FALSE,
           title = "scRNA-seq over CNV segments",
           font_size = 20) 
  {
    counts <- cnv_sim$data$counts
    
    M <- long_counts(counts)
    
    cell_n_pl = plot_singlecell_gwide(
      M ,
      legend = "counts",
      ordering = order,
      by_seg = by_seg,
      title = title,
      font_size = font_size
    )
    
    if (plot) {
      plot(cell_n_pl)
    }
    
    clus_inf <-  reshape2::melt(as.matrix(cnv_sim$data$clust_ids))
    
    
    clust_plot <- ggplot(data = clus_inf) +
      geom_tile(aes(
        x = Var2,
        y = Var1,
        fill = paste(value)
      ))  + 
      scale_fill_discrete("") +
      scale_y_discrete("", limits = gtools::mixedsort(unique(clus_inf$cell_id))) +
      theme(
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none"
      ) +
      ggtitle("") + 
      xlab("")
    
    cell_counts_and_clust_plot <-
      cowplot::plot_grid(
        clust_plot,
        cell_n_pl,
        ncol = 2,
        rel_widths = c(1, 18),
        axis = "tb",
        align = "hv"
      )
    
    return(cell_counts_and_clust_plot)
    
  }


plot_genome.CNVSimulation <-
  function(cnv_sim,
           plot = T,
           title = "",
           fsize = 10) {
    cnv <- cnv_sim$cnv
    
    
    bulk_pl = plot_bulk(cnv, title = title) + theme(text = element_text(size =
                                                                          fsize))
    if (plot) {
      plot(bulk_pl)
    }
    return(bulk_pl)
    
  }



plot_report.CNVSimulation <- function(cnv_sim, file = "report.pdf") {
  res <-  cowplot::plot_grid(
    plot(cnv_sim, FALSE),
    plot_counts(cnv_sim, FALSE),
    plot_genome(cnv_sim, FALSE),
    ncol = 1,
    align = 'vh',
    axis = "tblr",
    rel_heights = c(3, 3, 1.5)
  ) %>% ggsave(filename = file,
               width = 15,
               height = 12)
  
  
  return(res)
}