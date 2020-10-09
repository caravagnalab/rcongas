plot_input_rna = function(x)
{
  # Try to get data
  M = get_input_rna_data(x)

  if(all(is.null(M))) return(CNAqc:::eplot())

  # Actual stuff
  dim(M)

  RNA = reshape2::melt(M)
  colnames(RNA) = c("cell", "segment", "value")

  # Input segmentation - get size in Megabases (Mb)
  input_segments = get_input_segmentation(x) %>%
    dplyr::mutate(size = ceiling((to - from)/10^6), label_chr = paste0(chr, " (", size, 'Mb)')) %>%
    idify()

  # Pheatmap
  require(pheatmap)

  pheatmap(
    M,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    border_color = NA,
    show_rownames = FALSE,
    fontsize_col = 6,
    fontsize_number = 8,
    # gaps_row = cuts_classes,
    # color = rna_colors,
    # annotation_row = clusters_both_run %>% select(-cell),
    # annotation_colors = list(
    #   `run_1` =
    #     pio:::nmfy(
    #       unique(clusters_both_run$run_1),
    #       rev(ggsci::pal_npg()(unique(clusters_both_run$run_1) %>% length))
    #     ),
    #   `run_2` = pio:::nmfy(
    #     unique(clusters_both_run$run_2),
    #     ggsci::pal_nejm()(unique(clusters_both_run$run_2) %>% length)
    #   )
    )
  )


}