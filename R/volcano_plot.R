volcano_plot <- function(X, pCutoff = 10e-5, FCcutoff = 1, colorChr = NULL, ...){
  table <-  get_DE_table(X, cut_pvalue = 1, cut_lfc = 0)

  if(!is.null(colorChr)){
    nChr <- length(colorChr)
    cols <- RColorBrewer::brewer.pal(8, name = "Dark2")
    cols <- c(cols[1:nChr], "grey")
    names(cols) <- c(colorChr, "other chrs")
    colorChr <- cols[table$chr]
    colorChr[is.na(colorChr)] <-  "grey"
    for(j in seq(nChr + 1)) {
      names(colorChr)[colorChr == cols[j]] <- names(cols)[j]
    }
  }
  group1 <- X$DE$params$clone1
  group2 <- ifelse(is.null(X$DE$params$clone2), "all", X$DE$params$clone2)
  test_method <-  X$DE$params$method
  EnhancedVolcano::EnhancedVolcano(table,
                  lab = table$gene,
                  x = 'avg_log2FC',
                  y = 'p_val_adj',
                  ylab = "P-value adjusted",
                  pCutoff = pCutoff,
                  FCcutoff = FCcutoff,
                  title = paste0(group1, " vs ", group2),
                  subtitle = paste0("Test used: ", test_method),
                  colCustom = colorChr,
                  ...
                  )

}