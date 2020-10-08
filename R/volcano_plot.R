volcano_plot <- function(x, pCutoff = 10e-5, FCcutoff = 1, ...){

  group1 <- x$DE$params$clone1
  group2 <- ifelse(is.null(x$DE$params$clone2), "all", x$DE$params$clone2)
  test_method <-  x$DE$params$method

  EnhancedVolcano::EnhancedVolcano(x$DE$table,
                  lab = DE$gene,
                  x = 'avg_log2FC',
                  y = 'p_val_adj',
                  ylab = "P-value adjusted",
                  pCutoff = pCutoff,
                  FCcutoff = FCcutoff,
                  title = paste0(group1, " vs ", group2),
                  subtitle = paste0("Test used ", test_method),
                  ...
                  )

}