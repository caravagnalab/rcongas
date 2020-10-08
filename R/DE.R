#
# x = calculate_DGE(get_best_model(x), mat_pre2 %>% t, 1, 2, method = 'wilcox')

calculate_DE <-
  function(x,
           input,
           clone1,
           clone2,
           method = "wilcox",
           normalize = T,
           logfc.threshold = 0.25)
{
    # Test input corretto - data matrix e oggetto rcongas

    # Work with the best model fit
    best_m = get_best_model(x)

    # Create Seurat object
    so <- CreateSeuratObject(input)

    if (normalize) so <- NormalizeData(so)

    # TODO - estensione a gruppi
    so@meta.data$membership <-
      factor(paste0(best_m$parameters$assignement))

    so <- SetIdent(object = so, value = "membership")

    # DE analysis
    df_genes <-
      Seurat::FindMarkers(
        so,
        ident.1 = clone1,
        ident.2 = clone2,
        test.use = method,
        logfc.threshold
        = logfc.threshold
      )

    colnames(df_genes) <-
      gsub(pattern = "\\.",
           replacement = "_",
           colnames(df_genes))

    # Extend results with gene mapping
    gene_ann = get_gene_annotations(x)

    df_genes$gene = rownames(df_genes)

    DE = gene_ann %>%
      dplyr::right_join(df_genes, by = 'gene') %>%
      dplyr::mutate(avg_log2FC = log2(exp(avg_logFC))) %>%
      dplyr::rename(pct_detec_1 = pct_1, ptc_detec_2 = pct_2) %>%
      dplyr::select(-avg_logFC) %>%
      dplyr::as_tibble()

    # Store DE setup and results
    x$DE = list(
      params = list(
        clone1 = clone1,
        clone2 = clone2,
        method = method,
        normalize = normalize,
        logfc.threshold = logfc.threshold
      ),
      table = DE
    )

    x$DE <- df_genes

    return(x)
}
