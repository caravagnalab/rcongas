#
# X= calculate_DE(get_best_model(X), mat_pre2 %>% t, 1, 2, method = 'wilcox')

# input = get_input_raw_data(x)
# xx = calculate_DE(x, input, 1, 2)

#' Title
#'
#' @param X
#' @param input
#' @param clone1
#' @param clone2
#' @param method
#' @param normalize
#' @param logfc.threshold
#'
#' @return
#' @export
#'
#' @examples
calculate_DE <-
  function(x,
           input,
           clone1,
           clone2,
           method = "wilcox",
           normalize = TRUE,
           logfc.threshold = 0.25)
  {
    # Test input corretto - data matrix e oggetto rcongas
    stopifnot(inherits(x, "rcongas"))
    stopifnot(!is.null(input))

    # Screen intro
    cli::cli_h1("Differential Expression analysis with Seurat")
    cat("\n")

    cat(cli::boxx(
      paste0("Groups: ", clone1, " vs ", clone2, "; using ", method),
      padding = 1,
      float = "center",
      background_col = "orange",
      col = 'black'
    ))

    cat("\n\n")

    cli::cli_text(
      '- With normalisation: {.field {normalize}}; LFC threshold {.field {logfc.threshold}}. '
    )

    cat("\n")

    # Work with the best model fit
    best_m = Rcongas:::get_best_model(x)

    # Create Seurat object
    require(Seurat)

    so <- suppressWarnings(CreateSeuratObject(input))

    if (normalize)
      so <- NormalizeData(so)

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
    gene_ann = Rcongas::get_gene_annotations(x)

    df_genes$gene = rownames(df_genes)

    DE = gene_ann %>%
      dplyr::right_join(df_genes, by = 'gene') %>%
      dplyr::rename(pct_detec_1 = pct_1, ptc_detec_2 = pct_2) %>%
      dplyr::select(-avg_log2FC, -p_val) %>%
      dplyr::as_tibble() %>% dplyr::arrange(p_val_adj)

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

    n_p_sign = sum(x$DE$table$p_val_adj < 0.05, na.rm = T)

    if(n_p_sign == 0)
      cli::cli_alert_warning("No DEG at alpha level 0.05!")
    else
      cli::cli_alert_success("Found {.field {n_p_sign}} EG at alpha level 0.05.")

    return(x)
  }
