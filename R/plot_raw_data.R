#' Title
#'
#' @param genes 
#' @param lognormalise 
#' @param description 
#' @param clusters 
#' @param prompt 
#' @param ... 
#' @param x
#'
#' @return
#' @export
#'
#' @examples
#'
#' x = Rcongas::congas_example
#'
#' print(x)
#'
#' # Extract raw data
#' raw_x = get_input_raw_data(x)
#' 
#' # Get some genes to plot from those with DE
#' genes = get_DE_table(x) %>% pull(gene)
#'
#' # Get clusters
#' clusters = get_clusters(x)
#' 
#' # Default plot (without prompt for automatic documentation generation)
#' plot_raw_data(raw_x, genes, prompt = FALSE)
plot_raw_data = function(x,
                         genes,
                         lognormalise = TRUE,
                         description = "My CONGAS model",
                         clusters = NULL,
                         prompt = TRUE,
                         ...)
{
  # xdim(x)
  # genes = get_gene_annotations(rcongas_example)$gene[1:1000]
  # x=get_input_raw_data(rcongas_example)
  
  # Subset by gene
  x = x[rownames(x) %in% genes, , drop = FALSE]
  
  ng = nrow(x)
  nc = ncol(x)
  
  if (length(ng) == 0) {
    cli::cli_alert_warning("No genes in the input list, returning empty plot...")
    return(CNAqc:::eplot())
  }
  
  if (nrow(x) == 0 | ncol == 0) {
    cli::cli_alert_warning("0 rows or column in the input data, returning empty plot...")
    return(CNAqc:::eplot())
  }
  
  # Prompt
  if (prompt)
  {
    cli::cli_alert_warning(
      "You are about to plot {.field {ng}} genes and {.field {nc}} cells with pheatmap.
                           This may take a while if the values are large..."
    )
    
    ask = menu(c("Yes", "No"), title = "Do you want this?")
    
    if (ask == 2)
      return(CNAqc:::eplot())
  }
  
  # This plot uses pheatmap
  require(pheatmap)
  
  # colnames/rownames smart show, same for clustering
  show_rownames = show_colnames = FALSE
  cluster_cols = TRUE
  
  if (ng < 100)
    show_rownames = TRUE
  if (nc < 100)
    show_colnames = TRUE
  
  if (!all(is.null(clusters)))
    cluster_cols = FALSE
  
  if (!show_rownames)
    cli::cli_alert_info("With more than 100 genes, gene names are hidden.")
  if (!show_colnames)
    cli::cli_alert_info("With more than 100 cells, cell ids are hidden.")
  
  # Create annotations, only for clusters
  annotation_columns = NULL
  if (!all(is.null(clusters)))
  {
    # set ordering based on annotations
    clusters = clusters %>% dplyr::arrange(cluster)
    x = x[, clusters$cell, drop = FALSE]
    
    annotation_columns = data.frame(Cluster = clusters$cluster, stringsAsFactors = FALSE)
    rownames(annotation_columns) = clusters$cell
  }
  
  # Call
  data_transform = x
  main = description
  
  if (lognormalise) {
    x = log(x + 1)
    main = paste0(description, " log-transformed counts")
  }
  
  pheatmap(
    data_transform,
    show_rownames = show_rownames,
    show_colnames = show_colnames,
    cluster_cols = cluster_cols,
    color = colorRampPalette(rev(
      RColorBrewer::brewer.pal(n = 9, "Spectral")
    ))(100),
    annotation_col = annotation_columns,
    annotation_colors = list(
      Cluster = Rcongas:::get_clusters_colors(annotation_columns$Cluster)
    ),
    main = main,
    ...
  )
}
#
#
# plot_raw_data(
#   get_input_raw_data(x),
#   genes =  c('APC', 'KRAS', "TP53")
# )
