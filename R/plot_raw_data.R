#' Title
#'
#' @param x
#' @param normalised
#'
#' @return
#' @export
#'
#' @examples
#'
#'   x = get_input_raw_data(x)

plot_raw_data = function(x, genes, clusters = NULL, prompt = TRUE)
{
  # xdim(x)
  # genes = get_gene_annotations(rcongas_example)$gene[1:1000]
  # x=get_input_raw_data(rcongas_example)

  # Subset by gene
  x = x[rownames(x) %in% genes, , drop = FALSE]

  ng = nrow(x)
  nc = ncol(x)

  if(length(ng) == 0) {
    cli::cli_alert_warning("No genes in the input data, returning empty plot...")
    return(CNAqc:::eplot())
  }


  if(prompt)
  {
    cli::cli_alert_warning("You are about to plot {.field {ng}} genes and {.field {nc}} cells with pheatmap.
                           This may take a while if the values are large...")

    ask = menu(c("Yes", "No"), title="Do you want this?")

    if(ask == 2) return(CNAqc:::eplot())
  }

  # This plot uses pheatmap
  require(pheatmap)

  # colnames/rownames smart show
  show_rownames = show_colnames = FALSE

  if(ng < 100) show_rownames = TRUE
  if(nc < 100) show_colnames = TRUE

  if(!show_rownames) cli::cli_alert_info("With more than 100 genes, gene names are hidden.")
  if(!show_colnames) cli::cli_alert_info("With more than 100 cells, cell ids are hidden.")

  # Cluster annotation
  # clusters = get_clusters(rcongas_example)

  annotation_columns = NULL
  if(!all(is.null(clusters)))
  {
    annotation_columns = data.frame(Cluster = clusters$cluster, stringsAsFactors = FALSE)
    rownames(annotation_columns) = clusters$cell
  }

  # Call
  pheatmap(
    x,
    show_rownames = show_rownames,
    show_colnames = show_colnames,
    color = colorRampPalette(rev(
      RColorBrewer::brewer.pal(n = 9, "Spectral")
    ))(100),
    annotation_col = annotation_columns,
    annotation_colors = list(
      Cluster = get_clusters_colors(annotation_columns$Cluster)
    )
  )
}
#
#
# plot_raw_data(
#   get_input_raw_data(x),
#   genes =  c('APC', 'KRAS', "TP53")
# )
