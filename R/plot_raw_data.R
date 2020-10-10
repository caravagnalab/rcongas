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

plot_raw_data = function(x, genes, prompt = TRUE)
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
    cli::cli_alert_warning("You are about to plot {.field {ng}} genes and {.field {nc}} cells using pheatmap, this may take a while...")

    ask = menu(c("Yes", "No"), title="Do you want this?")

    if(ask == 2) return(CNAqc:::eplot())
  }

  require(pheatmap)

  # colnames/rownames smart show
  show_rownames = show_colnames = FALSE
  fontsize_row = fontsize_col = 2

  if(ng < 100) show_rownames = TRUE
  if(nc < 100) show_colnames = TRUE

  if(!show_rownames) cli::cli_alert_info("With more than 100 genes, gene names are hidden.")
  if(!show_colnames) cli::cli_alert_info("With more than 100 cells, cell ids are hidden.")

  pheatmap(
    x,
    show_rownames = show_rownames,
    show_colnames = show_colnames,
    # fontsize_row = 2,
    # fontsize_col = 2,
    color = colorRampPalette(rev(
      RColorBrewer::brewer.pal(n = 9, "Spectral")
    ))(100)
  )
}
#
#
# plot_raw_data(
#   get_input_raw_data(x),
#   genes =  c('APC', 'KRAS', "TP53")
# )
