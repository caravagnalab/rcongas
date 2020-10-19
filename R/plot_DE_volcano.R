volcano_plot <-
  function(X,
           pCutoff = 10e-5,
           FCcutoff = 1,
           colorChr = NULL,
           ...) {
    table <-  get_DE_table(X, cut_pvalue = 1, cut_lfc = 0)

    if (!is.null(colorChr)) {
      nChr <- length(colorChr)
      cols <- RColorBrewer::brewer.pal(8, name = "Dark2")
      cols <- c(cols[1:nChr], "grey")
      names(cols) <- c(colorChr, "other chrs")
      colorChr <- cols[table$chr]
      colorChr[is.na(colorChr)] <-  "grey"
      for (j in seq(nChr + 1)) {
        names(colorChr)[colorChr == cols[j]] <- names(cols)[j]
      }
    }
    group1 <- X$DE$params$clone1
    group2 <-
      ifelse(is.null(X$DE$params$clone2), "all", X$DE$params$clone2)
    test_method <-  X$DE$params$method
    EnhancedVolcano::EnhancedVolcano(
      table,
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


#' Title
#'
#' @param x
#' @param annotate_top
#' @param cut_pvalue
#' @param cut_lfc
#'
#' @return
#' @export
#'
#' @examples
plot_DE_volcano <-
  function(x,
           chromosomes = paste0("chr", c(1:22, "X", "Y")),
           annotate_top = 5,
           cut_pvalue = 0.001,
           cut_lfc = 0.25)
  {
    # Get all analysis
    DE_table = get_DE_table(x,
                            chromosomes = chromosomes,
                            cut_pvalue = 1,
                            cut_lfc = 0)

    range_lfc = max(abs(min(DE_table$avg_log2FC)), max(DE_table$avg_log2FC))

    # I would not do this
    # DE_table = DE_table %>%
    #   dplyr::mutate(
    #     p = p_val_adj < cut_pvalue,
    #     l = abs(avg_log2FC) > cut_lfc,
    #     class = case_when(
    #       p & l ~ "p-value & LFC",
    #       p & !l ~ "p-value",
    #       !p & l  ~ "LFC",
    #       !p & !l ~ "ns"
    #     )
    #   )

    # I would do this
    DE_table = DE_table %>%
      dplyr::mutate(
        p = p_val_adj < cut_pvalue,
        lu = avg_log2FC > cut_lfc,
        ld = avg_log2FC < -cut_lfc,
        class = case_when(p & lu ~ "Upregulated",
                          p & ld ~ "Downregulated",
                          TRUE ~ "n.s.")
      )

    # Colors
    colors = c("gainsboro", "darkorchid4", 'goldenrod1')
    names(colors) = c('n.s.', 'Upregulated', 'Downregulated')

    # Top
    top = DE_table %>%
      dplyr::filter(class != 'n.s.') %>%
      dplyr::arrange(class, (p_val_adj)) %>%
      dplyr::group_by(class) %>%
      dplyr::filter(row_number() <= annotate_top)


    # Chromosomes colours
    chromosomes = paste0("chr", c(1:22, "X", "Y"))

    # chr_colors = c(
    #   RColorBrewer::brewer.pal(n = 9, name = "Set1"),
    #   RColorBrewer::brewer.pal(n = 8, name = "Set2"),
    #   RColorBrewer::brewer.pal(n = 7, name = "Dark2")
    # )

    # chr_colors = c(
    #   ggsci::pal_futurama()(6),
    #   ggsci::pal_jco()(6),
    #   ggsci::pal_jama()(6),
    #   ggsci::pal_tron()(6)
    # )
    # chr_colors = sample(chr_colors, length(chr_colors), replace = T)

    chr_colors =
      c(
        "#F2F3F4",
        ggplot2::alpha("#222222", .4),
        "#F3C300",
        "#875692",
        "#F38400",
        "#A1CAF1",
        "#BE0032",
        "#C2B280",
        "#848482",
        "#008856",
        "#E68FAC",
        "#0067A5",
        "#F99379",
        "#604E97",
        "#F6A600",
        "#B3446C",
        "#DCD300",
        "#882D17",
        "#8DB600",
        "#654522",
        "#E25822",
        "#2B3D26",
        "indianred3",
        "steelblue"
      )

    names(chr_colors) = chromosomes

    group1 <- x$DE$params$clone1
    group2 <-
      ifelse(is.null(x$DE$params$clone2), "all", x$DE$params$clone2)

    test_method <-
      case_when(x$DE$params$method == "Wilcox" ~ "Wilcoxon test",
                TRUE ~ x$DE$params$method)

    # Plot
    ggplot(DE_table,
           aes(x = avg_log2FC, y = -log(p_val_adj))) +
      CNAqc:::my_ggplot_theme() +
      scale_color_manual(values = colors) +
      geom_vline(
        xintercept = c(-cut_lfc, cut_lfc),
        size = .3,
        linetype = 'dashed'
      ) +
      geom_hline(yintercept = -log(c(cut_pvalue)),
                 size = .3,
                 linetype = 'dashed') +
      geom_point(aes(color = class)) +
      geom_point(data = top, aes(color = class, fill = chr)) +
      guides(color = guide_legend("", ncol = 1)) +
      labs(y = "-log(p)",
           x = "LFC") +
      ggrepel::geom_label_repel(
        data = top,
        aes(label = gene, fill = chr),
        segment.size = 0.2,
        size = 2,
        nudge_x = .5,
        nudge_y = .5,
        color = 'black'
      ) +
      scale_fill_manual(values = chr_colors) +
      guides(fill = guide_legend("", override.aes = list(color = NA))) +
      labs(
        title = paste0("DE (", group1, ' vs ', group2, " with ", test_method, ")"),
        caption = paste0(
          "Adjusted p-value cut ",
          cut_pvalue,
          ", LFC cut ",
          cut_lfc,
          '.'
        )
      ) +
      xlim(-range_lfc, range_lfc)

  }

#
# volcano_plot(rcongas_example)

# plot_DE_volcano(rcongas_example, annotate_top = 10, cut_pvalue = 1e-6)
