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


# Get all analysis
DE_table = get_DE_table(x, cut_pvalue = 1, cut_lfc = 0)

cut_pvalue = 0.001
cut_lfc = 0.25
annotate_top = 10

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
    class = case_when(
      p & lu ~ "Upregulated",
      p & ld ~ "Downregulated",
      TRUE ~ "ns"
    )
  )

# Colors
colors = c("gainsboro", "indianred", 'steelblue')
names(colors) = c('ns', 'Upregulated', 'Downregulated')

# Top
top = DE_table %>%
  dplyr::filter(class != 'ns') %>%
  dplyr::arrange(class, (p_val_adj)) %>%
  dplyr::group_by(class) %>%
  dplyr::filter(row_number() <= annotate_top)


# Chromosomes colours
chromosomes = paste0("chr", c(1:22, "X", "Y"))

chr_colors = c(
  RColorBrewer::brewer.pal(n = 9, name = "Set1"),
  RColorBrewer::brewer.pal(n = 8, name = "Set2"),
  RColorBrewer::brewer.pal(n = 7, name = "Dark2")
)
names(chr_colors) = chromosomes
chr_colors = alpha(chr_colors, .5)

# Plot
ggplot(DE_table,
       aes(x = avg_log2FC, y = -log(p_val_adj))) +
  CNAqc:::my_ggplot_theme() +
  scale_color_manual(values = colors) +
  geom_vline(
    xintercept = c(-cut_lfc, cut_lfc), size = .3, linetype = 'dashed'
  ) +
  geom_hline(
    yintercept = c(-cut_lfc, cut_lfc), size = .3, linetype = 'dashed'
  ) +
  geom_point(aes(color = class)) +
  guides(color = guide_legend("Gene")) +
  labs(
    y = "-log(p)",
    x = "LFC"
  ) +
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
  guides(fill = guide_legend('Chr'))


volcano_plot(rcongas_example)

