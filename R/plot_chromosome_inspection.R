# x = readRDS("~/Downloads/Ric/I09_1_003_1014PDB_S3/I09_1_003_1014PDB_S3_fit.rds")
# plot_gene_counts_on_genome(x, chromosomes = 'chr12')
#
# chr = 'chr12'


#' Title
#'
#' @param x
#' @param chr
#' @param normalise
#' @param top
#' @param n_sorrounding
#'
#' @return
#' @export
#'
#' @examples
plot_chromosome_inspection = function(x,
                                      chr,
                                      normalise = 'ls',
                                      top = 5,
                                      n_sorrounding = 5,
                                      ...)
{
  get_ref_off = function(x) {
    # Scale by reference
    ref = get_reference_genome(x)
    
    off = NULL
    
    if (ref %in% c("hg19", "GRCh37"))
      off = CNAqc::chr_coordinates_hg19 %>% filter(chr == !!chr) %>% pull(from)
    
    if (ref %in% c("hg38", "GRCh38"))
      off = CNAqc::chr_coordinates_GRCh38 %>% filter(chr == !!chr) %>% pull(from)
    
    if (is.null(off))
      stop("Unknown reference?")
    
    off
  }
  
  # Required data
  if (has_inference(x))
  {
    input_RNA = get_input_raw_data(x, add_locations = TRUE, add_clusters = TRUE,  ...)
  }
  else
  {
    input_RNA = get_input_raw_data(x, add_locations = TRUE, add_clusters = FALSE,  ...)
    
    input_RNA$cluster = "Not available"
  }
  
  # Normalise by library size (ls)
  if ('ls' %in% normalise) {
    ls_factors = input_RNA %>% group_by(cell) %>% summarise(theta = sum(n))
    ls_factors = pio:::nmfy(ls_factors$cell, ls_factors$theta)
    
    input_RNA = input_RNA %>% mutate(n = n / ls_factors[cell])
  }
  
  input_RNA = input_RNA %>% filter(chr %in% !!chr)
  
  
  # Top variability across genes
  rank = input_RNA %>%
    group_by(gene, cluster) %>%
    summarise(mean = mean(n), stdev = sd(n)) %>%
    arrange(desc(stdev * mean)) %>%
    ungroup()
  
  rank_top = rank %>% distinct(gene) %>% filter(row_number() <= top) %>% pull(gene)
  
  plot_toprank = ggplot() +
    geom_point(
      data = rank %>% filter(!(gene %in% rank_top)),
      aes(x = mean, y = stdev),
      color = 'gray',
      size = .5
    ) +
    geom_point(
      data = rank %>% filter((gene %in% rank_top)),
      aes(x = mean, y = stdev, color = cluster),
      size = 1
    ) +
    ggrepel::geom_text_repel(data = rank %>% filter((gene %in% rank_top)),
                             aes(x = mean, y = stdev, label = gene),
                             size = 2) +
    scale_color_brewer(palette = 'Set1') +
    CNAqc:::my_ggplot_theme() +
    labs(
      x = "Mean",
      y = "Standard deviation",
      title = paste0("Top-", top, " on chromosome ", chr)
    ) +
    guides(color = guide_legend("Cluster", nrow = 1))
  
  # input_RNA %>% filter(gene %in% rank) %>%   pull(from) %>% unique() %>% sort
  #
  # wh = rank[1]
  # wh_from = (input_RNA %>% filter(gene == wh) %>% pull(from))[1]
  # wh_to = (input_RNA %>% filter(gene == wh) %>% pull(to))[1]
  #
  # with_MB_function = function(x, f, t) {
  #   x %>% mutate(within = (from > f - 8e6) & (to < t + 8e6))
  # }
  # plot_top_rank_wide = aux_gene_spatial_plot(x,
  #                                 chr,
  #                                 input_RNA,
  #                                 0,
  #                                 normalise = 'ls',
  #                                 second_plot = F) +
  #   geom_point(
  #     data = input_RNA %>% filter(gene %in% rank) %>% distinct(gene, from),
  #     aes(x = from + get_ref_off(x)),
  #     y = Inf,
  #     size = 2,
  #     pch = 18
  #   )
  
  annot = input_RNA %>%
    filter(gene %in% rank_top) %>%
    group_by(gene) %>%
    arrange(desc(n)) %>%
    filter(row_number() == 1)
  
  plot_top_rank_wide = aux_gene_spatial_plot(
    x,
    chr,
    input_RNA %>% mutate(cluster = "Outliers"),
    0,
    normalise = 'ls',
    second_plot = F
  ) +
    ggrepel::geom_text_repel(
      data = annot %>% mutate(cluster = "Outliers"),
      aes(
        x = from + get_ref_off(x),
        y = n,
        label = gene
      ),
      size = 2,
      inherit.aes = FALSE
    )
  
  
  
  with_ngenes_function = function(x, g, n) {
    rg = input_RNA %>% arrange(from) %>% distinct(gene, from)
    iOf = which(rg$gene == g)
    wg = rg$gene[(iOf - n):(iOf + n)]
    
    x %>% mutate(within = gene %in% wg)
  }
  
  
  # input_RNA = input_RNA %>% with_MB_function(f = wh_from, t = wh_to)
  
  # i_wh = input_RNA %>% with_ngenes_function(g = wh, n = 10) %>%  filter(within)
  
  special_points = Reduce(bind_rows,
                          lapply(rank_top, function(g)
                          {
                            w = input_RNA %>%
                              with_ngenes_function(g = g, n = n_sorrounding) %>%
                              filter(within) %>%
                              select(n, gene, cluster)
                            
                            nh = w$gene %>% unique %>% length()
                            
                            w %>%
                              mutate(G = paste0(g, " +/- ",!!n_sorrounding))
                          }))
  
  # special_points$gene %>% unique %>% length
  
  nonspecial_points = input_RNA %>%
    filter(!(gene %in% special_points$gene)) %>%
    select(n, gene, cluster)
  
  # nonspecial_points$gene %>% unique %>% length
  
  nonspecial_points = nonspecial_points %>%
    mutate(G = paste0("Others (n = ", nonspecial_points$gene %>% unique %>% length,
                      ')'))
  
  
  # ggplot(special_points %>% filter(gene == 'FTL'),
  #        aes(n, fill = cluster)) +
  #   geom_histogram(bins = 30) +
  #   # scale_y_log10() +
  #   xlim(0, .01 + max(special_points %>% filter(gene == 'FTL') %>% pull(n)))
  
  plot_counts = ggplot(special_points %>% bind_rows(nonspecial_points),
                       aes(n, fill = cluster)) +
    geom_histogram(bins = 30) +
    scale_fill_brewer(palette = 'Set1') +
    CNAqc:::my_ggplot_theme() +
    scale_y_log10() +
    facet_wrap(~ G, scales = 'free_x') +
    labs(
      subtitle = paste0(
        "Chromosome ",
        chr,
        ': ',
        input_RNA %>%  distinct(gene) %>% pull(gene) %>% length,
        " genes"
      ),
      x = paste0('Counts (', normalise, ')'),
      y = 'Observations'
    ) +
    guides(fill = guide_legend("Cluster", nrow = 1)) +
    coord_cartesian(clip = 'off') 
  
  
  cowplot::plot_grid(
    plot_toprank,
    cowplot::plot_grid(
      plot_top_rank_wide,
      plot_counts,
      nrow = 2,
      ncol = 1,
      axis = 'lr',
      align = 'v',
      rel_heights = c(1, 2)
    ),
    nrow = 1,
    ncol = 2
  )
  
}




