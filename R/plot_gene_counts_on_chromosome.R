#' plot_gene_counts_on_genome
#' 
#' @description 
#' 
#' Plots gene per cell counts on the genome.
#'
#' @param x Input object.
#' @param cutoff_quantile_sd Shows only counts with standard deviation above this quantile (default >95%); 
#' it should avoid plotting a number of low-frequency data points.
#' @param normalise Normalise by library size (\code{"ls"}) and/or by input copy number (\code{"ls"}).
#' @param ... Parameters forwarded to \code{get_genes_mapped_to_segments} and \code{get_input_raw_data}.
#'
#' @return A ggplot figure with multiple plots.
#' @export
#'
#' @examples
#' 
#' plot_gene_counts_on_genome(Rcongas::congas_example)
plot_gene_counts_on_genome = function(x,
                                      cutoff_quantile_sd = .95,
                                      normalise = c("ls", "cn"),
                                      chromosomes = paste0("chr", c(1:22, "X", "Y")),
                                      assembly = TRUE,
                                      ...)
{

  if (has_inference(x))
  {
    input_RNA = get_input_raw_data(x, add_locations = TRUE, add_clusters = TRUE,  ...)
  }
  else
  {
    input_RNA = get_input_raw_data(x, add_locations = TRUE, add_clusters = FALSE,  ...)
    
    input_RNA$cluster = "Not available"
  }
  
  if (any(is.na(input_RNA$chr)))
    warning("Some genes are not mapped to the, check raw data.")
  
  # Normalise by library size (ls)
  if ('ls' %in% normalise) {
    ls_factors = input_RNA %>% group_by(cell) %>% summarise(theta = sum(n))
    ls_factors = pio:::nmfy(ls_factors$cell, ls_factors$theta)
    
    input_RNA = input_RNA %>% mutate(n = n / ls_factors[cell])
  }
  
  # Normalise by Copy Number (cn)
  if ('cn' %in% normalise) {

    is = get_input_segmentation(x)
    ls_factors = input_RNA %>% group_by(cell) %>% summarise(theta = sum(n))
    ls_factors = pio:::nmfy(ls_factors$cell, ls_factors$theta)
    
    input_RNA = input_RNA %>% mutate(n = n / ls_factors[cell])
  }
  
  # input_RNA$chr %>% table
  
  chrs = input_RNA %>% pull(chr) %>% unique(na.rm = T)
  chrs = chrs[!is.na(chrs)]

  plots = suppressWarnings(
    lapply(
      chrs[chrs %in% chromosomes],
      aux_gene_spatial_plot,
      x = x,
      input_RNA = input_RNA,
      cutoff_quantile_sd = cutoff_quantile_sd,
      normalise = normalise
    )
  )
  
  if(!assembly) return(plots)
  
  n = plots %>% length %>% sqrt %>% ceiling()
  
  fig = ggarrange(
    plotlist = plots,
    ncol = n,
    nrow = ifelse(n * (n - 1) < plots %>% length, n, n - 1),
    common.legend = TRUE,
    legend = 'bottom'
  )
  
  fig
}

aux_gene_spatial_plot = function(x,
                                 chr,
                                 input_RNA,
                                 cutoff_quantile_sd,
                                 normalise,
                                 second_plot = TRUE)
{
  onchr = input_RNA %>%
    filter(chr == !!chr)
  
  input_RNA%>% group_by(gene, chr) 
  
  onchr_sd = onchr %>% group_by(gene) %>% summarise(sd = sd(n))
  
  qc = quantile(onchr_sd$sd, cutoff_quantile_sd, na.rm = T)
  
  # Scale by reference
  ref = get_reference_genome(x)
  
  off = NULL
  
  if (ref %in% c("hg19", "GRCh37"))
    off = CNAqc::chr_coordinates_hg19 %>% filter(chr == !!chr) %>% pull(from)
  
  if (ref %in% c("hg38", "GRCh38"))
    off = CNAqc::chr_coordinates_GRCh38 %>% filter(chr == !!chr) %>% pull(from)
  
  if (is.null(off))
    stop("Unknown reference?")
  
  # onchr$from =  onchr$from + off
  
  ranges_y = c(min(onchr %>% filter(n > qc) %>% pull(n)),
               max(onchr %>% filter(n > qc) %>% pull(n)))
  
  pl = CNAqc:::blank_genome(chromosomes = chr) +
    geom_point(
      data = onchr %>% filter(n > qc),
      aes(x = from + off, y = n, color = cluster),
      size = .1
    )  +
    scale_color_brewer(palette = "Set1") +
    labs(y = paste('Counts (', paste(normalise, collapse = ','), ')')) +
    ylim(
      onchr %>% filter(n > qc) %>% pull(n) %>%  min,
      onchr %>% filter(n > qc) %>% pull(n) %>%  max
    ) +
    facet_wrap( ~ cluster, ncol = 1) +
    guides(color = FALSE) +
    theme(axis.text.y = element_text(angle = 0))
  
  if(!second_plot) return(pl)
  
  to_high = highlights(x) %>% filter(highlight)
  
  tp = CNAqc:::blank_genome(chromosomes = chr) +
    geom_segment(
      data = get_input_segmentation(x, chromosomes = chr) %>% 
        mutate(highlight = segment_id %in% to_high$segment_id),
      aes(
        y = ploidy_real,
        x = from + off,
        yend = ploidy_real,
        xend = to + off,
        color = highlight
      ),
      size = 1.5
    ) +
    theme_linedraw() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank()) +
    labs(x = NULL, y = NULL) +
    ylim(
      -0.1 + get_input_segmentation(x, chromosomes = chr)$ploidy_real %>% min,
      0.1 + get_input_segmentation(x, chromosomes = chr)$ploidy_real %>% max
    ) +
    # labs(title = paste0(chr, " (quantile stdev >", cutoff_quantile_sd, ')'))  +
    facet_wrap( ~ paste(chr, " (quantile stdev >", cutoff_quantile_sd, ')')) +
    guides(color = FALSE) +
    scale_color_manual(values = c(`TRUE` = 'indianred3', `FALSE` = 'black'))
  
  cowplot::plot_grid(
    tp,
    pl,
    rel_heights = c(1, 3),
    axis = 'lr',
    align = 'v',
    ncol = 1
  )
}
