#' Filter genes from an RNA tibble.
#'
#' @param x The tibble for RNA.
#' @param what A code string. If this contains "r" ribosomal genes are removed.
#' If this contains "m" mitocondrial genes are removed. Both are matched by
#' a regular expression.
#' @param specials A list of genes to remove a priori.
#'
#' @return The RNA tibble without the required gene entries.
#' @export
#'
#' @examples
#' data('example_input')
#' filter_known_genes(example_input$rna)
filter_known_genes = function(x, what = 'rm', specials = 'MALAT1')
{
  if(!('gene' %in% names(x)))
  {
    cli::cli_alert_danger("The {.field {'gene'}} coulumn is missing from input, returning input data.")
    return(x)
  }
    
  all_genes = x$gene %>% unique
  codes = strsplit(what, split = '')[[1]]
  
  if ('r' %in% codes)
  {
    # Ribosomal genes
    cli::cli_h3("Ribosomal")
    
    ribo = all_genes %>% grepl(pattern = 'RP[SL]')
    nribo = ribo %>% sum
    
    cli::cli_alert_info("{.field n = {nribo}} genes found.")
    
    if (nribo > 0) {
      cat(all_genes[ribo] %>% head %>% paste(collapse = ', '), ', ...\n')
      all_genes = setdiff(all_genes, all_genes[ribo])
    }
  }
  
  if ('m' %in% codes)
  {
    # Mitochondrial genes
    cli::cli_h3("Mitochondrial")
    
    mito = all_genes %>% grepl(pattern = 'MT')
    nmito = mito %>% sum
    
    cli::cli_alert_info("{.field {nmito}} genes found.")
    
    if (nmito > 0)
    {
      cat(all_genes[mito] %>% head %>% paste(collapse = ', '), ', ...\n')
      all_genes = setdiff(all_genes, all_genes[mito])
    }
  }
  
  if (length(specials) >= 1)
  {
    # Ribosomal genes
    cli::cli_h3("Specials")
    
    spec = all_genes %>% intersect(specials)
    nspec = spec %>% length
    
    cli::cli_alert_info("{.field n = {nspec}} genes found.\n")
    
    if(nspec > 0) all_genes = setdiff(all_genes, spec)
  }
  
  return(x %>% filter(gene %in% all_genes))
}

#' Filter observed counts by quantile.
#'
#' @param x An input RNA/ATAC dataset where entries are indexeable by genomic 
#' coordinate: "chr", "from" and "to".
#' @param upper_quantile The maximum quantile to determine cuts. If a value
#' is above the quantile it is removed
#'
#' @return The input data with removed entries.
#' @export
#'
#' @examples
#' data('example_input')
#' filter_counts_by_quantile(example_input$rna, upper_quantile = .98)
#' filter_counts_by_quantile(example_input$atac, upper_quantile = .98)
filter_counts_by_quantile = function(x, upper_quantile = .98)
{
  map_quantiles = x %>% 
    group_by(chr, from, to) %>% 
    summarise(
      q_max = quantile(value, upper_quantile),
      # sd = sd(value),
      # median = median(value),
      .groups = 'drop'
    )
  
  x = x %>% 
    left_join(map_quantiles, by = c('chr', 'from', 'to')) %>% 
    mutate(del = value > q_max)
  
  # x %>% arrange(desc(sd))
  # x %>% arrange(desc(value))
  # hist(x %>% filter(chr == 'chr9', to == '140135734', from== '140135234') %>% pull(value), breaks = 100) 
  # 
  # m_atac = quantile(x$value, .98)
  # x %>% filter
  
  # hist(norm_atac$normalisation_factor, breaks = 100)
  
  
  # hist(log(x$sd), breaks = 100) 
  
  # ggplot(x %>% filter(log(x$sd) > 2), aes(x = median, y = sd, color = del)) +
  #   geom_point()
  
  r_cap = x %>% filter(del)
  
  if(nrow(r_cap) > 0)
  {
    cli::cli_h3("Upper quantile {.field {upper_quantile}}")
    
    cli::cli_alert_info("{.field n = {nrow(r_cap)}} entries to remove")
    
    cat("\n")
    r_cap  %>% print
    
    x = x %>% 
      filter(!del)
  }
  
  return(x %>% select(-q_max, -del))
}

# Return list of outlier genes/peaks in each segment.
#' @param segmentation
#' @param modality Either 'RNA' or 'ATAC'
#' @param data Tibble with chr from to value
#' @param normalisation_factors tibble that for each cell contained in \code{data} contains information about normalization facctors. COmputed with Rcongas function
#' \code{auto_normalisation_factor}.
clean_outlers_persegment <- function(segmentation, modality, data, normalisation_factors) {
  if (modality == 'ATAC') {
    data = data %>% mutate(peakId = paste(chr, from, to, sep =':'))
    featureVar = 'peakId'
  } else if (modality == 'RNA') {
    featureVar = 'gene'
  } else {
    stop('Please provide a valid modality, must be either RNA or ATAC')
  }
  
  segmentation = segmentation %>% Rcongas:::idify()
  
  evt_lbs = paste0(modality, '_nonzerovals')
  loc_lbs = ifelse(
    modality == 'RNA',
    paste0(modality, '_genes'),
    paste0(modality, '_peaks')
  )
  cells_lbs = paste0(modality, '_nonzerocells')
  segmentation[[evt_lbs]] = 0
  segmentation[[loc_lbs]] = 0
  segmentation[[cells_lbs]] = 0
  
  data$segment_id = NA
  
  pb <- progress::progress_bar$new(total = nrow(segmentation))
  
  pb$tick(0)
  
  for(i in 1:nrow(segmentation))
  {
    pb$tick()
    
    what_maps = which(
      data$chr == segmentation$chr[i] &
        data$from >= segmentation$from[i] &
        data$to <= segmentation$to[i]
    )
    
    if(length(what_maps) == 0) next;
    
    data$segment_id[what_maps] = segmentation$segment_id[i]
    
    segmentation[[evt_lbs]][i] = what_maps %>% length
    segmentation[[loc_lbs]][i] = data[what_maps, ] %>%
      distinct(chr, from, to) %>%
      nrow()
    
    segmentation[[cells_lbs]][i] = data[what_maps, ] %>%
      pull(cell) %>% unique %>% length
  }
  
  n_na = is.na(data$segment_id) %>% sum()
  nn_na = (data %>% nrow) - n_na
  
  cli::cli_alert("Entries mapped: {.field {nn_na}}, with {.field {n_na}} outside input segments that will be discarded.")
  if(n_na > 0) data = data %>% filter(!is.na(segment_id))
  likelihood = 'NB'
  cli::cli_alert("Using likelihood: {.field {likelihood}}.")
  
  
  data_norm = Rcongas:::normalise_modality(data %>% mutate(modality = !!modality), normalisation_factors)
  
  plist = sapply(unique(data_norm$segment_id),  simplify = F, USE.NAMES = T, function(s) {
    tmp = data_norm %>% filter(segment_id == !!s) %>% group_by(chr, from, to, !!featureVar) %>% mutate(mean_expr = mean(value))
    qmax = quantile(tmp$value, 0.97)
    tmp = tmp %>% mutate(isOutlier = mean_expr > !!qmax)
    
    p = ggplot(tmp, 
               aes(x = from, y = mean_expr, color = isOutlier)) + 
      geom_point(size = 0.5) +#+ ylim(2,150) +
      scale_color_manual(values = c('#1D1D1D', '#B50717')) 
    if (modality == 'RNA') {
      p = p + geom_text(aes(label=ifelse(mean_expr>!!qmax,as.character(!!featureVar),'')),
                hjust=-0.1,vjust=0, size = 2)+ theme_bw() + guides(color = FALSE, size = FALSE, label = F)    
    }
    return(p)
  })
  
  outliers = unlist(sapply(plist, function(p) {
    return(unique(p$data[[featureVar]][p$data$isOutlier]))
  }))
  # View(plist[[1]])
  # outdir = 'figures_congasp/outlier_genes/'
  # p = cowplot::plot_grid(plotlist = plist)
  # ggsave(paste0(outdir, 'segments_ouliers.pdf'), p, width = 15, height = 15)
  data = data %>% filter(!(!!featureVar %in% outliers))
  return(list(persegment_plot = plist, data_cleaned = data))
}

#' #' Cap observed values by quantile.
#' #'
#' #' @param x An input RNA/ATAC dataset where entries are indexeable by genomic 
#' #' coordinate: "chr", "from" and "to".
#' #' @param upper_quantile The maximum quantile to determine cuts. If a value
#' #' is above the quantile it is capped.
#' #'
#' #' @return The input data with capped entries.
#' #' @export
#' #'
#' #' @examples
#' #' data('example_input')
#' #' cap_values_by_quantile(example_input$x_rna, upper_quantile = .98)
#' #' cap_values_by_quantile(example_input$x_atac, upper_quantile = .98)
#' cap_values_by_quantile = function(x, upper_quantile = .98)
#' {
#'   map_quantiles = x %>% 
#'     group_by(chr, from, to) %>% 
#'     summarise(
#'       q_max = quantile(value, upper_quantile),
#'       .groups = 'drop'
#'     )
#'   
#'   x = x %>% 
#'     left_join(map_quantiles, by = c('chr', 'from', 'to')) %>% 
#'     mutate(cap = value > q_max)
#'   
#'   r_cap = x %>% filter(cap)
#'   
#'   if(nrow(r_cap) > 0)
#'   {
#'     cli::cli_h3("Upper quantile {.field {upper_quantile}}")
#'     
#'     cli::cli_alert_info("{.field n = {nrow(r_cap)}} entries to cap")
#'     
#'     cat("\n")
#'     r_cap  %>% print
#'     
#'     x = x %>% 
#'       mutate(value = ifelse(cap, q_max, value))
#'   }
#'   
#'   return(x %>% select(-q_max, -cap))
#' }