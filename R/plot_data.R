#' Plot data.
#'
#' @description General plotting function for data. This function uses a \code{what} parameter
#' to dispatch visualization to a number of internal functions. For the input data one can visualize:
#' 
#' * (\code{what = "histogram"}) a histogram of input values per segment, per modality;
#' * (\code{what = "lineplot"}) a lineplot showing, for all cells and modality, all values drawn as segments on a genome-wide plot.
#' In this case also the input segmentation is visualized;
#' * (\code{what = "heatmap"}) a heatmap showing, for all cells and modality, all values per segment;
#' * (\code{what = "mapping"}) a tile plot reporting the number of RNA genes or ATAC peaks associated to each segment.
#' 
#' Where appropriate, input values are normalized by input factors. In some cases (lineplot), they are also
#' scaled by the number of events mapped to each segment - this makes them comparable across segments.
#' 
#' In most cases all the plot functions return one \code{ggplot} figure. An exception is made for the
#' heatmap visualisation when there are more than one modalities: in that case all the generated figures
#' are assembled into a \code{cowplot} figure.
#'
#' The internal plotting functions can have some parameters in input. Passage of parameters to those functions
#' is done by a top-level ellipsis. The formals of the internal parameters are described in the Plotting
#' vignette of the package, where example runs are shown. Please refer to that to see how to customise input 
#' plots.
#'
#' @param x An object of class \code{rcongasplus}.
#' @param what Any of \code{"histogram"}, \code{"lineplot"}, \code{"heatmap"} or
#' \code{"mapping"}.
#' @param ... Parameters forwarded to the internal plotting functions.
#'
#' @return A \code{ggplot} or \code{cowplot} figure, depending on the number of
#' modalities and plot required.
#'
#' @export
#'
#' @examples
#' 
#' # Formals of all internal functions (see package Plotting vignette)
#'formals(Rcongas:::plot_data_histogram) %>% names()
#'formals(Rcongas:::plot_data_lineplot) %>% names()
#'formals(Rcongas:::plot_data_heatmap) %>% names()
#'formals(Rcongas:::plot_data_mapping) %>% names()
#'
#' data("example_object")
#' 
#' # Data histogram plot (default all segments)  
#' plot_data(example_object, what = 'histogram')
#' 
#' # Subset what to plot
# which_segments = get_input(example_object, what = 'segmentation') %>%
#    filter(row_number() <= 3) %>%
#    pull(segment_id)
# 
# # Pass formal "segments"
# plot_data(example_object, what = 'histogram', segments = which_segments)
#'  
#' # Lineplot segments (default)  
#' plot_data(example_object, what = 'lineplot')
#' 
#' # Data heatmap
#' plot_data(example_object, what = 'heatmap')
#' 
#' # Events mapping per segment
#' plot_data(example_object, what = 'mapping')
plot_data = function(x,
                     what = 'histogram',
                     ...)
{
  x %>% sanitize_obj()
  
  if (what == 'histogram')
    return(x %>% plot_data_histogram(...))
  
  if (what == 'lineplot')
    return(x %>% plot_data_lineplot(...))
  
  if (what == 'heatmap'){
    # This returns a more complex assembly
    all_plots = x %>% plot_data_heatmap(...)
    return(all_plots$figure)
  }
  
  if (what == 'mapping')
    return(x %>% plot_data_mapping())
  
  stop("Unrecognised 'what': use any of 'histogram', 'lineplot', 'heatmap' or 'mapping'.")
}

plot_data_histogram = function(x,
                               segments = get_input(x, what = 'segmentation') %>% pull(segment_id),
                               whichfacet = ggplot2::facet_wrap
                               )
{
  stats_data = stat(x, what = 'data')
  
  # Data stats
  subtitle = paste0(
    "RNA (",
    stats_data$ncells_RNA,
    " cells, ",
    stats_data$rna_genes,
    " genes) ATAC (",
    stats_data$ncells_ATAC,
    " cells, ",
    stats_data$atac_peaks,
    " peaks)"
  )
  
  # What to print, enhanced facets
  what_rna_lik = case_when(
    stats_data$rna_dtype == "NB" ~ "Negative Binomial",
    stats_data$rna_dtype == "P" ~ "Poisson",
    stats_data$rna_dtype == "G" ~ "Gaussian"
  )
  
  what_atac_lik = case_when(
    stats_data$atac_dtype == "NB" ~ "Negative Binomial",
    stats_data$atac_dtype == "P" ~ "Poisson",
    stats_data$atac_dtype == "G" ~ "Gaussian"
  )
  
  # Normalisation factors
  norm_factors = get_input(x, what = 'normalisation')
  
  # RNA_data values
  what_RNA = get_input(x, what = 'data') %>%
    filter(modality == "RNA")
  
  if (nrow(what_RNA) > 0  && stats_data$rna_dtype != "G")
    what_RNA = normalise_modality(what_RNA, norm_factors %>% filter(modality == "RNA"))
  
  # ATAC_data values
  what_ATAC = get_input(x, what = 'data') %>%
    filter(modality == "ATAC")
  
  if (nrow(what_ATAC) > 0 && stats_data$atac_dtype != "G")
    what_ATAC = normalise_modality(what_ATAC, norm_factors %>% filter(modality == "ATAC"))
  
  what = bind_rows(what_RNA, what_ATAC) %>%
    filter(segment_id %in% segments) %>%
    mutate(modality = case_when(
      modality == "RNA" ~ paste(modality, '(', what_rna_lik, ')'),
      modality == "ATAC" ~ paste(modality, '(', what_atac_lik, ')')
    )) 

  # Call
  what %>%
    # deidify() %>% 
    # mutate(segment_id = factor(chr, levels = gtools::mixedsort(chr) %>% unique)) %>% 
    mutate(segment_id = factor(segment_id, levels = gtools::mixedsort(segment_id) %>% unique)) %>%
    ggplot(aes(value, fill = modality)) +
    geom_histogram(bins = 70) +
    whichfacet(segment_id ~ modality, scales = 'free') +
    # facet_wrap(chr ~ modality, scales = 'free') +
    labs(title = x$description,
         subtitle = subtitle) +
    guides(fill = FALSE) +
    theme_linedraw(base_size = 9) +
    scale_fill_manual(values = modality_colors(what$modality %>% unique)) +
    labs(x = "Input",
         y = 'Observations') +
    theme(strip.text.y.right = element_text(angle = 0)) 
}

plot_data_lineplot = function(x, 
                              segments = get_input(x, what = 'segmentation') %>% pull(segment_id), 
                              alpha = 1)
{
  stats_data = stat(x, what = 'data')
  
  # Data stats
  subtitle = paste0(
    "RNA (",
    stats_data$ncells_RNA,
    " cells, ",
    stats_data$rna_genes,
    " genes) ATAC (",
    stats_data$ncells_ATAC,
    " cells, ",
    stats_data$atac_peaks,
    " peaks)"
  )
  
  # What to print, enhanced facets
  what_rna_lik = case_when(
    stats_data$rna_dtype == "NB" ~ "Negative Binomial",
    stats_data$rna_dtype == "P" ~ "Poisson",
    stats_data$rna_dtype == "G" ~ "Gaussian"
  )
  
  what_atac_lik = case_when(
    stats_data$atac_dtype == "NB" ~ "Negative Binomial",
    stats_data$atac_dtype == "P" ~ "Poisson",
    stats_data$atac_dtype == "G" ~ "Gaussian"
  )
  
  # Normalisation factors
  norm_factors = get_input(x, what = 'normalisation')
  
  # RNA_data values
  what_RNA = get_input(x, what = 'data') %>%
    filter(modality == "RNA") %>%
    mutate(alpha = alpha,
           size = .5)
  
  if (nrow(what_RNA) > 0 && stats_data$rna_dtype != "G")
    what_RNA = normalise_modality(what_RNA, norm_factors %>% filter(modality == "RNA"))
  
  # ATAC_data values
  what_ATAC = get_input(x, what = 'data') %>%
    filter(modality == "ATAC") %>%
    mutate(alpha = alpha,
           size = .5)
  
  if (nrow(what_ATAC) > 0 && stats_data$atac_dtype != "G")
    what_ATAC = normalise_modality(what_ATAC, norm_factors %>% filter(modality == "ATAC"))
  
  # Normalise observed values by number of events per segment to make them comparable
  if (nrow(what_RNA) > 0)
  {
    cli::cli_alert("Scaling RNA observed values by number of genes mapped per segment.")
    
    what_RNA = what_RNA %>%
      left_join(x %>% get_input(what = 'segmentation') %>% select(segment_id, RNA_genes),
                by = 'segment_id') %>%
      mutate(value = value / RNA_genes)
  }
  
  if (nrow(what_ATAC) > 0)
  {
    cli::cli_alert("Scaling ATAC observed values by number of peaks mapped per segment.")
    
    what_ATAC = what_ATAC %>%
      left_join(x %>% get_input(what = 'segmentation') %>% select(segment_id, ATAC_peaks),
                by = 'segment_id') %>%
      mutate(value = value / ATAC_peaks)
  }
  
  # Pool all modalities plus the segmentation, apply filters
  what = bind_rows(
    what_RNA %>% deidify(),
    what_ATAC %>% deidify(),
    get_input(x, what = 'segmentation') %>%
      mutate(modality = " Input segmentation") %>% # add space for factors ordering
      rename(value = copies) %>%
      mutate(alpha = 1,
             size = 2)
  ) %>%
    idify() %>%
    filter(segment_id %in% segments) %>%
    mutate(modality = case_when(
      modality == "RNA" ~ paste(modality, '(', what_rna_lik, ')'),
      modality == "ATAC" ~ paste(modality, '(', what_atac_lik, ')'),
      TRUE ~ modality
    ))
  
  # Shift CNA coordinates to absolute
  what = CNAqc:::relative_to_absolute_coordinates(x,
                                                  what)
  
  # Genome plot call
  blank_genome(ref = x$reference_genome,
               chromosomes = what$chr %>% unique()) +
    geom_segment(
      data = what,
      aes(
        x = from,
        xend = to,
        y = value,
        yend =  value,
        color = modality,
        alpha = alpha,
        size = size
      ),
      inherit.aes = FALSE
    ) +
    facet_wrap(~ modality, ncol = 1, scales = 'free_y') +
    guides(color = FALSE) +
    scale_color_manual(values = modality_colors(what$modality %>% unique)) +
    labs(title = x$description,
         subtitle = subtitle) +
    theme_linedraw(base_size = 9) +
    labs(x = "Chromosomes",
         y = 'Input',
         caption = "Per segment values are normalised by number of mapped RNA genes/ATAC peaks.") +
    guides(alpha = FALSE,
           size = FALSE)
}

plot_data_heatmap = function(x, segments = get_input(x, what = 'segmentation') %>% pull(segment_id))
{
  stats_data = stat(x, what = 'data')
  
  # Data stats
  subtitle_RNA = paste0("RNA (",
                        stats_data$ncells_RNA,
                        " cells, ",
                        stats_data$rna_genes,
                        " genes)")
  
  subtitle_ATAC = paste0("ATAC (",
                         stats_data$ncells_ATAC,
                         " cells, ",
                         stats_data$atac_peaks,
                         " peaks)")
  
  # What to print, enhanced facets
  what_rna_lik = case_when(
    stats_data$rna_dtype == "NB" ~ "Negative Binomial",
    stats_data$rna_dtype == "P" ~ "Poisson",
    stats_data$rna_dtype == "G" ~ "Gaussian"
  )
  
  what_atac_lik = case_when(
    stats_data$atac_dtype == "NB" ~ "Negative Binomial",
    stats_data$atac_dtype == "P" ~ "Poisson",
    stats_data$atac_dtype == "G" ~ "Gaussian"
  )
  
  # Normalisation factors
  norm_factors = get_input(x, what = 'normalisation')
  
  # RNA_data values
  what_RNA = get_input(x, what = 'data') %>%
    filter(modality == "RNA")
  
  if (nrow(what_RNA) > 0 && stats_data$rna_dtype != "G")
    what_RNA = normalise_modality(what_RNA, norm_factors %>% filter(modality == "RNA"))
  
  # ATAC_data values
  what_ATAC = get_input(x, what = 'data') %>%
    filter(modality == "ATAC")
  
  if (nrow(what_ATAC) > 0 && stats_data$atac_dtype != "G")
    what_ATAC = normalise_modality(what_ATAC, norm_factors %>% filter(modality == "ATAC"))
  
  # Pool all modalities plus the segmentation, apply filters
  what = bind_rows(what_RNA %>% deidify(),
                   what_ATAC %>% deidify()) %>%
    idify() %>%
    filter(segment_id %in% segments) %>%
    mutate(modality = case_when(
      modality == "RNA" ~ paste(modality, '(', what_rna_lik, ')'),
      modality == "ATAC" ~ paste(modality, '(', what_atac_lik, ')'),
      TRUE ~ modality
    ))
  
  # Genome plot call
  RNA_plot = ATAC_plot = NULL
  
  if (nrow(what_RNA) > 0)
    # RNA
  {
    RNA_plot = ggplot(what_RNA %>% filter(segment_id %in% segments)) +
      geom_tile(aes(x = segment_id, y = cell, fill = value)) +
      theme_linedraw(base_size = 9) +
      labs(x = "Segments",
           y = subtitle_RNA) +
      theme(axis.text.y = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),) +
      guides(fill = guide_colorbar(paste('RNA (', what_rna_lik, ')'),
                                   barheight = unit(3, 'cm')))
    
    if (stats_data$rna_dtype == "G")
      RNA_plot = RNA_plot +
        scale_fill_gradient2(low = "steelblue", high = 'indianred3')
    else
      RNA_plot = RNA_plot +
        scale_fill_distiller(palette = 'GnBu', direction = 1)
  }
  
  if (nrow(what_ATAC) > 0)
    # RNA
  {
    ATAC_plot = ggplot(what_ATAC %>% filter(segment_id %in% segments)) +
      geom_tile(aes(x = segment_id, y = cell, fill = value)) +
      theme_linedraw(base_size = 9) +
      labs(x = "Segments",
           y = subtitle_ATAC) +
      theme(axis.text.y = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),) +
      guides(fill = guide_colorbar(paste('ATAC (', what_atac_lik, ')'),
                                   barheight = unit(3, 'cm')))
    
    if (stats_data$atac_dtype == "G")
      ATAC_plot = ATAC_plot +
        scale_fill_gradient2(low = "steelblue", high = 'indianred3')
    else
      ATAC_plot = ATAC_plot + scale_fill_distiller(palette = 'GnBu', direction = 1)
  }
  
  # Figure assembly
  if (stats_data$nmodalities == 1)
  {
    if ("RNA" %in% stats_data$modalities)
      return(list(figure = RNA_plot, plots = list(RNA_plot, ATAC_plot)))
    if ("ATAC" %in% stats_data$modalities)
      return(list(figure = ATAC_plot, plots = list(RNA_plot, ATAC_plot)))
  }
  else
  {
    # Proportional to number of cells per assay
    rel_height = c(stats_data$ncells_ATAC, stats_data$ncells_RNA) /
      (stats_data$ncells_RNA + stats_data$ncells_ATAC)
    
    # Place ATAC on top of RNA, remove ATAC segment ids, move labels
    figure =
      cowplot::plot_grid(
        ATAC_plot + theme(axis.text.x = element_blank()) + labs(x = NULL),
        RNA_plot,
        rel_heights = rel_height,
        ncol = 1,
        align = 'v',
        axis = 'lr'
      )
    
    return(list(figure = figure, plots = list(RNA_plot, ATAC_plot)))
  }
}

plot_data_mapping = function(x)
{
  mapping = x %>%
    get_input(what = 'segmentation')
  
  w = c("chr", "segment_id")
  if ("RNA" %in% stat(x)$modalities)
    w = c(w, "RNA_genes")
  if ("ATAC" %in% stat(x)$modalities)
    w = c(w, "ATAC_peaks")
  
  order_segments = mapping %>% 
    mutate(E = ATAC_peaks + RNA_genes) %>% 
    arrange(chr, (E)) %>% 
    pull(segment_id)
  
  reshape2::melt(mapping[w], id = c("segment_id", "chr")) %>%
    # deidify() %>%
    # mutate(segment_id = paste(from, ':', to)) %>%
    mutate(chr = factor(chr, levels = gtools::mixedsort(chr) %>% unique)) %>% 
    mutate(value = ifelse(value == 0, NA, value)) %>%
    ggplot(aes(
      y = factor(segment_id, levels = order_segments),
      x = variable,
      fill = value,
      label = value
    )) +
    geom_tile() +
    geom_text(size = 2) +
    theme_linedraw(base_size = 9) +
    theme(strip.text.y.right = element_text(angle = 0)) +
    scale_fill_distiller(palette = 'Spectral',
                         direction = -1,
                         na.value = 'gray') +
    facet_grid(chr ~ "Modality", scales = "free_y", space = "free_y") +
    labs(x = "Modality", y = "Segments") 
    # scale_y_discrete(limits = order_segments)
}



blank_genome = function(ref = "GRCh38",
                        chromosomes = paste0("chr", c(1:22, "X", "Y")),
                        label_chr = -0.5,
                        cex = 1)
{
  reference_coordinates = CNAqc:::get_reference(ref) %>%
    filter(chr %in% chromosomes)
  
  low = min(reference_coordinates$from)
  upp = max(reference_coordinates$to)
  
  pl = ggplot(reference_coordinates) +
    geom_rect(
      aes(
        xmin = centromerStart,
        xmax = centromerEnd,
        ymin = 0,
        ymax = Inf
      ),
      alpha = 0.3,
      colour = "gainsboro"
    ) +
    geom_segment(
      data = reference_coordinates,
      aes(
        x = from,
        xend = from,
        y = 0,
        yend = Inf
      ),
      size = 0.1,
      color = "black",
      linetype = 8
    ) +
    labs(x = "Chromosome",
         y = "Value") + 
    theme(
      axis.text.x = element_text(angle = 90)
    ) +
    scale_x_continuous(
      breaks = c(0, reference_coordinates$from,
                 upp),
      labels = c(
        "",
        gsub(
          pattern = "chr",
          replacement = "",
          reference_coordinates$chr
        ),
        ""
      )
    )
  return(pl)
}
