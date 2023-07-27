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
                              to_plot = NULL, 
                              segments = get_input(x, what = 'segmentation') %>% pull(segment_id),
                              whichfacet = ggplot2::facet_wrap, bins = 60, position = 'stack',
                              highlights = FALSE,
                              single_segment_mode = FALSE, 
                              palette_name = NULL, 
                              colors = NULL)
{
  if (is.null(to_plot)) {
    to_plot = 'modality'
  } else {
    if (sum(! to_plot %in% colnames(x$input$metadata) & to_plot != 'clusters') > 0 ){
      stop('Error, to_plot must be either the name of a metadata column or "clusters". Exiting')
    }
  }
  stats_data = Rcongas:::stat(x, what = 'data')
  
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
    what_RNA = Rcongas:::normalise_modality(what_RNA, norm_factors %>% filter(modality == "RNA"))
  
  # ATAC_data values
  what_ATAC = get_input(x, what = 'data') %>%
    filter(modality == "ATAC")
  
  if (nrow(what_ATAC) > 0 && stats_data$atac_dtype != "G")
    what_ATAC = Rcongas:::normalise_modality(what_ATAC, norm_factors %>% filter(modality == "ATAC"))
  
  # Rescale against average normalization factors
  what_ATAC = what_ATAC %>%  mutate(value = value * mean(norm_factors %>% filter(modality == "ATAC") %>%  pull(normalisation_factor)))
  
  what_RNA = what_RNA %>%  mutate(value = value * mean(norm_factors %>% filter(modality == "RNA") %>%  pull(normalisation_factor)))
  
  what = bind_rows(what_RNA, what_ATAC) %>%
    filter(segment_id %in% segments) %>%
    mutate(modality = case_when(
      modality == "RNA" ~ paste(modality, '(', what_rna_lik, ')'),
      modality == "ATAC" ~ paste(modality, '(', what_atac_lik, ')')
    ))
  
  what$modality <- sapply(what$modality %>% strsplit(" "), 
                          function(y) y[1])
  
  # Set some properties of the plot
  if (position != 'stack') alpha = 0.8 else alpha = 1
  
  # Get segments to plot
  if (highlights) {
    CNAs = get_fit(x, "CNA")
    nclusters = CNAs$cluster %>% unique %>% length()
    CNAs = CNAs %>% group_by(segment_id, value) %>% mutate(grp_size = n()) %>% 
      filter(grp_size != nclusters) %>% pull(segment_id) %>% 
      unique
    cli::cli_alert("Plotting segments where different CNAs are present: {.field {CNAs}}.")
  } else {
    CNAs = x$input$segmentation %>% pull(segment_id) %>% unique
    cli::cli_alert("Showing all segments (this plot can be large).")
  }
  
  # Possible behaviours of this function: 
  # 1. to_plot == 'clusters'. In this case you need to get the cluser assignments from the object and then you return a list with one plot per element, colored according to cluster assignments.
  if (to_plot == 'clusters') {
    what = dplyr::left_join(what, get_fit(x, what = "cluster_assignments") %>% dplyr::select(-modality))
    if (length(CNAs) == 0) {
      return(ggplot() + geom_blank())
    }
    
    CNA <- Rcongas:::get_CNA(x)
    colnames(CNA)[3] <- "CNA"
    what <- dplyr::left_join(what, CNA)
    what$clusterCNA = paste0(what$cluster, " (CN = ", what$CNA, ")")
    
    if (!is.null(colors)) {
      cols = colors
    } else {
      palette_name = if (!is.null(palette_name)) palette_name else 'Set1'
      nclusters = CNA$cluster %>% unique %>% length()
        if (nclusters > 9) {
          getPalette = grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, palette_name))
          cols = getPalette(nclusters)
        } else {
          cols <- RColorBrewer::brewer.pal(nclusters, palette_name)
        }
    }

    CNAs = gtools::mixedsort(CNAs)
    ret <- sapply(CNAs,  simplify = F, USE.NAMES = T, function(s) {
      ggplot() + 
        geom_histogram(aes_string("value", fill = 'clusterCNA'), 
                       bins = bins, 
                       data = what %>% filter(segment_id == s), 
                       position=position,
                       alpha = alpha) + 
        #scale_color_manual("Cluster", values = cols, drop = FALSE) +
        facet_wrap(segment_id ~ modality, scales = "free") + 
        guides() + 
        theme_light(base_size = 9) + 
        theme(strip.text = element_text(colour = 'black')) +
        scale_fill_manual("Cluster ", values = cols) + 
        labs(x = "Input", y = "Observations") + 
        theme(strip.text.y.right = element_text(angle = 0), 
              # axis.text.y = element_blank(), 
              # axis.ticks.y = element_blank(), 
              # axis.ticks = element_blank(), 
              legend.position="top",
              legend.direction='vertical') 
    })
    return(ret)
  } 
  
  # 2. to_plot == any metadata column. IN this case you need to get the values to plot from the metadata column. 
  # 3. to_plot == 'modality'. In this case, we don't need to extract anything from the metadata, the variable what 
  # already contains the data we need for plotting.
  
  # Here we are putting inside the variable what the data we need for plotting (so we are finishing preparing 'what')
  # In case to_plot is not null we take its value from the metadata dataframe
  # In case to_plot is 'modality', then there is no need for the object to
  # have a metadata field.
  if (to_plot != 'modality') {
    what = dplyr::left_join(what, x$input$metadata %>% dplyr::select(cell, to_plot))
    what[[to_plot]] = factor(what[[to_plot]], levels = gtools::mixedsort(what[[to_plot]] %>% unique))
  } 
  
  what = what %>%
    mutate(segment_id = factor(segment_id, levels = gtools::mixedsort(segment_id %>% unique))) 
  
  # Two possible behaviours: if single_segment_mode == TRUE return a list of plots, otherwise return only one plot faceted.
  
  plot_list = list()
  if (single_segment_mode) {
    plot_list = sapply(CNAs, simplify = F, USE.NAMES = T, function(s) plot_data_histogram_aux(what, s, to_plot, alpha, position, palette_name, colors))

  } else {
    plot_list[[1]] = plot_data_histogram_aux(what, s = NULL, to_plot = to_plot, alpha = alpha, position = position, palette_name, colors) + 
      labs(title = x$description, subtitle = subtitle)
  }
  
  if (length(plot_list) == 1) {
    return(plot_list[[1]])
  } else 
    return(plot_list)
}

plot_data_histogram_aux = function(what, s, to_plot, alpha, position, palette_name, colors = NULL) {
  
  if (!is.null(s)) {
    what = what %>% filter(segment_id == s)
  }
  
  p = ggplot(what) + geom_histogram(aes_string("value", fill = to_plot), alpha=alpha, bins = 50, position=position) +
    facet_wrap(segment_id ~ modality, scales = 'free') +
    theme_light(base_size = 9) + 
    labs(x = "Input",
         y = 'Observations') +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 4, min.n = 2)) +
		scale_y_continuous(breaks = scales::pretty_breaks(n = 4, min.n = 2)) +
    theme(strip.text.y.right = element_text(angle = 0), 
          legend.position="top", legend.direction='horizontal',
          legend.title = element_blank()) +
    theme(strip.text = element_text(colour = 'black')) 

  if (!is.null(palette_name)) {
    cols <- RColorBrewer::brewer.pal(length(unique(what[,to_plot])), palette_name)
    p = p + scale_fill_manual(values = cols) 
  } else if(!is.null(colors)) {
    p = p + scale_fill_manual(values = colors) 
  }
  
  if (to_plot == 'modality') {
    cols = Rcongas:::modality_colors(what$modality %>% unique)
    p = p + scale_fill_manual(values = cols) 
  } 
  return(p)
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

plot_data_heatmap = function(x, segments = get_input(x, what = 'segmentation') %>% pull(segment_id), scale = FALSE, scale_min_lim = -Inf, scale_max_lim = Inf)
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
    if(scale) what_RNA <- what_RNA %>% group_by(segment_id) %>% mutate(value = (value - mean(value)) / sd(value)) %>% ungroup() %>%
        mutate(value = ifelse(value > scale_min_lim, value, scale_min_lim)) %>% mutate(value = ifelse(value > scale_max_lim, scale_max_lim, value)) 
    
    RNA_plot = ggplot(what_RNA %>% filter(segment_id %in% segments)) +
      geom_tile(aes(x = segment_id, y = cell, fill = value)) +
      theme_linedraw(base_size = 9) +
      labs(x = "Segments",
           y = subtitle_RNA) +
      theme(axis.text.y = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),) +
      guides(fill = guide_colorbar(paste('RNA (', what_rna_lik, ')'),
                                   barheight = unit(3, 'cm')))

    if (stats_data$rna_dtype == "G" | scale)
      RNA_plot = RNA_plot +
        scale_fill_gradient2(low = "steelblue", high = 'indianred3')
    else
      RNA_plot = RNA_plot +
        scale_fill_distiller(palette = 'GnBu', direction = 1)
  }

  if (nrow(what_ATAC) > 0)
    # RNA
  {
    if(scale) what_ATAC <- what_ATAC %>% group_by(segment_id) %>% mutate(value = (value - mean(value)) / sd(value)) %>% ungroup() %>%
        mutate(value = ifelse(value > scale_min_lim, value, scale_min_lim)) %>% mutate(value = ifelse(value > scale_max_lim, scale_max_lim, value)) 
    
    ATAC_plot = ggplot(what_ATAC %>% filter(segment_id %in% segments)) +
      geom_tile(aes(x = segment_id, y = cell, fill = value)) +
      theme_linedraw(base_size = 9) +
      labs(x = "Segments",
           y = subtitle_ATAC) +
      theme(axis.text.y = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),) +
      guides(fill = guide_colorbar(paste('ATAC (', what_atac_lik, ')'),
                                   barheight = unit(3, 'cm')))

    if (stats_data$atac_dtype == "G" | scale)
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
