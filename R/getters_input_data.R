

#' Title
#'
#' @param inf_obj
#' @param segments_input
#' @param chromosomes
#' @param norm
#'
#' @return
#' @export
#'
#' @examples
get_counts <-
  function(x,
           segments_input,
           chromosomes = paste0("chr", c(1:22, "X", "Y")),
           norm = T)
  {
    best_model <- get_best_model(x)

    # Era sbagliato entrare dentro counts visto che input E' la matrice dei counts.
    # poi la divisione da nonconformable arrays, ho fatto un for expplicito
    # if (norm)
    # input$counts <-
    #   input$counts / best_model$parameters$norm_factor

    if (norm)
      for (i in 1:nrow(segments_input))
        segments_input[i, ] = segments_input[i, ] / best_model$parameters$norm_factor[i]

      M <-
        long_counts(segments_input) %>%  select(chr, from, to, cell, n) %>%
        dplyr::arrange(chr, from, to)

      clts <- as.data.frame(best_model$parameters$assignement)
      colnames(clts) <- "cluster"
      clts$cell <- rownames(clts)

      ret <-  dplyr::inner_join(M, clts, by = "cell") %>%
        dplyr::select(chr, from, to, cell, n, cluster) %>%
        dplyr::mutate(cluster = paste(cluster))

      if (!grepl('chr', ret$chr[1]))
      {
        cli::cli_alert_warning("Missing `chr` prefix in chromosomes labels, added now.")

        ret <-  ret %>% dplyr::mutate(chr = paste0("chr", chr))
      }


      return(ret %>% filter(chr %in% chromosomes))

  }


#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
get_counts_matrix = function(x)
{
  if (all(is.null(x$data$counts))) {
    cli::cli_alert_warning("Input data has not been stored in the object, re-run the analysis with XXX = TRUE ...")
    return(NULL)
  }

  return(x$data$counts)
}



#' Title
#'
#' @param x
#' @param chromosomes
#'
#' @return
#' @export
#'
#' @examples
get_input_segmentation = function(x,
                                  chromosomes = paste0("chr", c(1:22, "X", "Y")))
{
  # Break the internal naming system
  df_segments = lapply(strsplit(x$inference$models[[1]]$dim_names$seg_names, split = ':'),
                       function(x) {
                         data.frame(
                           chr = x[1],
                           from = x[2],
                           to = x[3],
                           stringsAsFactors = FALSE
                         )
                       })

  df_segments = Reduce(dplyr::bind_rows, df_segments) %>% dplyr::as_tibble()
  df_segments$from = as.numeric(df_segments$from)
  df_segments$to = as.numeric(df_segments$to)

  # Fix missing chr label
  if (!grepl('chr', df_segments$chr[1]))
  {
    cli::cli_alert_warning("Missing `chr` prefix in chromosomes labels, added now.")

    df_segments = df_segments %>% dplyr::mutate(chr = paste0("chr", chr))
  }

  # Add other information(s)
  df_segments = df_segments %>%
    dplyr::left_join(x$data$cnv, by = c('chr', 'from', 'to')) %>%
    dplyr::mutate(size = to - from) %>%
    dplyr::select(chr, from, to, size, dplyr::everything())


  return(df_segments %>% dplyr::filter(chr %in% chromosomes))
}

#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
get_input_raw_data = function(x)
{
  if (all(is.null(x$data$gene_counts))) {
    cli::cli_alert_warning("Input data has not been stored in the object, re-run the analysis with XXX = TRUE ...")
    return(NULL)
  }

  return(x$data$gene_counts)
}

#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
get_reference_genome = function(x)
{
  return(x$reference_genome)
}

#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
get_gene_annotations <-  function(x)
{
  # TODO cambiare a questo appena riesco ad installare il package
  ld = function(x)
  {
    if('Rcongas' %in% (installed.packages() %>% rownames()))
      data(x, package = 'Rcongas')
    else
      load(paste0("~/Documents/GitHub/rcongas/data/", x, '.rda'))
  }

  if (x$reference_genome %in% c('hg19', 'GRCh37'))
  {
    ld('hg19_gene_coordinates')
    return(hg19_gene_coordinates)
  }

  if (x$reference_genome %in% c('hg38', 'GRCh38'))
  {
    ld('hg38_gene_coordinates')
    return(hg38_gene_coordinates)
  }

  stop("reference unknown, you should use a reference supported by CONGAS!")
}

