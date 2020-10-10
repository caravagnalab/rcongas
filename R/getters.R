#' Title
#'
#' @param x
#' @param chromosomes
#'
#' @return
#' @export
#'
#' @examples
get_clones_ploidy <- function(x,
                              chromosomes = paste0("chr", c(1:22, "X", "Y")),
                              clusters = NULL)
{
  best_model <- get_best_model(x)
  res <-
    data.frame(t(best_model$parameters$cnv_probs), stringsAsFactors = FALSE)
  res <-
    res %>% dplyr::mutate(tmp = rownames(res)) %>% tidyr::separate(col = tmp,
                                                            into = c("chr", "from", "to"),
                                                            sep = ":")
  colnames(res)[seq_along(best_model$parameters$mixture_weights)] <-
    paste(seq_along(best_model$parameters$mixture_weights))
  res <-
    reshape2::melt(
      res,
      id.vars = c("chr", "from", "to"),
      variable.name = "cluster",
      value.name = "CN"
    ) %>%
    dplyr::mutate(
      cluster = as.character(cluster),
      from = as.numeric(from),
      to = as.numeric(to)
    ) %>%
    dplyr::as_tibble()

  # Normalise CNA values for comparisons
  means = res %>%
    dplyr::group_by(chr, from, to) %>%
    dplyr::summarise(segment_mean = mean(CN), .groups = 'keep') %>%
    dplyr::ungroup()

  joined_res = res %>%
    dplyr::full_join(means, by = c('chr', 'from', 'to')) %>%
    dplyr::mutate(lognorm_CN = CN,
                  CN = CN - segment_mean)

  if (!grepl('chr', joined_res$chr[1]))
  {
    cli::cli_alert_warning("Missing `chr` prefix in chromosomes labels, added now.")

    joined_res = joined_res %>% dplyr::mutate(chr = paste0("chr", chr))
  }

  # Apply filters
  joined_res = joined_res %>% dplyr::filter(chr %in% chromosomes)

  if (!is.null(clusters))
    joined_res = joined_res %>% dplyr::filter(cluster %in% clusters)


  return(joined_res)
}

#' Title
#'
#' @param x
#' @param normalised
#'
#' @return
#' @export
#'
#' @examples
get_clusters_size = function(x, normalised = FALSE)
{
  best_model = get_best_model(x)

  # Counts from the assignments
  n = table(best_model$parameters$assignement)

  # Make it a named vector
  v = as.vector(n)
  names(v) = names(n)

  if(normalised)
    v = v/sum(v)

  return(v)
}


#' Title
#'
#' @param x
#' @param cluster_label
#'
#' @return
#' @export
#'
#' @examples
get_clusters = function(x, clusters = NULL)
{
  # TODO - take all assingments with z_nk > c, c >= 0
  best_model = get_best_model(x)

  clusters_table = data.frame(
    cell = names(best_model$parameters$assignement),
    cluster = paste(best_model$parameters$assignement),
    stringsAsFactors = FALSE
  ) %>%
    as_tibble()

  if (!is.null(clusters))
    clusters_table = clusters_table %>% dplyr::filter(cluster %in% clusters)

  return(clusters_table)
}

#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
get_k = function(x)
{
  best = get_best_model(x)
  best$parameters$assignement %>%
    unique() %>%
    length()
}


#' Title
#'
#' @param inf_obj
#' @param input
#' @param chromosomes
#' @param norm
#'
#' @return
#' @export
#'
#' @examples
get_counts <-
  function(x,
           input,
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
      for (i in 1:nrow(input))
        input[i, ] = input[i, ] / best_model$parameters$norm_factor[i]

      M <-
        long_counts(input) %>%  select(chr, from, to, cell, n) %>%
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
#' @param chromosomes
#' @param cut_pvalue
#' @param cut_lfc
#'
#' @return
#' @export
#'
#' @examples
get_DE_table <- function(x,
                         chromosomes = paste0("chr", c(1:22, "X", "Y")),
                         cut_pvalue = 0.001,
                         cut_lfc = 0.25)
{
  if (!has_DE(x))
  {
    cli::cli_alert_danger("No DE table in this object, compute it first!")
    return(NULL)
  }

  table_data = x$DE$table %>%
    dplyr::filter(p_val_adj < cut_pvalue,
                  abs(avg_log2FC) > cut_lfc,
                  chr %in% chromosomes)

  return(table_data)
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
  if (all(is.null(x))) {
    cli::cli_alert_warning("Input data has not been stored in the object, re-run the analysis with XXX = TRUE ...")
    return(NULL)
  }

  return(x$data$counts)
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


# PRIVATE GETTERS

get_karyotype <-  function(x) {
  if (x$reference_genome %in% c('hg19', 'GRCh37'))
    data('hg19_karyo')
  if (x$reference_genome %in% c('hg38', 'GRCh38'))
    data('hg38_karyo')
}

get_best_model <- function(X) {
  X$inference$models[[X$inference$model_selection$best_K]]
}

has_DE <-  function(X) {
  return(!is.null(X$DE))

}


# Key creation and decrypt
idify = function(y) {
  y %>% dplyr::mutate(segment_id = paste(chr, as.integer(from), as.integer(to), sep = ":"))
}

deidify = function(y) {
  y %>% tidyr::separate(segment_id,
                        into = c('chr', 'from', 'to'),
                        sep = ":")
}
