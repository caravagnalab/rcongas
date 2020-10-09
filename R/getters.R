#' Title
#'
#' @param x
#' @param chromosomes
#'
#' @return
#' @export
#'
#' @examples
get_clones_ploidy <-
  function(x, chromosomes = paste0("chr", c(1:22, "X", "Y"))) {
    best_model <- get_best_model(x)
    res <-
      data.frame(t(best_model$parameters$cnv_probs), stringsAsFactors = FALSE)
    res <-
      res %>% mutate(tmp = rownames(res)) %>% tidyr::separate(col = tmp,
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
      dplyr::summarise(segment_mean = mean(CN))

    joined_res = res %>%
      dplyr::full_join(means, by = c('chr', 'from', 'to')) %>%
      dplyr::mutate(lognorm_CN = CN,
                    CN = CN - segment_mean)

    if (!grepl('chr', joined_res$chr[1]))
    {
      cli::cli_alert_warning("Missing `chr` prefix in chromosomes labels, added now.")

      joined_res = joined_res %>% dplyr::mutate(chr = paste0("chr", chr))
    }

    return(joined_res %>% dplyr::filter(chr %in% chromosomes))
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
  function(inf_obj,
           input,
           chromosomes = paste0("chr", c(1:22, "X", "Y")),
           norm = T)
  {
    best_model <- get_best_model(inf_obj)
    if (norm)
      input$counts <-
        input$counts / best_model$parameters$norm_factor

    M <-
      long_counts(input$counts) %>%  select(chr, from, to, cell, n) %>%
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
    dplyr::filter(
      p_val < cut_pvalue,
      abs(avg_log2FC) > cut_lfc,
      chr %in% chromosomes
      )

  return(table_data)
}



#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
get_input_segmentation = function(x)
{
  # Break the internal naming system
  df_segments = lapply(strsplit(x$models[[1]]$dim_names$seg_names, split = ':'),
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

  return(df_segments)
}


# PRIVATE GETTERS

get_best_model <- function(inf_obj) {
  inf_obj$models[[inf_obj$best_K]]
}



get_gene_annotations = function(x)
{
  if (x$reference_genome %in% c('hg19', 'GRCh37'))
  {
    data('hg19_gene_coordinates')
    return(hg19_gene_coordinates)
  }

  # if(x$reference_genome %in% c('hg38', 'GRCh38')) data('hg38_gene_coordinates')

  stop("reference unknown")
}
