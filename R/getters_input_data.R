#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
get_model_description = function(x)
{
  if(all(is.null(x$description))) return("My CONGAS model")

  return(x$description)
}


#' Title
#'
#' @param inf_obj
#' @param segments_input
#' @param chromosomes
#' @param normalise
#'
#' @return
#' @export
#'
#' @examples
get_counts <-
  function(x,
           chromosomes = paste0("chr", c(1:22, "X", "Y")),
           normalise = TRUE,
           z_score = FALSE, sum_denominator = TRUE)
  {
    best_model <- Rcongas:::get_best_model(x)

    data_matrix = x$data$counts


    if (normalise){

      if(is.null( best_model$parameters$norm_factor))  best_model$parameters$norm_factor <-  rep(1, length( best_model$parameters$assignement))
      normalisation_factors = best_model$parameters$norm_factor

      assignments = get_cluster_assignments(x) %>% gsub(pattern = "C", replacement = "", ignore.case = T)  %>% as.numeric
      mu = get_input_segmentation(x)$mu
      total_cn =  rowSums(best_model$parameters$cnv_probs * mu)


      # Handle this special case which happens for already normalised data
      if(is.null(normalisation_factors)) {
        warning("normalisation_factors are NULL - replacing them with all 1s. Maybe you're using normalised data.")

        normalisation_factors = rep(1, ncol(data_matrix))
      }

      for (i in 1:ncol(data_matrix))
        if(sum_denominator)
          data_matrix[,i] = (data_matrix[,i] / normalisation_factors)  * (total_cn[assignments] / sum(mu))
        else
          data_matrix[,i] = (data_matrix[,i] / normalisation_factors)
    }

    if(z_score){
      rnames <-  rownames(data_matrix)
      data_matrix <-  apply(data_matrix, 2,scale, center = T, scale = T)
      rownames(data_matrix) <-  rnames
    }

      M <-
        long_counts(data_matrix) %>%  select(chr, from, to, cell, n) %>%
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


# get_counts_matrix = function(x)
# {
#   if (all(is.null(x$data$counts))) {
#     cli::cli_alert_warning("Input data has not been stored in the object, re-run the analysis with XXX = TRUE ...")
#     return(NULL)
#   }
#
#   return(x$data$counts)
# }



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
    dplyr::select(chr, from, to, size, dplyr::everything()) %>%
    deidify()


  return(df_segments %>% dplyr::filter(chr %in% chromosomes))
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
get_mapped_genes = function(x,
                            chromosomes = paste0("chr", c(1:22, "X", "Y"))
  )
{
  x$data$gene_locations %>%
    dplyr::select(-segment_id) %>%
    dplyr::filter(chr %in% chromosomes)
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
  # # TODO cambiare a questo appena riesco ad installare il package
  # ld = function(x)
  # {
  #   if('Rcongas' %in% (installed.packages() %>% rownames()))
  #     data(x, package = 'Rcongas')
  #   else
  #     load(paste0("~/Documents/GitHub/rcongas/data/", x, '.rda'))
  # }

  if (x$reference_genome %in% c('hg19', 'GRCh37'))
  {
    # ld('hg19_gene_coordinates')

    return(Rcongas::hg19_gene_coordinates)
  }

  if (x$reference_genome %in% c('hg38', 'GRCh38'))
  {
    # ld('hg38_gene_coordinates')
    return(Rcongas::hg38_gene_coordinates)
  }

  stop("reference unknown, you should use a reference supported by CONGAS!")
}

#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
get_dataset_stats = function(x)
{
  # number of cells, segments
  ncells = x$data$counts %>% nrow()
  nsegments = x$data$counts %>% ncol()

  # clustering stats
  nclusters = Rcongas::get_k(x)
  sizes_cl = Rcongas::get_clusters_size(x)
  psizes_cl = Rcongas::get_clusters_size(x, normalised = TRUE)

  # scores
  x$inference$model_selection


  return(
    list(
      ncells = ncells,
      nsegments = nsegments,
      clusters_k = nclusters,
      clusters_n = sizes_cl,
      clusters_pi = psizes_cl,
      score_type = x$inference$model_selection$IC_type,
      score = x$inference$model_selection$IC[x$inference$model_selection$best_K]
      )
  )
}


