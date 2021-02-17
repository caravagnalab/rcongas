#' Title
#'
#' @param data 
#' @param cnv_data 
#' @param chromosomes 
#' @param fun 
#' @param online 
#' @param reference_genome 
#' @param correct_bins 
#' @param gene_format 
#' @param description 
#'
#' @return
#' @export
#'
#' @examples
init <-
  function(data,
           cnv_data,
           chromosomes = paste0("chr", (1:22)),
           fun = sum,
           online = FALSE,
           reference_genome = "hg38",
           correct_bins = TRUE,
           gene_format = "hgnc_symbol",
           description = "My CONGAS dataset"
           )
    
  {
    cli::cli_alert_info("Extracting gene data from reference.")
    
    sanitise_input(data, cnv_data)
    
    # Salvatore's mapping
    df_list <-
      getChromosomeDF(
        data,
        genome = reference_genome,
        online = online,
        chrs = chromosomes,
        filters = gene_format
      )
    
    data_splitted <- df_list[[1]]
    chrs_cord <- df_list[[2]]
    
    
    # Process data
    cli::cli_alert_info("Processing input counts.")
    
    cp <-
      cnv_data %>% mutate(dist = to - from,
                          to = as.integer(to),
                          from = as.integer(from)) %>% filter(chr %in% chromosomes)
    
    if (nrow(cp) == 0)
      stop("Nothing mapped? Check your inputs..")
    
    if (correct_bins)
      cp <- correct_bins(cp, filt = 1.0e+07, hard = TRUE)
    
    # Data split
    bins_splitted <- split(cp , cp$chr)
    bins_splitted <-
      bins_splitted[gtools::mixedsort(names(bins_splitted))]
    
    chrs_cord <-
      chrs_cord[names(chrs_cord) %in% names(bins_splitted)]
    
    data_splitted <-
      data_splitted[names(data_splitted) %in% names(bins_splitted)]
    
    mask <- mapply(
      bins_splitted,
      chrs_cord,
      FUN = function(x, y) {
        res <- vector(length = nrow(x))
        for (i in seq_len(nrow(x))) {
          res[i] <-
            any(x$from[i] <  y$start_position &
                  x$to[i] > y$end_position)
        }
        return(res)
      }
    )
    
    data_binned <-
      mapply(
        data_splitted,
        chrs_cord,
        bins_splitted ,
        FUN = function(x, y, z)
          custom_fixed_binned_apply(x, y, z, fun)
      )
    
    # Assembly final object
    cli::cli_alert_info("Assembly Rcongas object.")
    
    result <- data_binned[seq(1, length(data_binned), 4)]
    result_bindim <- data_binned[seq(2, length(data_binned), 4)]
    fixed_bindim <- data_binned[seq(3, length(data_binned), 4)]
    genes_in_bins <- data_binned[seq(4, length(data_binned), 4)]
    
    result <- t(do.call(rbind, result))
    result_bindim <- t(do.call(rbind, result_bindim))
    fixed_bindim <- do.call(c, fixed_bindim)
    colnames(result_bindim) <- colnames(result)
    cp <- do.call(rbind, bins_splitted)
    
    cp$mu <- apply(result_bindim, 2, median)
    cp$fixed_mu <- fixed_bindim
    colnames(cp)[4] <- "ploidy_real"
    
    gene_locations <-
      as.matrix(do.call(rbind, genes_in_bins)) %>%
      dplyr::as_tibble()
    
    locations = Reduce(bind_rows, chrs_cord) %>% distinct()
    colnames(locations) = c('gene', 'chr', 'from', 'to')
    
    gene_locations = gene_locations %>%
      left_join(locations, by = 'gene') %>%
      select(gene, chr, from, to, segment_id)
    
    # Retain only what is mappable for sure
    nr = gene_locations %>% nrow
    gene_locations = gene_locations[complete.cases(gene_locations), ]
    nr2 = gene_locations %>% nrow

    if(nr - nr2 > 0) 
    {
      cli::cli_alert_warning("Mapping inconsistent for {.field {nr-nr2}} genes out {.field {nr}}, removing those from the raw data table.")
    }
    
    # Tibble data
    data_tb = data %>% as_tibble()
    data_tb$gene = rownames(data)
    
    data_tb = data_tb %>%
      reshape2::melt(id = 'gene', value.name = 'n') %>%
      mutate(variable = paste(variable)) %>% 
      as_tibble() %>%
      rename(cell = variable) %>%
      filter(n > 0)
    
    ds = format(object.size(data_tb), units = 'Mb') %>% crayon::blue()
    dsm = format(object.size(data), units = 'Mb') %>% crayon::red()
    
    cli::cli_alert_info(
      "Retaining {.value {ds}} long-format tibble data with {.field {data_tb %>% nrow}} points, matrix was {.value {dsm}}."
    )
    
    data_list <-
      list(
        counts = as.matrix(result),
        bindims = result_bindim,
        cnv = cp %>%  dplyr::as_tibble() %>% idify(),
        gene_locations = gene_locations,
        raw = data_tb
      )
    
    ret <-
      structure(
        list(
          data = data_list,
          reference_genome = reference_genome,
          description = description
        ),
        class = "rcongas"
      )
    
    ds = format(object.size(ret), units = 'Mb') %>% crayon::blue()
    
    cli::cli_alert_info("Object size in memory: {.value {ds}}.")
    cat("\n")
    print(ret)
    
    return(ret)
    
  }

# Sanitise input requrements
sanitise_input = function(data, cnv_data)
{
  # Duplicate gene or cell ids
  cell_names = colnames(data) %>% table
  gene_names = rownames(data) %>% table
  
  if(cell_names[cell_names>1] %>% length > 0)
    stop("Cell ids cannot be duplicated: ", cell_names[cell_names>1] %>% names %>% paste(collapse = ', '))

  if(gene_names[gene_names>1] %>% length > 0)
    stop("Gene ids cannot be duplicated: ", gene_names[gene_names>1] %>% names %>% paste(collapse = ', '))
  
  # CNV data format check
  if (!all(c("from", "to", "tot") %in% colnames(cnv_data))) {
    stop("Chromsomome should have start/end columns.")
  }
  
  if (!all(grepl('chr', cnv_data$chr))) {
    stop("Chromsomome format has to be chr1, chr2, ....")
  }
  
  cli::cli_alert_success(
    "Validated input(s): {.field {cell_names %>% length}} cells, {.field {gene_names %>% length}} genes and {.field {cnv_data %>% nrow}} CNA segments."
  )
}
