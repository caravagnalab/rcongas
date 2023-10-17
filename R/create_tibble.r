#' Format the input counts before creating the CONGAS+ object
#' 
#' @description
#'   
#' This function takes in input a SummarizedExperiment object with rowranges and  and creates a tibble with counts needed to create the Rcongas object.
#' @param counts (Required). Count matrix with cells in columns and genes/bins in rows. 
#' @param features (Required) Dataframe containing columns chr, from, to and (only if RNA) gene. 
#' @param modality. Default 'ATAC' It has to be either ATAC or RNA based on the content of the count matrix.
#' @param output_dir. Optional. Default = NULL. Directory where to save the tibble as a \code{tsv} file. If NULL, the tibble object is only returned. 
#' @param output_file. Optional, default is 'counts_final.tsv'. Name to give to the counts file. Note that this has not effect is \code{output_dir} is NULL. 
#' @import Matrix
#' @import dplyr
#' @importFrom readr read_delim write_tsv
#' @export 
create_congas_tibble = function(counts, 
                                features = NULL,
                                modality = 'ATAC', 
                                save_dir = NULL,
                                output_file = 'counts_final.tsv') {

  cli::cli_alert_info("Saving temporary counts file\n")
  writeMM(as(counts, "sparseMatrix"), "counts_tmp.mtx")
  cells = data.frame(cell = colnames(counts)) %>% mutate(cell_index = seq(1, nrow(.)))
  # Now create tibble cell chr from to value  
  counts = readr::read_delim("counts_tmp.mtx",
                             comment = '%', delim = ' ') 
  unlink("counts_tmp.mtx")
  
  colnames(counts) =  c('feature_index', 'cell_index', 'value')
  
  # if (!is.null(features) && modality == 'RNA'){
  #   gene_names = rownames(counts)

  #   features = tibble(gene = gene_names) %>% left_join(features) 
  # } 
  if (!is.null(features) ){
    features = features %>% 
                mutate(feature_index = seq(1, nrow(.)))
  } 
  
  if (modality == 'ATAC' & is.null(features)) {
    stop('Please provide a dataframe with c("chr", "from", "to") through the parameter features')

  } else if (modality == 'RNA' & is.null(features)) { 
    stop('Please provide a dataframe with c("chr", "from", "to", "gene") through the parameter features')
  }
  
  
  counts = counts %>% dplyr::left_join(features) %>% dplyr::left_join(cells)
  
  if (modality == 'ATAC') {
    counts = counts %>% dplyr::select(c("cell", "value",  'chr', 'from', 'to')) %>%
      mutate(from = as.integer(from), to = as.integer(to), chr = as.character(chr))
  } else {
    counts = counts %>% dplyr::select(c("cell", "value",  'gene', 'chr', 'from', 'to')) %>%
      mutate(from = as.integer(from), to = as.integer(to))
  }
  
  if (!is.null(save_dir)) {
    readr::write_tsv(counts, file = paste0(save_dir, '/', output_file))
  } 
  
  # saveRDS(counts, paste0(save_dir, '/counts_final.rds'))
  return(counts)
  
}