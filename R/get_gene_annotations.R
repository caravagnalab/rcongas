#' Genome coordinates available in the package.
#' 
#' @description From a CONGAS object, it teturns the  
#' tibbles for the genomic coordinates available in the package
#' that have been used to assemble the object. 
#'
#' @param x An input CONGAS object.
#'
#' @return The required tibble.
#' @export
#'
#' @examples
#' 
#' x = Rcongas::congas_example
#'
#' get_gene_annotations(x)
get_gene_annotations <-  function(x)
{
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
  
  if (x$reference_genome %in% c('mm10', 'GRCm38'))
    return(Rcongas::mm10_gene_coordinates)
  
  stop("reference unknown, you should use a reference supported by CONGAS!")
}