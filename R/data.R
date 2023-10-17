#' Example tibbles with input ATAC data in the format for (R)CONGAS+
#'
#' @description Example tibbles with input data in the format for (R)CONGAS+
#'
#' @docType data
#'
#' @usage data(atac_counts)
#'
#' @format Example tibbles with input ATAC data.
#'
#' @keywords datasets
#'
#' @examples
#' data(atac_counts)
#' print(atac_counts)
"atac_counts"

#' Example tibbles with input RNA data in the format for (R)CONGAS+
#'
#' @description Example tibbles with input data in the format for (R)CONGAS+
#'
#' @docType data
#'
#' @usage data(rna_counts)
#'
#' @format Example tibble with input RNA data.
#'
#' @keywords datasets
#'
#' @examples
#' data(rna_counts)
#' print(rna_counts)
"rna_counts"

#' Dataframe containing cell type annotations for the 10x multiome lymphoma dataset.
#'
#' @description Dataframe with cell type annotations for the 10x multiome lymphoma dataset.
#'
#' @docType data
#'
#' @usage data(metadata)
#'
#' @format Dataframe.
#'
#' @keywords datasets
#'
#' @examples
#' data(metadata)
#' head(metadata)
"metadata"


#' Breast cancer patient-derived xenograft (10x RNA, counts-based).
#' 
#' @description (R)CONGAS+ object for the breast cancer patient-derived xenograft,
#' first used in the original CONGAS paper. The original paper describing
#' this data is Campbell, K.R., Steif, A., Laks, E. et al. clonealign: 
#' statistical integration of independent single-cell RNA and DNA sequencing 
#' data from human cancers. Genome Biol 20, 54 (2019). 
#' https://doi.org/10.1186/s13059-019-1645-z
#'
#' @docType data
#'
#' @usage data(campbell_bcpdx)
#'
#' @format (R)CONGAS+ object.
#'
#' @keywords internal
#'
#' @examples
#' data(campbell_bcpdx)
#' print(campbell_bcpdx)
"campbell_bcpdx"

#' Glioblastoma tumour and normal (SmartSeq RNA, z-score)
#'
#'  @description (R)CONGAS+ object for the glioblastoma dataset,
#' first used in the original CONGAS paper. The paper providing these data
#' is Patel, Anoop P., et al. "Single-cell RNA-seq highlights intratumoral 
#' heterogeneity in primary glioblastoma." Science 344.6190 (2014): 1396-1401.
#'
#'
#' @docType data
#'
#' @usage data(patel_gbmtn)
#'
#' @format (R)CONGAS+ object.
#'
#' @keywords internal
#'
#' @examples
#' data(patel_gbmtn)
#' print(patel_gbmtn)
"patel_gbmtn"

#' Hematopoietic bone marrow failure (10x RNA, counts-based)
#'
#' @description (R)CONGAS+ object for the hematopoietic dataset,
#' first used in the original CONGAS paper. The paper providing these data
#' is Zhao, Xin, et al. "Single-cell RNA-seq reveals a distinct transcriptome
#' signature of aneuploid hematopoietic cells." Blood 130.25 (2017): 2762-2773.
#'
#' @docType data
#'
#' @usage data(zaho_hemato)
#'
#' @format (R)CONGAS+ object.
#'
#' @keywords internal
#'
#' @examples
#' data(zaho_hemato)
#' print(zaho_hemato)
"zaho_hemato"

#' Karyotype hg38
#'
#' @description Data about the coordinates of hg38
#'
#' @docType data
#'
#' @usage data(hg38_karyo)
#'
#' @format dataframe
#'
#' @keywords datasets
#'
#' @examples
#' data(hg38_karyo)
#' print(hg38_karyo)
"hg38_karyo"

#' Multiome congas object
#'
#' @description Example of a CONGAS+ object containg a subset of the multiome Lymphoma data.
#'
#' @docType data
#'
#' @usage data(multiome_congas_object)
#'
#' @format CONGAS+ object
#'
#' @keywords datasets
#'
#' @examples
#' data(multiome_congas_object)
#' print(multiome_congas_object)
"multiome_congas_object"
