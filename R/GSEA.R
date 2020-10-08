

calculate_GSEA <- function(inference,fc_df,clone1, clone2) {

  if(require("org.Hs.eg.db")){
    organism <- org.Hs.eg.db::org.Hs.eg.db
  } else {
    BiocManager::install("org.Hs.eg.db", character.only = TRUE)
    organism <- org.Hs.eg.db::org.Hs.eg.db
  }

  change <- log2(unlist(gtools::foldchange(fc_df[clone1,], fc_df[clone2,])))
  names(change) <-  colnames(fc_df)
  gene_list <- sort(change, decreasing = TRUE)
  gse <- gseGO(geneList=gene_list,
               ont ="ALL",
               keyType = "SYMBOL",
               minGSSize = 3,
               maxGSSize = 800,
               pvalueCutoff = 0.05,
               verbose = TRUE,
               OrgDb = organism,
               pAdjustMethod = "none", eps = 0)


  ids<-bitr(names(change), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)
  dedup_ids <-  ids[!duplicated(ids[c("SYMBOL")]),]
  change2 <-  change[names(change) %in% dedup_ids$SYMBOL]
  names(change2) <-  dedup_ids$ENTREZID
  kegg_gene_list <- change2
  kegg_gene_list<-na.omit(kegg_gene_list)
  kegg_gene_list <-  sort(kegg_gene_list, decreasing = TRUE)

  kegg_organism <-  "hsa"
  kk2 <- gseKEGG(geneList     = kegg_gene_list,
                 organism     = kegg_organism,
                 minGSSize    = 3,
                 maxGSSize    = 800,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "none",
                 keyType       = "ncbi-geneid", eps = 0)

  return(list(kegg = kk2, go = gse))

}
