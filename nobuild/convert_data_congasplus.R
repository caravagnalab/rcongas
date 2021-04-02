load("~/Downloads/congas_example.rda")
devtools::load_all("~/Documents/GitHub/rcongas/")

require(dplyr)

cells = congas_example$data$counts %>% rownames()
genes = congas_example$data$gene_locations %>% pull(gene)

df_counts = congas_example$data$gene_counts %>% 
  reshape2::melt() %>% 
  rename(cell = Var1, gene = Var2) %>% 
  as_tibble() %>% 
  filter(value > 0)

df_counts = df_counts %>% filter(cell %in% !!cells)
df_counts = df_counts %>% filter(gene %in% genes)
# !(df_counts$gene %in% genes)


load("~/Downloads/hg19_gene_coordinates.rda")
df_counts = df_counts %>% 
  left_join(hg19_gene_coordinates)

df_counts = df_counts[complete.cases(df_counts), ]

df_counts$cell = paste(df_counts$cell)
df_counts$value = as.integer(df_counts$value)

segments = congas_example$data$cnv %>% 
  select(chr, from, to, tot) %>% 
  rename(copies = tot)

rna_f = df_counts %>% auto_normalisation_factor()

x = init(
  rna = df_counts,
  atac = NULL,
  segmentation = segments,
  rna_likelihood = 'NB',
  reference_genome = 'hg19',
  description = "Campbell et al. Breast cancer PDX used in CONGAS",
  smooth = FALSE
)

x$input$segmentation %>% idify() %>% pull(segment_id)
plot_data(x, segments = c("chr15:67050001:102600000" , "chr16:1:3750000" , "chr18:32400001:55950000" ))
plot_data(x, segments = c("chr15:1:102600000" , "chr16:1:3750000" , "chr18:7950001:55950000" ))

campbell_bcpdx = x
campbell_bcpdx = filter_missing_data(campbell_bcpdx, proportion_RNA = 0)
campbell_bcpdx = filter_outliers(campbell_bcpdx)
campbell_bcpdx

plot_data(campbell_bcpdx, segments = c("chr15:67050001:102600000" , "chr16:1:3750000" , "chr18:32400001:55950000" ))
plot_data(x, segments = c("chr15:1:102600000" , "chr16:1:3750000" , "chr18:7950001:55950000" ))


setwd('~/Documents/GitHub/rcongas/')
usethis::use_data(campbell_bcpdx, overwrite = TRUE)

plot_data(campbell_bcpdx, segments = x$input$segmentation$segment_id[1:10])

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

load("~/Downloads/GMB/GBM_smartseq_normal_vs_tumor.rda")
devtools::load_all("~/Documents/GitHub/rcongas/")

require(dplyr)

cells = GBM_smartseq_normal_vs_tumor$data$counts %>% rownames()
genes = GBM_smartseq_normal_vs_tumor$data$gene_locations %>% pull(gene)

df_counts = GBM_smartseq_normal_vs_tumor$data$gene_counts %>% 
  reshape2::melt() %>% 
  rename(cell = Var2, gene = Var1) %>% 
  as_tibble() %>% 
  filter(value > 0)

df_counts = df_counts %>% filter(cell %in% !!cells)
df_counts = df_counts %>% filter(gene %in% genes)
# !(df_counts$gene %in% genes)

load("~/Downloads/hg19_gene_coordinates.rda")
df_counts = df_counts %>% 
  left_join(hg19_gene_coordinates)

df_counts = df_counts[complete.cases(df_counts), ]

df_counts$cell = paste(df_counts$cell)
df_counts$value = as.integer(df_counts$value)

segments = GBM_smartseq_normal_vs_tumor$data$cnv %>% 
  select(chr, from, to, ploidy_real) %>% 
  rename(copies = ploidy_real) %>% 
  mutate(copies = as.integer(2))

x = init(
  rna = df_counts,
  atac = NULL,
  segmentation = segments,
  rna_likelihood = 'G',
  reference_genome = 'GRCh38',
  description = "Patel et al. GBM tumour/normal used in CONGAS (SmartSeq)",
  smooth = FALSE
)

plot_data(patel_gbmtn)
patel_gbmtn = x

setwd('~/Documents/GitHub/rcongas/')
usethis::use_data(patel_gbmtn, overwrite = TRUE)

plot_data(campbell_bcpdx, segments = x$input$segmentation$segment_id[1:10])

############# 

load("~/Projects/2020. CONGAS+/monosomy_res.rda")
devtools::load_all("~/Documents/GitHub/rcongas/")

require(dplyr)

cells = monosomy_res$data$counts %>% rownames()
genes = monosomy_res$data$gene_locations %>% pull(gene)
locs = monosomy_res$data$gene_locations 

df_counts = monosomy_res$data$gene_counts %>% 
  reshape2::melt() %>% 
  rename(cell = Var2, gene = Var1) %>% 
  as_tibble() %>% 
  filter(value > 0)

df_counts = df_counts %>% filter(cell %in% !!cells)
df_counts = df_counts %>% filter(gene %in% genes)
# !(df_counts$gene %in% genes)


load("~/Projects/2020. CONGAS+/hg38_gene_coordinates.rda")
df_counts = df_counts %>% 
  left_join(hg38_gene_coordinates)

df_counts = df_counts[complete.cases(df_counts), ]

df_counts$cell = paste(df_counts$cell)
df_counts$value = as.integer(df_counts$value)

segments = monosomy_res$data$cnv %>% 
  dplyr::select(chr, from, to, ploidy_real) %>% 
  dplyr::rename(copies = ploidy_real) %>% 
  mutate(copies = as.integer(2))

x = init(
  rna = df_counts,
  atac = NULL,
  segmentation = segments,
  rna_likelihood = 'G',
  reference_genome = 'GRCh38',
  description = "Zaho et al. MDS hematopoietical used in CONGAS (HCA)",
  smooth = FALSE
)

plot_data(zaho_hemato)
zaho_hemato = x

zaho_hemato = filter_outliers(zaho_hemato)

setwd('~/Documents/GitHub/rcongas/')
usethis::use_data(zaho_hemato, overwrite = TRUE)

plot_data(campbell_bcpdx, segments = x$input$segmentation$segment_id[1:10])





library(biomaRt)
mart <- useMart("ensembl")
mart <- useDataset("hsapiens_gene_ensembl", mart)
attributes <- c("ensembl_gene_id","start_position","end_position","hgnc_symbol","chromosome_name")
filters <- c("chromosome_name","start","end")

coo = CNAqc::chr_coordinates_GRCh38
all_genes = NULL
for(i in 1:nrow(coo)){
  
  values <- list(chromosome=gsub("chr", "", coo$chr[i]),start="0",end=paste(coo$to[i]))
  all.genes <- getBM(attributes=attributes, filters=filters, values=values, mart=mart) %>% filter(hgnc_symbol!='') %>% 
    group_by(hgnc_symbol) %>% 
    filter(row_number() == 1) %>% 
    ungroup()
  
  colnames(all.genes)[c(2,3,4,5)] = c('from', 'to', 'gene', 'chr')

  all.genes = all.genes %>% 
    mutate(chr = paste0('chr', chr)) %>% 
    dplyr::select(chr, from, to, gene, ensembl_gene_id)
  
  all_genes = bind_rows(all_genes, all.genes)

}

all_genes

