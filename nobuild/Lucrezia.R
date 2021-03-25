devtools::load_all()

# ATAC
atac = readr::read_csv('~/Downloads/tmp_scatac.csv')
atac = atac %>% filter(value > 0)
colnames(atac) = c('cell', 'segment_id', 'value')

norm_atac = atac %>% 
  group_by(cell) %>% 
  summarise(normalisation_factor = sum(value)) %>% 
  mutate(modality = 'ATAC')

atac$value = atac$value %>% as.integer()
atac = atac %>% deidify()

# RNA
rna = readr::read_csv('~/Downloads/tmp_scrnac.csv')
rna = rna %>% filter(value > 0)
colnames(rna) = c( 'gene', 'cell', 'value')

norm_rna = rna %>% 
  group_by(cell) %>% 
  summarise(normalisation_factor = sum(value)) %>% 
  mutate(modality = 'RNA')

load("~/Downloads/hg19_gene_coordinates.rda")

rna = rna %>% 
  left_join(hg19_gene_coordinates, by = 'gene')
rna = rna[complete.cases(rna), ]

rna$value = rna$value %>% as.integer()

# Segmentation
segmentation = readr::read_csv("/Volumes/Data/Dropbox/Research Projects/2021. Organoidi HSR/WEX/CRC17/CRC17_smoothed_segments.csv")

segmentation = segmentation %>% rename(copies = CNt)

segmentation$from = segmentation$from %>% as.integer()
segmentation$to = segmentation$to %>% as.integer()
segmentation$copies = segmentation$copies %>% as.integer()

x = init(
  rna = rna %>% filter_known_genes(),
  atac = atac,
  segmentation = segmentation,
  normalisation_factors = rbind(norm_rna, norm_atac),
  rna_likelihood = "G", 
  atac_likelihood = 'NB',
  description = 'CRC17 HSR Organoid')

# Given mapping, keep only segments with joint signal
segmentation = get_input(x, 'segmentation') %>% 
  filter(RNA_genes > 0, ATAC_peaks > 0)

x = init(
  rna = rna %>% filter_known_genes(),
  atac = atac,
  segmentation = segmentation,
  normalisation_factors = rbind(norm_rna, norm_atac),
  rna_likelihood = "G", 
  atac_likelihood = 'NB',
  description = 'CRC17 HSR Organoid - segments with both RNA/ATAC')

example_object = x
usethis::use_data(example_object)

df = data.frame(
  segment_id = segments$segment_id,
  parameter = 'map',
  value = sample(1:5, nrow(segments), replace = T),
  cluster = 'C1'
)

df2 = df
df2$cluster='C2'
df2$value[c(3, 4, 12)] = 2


x$best_fit$CNA = rbind(df,df2)

cells = x$input$dataset$cell %>% unique
rna_cells =  x$input$dataset %>% filter(modality == 'RNA') %>% pull(cell)
a_cells =  x$input$dataset %>% filter(modality != 'RNA') %>% pull(cell)

df = data.frame(
  cell = cells,
  cluster =  sample(c("C1", "C2"), length(cells), replace = T)
) %>% 
  mutate(  
    modality = ifelse(cell %in% rna_cells, "RNA", "ATAC")
)
x$best_fit$cluster_assignments = df


lucrez = readr::read_tsv("~/Downloads/counts_final.tsv")
atac = lucrez
colnames(atac) = c("cell", "value",  'chr', 'from', 'to') 
atac$value = atac$value %>% as.integer()
atac$from = atac$from %>% as.integer()
atac$to = atac$to %>% as.integer()

rna = readr::read_tsv("~/Downloads/RNA_counts_final.tsv")
colnames(rna) = c("gene", "chr", 'from', 'to', 'cell', 'value') 
rna$value = rna$value %>% as.integer()
rna$from = rna$from %>% as.integer()
rna$to = rna$to %>% as.integer()

# Apply some basic filters - cap outliers and remove known genes
atac = cap_values_by_quantile(atac)
rna = cap_values_by_quantile(rna)
rna = rna %>% filter_known_genes()

# Subset cell types
table_celltypes = readr::read_tsv("~/Downloads/atac_rna_celltypes.tsv", col_names = FALSE)
colnames(table_celltypes) = c("cell", 'type')

table_celltypes$type %>% table

tum_cells = table_celltypes %>% filter(type == "Tumor")
imm_cells = table_celltypes %>% filter(type == "Immune") %>% sample_n(1256)
tcl_cells = table_celltypes %>% filter(type == "T-cells") %>% sample_n(1256)

subset_cells = c(tum_cells$cell, imm_cells$cell, tcl_cells$cell)

atac = atac %>% filter(cell %in% subset_cells)
rna = rna %>% filter(cell %in% subset_cells)

# Empirical normalisation factors
norm_rna = rna %>% 
  group_by(cell) %>% 
  summarise(normalisation_factor = sum(value)) %>% 
  mutate(modality = 'RNA')

norm_atac = atac %>% 
  group_by(cell) %>% 
  summarise(normalisation_factor = sum(value)) %>% 
  mutate(modality = 'ATAC')

# Segment files
segments = readr::read_tsv('~/Downloads/segments (1).tsv')
colnames(segments)[3:5] = c('chr', 'from', 'to')
segments$copies = segments$Corrected_Copy_Number %>% as.integer()
segments$from = segments$from %>% as.integer()
segments$to = segments$to %>% as.integer()
segments = segments %>% select(chr, from, to, copies) %>% filter(copies>0)

# Check types
rna$value = rna$value %>% as.integer()
atac$value = atac$value %>% as.integer()

x = init(
  rna = rna,
  atac = atac,
  segmentation = segments,
  normalisation_factors = rbind(norm_rna, norm_atac),
  rna_likelihood = "NB", 
  atac_likelihood = 'NB',
  description = 'Lucrezia RNA/ATAC')

segmentation = get_input(x, 'segmentation') %>% 
  filter(RNA_genes > 50, ATAC_peaks > 250)

x = init(
  rna = rna,
  atac = atac,
  segmentation = segmentation,
  normalisation_factors = rbind(norm_rna, norm_atac),
  rna_likelihood = "NB", 
  atac_likelihood = 'NB',
  description = 'Lucrezia RNA/ATAC')

# Errors in the fit (NA!)
hyperparams <- auto_config_run(x, 1:4)
fit <- fit_congas(x, K = 1:4, learning_rate = 0.05, steps = 500, model_parameters = hyperparams)

# Select segments
pdf("segmentation.pdf", width = 8, height = 4)
for(segment in x$input$segmentation$segment_id)
  plot_data(x, what = 'histogram', segments = segment) %>% print
dev.off()

# Inspected
segments_retain =c("chr3:129807186:195300676", 
                   "chr4:40332:190986668",
                   "chr6:142840:29844875",
                   "chr6:32604274:57301860",
                   "chr8:102060353:146295114",
                   "chr17:41382596:79601279",
                   "chr18:18917714:77990447",
                   "chr20:1895949:25700309",
                   "chr20:30437522:50255762",
                   "chr21:11186714:48117308"
                   )

segmentation = segmentation %>% filter(segment_id %in% segments_retain)

x = init(
  rna = rna,
  atac = atac,
  segmentation = segmentation,
  normalisation_factors = rbind(norm_rna, norm_atac),
  rna_likelihood = "NB", 
  atac_likelihood = 'NB',
  description = 'Lucrezia RNA/ATAC, 10 selected segments, balanced cells')

# Still explodes!
hyperparams <- auto_config_run(x, 1:4)
fit <- fit_congas(x, K = 1:4, learning_rate = 0.05, steps = 500, model_parameters = hyperparams)
plot_fit_scores(fit)

plot_fit(fit, 'density', highlights = FALSE)

# Data plots
plot_data(x, what = 'histogram')
plot_data(x, what = 'lineplot')
plot_data(x, what = 'heatmap')

save(x, file = 's.RData')

# c(extraDistr::rlaplace(1000),extraDistr::rlaplace(400, mu = 2)) %>% hist(breaks = 100)

x$input$dataset %>% filter()

summary(x$input$segmentation$ATAC_peaks)

ggplot(x$input$segmentation, aes(x = ATAC_peaks, y = RNA_genes)) +geom_point()

plot_data(x, 'mapping')
ggsave("a.png", width = 5, height = 13)

hist(x$input$segmentation$ATAC_peaks, breaks = 50)
hist(x$input$segmentation$RNA_genes, breaks = 50)

x = init(
  rna = rna %>% filter_known_genes(),
  atac = atac,
  segmentation = segmentation,
  normalisation_factors = rbind(norm_rna, norm_atac),
  rna_likelihood = "G", 
  atac_likelihood = 'NB',
  description = 'Lucrezia RNA/ATAC (segs. with > 100 ATAC peaks and 20 genes)')

plot_data(x)

# Generate report
plot_data(x, what = 'lineplot') %>% ggsave(filename = 'lineplot.png', width = 12, height = 8)

pdf("segmentation.pdf", width = 8, height = 4)
for(segment in x$input$segmentation$segment_id)
  plot_data(x, what = 'histogram', segments = segment) %>% print
dev.off()

plot_data(x, what = 'heatmap') %>% ggsave(filename = 'hmap.png', width = 12, height = 18)

# 
table_celltypes = readr::read_tsv("~/Downloads/atac_rna_celltypes.tsv", col_names = FALSE)
colnames(table_celltypes) = c("cell", 'type')

table_celltypes$type %>% table

tum_cells = table_celltypes %>% filter(type == "Tumor")
imm_cells = table_celltypes %>% filter(type == "Immune") %>% sample_n(1256)
tcl_cells = table_celltypes %>% filter(type == "T-cells") %>% sample_n(1256)

subset_cells = c(tum_cells$cell, imm_cells$cell, tcl_cells$cell)


x = init(
  rna = rna %>% filter_known_genes() %>% filter(cell %in% subset_cells),
  atac = atac %>% filter(cell %in% subset_cells),
  segmentation = segmentation,
  normalisation_factors = rbind(norm_rna, norm_atac)  %>% filter(cell %in% subset_cells),
  rna_likelihood = "NB", 
  atac_likelihood = 'NB',
  description = 'Lucrezia RNA/ATAC (segs. with > 100 ATAC peaks and 20 genes), randomised 1256 cells per class')

# Generate report
plot_data(x, what = 'lineplot') %>% ggsave(filename = 'lineplot.png', width = 12, height = 8)

pdf("segmentation.pdf", width = 8, height = 4)
for(segment in x$input$segmentation$segment_id)
  plot_data(x, what = 'histogram', segments = segment) %>% print
dev.off()

plot_data(x, what = 'heatmap') %>% ggsave(filename = 'hmap.png', width = 12, height = 18)




