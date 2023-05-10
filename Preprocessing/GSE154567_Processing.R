# This script is used to process the GSE154567 dataset 
# for the scANNA manuscript.

# # Download data -----------------------------------------------------------
# # download data from GEO

# curl "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE154567&format=file" -o "GSE154567_RAW.tar"

# # download metadata from GEO
# curl "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154567/suppl/GSE154567%5Fannotation%2Ecsv%2Egz" -o "GSE154567_annotation.csv.gz"

# # make directory for data
# mkdir GSE154567

# # move data to directory
# mv GSE154567_RAW.tar GSE154567/GSE154567_RAW.tar
# mv GSE154567_annotation.csv.gz GSE154567/GSE154567_annotation.csv.gz


# Loading packages and setting up directories -----------------------------
library(tidyverse)
library(ggthemes)
library(Seurat)
library(harmony)


sample_no <- 'GSE154567' # change accordingly

# directories
workingdir <- paste0(getwd(), '/',sample_no)
plotdir <- 'plots'
datadir <- 'data'
resultsdir <- 'results'

# check if sample directory exists
ifelse(!dir.exists(file.path(workingdir)),
       dir.create(file.path(workingdir)), FALSE)

# check for plot directory
ifelse(!dir.exists(file.path(workingdir, plotdir)),
       dir.create(file.path(workingdir, plotdir)), FALSE)

# check for data directory
ifelse(!dir.exists(file.path(workingdir, datadir)),
       dir.create(file.path(workingdir, datadir)), FALSE)

# check for results directory
ifelse(!dir.exists(file.path(workingdir, resultsdir)),
       dir.create(file.path(workingdir, resultsdir)), FALSE)


# Read and QC data --------------------------------------------------------

# metadata
metadata <- read_csv(file = '/Volumes/BigDrive/GradSchool/NACT_NewDatasets/GSE154567/GSE154567_annotation.csv')
metadata <- metadata %>%
  rename(barcodes = `...1`,
         celltypes = celltype5)


# read in data CM1 and CM2
CM1_2 <- Read10X_h5(file = "/Volumes/BigDrive/GradSchool/NACT_NewDatasets/GSE154567/GSE154567_RAW/raw_h5/GSM4674661_CM1_CM2_filtered_feature_bc_matrix.h5")
CM1_2 <- CreateSeuratObject(counts = CM1_2, project = "CM1_2", min.cells = 3, min.features = 200)
CM1_2@meta.data$barcodes <- rownames(CM1_2@meta.data)
# temp subseted metadata df
tmp_MD <- metadata %>%
  filter(orig.ident == 'CM1_2') %>%
  mutate(barcodes = str_replace_all(barcodes, '_[0-9]', ''))

# merge metadata to seurat object
CM1_2@meta.data <- left_join(CM1_2@meta.data, tmp_MD, by='barcodes') %>% # append metadata 
  select(-orig.ident.y) %>%
  rename(orig.ident = orig.ident.x)
rownames(CM1_2@meta.data) <- CM1_2@meta.data$barcodes # re-append barcodes as rownames 
rm(tmp_MD)

# read in data CM3 and CM4
CM3_4 <- Read10X_h5(file = "/Volumes/BigDrive/GradSchool/NACT_NewDatasets/GSE154567/GSE154567_RAW/raw_h5/GSM4674662_CM3_CM4_filtered_feature_bc_matrix.h5")
CM3_4 <- CreateSeuratObject(counts = CM3_4, project = "CM3_4", min.cells = 3, min.features = 200)
CM3_4@meta.data$barcodes <- rownames(CM3_4@meta.data)
# temp subseted metadata df
tmp_MD <- metadata %>%
  filter(orig.ident == 'CM3_4') %>%
  mutate(barcodes = str_replace_all(barcodes, '_[0-9]', ''))

# merge metadata to seurat object
CM3_4@meta.data <- left_join(CM3_4@meta.data, tmp_MD, by='barcodes') %>% # append metadata 
  select(-orig.ident.y) %>%
  rename(orig.ident = orig.ident.x)
rownames(CM3_4@meta.data) <- CM3_4@meta.data$barcodes # re-append barcodes as rownames 
rm(tmp_MD)

# read in data CM5
CM5 <- Read10X_h5(file = "/Volumes/BigDrive/GradSchool/NACT_NewDatasets/GSE154567/GSE154567_RAW/raw_h5/GSM4674663_CM5_filtered_feature_bc_matrix.h5")
CM5 <- CreateSeuratObject(counts = CM5, project = "CM5", min.cells = 3, min.features = 200)
CM5@meta.data$barcodes <- rownames(CM5@meta.data)
# temp subseted metadata df
tmp_MD <- metadata %>%
  filter(orig.ident == 'CM5') %>%
  mutate(barcodes = str_replace_all(barcodes, '_[0-9]', ''))

# merge metadata to seurat object
CM5@meta.data <- left_join(CM5@meta.data, tmp_MD, by='barcodes') %>% # append metadata 
  select(-orig.ident.y) %>%
  rename(orig.ident = orig.ident.x)
rownames(CM5@meta.data) <- CM5@meta.data$barcodes # re-append barcodes as rownames 
rm(tmp_MD)

# read in data CS7 and CS8
CS7_8<- Read10X_h5(file = "/Volumes/BigDrive/GradSchool/NACT_NewDatasets/GSE154567/GSE154567_RAW/raw_h5/GSM4674664_CS7_CS8_filtered_feature_bc_matrix.h5")
CS7_8 <- CreateSeuratObject(counts = CS7_8, project = "SeuratProject", min.cells = 3, min.features = 200) # original authors named this project 'SeuratProject'
CS7_8@meta.data$barcodes <- rownames(CS7_8@meta.data)
# temp subseted metadata df
tmp_MD <- metadata %>%
  filter(orig.ident == 'CS7_8') %>%
  mutate(barcodes = str_replace_all(barcodes, '_[0-9]', ''))

# merge metadata to seurat object
CS7_8@meta.data <- left_join(CS7_8@meta.data, tmp_MD, by='barcodes') %>% # append metadata 
  select(-orig.ident.y) %>%
  rename(orig.ident = orig.ident.x)
rownames(CS7_8@meta.data) <- CS7_8@meta.data$barcodes # re-append barcodes as rownames 
rm(tmp_MD)

# read in data CS9 and CS10
CS9_10<- Read10X_h5(file = "/Volumes/BigDrive/GradSchool/NACT_NewDatasets/GSE154567/GSE154567_RAW/raw_h5/GSM4674665_CS9_CS10_filtered_feature_bc_matrix.h5")
CS9_10 <- CreateSeuratObject(counts = CS9_10, project = "CS9_10", min.cells = 3, min.features = 200)
CS9_10@meta.data$barcodes <- rownames(CS9_10@meta.data)
# temp subseted metadata df
tmp_MD <- metadata %>%
  filter(orig.ident == 'CS9_10') %>%
  mutate(barcodes = str_replace_all(barcodes, '_[0-9]', ''))

# merge metadata to seurat object
CS9_10@meta.data <- left_join(CS9_10@meta.data, tmp_MD, by='barcodes') %>% # append metadata 
  select(-orig.ident.y) %>%
  rename(orig.ident = orig.ident.x)
rownames(CS9_10@meta.data) <- CS9_10@meta.data$barcodes # re-append barcodes as rownames 
rm(tmp_MD)

# read in data CS11 and CS12
CS11_12<- Read10X_h5(file = "/Volumes/BigDrive/GradSchool/NACT_NewDatasets/GSE154567/GSE154567_RAW/raw_h5/GSM4674666_CS11_CS12_filtered_feature_bc_matrix.h5")
CS11_12 <- CreateSeuratObject(counts = CS11_12, project = "CS11_12", min.cells = 3, min.features = 200)
CS11_12@meta.data$barcodes <- rownames(CS11_12@meta.data)
# temp subseted metadata df
tmp_MD <- metadata %>%
  filter(orig.ident == 'CS11_12') %>%
  mutate(barcodes = str_replace_all(barcodes, '_[0-9]', ''))

# merge metadata to seurat object
CS11_12@meta.data <- left_join(CS11_12@meta.data, tmp_MD, by='barcodes') %>% # append metadata 
  select(-orig.ident.y) %>%
  rename(orig.ident = orig.ident.x)
rownames(CS11_12@meta.data) <- CS11_12@meta.data$barcodes # re-append barcodes as rownames 
rm(tmp_MD)


# read in data CR13 and CR14
CR13_14<- Read10X_h5(file = "/Volumes/BigDrive/GradSchool/NACT_NewDatasets/GSE154567/GSE154567_RAW/raw_h5/GSM4674667_CR13_CR14_filtered_feature_bc_matrix.h5")
CR13_14 <- CreateSeuratObject(counts = CR13_14, project = "CR13_14", min.cells = 3, min.features = 200)
CR13_14@meta.data$barcodes <- rownames(CR13_14@meta.data)
# temp subseted metadata df
tmp_MD <- metadata %>%
  filter(orig.ident == 'CR13_14') %>%
  mutate(barcodes = str_replace_all(barcodes, '_[0-9]', ''))

# merge metadata to seurat object
CR13_14@meta.data <- left_join(CR13_14@meta.data, tmp_MD, by='barcodes') %>% # append metadata 
  select(-orig.ident.y) %>%
  rename(orig.ident = orig.ident.x)
rownames(CR13_14@meta.data) <- CR13_14@meta.data$barcodes # re-append barcodes as rownames 
rm(tmp_MD)

# read in data CR15 and CR16
CR15_16<- Read10X_h5(file = "/Volumes/BigDrive/GradSchool/NACT_NewDatasets/GSE154567/GSE154567_RAW/raw_h5/GSM4674668_CR15_CR16_filtered_feature_bc_matrix.h5")
CR15_16 <- CreateSeuratObject(counts = CR15_16, project = "CR15_16", min.cells = 3, min.features = 200)
CR15_16@meta.data$barcodes <- rownames(CR15_16@meta.data)
# temp subseted metadata df
tmp_MD <- metadata %>%
  filter(orig.ident == 'CR15_16') %>%
  mutate(barcodes = str_replace_all(barcodes, '_[0-9]', ''))

# merge metadata to seurat object
CR15_16@meta.data <- left_join(CR15_16@meta.data, tmp_MD, by='barcodes') %>% # append metadata 
  select(-orig.ident.y) %>%
  rename(orig.ident = orig.ident.x)
rownames(CR15_16@meta.data) <- CR15_16@meta.data$barcodes # re-append barcodes as rownames 
rm(tmp_MD)

# read in data CR17 and CR18
CR17_18<- Read10X_h5(file = "/Volumes/BigDrive/GradSchool/NACT_NewDatasets/GSE154567/GSE154567_RAW/raw_h5/GSM4674669_CR17_CR18_filtered_feature_bc_matrix.h5")
CR17_18 <- CreateSeuratObject(counts = CR17_18, project = "CR17_18", min.cells = 3, min.features = 200)
CR17_18@meta.data$barcodes <- rownames(CR17_18@meta.data)
# temp subseted metadata df
tmp_MD <- metadata %>%
  filter(orig.ident == 'CR17_18') %>%
  mutate(barcodes = str_replace_all(barcodes, '_[0-9]', ''))

# merge metadata to seurat object
CR17_18@meta.data <- left_join(CR17_18@meta.data, tmp_MD, by='barcodes') %>% # append metadata 
  select(-orig.ident.y) %>%
  rename(orig.ident = orig.ident.x)
rownames(CR17_18@meta.data) <- CR17_18@meta.data$barcodes # re-append barcodes as rownames 
rm(tmp_MD)

# merge all objects together
pbmc_covid <- merge(CM1_2, c(CM3_4, CM5, CS7_8, CS9_10,
                             CS11_12, CR13_14, CR15_16, CR17_18))

# remove objects to save space
rm(CM1_2, CM3_4, CM5, CS7_8, CS9_10, 
   CS11_12, CR13_14, CR15_16, CR17_18)

# clean up memory
gc()


# check if data has been processed or if it is raw
pbmc_covid@assays$RNA@counts[c(1:1000),c(1:10)] %>% colSums()

# check NA quantity
sum(is.na(pbmc_covid@meta.data$celltypes))

# rename NA's as 'DROP'
pbmc_covid@meta.data$celltypes <- if_else(is.na(pbmc_covid@meta.data$celltypes), 'DROP', pbmc_covid@meta.data$celltypes)

# subset out NA's labeled as 'DROP'
pbmc_covid <- subset(pbmc_covid, subset = celltypes != 'DROP')

# subset out platelets and pDCs since less than 100 cells per celltype
pbmc_covid <- subset(pbmc_covid, subset = celltypes != 'platelets')
pbmc_covid <- subset(pbmc_covid, subset = celltypes != 'pDCs')


# save raw object with metadata included
saveRDS(pbmc_covid,
        file = paste0(workingdir,
                      '/',
                      datadir,
                      '/',
                      sample_no,
                      "_raw.rds"))

# view cell type proportions
propsmd <- pbmc_covid@meta.data %>%
  group_by(celltypes) %>%
  summarise(n = n())

p_cellquants <- ggplot(propsmd, aes(celltypes, n, fill = celltypes)) + 
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = qual_col_pals$colors) +
  ylab('n') + 
  theme_bw(base_size = 14) + 
  theme(axis.text.x = element_blank(),
        axis.text = element_text(color='black'))

p_cellquants_log10 <- ggplot(propsmd, aes(celltypes, n, fill = celltypes)) + 
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = qual_col_pals$colors) +
  scale_y_log10() + 
  ylab('n log10 scale') + 
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_blank(),
        axis.text = element_text(color='black'))

cowplot::plot_grid(p_cellquants, p_cellquants_log10)
ggsave(
  paste0(
    workingdir,
    "/",
    plotdir,
    "/",
    sample_no,
    "_cell_props.pdf"
  ),
  height = 6,
  width = 12,
  units = "in",
  dpi = 300
)

# calculate percentages of mito and ribo genes
pbmc_covid <- PercentageFeatureSet(pbmc_covid, pattern = "^MT-", col.name = "percent_mito")

# Visualize QC metrics as a violin plot
feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito")
VlnPlot(pbmc_covid, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 3) + 
  NoLegend()

ggsave(
  paste0(
    workingdir,
    "/",
    plotdir,
    "/",
    sample_no,
    "_qc_violin.pdf"
  ),
  height = 8,
  width = 12,
  units = "in",
  dpi = 300
)


# Visualize QC metrics as a scatter plot
p_ctsmt <- FeatureScatter(pbmc_covid, feature1 = "nCount_RNA", feature2 = "percent_mito") + geom_hline(yintercept = 5, linetype='dashed') + 
  scale_color_manual(values = qual_col_pals$colors)
p_ctsfts <- FeatureScatter(pbmc_covid, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  scale_color_manual(values = qual_col_pals$colors)

cowplot::plot_grid(p_ctsfts, p_ctsmt)
ggsave(
  paste0(
    workingdir,
    "/",
    plotdir,
    "/",
    sample_no,
    "_qc_scatter.pdf"
  ),
  height = 6,
  width = 12,
  units = "in",
  dpi = 300
)

# subset out cells with less than 200 genes and more than 2500 genes and more than 5 percent mito
pbmc_covid <- subset(pbmc_covid, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent_mito < 5)

# normalize data
pbmc_covid <- NormalizeData(pbmc_covid)

# Identify variable features (5K features for use with scANNA)
pbmc_covid <- FindVariableFeatures(pbmc_covid, selection.method = "vst", nfeatures = 5000)

# Identify the 20 most highly variable genes
top20 <- head(VariableFeatures(pbmc_covid), 20)

# plot top20 variable features
vplot <- VariableFeaturePlot(pbmc_covid)
vplot <- LabelPoints(plot = vplot, points = top20, repel = TRUE)
vplot + ggtitle('Top 20 Variable Features', subtitle = paste0('Sample:',sample_no))

ggsave(
  paste0(
    workingdir,
    "/",
    plotdir,
    "/",
    sample_no,
    "_varfeatures.pdf"
  ),
  height = 8,
  width = 12,
  units = "in",
  dpi = 300
)


# save object with variable features
saveRDS(pbmc_covid,
        file = paste0(workingdir,
                      '/',
                      datadir,
                      '/',
                      sample_no,
                      "_qc_hvg.rds"))


# scale data
all.genes <- rownames(pbmc_covid)
pbmc_covid <- ScaleData(pbmc_covid, features = all.genes)

# add a bool list of variable genes to the metadata
variable_genes <- VariableFeatures(pbmc_covid)
hvgs <- all.genes %in% variable_genes
hvg_df <- data.frame(allgenes = all.genes,
                     hvgs = hvgs)


write_csv(hvg_df, 
          file=paste0(workingdir,
                      "/",
                      resultsdir,
                      "/",
                      sample_no,
                      "_hvg_genelist.csv"))

# Dimension Reduction/Clustering ------------------------------------------

# PCA
pbmc_covid <- RunPCA(pbmc_covid, features = VariableFeatures(object = pbmc_covid))

ElbowPlot(pbmc_covid, ndims = 50)

ggsave(
  paste0(
    workingdir,
    "/",
    plotdir,
    "/",
    sample_no,
    "_dr_elbowPCA.pdf"
  ),
  height = 8,
  width = 12,
  units = "in",
  dpi = 300
)


# Harmony batch correction
pbmc_covid <- RunHarmony(pbmc_covid, group.by.vars = "patient", reduction = "pca",
                      dims.use = 1:50, assay.use = "RNA", plot_convergence = TRUE)

# Clustering
use.pcs <- 1:50

pbmc_covid <- FindNeighbors(pbmc_covid, 
                         reduction="pca",
                         dims = use.pcs)

# run the clustering for the multiple resolutions
pbmc_covid <- FindClusters(object = pbmc_covid, 
                        resolution = seq(0.20, 1, 0.20) # change resolution range here
)

# Number of unique clusters generated per resolution
grep("res", colnames(pbmc_covid@meta.data), value = TRUE) %>%
  purrr::map_chr(~ paste(.x, "--> clusters generated:", length(unique(
    pbmc_covid@meta.data[, .x]
  ))))

Idents(pbmc_covid) <- 'RNA_snn_res.0.2'
table(pbmc_covid@active.ident)


# Plotting UMAP -----------------------------------------------------------

Idents(pbmc_covid) <- pbmc_covid@meta.data$celltypes

pbmc_covid <- RunUMAP(pbmc_covid, 
                   reduction = "harmony", 
                   dims = use.pcs,
                   umap.method = 'umap-learn')

p_conditions <- DimPlot(pbmc_covid,
                        reduction = "umap",
                        label = TRUE,
                        repel = TRUE,
                        group.by = 'group') +
  scale_color_manual(values = qual_col_pals$colors) +
  ggtitle('Conditions') + 
  theme(plot.title = element_text(hjust = 0.5)) +
  NoLegend()

p_clusters <- DimPlot(pbmc_covid,
                      reduction = "umap",
                      label = TRUE,
                      repel = TRUE,
                      group.by = 'RNA_snn_res.0.2') +
  scale_color_manual(values = qual_col_pals$colors) +
  ggtitle('Seurat_Clusters') + 
  theme(plot.title = element_text(hjust = 0.5)) +
  NoLegend()

p_celltypes <- DimPlot(pbmc_covid,
                       reduction = "umap",
                       label = TRUE,
                       repel = TRUE,
                       group.by = 'celltypes') +
  scale_color_manual(values = qual_col_pals$colors) +
  ggtitle('Celltypes') + 
  theme(plot.title = element_text(hjust = 0.5)) +
  NoLegend()


cowplot::plot_grid(p_clusters, p_celltypes, p_conditions)
ggsave(
  paste0(
    workingdir,
    "/",
    plotdir,
    "/",
    sample_no,
    "_dr_clustered_UMAP.pdf"
  ),
  height = 8,
  width = 12,
  units = "in",
  dpi = 300
)


# UMAP celltypes
DimPlot(pbmc_covid,
        reduction = "umap",
        label = TRUE,
        repel = TRUE,
        label.size = 8, 
        group.by = 'celltypes') + 
  ggthemes::scale_colour_tableau(palette = 'Tableau 20', direction = -1) +
  xlab('UMAP 1') +
  ylab('UMAP 2') +
  ggtitle('Celltypes') + 
  theme(plot.title = element_text(hjust = 0.5)) +
  NoLegend()


ggsave(
  paste0(
    workingdir,
    "/",
    plotdir,
    "/",
    sample_no,
    "_dr_clustered_UMAP_celltypes.pdf"
  ),
  height = 8,
  width = 12,
  units = "in",
  dpi = 300
)

# save the full object with the annotations
saveRDS(pbmc_covid,
        file = paste0(workingdir,
                      '/',
                      datadir,
                      '/',
                      sample_no,
                      "_qc_hvg_anno_full.rds"))


# remove scale.data so that raw counts gets saved to raw slot and not normalized counts
pbmc_covid <- DietSeurat(pbmc_covid, 
                         counts = TRUE, 
                         data = TRUE, 
                         scale.data = FALSE,
                         dimreducs = c('harmony', 'pca', 'umap'), 
                         graphs = c('RNA_nn', 'RNA_snn'))

# get the variable genes
variable_genes <- VariableFeatures(pbmc_covid)

# subset the object to only include the variable genes
pbmc_covid <- subset(pbmc_covid, features = variable_genes)

# save the object with the annotations
saveRDS(pbmc_covid,
        file = paste0(workingdir,
                      '/',
                      datadir,
                      '/',
                      sample_no,
                      "_qc_hvg_anno_5k.rds"))

# H5AD Object Creation ----------------------------------------------------

# saving the integrated file for analysis in scanpy & scANNA
library(SeuratDisk) 

SaveH5Seurat(pbmc_covid, filename = "~/GSE154567/GSE154567_qc_hvg_anno_new.h5Seurat")
Convert("~/GSE154567/GSE154567_qc_hvg_anno_new.h5Seurat", dest = "h5ad")



# Differential Expression -------------------------------------------------

DEWilcox <- FindAllMarkers(pbmc_covid, log2FC.threshold = 0.25, test.use = "wilcox",
                           min.pct = 0.1, min.diff.pct = 0.1, 
                           only.pos = TRUE, assay = "RNA")

write_csv(DEWilcox, 
          file=paste0(workingdir,
                      "/",
                      resultsdir,
                      "/",
                      sample_no,
                      "_diffexpression_wilcox.csv"))


top25 <- DEWilcox %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  top_n(n = 25, wt = avg_log2FC)
top25
table(top25$cluster)


top5 <- DEWilcox %>% 
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)



DotPlot(pbmc_covid, 
        features = rev(as.character(unique(top5$gene))),
        assay = "RNA", 
        cluster.idents = FALSE) +
  scale_color_viridis_c(option = 'viridis') + 
  xlab('Gene') + 
  ylab('') +
  coord_flip() +
  theme_classic(base_size = 18) + 
  theme(axis.text = element_text(colour = 'black'),
        legend.text = element_text(color = 'black'),
        legend.position = 'top')

ggsave(
  paste0(
    workingdir,
    "/",
    plotdir,
    "/",
    sample_no,
    "_top5_dotplot.pdf"
  ),
  height = 8,
  width = 12,
  units = "in",
  dpi = 300
)




