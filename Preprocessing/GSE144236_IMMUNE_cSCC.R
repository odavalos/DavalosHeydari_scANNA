# This script is for processing the GSE144236 dataset
# for the scANNA manuscript.

# # Download data -----------------------------------------------------------
# # download data from GEO

# curl "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144236/suppl/GSE144236%5FcSCC%5Fcounts%2Etxt%2Egz" -o "GSE144236_cSCC_counts.txt.gz"

# # download metadata from GEO
# curl "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144236/suppl/GSE144236%5Fpatient%5Fmetadata%5Fnew%2Etxt%2Egz" -o "GSE144236_patient_metadata_new.txt.gz"

# # make directory for data
# mkdir GSE144236

# # move data to directory
# mv GSE144236_cSCC_counts.txt.gz GSE144236/GSE144236_cSCC_counts.txt.gz
# mv GSE144236_patient_metadata_new.txt.gz GSE144236/GSE144236_patient_metadata_new.txt.gz

# # we used python to read in the very large file (challenging for R) 
# # and saved it as a .h5ad file 

# # this is the code we used to do that

# # import libraries
# import numpy as np
# import pandas as pd 
# import matplotlib.pyplot as plt
# import scanpy as sc

# # read in data
# adata = sc.read_text('~/GSE144236/GSE144236_cSCC_counts.txt', delimiter='\t', first_column_names=True)

# # transpose the dataset since it is read in flipped
# adata = adata.transpose()

# # remove junk names like 'Patient' and 'Tissue: 0=Normal, 1=Tumor'
# names_to_remove = ['Patient', 'Tissue: 0=Normal, 1=Tumor']
# non_junk_genes_list = [name for name in adata.var.index if name not in names_to_remove]
# adata = adata[:, non_junk_genes_list]

# # minimal filtering
# sc.pp.filter_cells(adata, min_genes=200)
# sc.pp.filter_genes(adata, min_cells=3)

# # add qc metrics
# adata.var["percent_mito"] = adata.var_names.str.startswith("MT-")
# sc.pp.calculate_qc_metrics(adata, qc_vars=["percent_mito"], percent_top=None, log1p=False, inplace=True)

# # save as h5ad file
# adata.write_h5ad('~/GSE144236/GSE144236_raw.h5ad')


# Loading packages and setting up directories -----------------------------
library(tidyverse)
library(ggthemes)
library(Seurat)
library(harmony)
library(SeuratDisk)


sample_no <- 'GSE144236' # change accordingly

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

# Load the Immune_cSCC dataset
Convert("~/GSE144236/GSE144236_raw.h5ad", dest = "h5seurat", overwrite = FALSE)
immune_cscc <- LoadH5Seurat("~/GSE144236/GSE144236_raw.h5seurat")
immune_cscc


# add barcodes to metadata for filtering
immune_cscc@meta.data$barcodes <- rownames(immune_cscc@meta.data)

# read in all metadata
metadata <- read.csv(file = '~/GSE144236/patient_metadata_new.txt', sep = '\t', row.names = 1)
metadata$barcodes <- rownames(metadata)

# join metadata to immune_cscc object
metadata <- metadata %>%
  rename(celltypes = level1_celltype, 
         celltypes_lvl2 = level2_celltype, 
         celltypes_lvl3 = level3_celltype) %>%
  select(-nCount_RNA, -nFeature_RNA)

immune_cscc@meta.data <- left_join(immune_cscc@meta.data, metadata, by = 'barcodes')

# converting some metadata columns to Seurat names 
immune_cscc@meta.data <- immune_cscc@meta.data %>%
  rename(nCount_RNA = total_counts,
         nFeature_RNA = n_genes) %>%
  select(-total_counts_percent_mito, 
         -pct_counts_percent_mito, 
         -n_genes_by_counts)

# adding barcodes back as rownames
rownames(immune_cscc@meta.data) <- immune_cscc@meta.data$barcodes

# check if data has been processed or if it is raw
immune_cscc@assays$RNA@counts[c(1:1000),c(1:10)] %>% colSums()

# check NA quantity
sum(is.na(immune_cscc@meta.data$celltypes))

# subset out Multiplet since this is junk
immune_cscc <- subset(immune_cscc, subset = celltypes != 'Multiplet')

# save raw object with metadata included
saveRDS(immune_cscc,
        file = paste0(workingdir,
                      '/',
                      datadir,
                      '/',
                      sample_no,
                      "_raw.rds"))

# view cell type proportions
propsmd <- immune_cscc@meta.data %>%
  group_by(celltypes) %>%
  summarise(n = n())

p_cellquants <- ggplot(propsmd, aes(reorder(celltypes, -n), n, fill = celltypes)) + 
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = qual_col_pals$colors) +
  xlab('') +
  ylab('Total') + 
  theme_bw(base_size = 14) + 
  theme(axis.text.x = element_blank(),
        axis.text = element_text(color='black'),
        axis.ticks.x.bottom = element_blank())

p_cellquants_log10 <- ggplot(propsmd, aes(reorder(celltypes, -n), n, fill = celltypes)) + 
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = qual_col_pals$colors) +
  scale_y_log10() + 
  xlab('') +
  ylab('log10(Total)') + 
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_blank(),
        axis.text = element_text(color='black'),
        axis.ticks.x.bottom = element_blank())

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


# calculate percentages of mitochondrial genes
immune_cscc <- PercentageFeatureSet(immune_cscc, pattern = "^MT-", col.name = "percent_mito")

# Visualize QC metrics as a violin plot
feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito")
VlnPlot(immune_cscc, group.by = "tum.norm", features = feats, pt.size = 0.1, ncol = 3) + 
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
p_ctsmt <- FeatureScatter(immune_cscc, feature1 = "nCount_RNA", feature2 = "percent_mito") + geom_hline(yintercept = 10, linetype='dashed') + 
  scale_color_manual(values = qual_col_pals$colors)
p_ctsfts <- FeatureScatter(immune_cscc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_hline(yintercept = 7000, linetype='dashed') + 
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

# Subset based on QC metrics
immune_cscc <- subset(immune_cscc, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent_mito < 10) # could change to 20 percent instead
gc()

# Normalize the data
immune_cscc <- NormalizeData(immune_cscc)
gc()

# Identify variable features (5K features for use with scANNA)
immune_cscc <- FindVariableFeatures(immune_cscc, selection.method = "vst", nfeatures = 5000)
gc()

# Identify the 20 most highly variable genes
top20 <- head(VariableFeatures(immune_cscc), 20)

# plot top20 variable features
vplot <- VariableFeaturePlot(immune_cscc)
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
  height = 6,
  width = 8,
  units = "in",
  dpi = 300
)
gc()


# save object with variable features
saveRDS(immune_cscc,
        file = paste0(workingdir,
                      '/',
                      datadir,
                      '/',
                      sample_no,
                      "_qc_hvg.rds"))

gc()

# Scale the data
all.genes <- rownames(immune_cscc)
immune_cscc <- ScaleData(immune_cscc, features = all.genes)
gc()

# add a bool list of variable genes to the metadata
variable_genes <- VariableFeatures(immune_cscc)
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
immune_cscc <- RunPCA(immune_cscc, features = VariableFeatures(object = immune_cscc))
gc()

ElbowPlot(immune_cscc, ndims = 50)
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

# harmonize the dataset and remove batch effects
immune_cscc <- RunHarmony(immune_cscc, group.by.vars = c("tum.norm", "patient"), reduction = "pca", dims.use = 1:50, assay.use = "RNA", plot_convergence = TRUE)
gc()

# Clustering
use.pcs <- 1:50

immune_cscc <- FindNeighbors(immune_cscc, 
                            reduction="harmony",
                            dims = use.pcs)

# run the clustering for the multiple resolutions
immune_cscc <- FindClusters(object = immune_cscc, 
                           resolution = seq(0.20, 1, 0.20) # change resolution range here
)

# Number of unique clusters generated per resolution
grep("res", colnames(immune_cscc@meta.data), value = TRUE) %>%
  purrr::map_chr(~ paste(.x, "--> clusters generated:", length(unique(
    immune_cscc@meta.data[, .x]
  ))))

Idents(immune_cscc) <- 'RNA_snn_res.0.2'
table(immune_cscc@active.ident)


# Plotting UMAP -----------------------------------------------------------

Idents(immune_cscc) <- immune_cscc@meta.data$celltypes

immune_cscc <- RunUMAP(immune_cscc, 
                      reduction = "harmony", 
                      dims = 1:20,
                      umap.method = 'umap-learn')

p_conditions <- DimPlot(immune_cscc,
                        reduction = "umap",
                        label = TRUE,
                        repel = TRUE,
                        #label.size = 8, 
                        group.by = 'tum.norm') +
  scale_color_manual(values = qual_col_pals$colors) +
  ggtitle('Conditions') + 
  theme(plot.title = element_text(hjust = 0.5)) +
  NoLegend()

p_clusters <- DimPlot(immune_cscc,
                      reduction = "umap",
                      label = TRUE,
                      repel = TRUE,
                      #label.size = 8, 
                      group.by = 'RNA_snn_res.0.2') +
  scale_color_manual(values = qual_col_pals$colors) +
  ggtitle('Seurat_Clusters') + 
  theme(plot.title = element_text(hjust = 0.5)) +
  NoLegend()

p_celltypes <- DimPlot(immune_cscc,
                       reduction = "umap",
                       label = TRUE,
                       repel = TRUE,
                       #label.size = 8, 
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


DimPlot(immune_cscc,
        reduction = "umap",
        label = TRUE,
        repel = TRUE,
        label.size = 8, 
        group.by = 'celltypes') + 
  ggthemes::scale_colour_tableau(palette = 'Tableau 20', direction = -1) +
  xlab('UMAP 1') +
  ylab('UMAP 2') +
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



# save the object with the annotations
saveRDS(immune_cscc,
        file = paste0(workingdir,
                      '/',
                      datadir,
                      '/',
                      sample_no,
                      "_qc_hvg_anno_full.rds"))


# H5AD Object Creation ----------------------------------------------------

# saving the integrated file for analysis in scanpy
library(SeuratDisk)

# remove scale.data so that raw counts gets saved to raw slot and not normalized counts
immune_cscc <- DietSeurat(immune_cscc, 
                         counts = TRUE, 
                         data = TRUE, 
                         scale.data = FALSE,
                         dimreducs = c('harmony', 'pca', 'umap'), 
                         graphs = c('RNA_nn', 'RNA_snn'))

# subset (5k variable genes)
variable_genes <- VariableFeatures(immune_cscc)

immune_cscc_5k <- subset(immune_cscc, features = variable_genes)

# save the object with the annotations
saveRDS(immune_cscc_5k,
        file = paste0(workingdir,
                      '/',
                      datadir,
                      '/',
                      sample_no,
                      "_qc_hvg_anno_5k.rds"))


SaveH5Seurat(immune_cscc_5k, filename = "~/GSE144236/GSE144236_qc_hvg_anno_5k.h5Seurat")
Convert("~/GSE144236/GSE144236_qc_hvg_anno_5k.h5Seurat", dest = "h5ad")



# Differential Expression -------------------------------------------------

DEWilcox <- FindAllMarkers(immune_cscc, log2FC.threshold = 0.25, test.use = "wilcox",
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



DotPlot(immune_cscc, 
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




