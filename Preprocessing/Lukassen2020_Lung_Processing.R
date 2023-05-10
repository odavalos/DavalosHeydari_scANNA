# This script processes the data from Lukassen et al. 2020
# for the scANNA manuscript.


# Download data -----------------------------------------------------------
# download data from figshare links provided by authors
# https://doi.org/10.6084/m9.figshare.11981034.v1

# curl "https://figshare.com/ndownloader/files/21999237" -o "Counts_lung_cells.csv"

# # download metadata from figshare
# curl "https://figshare.com/ndownloader/files/21999240" -o "Metadata_lung_cells.csv"

# # make directory for data
# mkdir Lukassen2020_Lung

# # move data to directory
# mv Counts_lung_cells.csv Lukassen2020_Lung/Counts_lung_cells.csv
# mv Metadata_lung_cells.csv Lukassen2020_Lung/Metadata_lung_cells.csv


# Loading packages and setting up directories -----------------------------
library(tidyverse)
library(ggthemes)
library(Seurat)
library(harmony)

sample_no <- 'Lukassen2020_Lung' # change accordingly

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

# Load the Immune_CSF dataset
lung <- read.csv(file = "~/Lukassen2020_Lung/Counts_lung_cells.csv")
lung <- CreateSeuratObject(counts = lung, project = "Lukassen2020_Lung", min.cells = 3, min.features = 200)

# add barcodes to metadata for filtering
lung@meta.data$barcodes <- rownames(lung@meta.data)
lung@meta.data$barcodes <- gsub('X','',lung@meta.data$barcodes)
lung@meta.data$barcodes <- gsub(pattern = '\\.', '\\-', lung@meta.data$barcodes)
lung@meta.data$orig.ident <- gsub('X','',lung@meta.data$orig.ident)
rownames(lung@meta.data) <- lung@meta.data$barcodes


lung@meta.data$barcodes

# Load the metadata
lung_metadata <- read.csv(file = "~/Lukassen2020_Lung/Metadata_lung_cells.csv", row.names = 1)
lung_metadata$barcodes <- rownames(lung_metadata)

lung_metadata$Cell.type <- gsub('Immuno_', '', lung_metadata$Cell.type)

# clean up barcodes
cell_names <- colnames(lung)
cell_names <- gsub('X','', cell_names)
cell_names <- gsub(pattern = '\\.', '\\-', cell_names)
colnames(lung) <- cell_names
lung <- RenameCells(lung, new.names = cell_names)


# add metadata to object
lung@meta.data <- cbind(lung@meta.data, lung_metadata)


# check if data has been processed or if it is raw
lung@assays$RNA@counts[c(1:1000),c(1:10)] %>% rowSums()

# save raw object with metadata included
saveRDS(lung,
        file = paste0(workingdir,
                      '/',
                      datadir,
                      '/',
                      sample_no,
                      "_raw.rds"))

# add celltype to metadata
lung@meta.data$celltypes <- lung@meta.data$Cell.type

# fix column names
lung@meta.data <- repair_names(lung@meta.data)

# view condition proportions
propsmd <- lung@meta.data %>%
  group_by(Smoking) %>%
  summarise(n = n())

propsct <- lung@meta.data %>%
  group_by(celltypes) %>%
  summarise(n = n())

p_cellquants <- ggplot(propsmd, 
                       aes(reorder(Smoking, -n), n, fill = Smoking)) + 
  geom_bar(stat = 'identity') +
  scale_fill_tableau(palette = 'Tableau 20') +
  xlab('') +
  ylab('Total') + 
  theme_bw(base_size = 14) + 
  theme(axis.text.x = element_blank(),
        axis.text = element_text(color='black'),
        axis.ticks.x.bottom = element_blank())

p_cellquants_log10 <- ggplot(propsct, 
                             aes(reorder(celltypes, -n), n, fill = celltypes)) + 
  geom_bar(stat = 'identity') +
  scale_fill_tableau(palette = 'Tableau 20') + 
  xlab('') +
  ylab('Total') + 
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



# convert to mitochondrial percentage
lung@meta.data$percent_mt <- lung@meta.data$MT.ratio *100

# Visualize QC metrics as a violin plot
feats <- c("nFeature_RNA", "nCount_RNA", "MT.ratio")
VlnPlot(lung, features = feats, pt.size = 0.1, ncol = 3) +
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
p_ctsfts <- FeatureScatter(lung, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_hline(yintercept = 300, linetype='dashed') +
  scale_color_manual(values = qual_col_pals$colors)

cowplot::plot_grid(p_ctsfts)
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

# Normalize the data
lung <- NormalizeData(lung)

# Identify variable features (5K features for use with scANNA)
lung <- FindVariableFeatures(lung, selection.method = "vst", nfeatures = 5000)

# Identify the 20 most highly variable genes
top20 <- head(VariableFeatures(lung), 20)

# plot top20 variable features
vplot <- VariableFeaturePlot(lung)
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


# save object with variable features
saveRDS(lung,
        file = paste0(workingdir,
                      '/',
                      datadir,
                      '/',
                      sample_no,
                      "_qc_hvg.rds"))


# Scale the data
all.genes <- rownames(lung)
lung <- ScaleData(lung, features = all.genes)

# add a bool list of variable genes to the metadata
variable_genes <- VariableFeatures(lung)
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
lung <- RunPCA(lung, features = VariableFeatures(object = lung))

pca_1 <- ElbowPlot(lung, ndims = 50)

pca_2 <- DimPlot(lung, group.by = 'ID')

DimPlot(lung, group.by = 'Smoking')

cowplot::plot_grid(pca_1, pca_2)
ggsave(
  paste0(
    workingdir,
    "/",
    plotdir,
    "/",
    sample_no,
    "_elbow_pca.pdf"
  ),
  height = 6,
  width = 12,
  units = "in",
  dpi = 300
)

# harmonize the dataset and remove batch effects
lung <- RunHarmony(lung, 
                     group.by.vars = c("ID", 'Smoking'), 
                     reduction = "pca", 
                     dims.use = 1:50, 
                     assay.use = "RNA", 
                     plot_convergence = TRUE)


# Clustering
use.pcs <- 1:50

lung <- FindNeighbors(lung, 
                      reduction="harmony",
                      dims = use.pcs, 
                      k.param = 60)

# run the clustering for the multiple resolutions
lung <- FindClusters(object = lung,
                     algorithm = 1,
                     verbose = TRUE,
                     resolution = seq(0.20, 1, 0.20) # change resolution range here
)

# Number of unique clusters generated per resolution
grep("res", colnames(lung@meta.data), value = TRUE) %>%
  purrr::map_chr(~ paste(.x, "--> clusters generated:", length(unique(
    lung@meta.data[, .x]
  ))))


Idents(lung) <- 'RNA_snn_res.0.2'
table(lung@active.ident)


# Plotting UMAP -----------------------------------------------------------

lung <- RunUMAP(lung, 
                  reduction = "harmony", 
                  dims = use.pcs,
                  umap.method = 'umap-learn')

p_conditions <- DimPlot(lung,
                        reduction = "umap",
                        label = TRUE,
                        repel = TRUE,
                        #label.size = 8, 
                        group.by = 'ID') +
  scale_color_tableau(palette = 'Tableau 20') +
  ggtitle('ID') + 
  theme(plot.title = element_text(hjust = 0.5)) +
  NoLegend()

p_clusters <- DimPlot(lung,
                      reduction = "umap",
                      label = TRUE,
                      repel = TRUE,
                      #label.size = 8, 
                      group.by = 'cluster') +
  scale_color_tableau(palette = 'Tableau 20') +
  ggtitle('Seurat_Clusters') + 
  theme(plot.title = element_text(hjust = 0.5)) +
  NoLegend()

p_celltypes <- DimPlot(lung,
                        reduction = "umap",
                        label = TRUE,
                        repel = TRUE,
                        #label.size = 8, 
                        group.by = 'celltypes') +
  scale_color_tableau(palette = 'Tableau 20') +
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

# UMAP with celltypes
DimPlot(lung,
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


# save the object with the annotations
saveRDS(lung,
        file = paste0(workingdir,
                      '/',
                      datadir,
                      '/',
                      sample_no,
                      "_qc_hvg_anno_full.rds"))


# H5AD Object Creation ----------------------------------------------------

# saving the integrated file for analysis in scanpy
library(SeuratDisk)

# subset (5k variable genes)
variable_genes <- VariableFeatures(lung)

lung_5k <- subset(lung, features = variable_genes)

# remove scale.data so that raw counts gets saved to raw
lung_5k <- DietSeurat(lung_5k, 
                        counts = TRUE, 
                        data = TRUE, 
                        scale.data = FALSE,
                        dimreducs = c('harmony', 'pca', 'umap'), 
                        graphs = c('RNA_nn', 'RNA_snn'))

# save the object with the annotations
saveRDS(lung_5k,
        file = paste0(workingdir,
                      '/',
                      datadir,
                      '/',
                      sample_no,
                      "_qc_hvg_anno_5k.rds"))

SaveH5Seurat(lung_5k, 
             filename = "~/Lukassen2020_Lung/Lukassen2020_Lung_qc_hvg_anno_5k.h5Seurat")

Convert("~/Lukassen2020_Lung/Lukassen2020_Lung_qc_hvg_anno_5k.h5Seurat", dest = "h5ad")


# Differential Expression -------------------------------------------------

DEWilcox <- FindAllMarkers(lung, log2FC.threshold = 0.25, test.use = "wilcox",
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



DotPlot(lung, 
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


