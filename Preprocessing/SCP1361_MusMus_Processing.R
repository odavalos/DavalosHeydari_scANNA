# This script is used to process the SCP1361 dataset
# for the scANNA manuscript.

# to download the data use the following link
# https://singlecell.broadinstitute.org/single_cell/study/SCP1361/single-cell-transcriptome-analysis-reveals-cellular-heterogeneity-in-the-ascending-aorta-of-normal-and-high-fat-diet-mice



# Loading packages and setting up directories -----------------------------
library(tidyverse)
library(ggthemes)
library(Seurat)
library(harmony)

sample_no <- 'scp1361' # change accordingly

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


# Load the mus_hfd dataset
mus <- Read10X(data.dir = "/Volumes/BigDrive/GradSchool/NACT_NewDatasets/SCP1361/expression/data/")

# Initialize the Seurat object with the raw (non-normalized data).
mus_hfd <- CreateSeuratObject(counts = mus, project = "scp1361", min.cells = 3, min.features = 200)


# check if data has been processed or if it is raw
mus_hfd@assays$RNA@counts[c(1:10000),c(1:10000)] %>% colSums()
mus_hfd@assays$RNA@counts[c(1:10000),c(1:10000)] %>% rowSums()


# Add metadata ------------------------------------------------------------

# metadata
metadata <- read_tsv(file = '/Volumes/BigDrive/GradSchool/NACT_NewDatasets/SCP1361/metadata/meta_data.txt')
metadata <- metadata[-1,] %>%
  rename(barcodes = NAME)

# clusters
clusters <- read_tsv(file = '/Volumes/BigDrive/GradSchool/NACT_NewDatasets/SCP1361/cluster/ordinations.txt')
clusters <- clusters[-1,] %>%
  select(NAME, Category) %>%
  rename(barcodes = NAME, celltypes = Category)

# extract barcodes from seurat dataset
mus_hfd@meta.data$barcodes <- rownames(mus_hfd@meta.data)

# merge metadata to seurat object
mus_hfd@meta.data <- left_join(mus_hfd@meta.data, metadata, by='barcodes') # append metadata
mus_hfd@meta.data <- left_join(mus_hfd@meta.data, clusters, by='barcodes') # append clusters
rownames(mus_hfd@meta.data) <- mus_hfd@meta.data$barcodes # re-append barcodes as rownames 

# calculate percentages of mito and ribo genes
mus_hfd <- PercentageFeatureSet(mus_hfd, pattern = "^mt-", col.name = "percent_mito")

# save raw object with metadata included
saveRDS(mus_hfd,
        file = paste0(workingdir,
                      '/',
                      datadir,
                      '/',
                      sample_no,
                      "_raw.rds"))


# Visualize QC metrics as a violin plot
feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo")
VlnPlot(mus_hfd, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 3) + 
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
p_ctsmt <- FeatureScatter(mus_hfd, feature1 = "nCount_RNA", feature2 = "percent.mt") #+ 
p_ctsfts <- FeatureScatter(mus_hfd, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") #+
cowplot::plot_grid(p_ctsfts, p_ctsmt, p_ctsribo)
ggsave(
  paste0(
    workingdir,
    "/",
    plotdir,
    "/",
    sample_no,
    "_qc_scatter.pdf"
  ),
  height = 8,
  width = 12,
  units = "in",
  dpi = 300
)


# Normalize data
mus_hfd <- NormalizeData(mus_hfd)

mus_hfd <- FindVariableFeatures(mus_hfd, selection.method = "vst", nfeatures = 5000)

# Identify the 20 most highly variable genes
top20 <- head(VariableFeatures(mus_hfd), 20)

# plot top20 variable features
vplot <- VariableFeaturePlot(mus_hfd)
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
saveRDS(mus_hfd,
        file = paste0(workingdir,
                      '/',
                      datadir,
                      '/',
                      sample_no,
                      "_qc_hvg.rds"))


# Scale data
all.genes <- rownames(mus_hfd)
mus_hfd <- ScaleData(mus_hfd, features = all.genes)

# add a bool list of variable genes to the metadata
variable_genes <- VariableFeatures(mus_hfd)
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
mus_hfd <- RunPCA(mus_hfd, features = VariableFeatures(object = mus_hfd))

ElbowPlot(mus_hfd, ndims = 50)
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

mus_hfd <- RunHarmony(mus_hfd, group.by.vars = "disease__ontology_label", reduction = "pca",
                              dims.use = 1:50, assay.use = "RNA")

# Clustering
use.pcs <- 1:50

mus_hfd <- FindNeighbors(mus_hfd, 
                         reduction="pca",
                         dims = 1:50,
                         k.param = 60)

# run the clustering for the multiple resolutions
mus_hfd <- FindClusters(object = mus_hfd, 
                          resolution = seq(0.20, 1, 0.20) # change resolution range here
)


# Number of unique clusters generated per resolution
grep("res", colnames(mus_hfd@meta.data), value = TRUE) %>%
  purrr::map_chr(~ paste(.x, "--> clusters generated:", length(unique(
    mus_hfd@meta.data[, .x]
  ))))

Idents(mus_hfd) <- 'RNA_snn_res.0.2'
table(mus_hfd@active.ident)


# Plotting UMAP -----------------------------------------------------------

# need to make a general annotation for celltypes --> 'celltypes' 
# more refined which is current celltypes --> celltypes_2
mus_hfd@meta.data$celltypes_2 <- mus_hfd@meta.data$celltypes
mus_hfd@meta.data$celltypes <- mus_hfd@meta.data$cell_type__ontology_label
mus_hfd@meta.data$clusters <- mus_hfd@meta.data$RNA_snn_res.0.2

Idents(mus_hfd) <- mus_hfd@meta.data$celltypes

mus_hfd <- RunUMAP(mus_hfd, 
                     reduction = "harmony", 
                     dims = use.pcs,
                     umap.method = 'umap-learn')

p_conditions <- DimPlot(mus_hfd,
                        reduction = "umap",
                        label = TRUE,
                        repel = TRUE,
                        group.by = 'disease__ontology_label') +
  ggthemes::scale_colour_tableau(palette = 'Tableau 20', direction = -1) +
  ggtitle('Conditions') + 
  theme(plot.title = element_text(hjust = 0.5)) +
  NoLegend()

p_clusters <- DimPlot(mus_hfd,
                      reduction = "umap",
                      label = TRUE,
                      repel = TRUE,
                      group.by = 'RNA_snn_res.0.2') +
  ggthemes::scale_colour_tableau(palette = 'Tableau 20', direction = -1) +
  ggtitle('Seurat_Clusters') + 
  theme(plot.title = element_text(hjust = 0.5)) +
  NoLegend()

p_celltypes <- DimPlot(mus_hfd,
                       reduction = "umap",
                       label = TRUE,
                       repel = TRUE, 
                       group.by = 'celltypes_2') +
  ggthemes::scale_colour_tableau(palette = 'Tableau 20', direction = -1) +
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
DimPlot(mus_hfd,
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
saveRDS(mus_hfd,
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
variable_genes <- VariableFeatures(mus_hfd)

mus_hfd_5k <- subset(mus_hfd, features = variable_genes)

# remove scale.data so that raw counts gets saved to raw
mus_hfd_5k <- DietSeurat(mus_hfd_5k, 
                      counts = TRUE, 
                      data = TRUE, 
                      scale.data = FALSE,
                      dimreducs = c('harmony', 'pca', 'umap'), 
                      graphs = c('RNA_nn', 'RNA_snn'))

# save the object with the annotations
saveRDS(mus_hfd_5k,
        file = paste0(workingdir,
                      '/',
                      datadir,
                      '/',
                      sample_no,
                      "_qc_hvg_anno_5k.rds"))


SaveH5Seurat(mus_hfd_5k, 
             filename = "~/SCP1361/scp1361_qc_hvg_anno_5k.h5Seurat")

Convert("~/SCP1361/scp1361_qc_hvg_anno_5k.h5Seurat", dest = "h5ad")



# Differential Expression -------------------------------------------------

DEWilcox <- FindAllMarkers(mus_hfd, log2FC.threshold = 0.25, test.use = "wilcox",
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



DotPlot(mus_hfd, 
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



