# This script is for processing the GSE163005 dataset
# for the scANNA manuscript.

# # Download data ------------------------------------------------------------
# # download data from GEO
# curl "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE163nnn/GSE163005/suppl/GSE163005%5Fbarcodes%2Etsv%2Egz" -o "barcodes.tsv.gz"
# curl "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE163nnn/GSE163005/suppl/GSE163005%5Ffeatures%2Etsv%2Egz" -o "genes.tsv.gz"
# curl "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE163nnn/GSE163005/suppl/GSE163005%5Fmatrix%2Emtx%2Egz" -o "matrix.mtx.gz"

# # download metadata from GEO
# curl "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE163nnn/GSE163005/suppl/GSE163005%5Fannotation%5Fcluster%2Ecsv%2Egz" -o "GSE163005_annotation_cluster.csv.gz"
# curl "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE163nnn/GSE163005/suppl/GSE163005%5Fannotation%5Fdx%2Ecsv%2Egz" -o "GSE163005_annotation_dx.csv.gz"
# curl "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE163nnn/GSE163005/suppl/GSE163005%5Fannotation%5Fpatients%2Ecsv%2Egz" -o "GSE163005_annotation_patients.csv.gz"

# # make directory for data
# mkdir GSE163005_IMMUNE_CSF

# # move data to directory
# mv barcodes.tsv.gz GSE163005_IMMUNE_CSF/barcodes.tsv.gz
# mv genes.tsv.gz GSE163005_IMMUNE_CSF/genes.tsv.gz
# mv matrix.mtx.gz GSE163005_IMMUNE_CSF/matrix.mtx.gz

# # move metadata to directory
# mv GSE163005_annotation_cluster.csv.gz GSE163005_IMMUNE_CSF/GSE163005_annotation_cluster.csv.gz
# mv GSE163005_annotation_dx.csv.gz GSE163005_IMMUNE_CSF/GSE163005_annotation_dx.csv.gz
# mv GSE163005_annotation_patients.csv.gz GSE163005_IMMUNE_CSF/GSE163005_annotation_patients.csv.gz


# Loading packages and setting up directories -----------------------------
library(tidyverse)
library(ggthemes)
library(Seurat)
library(harmony)


sample_no <- 'GSE163005_IMMUNE_CSF' # change accordingly

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
immune_csf <- Read10X(data.dir = "~/GSE163005_IMMUNE_CSF/data/")

immune_csf <- CreateSeuratObject(counts = immune_csf, project = "immune_csf", min.cells = 3, min.features = 200)
immune_csf

# add barcodes to metadata for filtering
immune_csf@meta.data$barcodes <- rownames(immune_csf@meta.data)

# read in all metadata
clusters <- read_csv(file = '~/GSE163005_IMMUNE_CSF/metadata/GSE163005_annotation_cluster.csv')
colnames(clusters) <- c('barcodes','celltype')
celltypes <- clusters$celltype

# subset the data to contain same barcodes as in the metadata
immune_csf <- subset(immune_csf, barcodes %in% clusters$barcodes)

# add metadata
immune_csf <- AddMetaData(immune_csf, celltypes, col.name = 'celltypes')

# read in condition metadata
conditions <- read_csv(file = '~/GSE163005_IMMUNE_CSF/metadata/GSE163005_annotation_dx.csv')
colnames(conditions) <- c('barcodes', 'condition')
conditions <- conditions$condition

# add metadata
immune_csf <- AddMetaData(immune_csf, conditions, col.name = 'condition')

# read in patient metadata
patients <- read_csv(file = '~/GSE163005_IMMUNE_CSF/metadata/GSE163005_annotation_patients.csv')
colnames(patients) <- c('barcodes','patient')
patients <- patients$patient

# add metadata
immune_csf <- AddMetaData(immune_csf, patients, col.name = 'patient')

immune_csf@meta.data$celltypes %>% table()
immune_csf@meta.data$condition %>% table()
immune_csf@meta.data$patient %>% table()

# check if data has been processed or if it is raw
immune_csf@assays$RNA@counts[c(1:1000),c(1:10)] %>% colSums()

# check NA quantity
sum(is.na(immune_csf@meta.data$celltypes))

# subset out matDC since less than 100 cells per celltype
immune_csf <- subset(immune_csf, subset = celltypes != 'matDC')


# save raw object with metadata included
saveRDS(immune_csf,
        file = paste0(workingdir,
                      '/',
                      datadir,
                      '/',
                      sample_no,
                      "_raw.rds"))

# view cell type proportions
propsmd <- immune_csf@meta.data %>%
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
immune_csf <- PercentageFeatureSet(immune_csf, pattern = "^MT-", col.name = "percent_mito")

# Visualize QC metrics as a violin plot
feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito")
VlnPlot(immune_csf, group.by = "condition", features = feats, pt.size = 0.1, ncol = 3) + 
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
p_ctsmt <- FeatureScatter(immune_csf, feature1 = "nCount_RNA", feature2 = "percent_mito") + geom_hline(yintercept = 5, linetype='dashed') + 
  scale_color_manual(values = qual_col_pals$colors)
p_ctsfts <- FeatureScatter(immune_csf, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_hline(yintercept = 5000, linetype='dashed') + 
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

# Subset data based on QC metrics
immune_csf <- subset(immune_csf, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent_mito < 5) # could change to 20 percent instead

# Normalize the data
immune_csf <- NormalizeData(immune_csf)

# Identify variable features (5K features for use with scANNA)
immune_csf <- FindVariableFeatures(immune_csf, selection.method = "vst", nfeatures = 5000)

# Identify the 20 most highly variable genes
top20 <- head(VariableFeatures(immune_csf), 20)

# plot top20 variable features
vplot <- VariableFeaturePlot(immune_csf)
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
saveRDS(immune_csf,
        file = paste0(workingdir,
                      '/',
                      datadir,
                      '/',
                      sample_no,
                      "_qc_hvg.rds"))

# Scale data
all.genes <- rownames(immune_csf)
immune_csf <- ScaleData(immune_csf, features = all.genes)

# add a bool list of variable genes to the metadata
variable_genes <- VariableFeatures(immune_csf)
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
immune_csf <- RunPCA(immune_csf, features = VariableFeatures(object = immune_csf))

ElbowPlot(immune_csf, ndims = 50)

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
immune_csf <- RunHarmony(immune_csf, group.by.vars = c("condition", "patient"), reduction = "pca", dims.use = 1:50, assay.use = "RNA", plot_convergence = TRUE)

# Clustering
use.pcs <- 1:50

immune_csf <- FindNeighbors(immune_csf, 
                            reduction="harmony",
                            dims = use.pcs)

# run the clustering for the multiple resolutions
immune_csf <- FindClusters(object = immune_csf, 
                           resolution = seq(0.20, 1, 0.20) # change resolution range here
)

# Number of unique clusters generated per resolution
grep("res", colnames(immune_csf@meta.data), value = TRUE) %>%
  purrr::map_chr(~ paste(.x, "--> clusters generated:", length(unique(
    immune_csf@meta.data[, .x]
  ))))

Idents(immune_csf) <- 'RNA_snn_res.0.2'
table(immune_csf@active.ident)


# Plotting UMAP -----------------------------------------------------------

Idents(immune_csf) <- immune_csf@meta.data$celltypes

immune_csf <- RunUMAP(immune_csf, 
                      reduction = "harmony", 
                      dims = use.pcs,
                      umap.method = 'umap-learn')

p_conditions <- DimPlot(immune_csf,
                        reduction = "umap",
                        label = TRUE,
                        repel = TRUE,
                        group.by = 'condition') +
  scale_color_manual(values = qual_col_pals$colors) +
  ggtitle('Conditions') + 
  theme(plot.title = element_text(hjust = 0.5)) +
  NoLegend()

p_clusters <- DimPlot(immune_csf,
                      reduction = "umap",
                      label = TRUE,
                      repel = TRUE,
                      group.by = 'RNA_snn_res.0.2') +
  scale_color_manual(values = qual_col_pals$colors) +
  ggtitle('Seurat_Clusters') + 
  theme(plot.title = element_text(hjust = 0.5)) +
  NoLegend()

p_celltypes <- DimPlot(immune_csf,
                       reduction = "umap",
                       label = TRUE,
                       repel = TRUE,
                       group.by = 'general_celltypes') +
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

# add a more general annotation called 'general_celltypes'
celltypes_general <- immune_csf@meta.data$celltypes
celltypes_general <- str_replace_all(celltypes_general, 'mono[0-9]', 'monocytes') # make all monocyte populations into one
celltypes_general <- str_replace_all(celltypes_general, 'granulo[0-9]', 'granulocytes') # make all granulocyte populations into one (Granulocytes: neutrophil, basophil, or eosinophil)
celltypes_general <- str_replace_all(celltypes_general, 'mDC[0-9]', 'mDC') # make all mature DC (mDC) populations into one
celltypes_general <- str_replace_all(celltypes_general, 'naiveBc', 'B cells') # change naiveBc to B cells
immune_csf <- AddMetaData(immune_csf, celltypes_general, col.name = 'general_celltypes')


# UMAP with celltypes
DimPlot(immune_csf,
        reduction = "umap",
        label = TRUE,
        repel = TRUE,
        label.size = 8, 
        group.by = 'general_celltypes') + 
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
saveRDS(immune_csf,
        file = paste0(workingdir,
                      '/',
                      datadir,
                      '/',
                      sample_no,
                      "_qc_hvg_anno_full.rds"))

# H5AD Object Creation ----------------------------------------------------

# saving the integrated file for analysis in scanpy & scANNA
library(SeuratDisk)

# remove scale.data so that raw counts gets saved to raw slot and not normalized counts
immune_csf <- DietSeurat(immune_csf, 
                         counts = TRUE, 
                         data = TRUE, 
                         scale.data = FALSE,
                         dimreducs = c('harmony', 'pca', 'umap'), 
                         graphs = c('RNA_nn', 'RNA_snn'))

# subset (5k variable genes)
variable_genes <- VariableFeatures(immune_csf)

immune_csf_5k <- subset(immune_csf, features = variable_genes)

# save the object with the annotations
saveRDS(immune_csf_5k,
        file = paste0(workingdir,
                      '/',
                      datadir,
                      '/',
                      sample_no,
                      "_qc_hvg_anno_5k.rds"))


SaveH5Seurat(immune_csf_5k, filename = "~/GSE163005_IMMUNE_CSF/GSE_ImmuneCSF_qc_hvg_anno_5k.h5Seurat")
Convert("~/GSE163005_IMMUNE_CSF/GSE_ImmuneCSF_qc_hvg_anno_5k.h5Seurat", dest = "h5ad")


# Differential Expression -------------------------------------------------

DEWilcox <- FindAllMarkers(immune_csf, log2FC.threshold = 0.25, test.use = "wilcox",
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



DotPlot(immune_csf, 
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



