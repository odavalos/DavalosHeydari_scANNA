# This script is used to generate the correlation plots for the attention 
# sums in the scANNA manuscript.

library(tidyverse)
library(ggthemes)


# Misc Functions ----------------------------------------------------------

# Adapted from this repo with ggplot tricks https://github.com/teunbrand/ggplot_tricks
contrast <- function(colour) {
  out   <- rep("black", length(colour))
  light <- farver::get_channel(colour, "l", space = "hcl")
  out[light < 50] <- "white"
  out
}

autocontrast <- aes(colour = after_scale(contrast(fill)))

# Functions for generating lower and upper triangle of correlation matrix
# plot code adapted from http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization

# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}



# GSE154567 Correlation Plots ---------------------------------------------

sample_no <- 'GSE154567' # change accordingly

# directories
workingdir <- paste0(getwd(), '/',sample_no)
plotdir <- 'plots'
datadir <- 'data'
resultsdir <- 'results'

# check for plot directory
ifelse(!dir.exists(file.path(workingdir, plotdir)), 
       dir.create(file.path(workingdir, plotdir)), FALSE)

# check for data directory
ifelse(!dir.exists(file.path(workingdir, datadir)), 
       dir.create(file.path(workingdir, datadir)), FALSE)

# check for results directory
ifelse(!dir.exists(file.path(workingdir, resultsdir)), 
       dir.create(file.path(workingdir, resultsdir)), FALSE)



# reading in the metadata
metadata <- read_csv(file = '~/Documents/NACT_ComparisonModels/DatasetSplits/Metadata_Splits/GSE154567_metadata_splits.csv')

# getting the celltypes
celltype_df <- metadata %>%
  select(celltypes, cluster) %>% 
  unique()

celltype_df$cluster_no <- paste0('Cluster_', celltype_df$cluster)


# loading in the attention data
attention <- read_csv(file = 'GSE154567/results/ICML_Files/GSE154567_attentionsums.csv')
attention <- attention %>%
  rename(genes = `...1`)

# rename the columns to contain celltypes
idx <- celltype_df$cluster_no %in% names(attention)

# rename the columns
attention <- attention %>%
  rename_with(~celltype_df$celltypes[idx], celltype_df$cluster_no[idx])

attention <- as.data.frame(attention)


# correlation matrix
corr <- round(cor(attention[,2:10], method = 'pearson'), 2)


# get upper triangle of correlation matrix
upper_tri <- get_upper_tri(corr)
upper_tri

# melt the correlation matrix
melted_cormat <- reshape2::melt(upper_tri, na.rm = TRUE)


# Heatmap
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(Var2, Var1, label = value, !!!autocontrast),
            size = 14) +
  scale_fill_viridis_c(option = 'viridis', 
                       direction = 1,
                       guide = guide_colorbar( 
                                              barwidth = 20,
                                              barheight = 4)) +
  labs(x = "",
       y = "",
       fill = "Pearson\nCorrelation") + 
  theme_minimal() + 
  theme(axis.text = element_text(color = 'black', 
                                 face = 'bold'), 
        axis.text.x = element_text(angle = 35, 
                                   vjust = 1,
                                   size = 18,
                                   hjust = 1),
        axis.text.y = element_text(size = 18),
        title = element_text(size = 14),
        legend.position = 'top',
        legend.text = element_text(size = 14, 
                                   #face = 'bold',
                                   color = 'black')) +
  coord_fixed()

# save plot
ggsave(filename = paste0(workingdir,'/plots/',sample_no,"attention_sums_correlations_supfig.pdf"), 
       units = 'in',
       width = 16,
       height = 14, 
       dpi = 300)




# GSE163005_IMMUNE_CSF Correlation Plots -----------------------------------

sample_no <- 'GSE163005_IMMUNE_CSF' # change accordingly

# directories
workingdir <- paste0(getwd(), '/',sample_no)
plotdir <- 'plots'
datadir <- 'data'
resultsdir <- 'results'

# check for plot directory
ifelse(!dir.exists(file.path(workingdir, plotdir)), 
       dir.create(file.path(workingdir, plotdir)), FALSE)

# check for data directory
ifelse(!dir.exists(file.path(workingdir, datadir)), 
       dir.create(file.path(workingdir, datadir)), FALSE)

# check for results directory
ifelse(!dir.exists(file.path(workingdir, resultsdir)), 
       dir.create(file.path(workingdir, resultsdir)), FALSE)


# reading in the metadata
metadata <- read_csv(file = '~/Documents/NACT_ComparisonModels/DatasetSplits/Metadata_Splits/GSE_ImmuneCSF_metadata_splits.csv')

# getting the celltypes
celltype_df <- metadata %>%
  select(general_celltypes, cluster) %>% 
  unique()

celltype_df$cluster_no <- paste0('Cluster_', celltype_df$cluster)


# loading in the attention data
attention <- read_csv(file = 'GSE163005_IMMUNE_CSF/results/GSE163005_IMMUNECSF_attentionsums_new.csv')
attention <- attention %>%
  rename(genes = `...1`)

# rename the columns to contain celltypes
idx <- celltype_df$cluster_no %in% names(attention)

# rename the columns
attention <- attention %>%
  rename_with(~celltype_df$general_celltypes[idx], celltype_df$cluster_no[idx])

attention <- as.data.frame(attention)
col_no <- ncol(attention)

# correlation matrix
corr <- round(cor(attention[,2:col_no], method = 'pearson'), 2)

# Get upper triangle of the correlation matrix
upper_tri <- get_upper_tri(corr)
upper_tri

# Melt the correlation matrix
melted_cormat <- reshape2::melt(upper_tri, na.rm = TRUE)

# Heatmap
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(Var2, Var1, label = value, !!!autocontrast),
            size = 12) +
  scale_fill_viridis_c(option = 'viridis', 
                       direction = 1,
                       guide = guide_colorbar( 
                         barwidth = 20,
                         barheight = 4)) +
  labs(x = "",
       y = "",
       fill = "Pearson\nCorrelation") + 
  theme_minimal() + 
  theme(axis.text = element_text(color = 'black', 
                                 face = 'bold'), 
        axis.text.x = element_text(angle = 35, 
                                   vjust = 1,
                                   size = 18,
                                   hjust = 1),
        axis.text.y = element_text(size = 18),
        title = element_text(size = 14),
        legend.position = 'top',
        legend.text = element_text(size = 14, 
                                   color = 'black')) +
  coord_fixed()

# save plot
ggsave(filename = paste0(workingdir,'/plots/',sample_no,"attention_sums_correlations_supfig.pdf"), 
       units = 'in',
       width = 16,
       height = 14, 
       dpi = 300)




# GSE144236 Correlation Plots ----------------------------------------------

sample_no <- 'GSE144236' # change accordingly

# directories
workingdir <- paste0(getwd(), '/',sample_no)
plotdir <- 'plots'
datadir <- 'data'
resultsdir <- 'results'

# check for plot directory
ifelse(!dir.exists(file.path(workingdir, plotdir)), 
       dir.create(file.path(workingdir, plotdir)), FALSE)

# check for data directory
ifelse(!dir.exists(file.path(workingdir, datadir)), 
       dir.create(file.path(workingdir, datadir)), FALSE)

# check for results directory
ifelse(!dir.exists(file.path(workingdir, resultsdir)), 
       dir.create(file.path(workingdir, resultsdir)), FALSE)


# reading in the metadata
metadata <- read_csv(file = '~/Documents/NACT_ComparisonModels/DatasetSplits/Metadata_Splits/GSE144236_metadata_splits.csv')

# getting the celltypes
celltype_df <- metadata %>%
  select(celltypes, cluster) %>% 
  unique()

celltype_df$cluster_no <- paste0('Cluster_', celltype_df$cluster)


# loading in the attention data
attention <- read_csv(file = 'GSE144236/results/GSE144236_attentionsums_new.csv')
attention <- attention %>%
  rename(genes = `...1`)

# rename the columns to contain celltypes
idx <- celltype_df$cluster_no %in% names(attention)

# rename the columns
attention <- attention %>%
  rename_with(~celltype_df$celltypes[idx], celltype_df$cluster_no[idx])

attention <- as.data.frame(attention)
col_no <- ncol(attention)

# correlation matrix
corr <- round(cor(attention[,2:col_no], method = 'pearson'), 2)

# Get upper triangle of the correlation matrix
upper_tri <- get_upper_tri(corr)
upper_tri

# Melt the correlation matrix
melted_cormat <- reshape2::melt(upper_tri, na.rm = TRUE)

# Heatmap
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(Var2, Var1, label = value, !!!autocontrast),
            size = 10) +
  scale_fill_viridis_c(option = 'viridis', 
                       direction = 1,
                       guide = guide_colorbar( 
                         barwidth = 20,
                         barheight = 4)) +
  labs(x = "",
       y = "",
       fill = "Pearson\nCorrelation") + 
  theme_minimal() + 
  theme(axis.text = element_text(color = 'black', 
                                 face = 'bold'), 
        axis.text.x = element_text(angle = 35, 
                                   vjust = 1,
                                   size = 18,
                                   hjust = 1),
        axis.text.y = element_text(size = 18),
        title = element_text(size = 14),
        legend.position = 'top',
        legend.text = element_text(size = 14, 
                                   color = 'black')) +
  coord_fixed()

# save plot
ggsave(filename = paste0(workingdir,'/plots/',sample_no,"attention_sums_correlations_supfig.pdf"), 
       units = 'in',
       width = 17,
       height = 15, 
       dpi = 300)



# Lukassen2020_Lung Correlation Plots ----------------------------------------------

sample_no <- 'Lukassen2020_Lung' # change accordingly

# directories
workingdir <- paste0(getwd(), '/',sample_no)
plotdir <- 'plots'
datadir <- 'data'
resultsdir <- 'results'

# check for plot directory
ifelse(!dir.exists(file.path(workingdir, plotdir)), 
       dir.create(file.path(workingdir, plotdir)), FALSE)

# check for data directory
ifelse(!dir.exists(file.path(workingdir, datadir)), 
       dir.create(file.path(workingdir, datadir)), FALSE)

# check for results directory
ifelse(!dir.exists(file.path(workingdir, resultsdir)), 
       dir.create(file.path(workingdir, resultsdir)), FALSE)


# reading in the metadata
metadata <- read_csv(file = '~/Documents/NACT_ComparisonModels/DatasetSplits/Metadata_Splits/Lukassen2020_Lung_metadata_splits.csv')

# getting the celltypes
celltype_df <- metadata %>%
  select(celltypes, cluster) %>% 
  unique()

celltype_df$cluster_no <- paste0('Cluster_', celltype_df$cluster)


# loading in the attention data
attention <- read_csv(file = 'Lukassen2020_Lung/results/Lukassen2020_Lung_attentionsums_new.csv')
attention <- attention %>%
  rename(genes = `...1`)

# rename the columns to contain celltypes
idx <- celltype_df$cluster_no %in% names(attention)

# rename the columns
attention <- attention %>%
  rename_with(~celltype_df$celltypes[idx], celltype_df$cluster_no[idx])

attention <- as.data.frame(attention)
col_no <- ncol(attention)

# correlation matrix
corr <- round(cor(attention[,2:col_no], method = 'pearson'), 2)

# Get upper triangle of the correlation matrix
upper_tri <- get_upper_tri(corr)
upper_tri

# Melt the correlation matrix
melted_cormat <- reshape2::melt(upper_tri, na.rm = TRUE)

# Heatmap
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(Var2, Var1, label = value, !!!autocontrast),
            size = 12) +
  scale_fill_viridis_c(option = 'viridis', 
                       direction = 1,
                       guide = guide_colorbar( 
                         barwidth = 20,
                         barheight = 4)) +
  labs(x = "",
       y = "",
       fill = "Pearson\nCorrelation") + 
  theme_minimal() + 
  theme(axis.text = element_text(color = 'black', 
                                 face = 'bold'), 
        axis.text.x = element_text(angle = 35, 
                                   vjust = 1,
                                   size = 18,
                                   hjust = 1),
        axis.text.y = element_text(size = 18),
        title = element_text(size = 14),
        legend.position = 'top',
        legend.text = element_text(size = 14, 
                                   color = 'black')) +
  coord_fixed()

# save plot
ggsave(filename = paste0(workingdir,'/plots/',sample_no,"attention_sums_correlations_supfig.pdf"), 
       units = 'in',
       width = 16,
       height = 14, 
       dpi = 300)



# SCP1361 Correlation Plots ----------------------------------------------

sample_no <- 'scp1361' # change accordingly

# directories
workingdir <- paste0(getwd(), '/',sample_no)
plotdir <- 'plots'
datadir <- 'data'
resultsdir <- 'results'

# check for plot directory
ifelse(!dir.exists(file.path(workingdir, plotdir)), 
       dir.create(file.path(workingdir, plotdir)), FALSE)

# check for data directory
ifelse(!dir.exists(file.path(workingdir, datadir)), 
       dir.create(file.path(workingdir, datadir)), FALSE)

# check for results directory
ifelse(!dir.exists(file.path(workingdir, resultsdir)), 
       dir.create(file.path(workingdir, resultsdir)), FALSE)


# reading in the metadata
metadata <- read_csv(file = '~/Documents/NACT_ComparisonModels/DatasetSplits/Metadata_Splits/SCP1361_metadata_splits.csv')

# getting the celltypes
celltype_df <- metadata %>%
  select(celltypes, cluster) %>% 
  unique()

celltype_df$cluster_no <- paste0('Cluster_', celltype_df$cluster)



# loading in the attention data
attention <- read_csv(file = 'scp1361/results/scp1361_attentionsums_new.csv')
attention <- attention %>%
  rename(genes = `...1`)

# rename the columns to contain celltypes
idx <- celltype_df$cluster_no %in% names(attention)

# rename the columns
attention <- attention %>%
  rename_with(~celltype_df$celltypes[idx], celltype_df$cluster_no[idx])

attention <- as.data.frame(attention)
col_no <- ncol(attention)

# correlation matrix
corr <- round(cor(attention[,2:col_no], method = 'pearson'), 2)

# Get upper triangle of the correlation matrix
upper_tri <- get_upper_tri(corr)
upper_tri

# Melt the correlation matrix
melted_cormat <- reshape2::melt(upper_tri, na.rm = TRUE)

# Heatmap
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(Var2, Var1, label = value, !!!autocontrast),
            size = 12) +
  scale_fill_viridis_c(option = 'viridis', 
                       direction = 1,
                       guide = guide_colorbar( 
                         barwidth = 20,
                         barheight = 4)) +
  labs(x = "",
       y = "",
       fill = "Pearson\nCorrelation") + 
  theme_minimal() + 
  theme(axis.text = element_text(color = 'black', 
                                 face = 'bold'), 
        axis.text.x = element_text(angle = 35, 
                                   vjust = 1,
                                   size = 18,
                                   hjust = 1),
        axis.text.y = element_text(size = 18),
        title = element_text(size = 14),
        legend.position = 'top',
        legend.text = element_text(size = 14, 
                                   color = 'black')) +
  coord_fixed()

# save plot
ggsave(filename = paste0(workingdir,'/plots/',sample_no,"attention_sums_correlations_supfig.pdf"), 
       units = 'in',
       width = 16,
       height = 14, 
       dpi = 300)
