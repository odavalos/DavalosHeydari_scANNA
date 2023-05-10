# This script is used to generate the plots for the feature selection
# evaluation of the scANNA manuscript. 


library(tidyverse)
library(rcartocolor)
library(patchwork)


# Directory Setup ---------------------------------------------------------

sample_no <- 'ModelEvaluation' # change accordingly

# directories
workingdir <- paste0(getwd(), '/',sample_no)
plotdir <- 'plots'
datadir <- 'data'
resultsdir <- 'results'

# check for sampledir
# check for plot directory
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



# Global Feature Selection ------------------------------------------------

# reading in the csv
featureselection_eval <- read_csv(file = paste0(
  file.path(file.path(workingdir, datadir), '/Model Evaluation - Supervised_GlobalGeneMarkerSelection_Final.csv')))

# little bit of prepro for plotting
featureselection_eval <- featureselection_eval %>%
  mutate(Dataset = case_when(Dataset == 'ImmuneCSF_PBMC' ~ 'ImmuneCSF',
                             Dataset == 'GSE154567_COVID PBMC' ~ 'GSE154567',
                             Dataset == 'GSE144236_Immune cSCC' ~ 'GSE144236',
                             Dataset == 'Lukassen2020_Lung' ~ 'Lukassen2020',
                             Dataset == 'Mouse HDF (SCP1361)' ~ 'SCP1361'))

# arranging for ploting
featureselection_eval$Feature_Selection <- factor(featureselection_eval$Feature_Selection, 
                                                  levels = c('scANNA', 
                                                             'SMaSH', 
                                                             'scGeneFit',
                                                             'Triku',
                                                             'Seurat', 
                                                             'CellRanger'))

# separate dataframe for plotting baselines
featureselection_eval_5k <- featureselection_eval %>%
  filter(Num_Markers == 5000) %>%
  group_by(Dataset) %>%
  summarise(baseline = mean(Weighted_F1))

# modify Rcolorbrewer's 'Dark2' color palette to emphasize scANNA
darkpal <- RColorBrewer::brewer.pal(n=8,'Dark2')
darkpal_v2 <- replace(darkpal, c(1,2,3,4), darkpal[c(4,2,1,3)])


# horizontal plots - grid
featureselection_eval %>%
  filter(Num_Markers != 5000) %>%
  ggplot(aes(Num_Markers, 
             Weighted_F1, 
             color = Feature_Selection, 
             group = Feature_Selection,
             alpha = Feature_Selection)) + 
  geom_line(linetype = 'dashed', size=1) + 
  geom_point(size = 4) +
  xlab('Number of Markers') + 
  ylab('Weighted F1') +
  labs(color = '',
       alpha = '') +
  scale_x_log10() + 
  scale_y_continuous(breaks = seq.int(0,1,0.1)) +
  scale_color_manual(values = darkpal_v2) +
  facet_grid(cols = vars(Dataset)) +
  geom_hline(data = featureselection_eval_5k, 
             aes(yintercept = baseline, 
                 linetype = 'Baseline\n5K Markers'),
             color = 'dodgerblue',
             size = 1,
             alpha = 0.5) +
  scale_linetype_manual(name = "", 
                        values = 'dashed') +
  scale_alpha_manual(values=c(1, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8)) +
  theme_bw(base_size = 18) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'bottom',
        axis.text = element_text(color = 'black'),
        legend.text = element_text(color = 'black'))

ggsave(filename = paste0(workingdir,'/plots/',sample_no,"global_featureselection.pdf"), 
       units = 'in',
       width = 12,
       height = 6, 
       dpi = 300)


# horizontal plots - wrap
p_others <- featureselection_eval %>%
  filter(Num_Markers != 5000) %>%
  filter(Dataset != "GSE144236") %>%
  ggplot(aes(Num_Markers, 
             Weighted_F1, 
             color = Feature_Selection, 
             group = Feature_Selection,
             alpha = Feature_Selection)) + 
  geom_line(linetype = 'dashed', size=1) + 
  geom_point(size = 4) +
  xlab('Number of Markers') + 
  ylab('Weighted F1') +
  labs(color = '',
       alpha = '') +
  scale_x_log10() + 
  scale_y_continuous(breaks = seq.int(0,1,0.1)) +
  scale_color_manual(values = darkpal_v2) +
  facet_wrap(~Dataset, scales = "free_x", ncol = 2) + 
  geom_hline(data = subset(featureselection_eval_5k, 
                           Dataset != "GSE144236"), 
             aes(yintercept = baseline, 
                 linetype = 'Baseline\n5K Markers'),
             color = 'dodgerblue',
             size = 1,
             alpha = 0.5) +
  scale_linetype_manual(name = "", 
                        values = 'dashed') +
  scale_alpha_manual(values=c(1, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8)) +
  theme_bw(base_size = 18) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'bottom',
        axis.text = element_text(color = 'black'),
        legend.text = element_text(color = 'black'))

p_GSE144236 <- featureselection_eval %>%
  filter(Num_Markers != 5000) %>%
  filter(Dataset == "GSE144236") %>%
  ggplot(aes(Num_Markers, 
             Weighted_F1, 
             color = Feature_Selection, 
             group = Feature_Selection,
             alpha = Feature_Selection)) + 
  geom_line(linetype = 'dashed', size=1) + 
  geom_point(size = 4) +
  xlab('Number of Markers') + 
  ylab('Weighted F1') +
  labs(color = '',
       alpha = '') +
  scale_x_log10() + 
  scale_y_continuous(breaks = seq.int(0,1,0.1)) +
  expand_limits(x = c(0,300),
                y = c(0,1)) +
  coord_cartesian(expand = TRUE) +
  scale_color_manual(values = darkpal_v2) +
  facet_wrap(~Dataset, scales = "free_x", ncol = 2) + 
  geom_hline(data = subset(featureselection_eval_5k, 
                           Dataset == "GSE144236"), 
             aes(yintercept = baseline, 
                 linetype = 'Baseline\n5K Markers'),
             color = 'dodgerblue',
             size = 1,
             alpha = 0.5) +
  scale_linetype_manual(name = "", 
                        values = 'dashed') +
  scale_alpha_manual(values=c(1, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8)) +
  theme_bw(base_size = 18) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'none',
        axis.text = element_text(color = 'black'),
        legend.text = element_text(color = 'black'))

(p_GSE144236| p_others) + plot_layout(widths = c(1, 2), 
                                      guides = "collect") & 
  theme(legend.position = 'bottom')


ggsave(filename = paste0(workingdir,'/plots/',sample_no,"global_featureselection_wrapped.pdf"), 
       units = 'in',
       width = 12,
       height = 8, 
       dpi = 300)



# Variance Feature Selection ----------------------------------------------

# reading in the data
ft_variance <- read_csv(file = paste0(
  file.path(file.path(workingdir, datadir), '/Model Evaluation - Variance_Supervised_GlobalGeneMarkerSelection_Final.csv')))


# little bit of prepro for plotting
ft_variance <- ft_variance %>%
  mutate(Dataset = case_when(Dataset == 'ImmuneCSF_PBMC' ~ 'ImmuneCSF',
                             Dataset == 'GSE154567_COVID PBMC' ~ 'GSE154567',
                             Dataset == 'GSE144236_Immune cSCC' ~ 'GSE144236',
                             Dataset == 'Lukassen2020_Lung' ~ 'Lukassen2020',
                             Dataset == 'Mouse HDF (SCP1361)' ~ 'SCP1361'))

# arranging for ploting
ft_variance$Feature_Selection <- factor(ft_variance$Feature_Selection, 
                                        levels = c('scANNA', 
                                                   'SMaSH', 
                                                   'scGeneFit',
                                                   'Triku',
                                                   'Seurat', 
                                                   'CellRanger'))

# separate dataframe for plotting baselines
ft_variance_5k <- ft_variance %>%
  filter(Num_Markers == 5000) %>%
  group_by(Dataset) %>%
  summarise(baseline = mean(Percentage))


max_y <- ft_variance %>%
  filter(Num_Markers != 5000) %>%
  slice_max(order_by = Percentage) %>%
  pull(Percentage) %>%
  round(., 1) 



# horizontal plots 
ft_variance %>%
  filter(Num_Markers != 5000) %>%
  ggplot(aes(Num_Markers, 
             Percentage, 
             color = Feature_Selection, 
             group = Feature_Selection, 
             alpha=Feature_Selection)) + 
  geom_line(linetype = 'dashed', 
            size = 1) + 
  geom_point(size = 4) +
  xlab('Number of Markers') + 
  ylab('Total Variance Explained %') +
  labs(color = '', 
       alpha = '') +
  ggtitle(label = 'Global Feature Selection') +
  scale_x_log10() + 
  scale_y_continuous(labels = scales::label_percent(scale = 1)) +
  scale_color_manual(values = darkpal_v2) +
  scale_alpha_manual(values=c(1, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8)) +
  facet_grid(cols = vars(Dataset)) +
  theme_bw(base_size = 18) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'bottom',
        axis.text = element_text(color = 'black'),
        legend.text = element_text(color = 'black'))

ggsave(filename = paste0(workingdir,'/plots/',sample_no,"variance_global_featureselection.pdf"), 
       units = 'in',
       width = 12,
       height = 6, 
       dpi = 300)


# horizontal plots - wrap
p_others <- ft_variance %>%
  filter(Num_Markers != 5000) %>%
  filter(Dataset != "GSE144236") %>%
  ggplot(aes(Num_Markers, 
             Percentage, 
             color = Feature_Selection, 
             group = Feature_Selection, 
             alpha=Feature_Selection)) + 
  geom_line(linetype = 'dashed', 
            size = 1) + 
  geom_point(size = 4) +
  xlab('Number of Markers') + 
  ylab('Total Variance Explained %') +
  labs(color = '', 
       alpha = '') +
  scale_x_log10() + 
  scale_y_continuous(labels = scales::label_percent(scale = 1)) +
  scale_color_manual(values = darkpal_v2) +
  scale_alpha_manual(values=c(1, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8)) +
  facet_wrap(~Dataset, scales = "free_x", ncol = 2) + 
  theme_bw(base_size = 18) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'bottom',
        axis.text = element_text(color = 'black'),
        legend.text = element_text(color = 'black'))

p_GSE144236 <- ft_variance %>%
  filter(Num_Markers != 5000) %>%
  filter(Dataset == "GSE144236") %>%
  ggplot(aes(Num_Markers, 
             Percentage, 
             color = Feature_Selection, 
             group = Feature_Selection, 
             alpha=Feature_Selection)) + 
  geom_line(linetype = 'dashed', 
            size = 1) + 
  geom_point(size = 4) +
  xlab('Number of Markers') + 
  ylab('Total Variance Explained %') +
  labs(color = '', 
       alpha = '') +
  scale_x_log10() + 
  scale_y_continuous(labels = scales::label_percent(scale = 1)) +
  scale_color_manual(values = darkpal_v2) +
  scale_alpha_manual(values=c(1, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8)) +
  facet_wrap(~Dataset, scales = "free_x", ncol = 2) + 
  theme_bw(base_size = 18) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'none',
        axis.text = element_text(color = 'black'),
        legend.text = element_text(color = 'black'))


(p_GSE144236| p_others) + plot_layout(widths = c(1, 2), 
                                      guides = "collect") & 
  theme(legend.position = 'bottom')


ggsave(filename = paste0(workingdir,'/plots/',sample_no,"variance_global_featureselection_wrapped.pdf"), 
       units = 'in',
       width = 12,
       height = 8, 
       dpi = 300)


