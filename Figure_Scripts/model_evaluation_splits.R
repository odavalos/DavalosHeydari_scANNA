# This script is used to generate the model evaluation plots for 
# supervised learning experiments in the scANNA manuscript.

library(tidyverse)
library(RColorBrewer)
library(rcartocolor)
library(patchwork)


# Directory setup ----------------------------------------------------------

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

# Load runtime data
model_eval <- read_csv(file = paste0(
  file.path(file.path(workingdir, datadir), 
            '/Model Evaluation - SupervisedModels_Splits.csv')))

# little bit of prepro for plotting
model_eval <- model_eval %>%
  mutate(Dataset = case_when(Dataset == 'Immune CSF' ~ 'ImmuneCSF',
                             Dataset == 'GSE154567_COVID PBMC' ~ 'GSE154567',
                             Dataset == 'GSE144236_Immune cSCC' ~ 'GSE144236',
                             Dataset == 'Lukassen2020_Lung' ~ 'Lukassen2020',
                             Dataset == 'Mouse HDF (SCP1361)' ~ 'SCP1361'))


# arranging for ploting
model_eval$Model <- factor(model_eval$Model, 
                           levels = c('scANNA', 
                                      'ACTINN',
                                      'scClassify',
                                      'SingleCellNet',
                                      'CHETAH',
                                      'Random Forest', 
                                      'scPred (SVM)', 
                                      'scPred (NN)'))


model_eval %>%
  group_by(Dataset, Model) %>%
  summarise(mean = mean(Accuracy),
            sd = sd(Accuracy))

# modify rcartocolor's 'Safe' color palette to emphasize scANNA
safepal <- rcartocolor::carto_pal(n=8, 'Safe')
safepal_v2 <- replace(safepal, c(1,5), safepal[c(5,1)])

# accuracy plot
model_eval %>%
  group_by(Dataset, Model) %>%
  summarise(mean = mean(Accuracy),
            sd = sd(Accuracy)) %>%
  ggplot(aes(Dataset, 
             y=mean, 
             fill = Model, 
             groups = Dataset,
             alpha = Model)) + 
  geom_col(width = 0.75,
           position =  position_dodge(0.75), 
           color = "white") +
  geom_errorbar(aes(ymin=mean-sd, 
                    ymax=mean+sd), 
                width=.4,
                position=position_dodge(.75)) + 
  ylab('Accuracy') + 
  xlab('') +
  coord_cartesian(ylim = c(0.5,1)) +
  scale_fill_manual(name = 'Model', 
                    values = safepal_v2) +
  scale_alpha_manual(values=c(1, 
                              0.75, 
                              0.75, 
                              0.75, 
                              0.75, 
                              0.75, 
                              0.75, 
                              0.75)) +
  theme_bw(base_size = 18) + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'right',
        axis.text = element_text(color = 'black'),
        legend.text = element_text(color = 'black'))

ggsave(filename = paste0(workingdir,'/plots/',sample_no,"accuracy_barplot.pdf"), 
       units = 'in',
       height = 8,
       width = 12,
       dpi = 300)


# accuracy plot - facet wrap
p_others <- model_eval %>%
  group_by(Dataset, Model) %>%
  summarise(mean = mean(Accuracy),
            sd = sd(Accuracy)) %>%
  filter(Dataset != "GSE144236") %>%
  ggplot(aes(Dataset, 
             y=mean, 
             fill = Model, 
             groups = Dataset,
             alpha = Model)) + 
  geom_col(width = 0.75,
           position =  position_dodge(0.75), 
           color = "white") +
  geom_errorbar(aes(ymin=mean-sd, 
                    ymax=mean+sd), 
                width=.4,
                position=position_dodge(.75)) + 
  ylab('Accuracy') + 
  xlab('') +
  coord_cartesian(ylim = c(0.5,1)) +
  scale_fill_manual(name = 'Model', 
                    values = safepal_v2) +
  scale_alpha_manual(values=c(1, 
                              0.75, 
                              0.75, 
                              0.75, 
                              0.75, 
                              0.75, 
                              0.75, 
                              0.75)) +
  facet_wrap(~Dataset, scales = "free_x", ncol = 2) + 
  theme_bw(base_size = 18) + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'bottom',
        axis.text = element_text(color = 'black'),
        legend.text = element_text(color = 'black'),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

p_GSE144236 <- model_eval %>%
  group_by(Dataset, Model) %>%
  summarise(mean = mean(Accuracy),
            sd = sd(Accuracy)) %>%
  filter(Dataset == "GSE144236") %>%
  ggplot(aes(Dataset, 
             y=mean, 
             fill = Model, 
             groups = Dataset,
             alpha = Model)) + 
  geom_col(width = 0.75,
           position =  position_dodge(0.75), 
           color = "white") +
  geom_errorbar(aes(ymin=mean-sd, 
                    ymax=mean+sd), 
                width=.4,
                position=position_dodge(.75)) + 
  ylab('Accuracy') + 
  xlab('') +
  coord_cartesian(ylim = c(0.5,1)) +
  scale_fill_manual(name = 'Model', 
                    values = safepal_v2) +
  scale_alpha_manual(values=c(1, 
                              0.75, 
                              0.75, 
                              0.75, 
                              0.75, 
                              0.75, 
                              0.75, 
                              0.75)) +
  facet_wrap(~Dataset, scales = "free_x", ncol = 2) + 
  theme_bw(base_size = 18) + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'none',
        axis.text = element_text(color = 'black'),
        legend.text = element_text(color = 'black'),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

(p_GSE144236| p_others) + plot_layout(widths = c(1, 2), 
                                      guides = "collect") & 
  theme(legend.position = 'bottom')

ggsave(filename = paste0(workingdir,'/plots/',sample_no,"accuracy_facetwrap_barplot.pdf"), 
       units = 'in',
       height = 8,
       width = 12,
       dpi = 300)


# accuracy reduced plot
model_eval %>%
  group_by(Dataset, Model) %>%
  summarise(mean = mean(Accuracy),
            sd = sd(Accuracy)) %>%
  ggplot(aes(Dataset, y=mean, fill = Model, groups = Dataset)) + 
  geom_col(position = 'dodge', color = 'white') +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9)) + 
  ylab('Accuracy') + 
  xlab('') +
  scale_fill_manual(name = 'Model', 
                    values = safepal_v2) +
  coord_cartesian(ylim = c(0.5, 1)) +
  theme_bw(base_size = 14) + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'right',
        axis.text = element_text(color = 'black'),
        legend.text = element_text(color = 'black'))

ggsave(filename = paste0(workingdir,'/plots/',sample_no,"accuracy_reducedaxis_barplot.pdf"), 
       units = 'in',
       height = 8,
       width = 12,
       dpi = 300)



# Macro F1 ----------------------------------------------------------------


# macro f1 plot
model_eval %>%
  group_by(Dataset, Model) %>%
  summarise(mean = mean(`Macro F1`),
            sd = sd(`Macro F1`)) %>%
  ggplot(aes(Dataset, 
             y=mean, 
             fill = Model, 
             groups = Dataset,
             alpha = Model)) + 
  geom_col(width = 0.75,
           position =  position_dodge(0.75), 
           color = "white") +
  geom_errorbar(aes(ymin=mean-sd, 
                    ymax=mean+sd), 
                width=.4,
                position=position_dodge(.75)) + 
  ylab('Macro F1') + 
  xlab('') +
  scale_fill_manual(name = 'Model', 
                    values = safepal_v2) +
  scale_alpha_manual(values=c(1, 
                              0.75, 
                              0.75, 
                              0.75, 
                              0.75, 
                              0.75, 
                              0.75, 
                              0.75)) +
  theme_bw(base_size = 18) + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'bottom',
        axis.text = element_text(color = 'black'),
        legend.text = element_text(color = 'black'))

ggsave(filename = paste0(workingdir,'/plots/',sample_no,"macroF1_barplot.pdf"), 
       units = 'in',
       height = 8,
       width = 12,
       dpi = 300)

# macro f1 plot - facet wrap
p_others <- model_eval %>%
  group_by(Dataset, Model) %>%
  summarise(mean = mean(`Macro F1`),
            sd = sd(`Macro F1`)) %>%
  filter(Dataset != "GSE144236") %>%
  ggplot(aes(Dataset, 
             y=mean, 
             fill = Model, 
             groups = Dataset,
             alpha = Model)) + 
  geom_col(width = 0.75,
           position =  position_dodge(0.75), 
           color = "white") +
  geom_errorbar(aes(ymin=mean-sd, 
                    ymax=mean+sd), 
                width=.4,
                position=position_dodge(.75)) + 
  ylab('Macro F1') + 
  xlab('') +
  scale_fill_manual(name = 'Model', 
                    values = safepal_v2) +
  scale_alpha_manual(values=c(1, 
                              0.75, 
                              0.75, 
                              0.75, 
                              0.75, 
                              0.75, 
                              0.75, 
                              0.75)) +
  facet_wrap(~Dataset, scales = "free_x", ncol = 2) + 
  theme_bw(base_size = 18) + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'bottom',
        axis.text = element_text(color = 'black'),
        legend.text = element_text(color = 'black'),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

p_GSE144236 <- model_eval %>%
  group_by(Dataset, Model) %>%
  summarise(mean = mean(`Macro F1`),
            sd = sd(`Macro F1`)) %>%
  filter(Dataset == "GSE144236") %>%
  ggplot(aes(Dataset, 
             y=mean, 
             fill = Model, 
             groups = Dataset,
             alpha = Model)) + 
  geom_col(width = 0.75,
           position =  position_dodge(0.75), 
           color = "white") +
  geom_errorbar(aes(ymin=mean-sd, 
                    ymax=mean+sd), 
                width=.4,
                position=position_dodge(.75)) + 
  ylab('Macro F1') + 
  xlab('') +
  scale_fill_manual(name = 'Model', 
                    values = safepal_v2) +
  scale_alpha_manual(values=c(1, 
                              0.75, 
                              0.75, 
                              0.75, 
                              0.75, 
                              0.75, 
                              0.75, 
                              0.75)) +
  facet_wrap(~Dataset, scales = "free_x", ncol = 2) + 
  theme_bw(base_size = 18) + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'none',
        axis.text = element_text(color = 'black'),
        legend.text = element_text(color = 'black'),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

(p_GSE144236| p_others) + plot_layout(widths = c(1, 2), 
                                      guides = "collect") & 
  theme(legend.position = 'bottom')

ggsave(filename = paste0(workingdir,'/plots/',sample_no,"macroF1_facetwrap_barplot.pdf"), 
       units = 'in',
       height = 8,
       width = 12,
       dpi = 300)


# macro f1 reduced plot
model_eval %>%
  group_by(Dataset, Model) %>%
  summarise(mean = mean(`Macro F1`),
            sd = sd(`Macro F1`)) %>%
  ggplot(aes(Dataset, y=mean, fill = Model, groups = Dataset)) + 
  geom_col(position = 'dodge', color = 'white') +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9)) + 
  ylab('Macro F1') + 
  xlab('') +
  scale_fill_manual(name = 'Model', 
                    values = safepal_v2) +
  coord_cartesian(ylim = c(0.2, 1)) +
  theme_bw(base_size = 14) + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'right',
        axis.text = element_text(color = 'black'),
        legend.text = element_text(color = 'black'))

ggsave(filename = paste0(workingdir,'/plots/',sample_no,"macroF1_reducedaxis_barplot.pdf"), 
       units = 'in',
       height = 8,
       width = 12,
       dpi = 300)



# Weighted F1 -------------------------------------------------------------

# weighted f1 plot
model_eval %>%
  group_by(Dataset, Model) %>%
  summarise(mean = mean(`Weighted F1`),
            sd = sd(`Weighted F1`)) %>%
  ggplot(aes(Dataset, 
             y=mean, 
             fill = Model, 
             groups = Dataset,
             alpha = Model)) + 
  geom_col(width = 0.75,
           position =  position_dodge(0.75), 
           color = "white") +
  geom_errorbar(aes(ymin=mean-sd, 
                    ymax=mean+sd), 
                width=.4,
                position=position_dodge(.75)) + 
  ylab('Weighted F1') + 
  xlab('') +
  coord_cartesian(ylim = c(0.5, 1)) +
  scale_fill_manual(name = 'Model', 
                    values = safepal_v2) +
  scale_alpha_manual(values=c(1, 
                              0.75, 
                              0.75, 
                              0.75, 
                              0.75, 
                              0.75, 
                              0.75, 
                              0.75)) +
  theme_bw(base_size = 18) + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'bottom',
        axis.text = element_text(color = 'black'),
        legend.text = element_text(color = 'black'))

ggsave(filename = paste0(workingdir,'/plots/',sample_no,"weightedF1_barplot.pdf"), 
       units = 'in',
       height = 8,
       width = 12,
       dpi = 300)


# weighted f1 plot - facet wrap
p_others <- model_eval %>%
  group_by(Dataset, Model) %>%
  summarise(mean = mean(`Weighted F1`),
            sd = sd(`Weighted F1`)) %>%
  filter(Dataset != "GSE144236") %>%
  ggplot(aes(Dataset, 
             y=mean, 
             fill = Model, 
             groups = Dataset,
             alpha = Model)) + 
  geom_col(width = 0.75,
           position =  position_dodge(0.75), 
           color = "white") +
  geom_errorbar(aes(ymin=mean-sd, 
                    ymax=mean+sd), 
                width=.4,
                position=position_dodge(.75)) + 
  ylab('Weighted F1') + 
  xlab('') +
  coord_cartesian(ylim = c(0.5,1)) + 
  scale_fill_manual(name = 'Model', 
                    values = safepal_v2) +
  scale_alpha_manual(values=c(1, 
                              0.75, 
                              0.75, 
                              0.75, 
                              0.75, 
                              0.75, 
                              0.75, 
                              0.75)) +
  facet_wrap(~Dataset, scales = "free_x", ncol = 2) + 
  theme_bw(base_size = 18) + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'right',
        axis.text = element_text(color = 'black'),
        legend.text = element_text(color = 'black'),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

p_GSE144236 <- model_eval %>%
  group_by(Dataset, Model) %>%
  summarise(mean = mean(`Weighted F1`),
            sd = sd(`Weighted F1`)) %>%
  filter(Dataset == "GSE144236") %>%
  ggplot(aes(Dataset, 
             y=mean, 
             fill = Model, 
             groups = Dataset,
             alpha = Model)) + 
  geom_col(width = 0.75,
           position =  position_dodge(0.75), 
           color = "white") +
  geom_errorbar(aes(ymin=mean-sd, 
                    ymax=mean+sd), 
                width=.4,
                position=position_dodge(.75)) + 
  ylab('Weighted F1') + 
  xlab('') +
  coord_cartesian(ylim = c(0.5,1)) +
  scale_fill_manual(name = 'Model', 
                    values = safepal_v2) +
  scale_alpha_manual(values=c(1, 
                              0.75, 
                              0.75, 
                              0.75, 
                              0.75, 
                              0.75, 
                              0.75, 
                              0.75)) +
  facet_wrap(~Dataset, scales = "free_x", ncol = 2) + 
  theme_bw(base_size = 18) + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'none',
        axis.text = element_text(color = 'black'),
        legend.text = element_text(color = 'black'),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

(p_GSE144236| p_others) + plot_layout(widths = c(1, 2), 
                                      guides = "collect") & 
  theme(legend.position = 'bottom')

ggsave(filename = paste0(workingdir,'/plots/',sample_no,"weightedF1_facetwrap_barplot.pdf"), 
       units = 'in',
       height = 8,
       width = 12,
       dpi = 300)



# weighted f1 reduced plot
model_eval %>%
  group_by(Dataset, Model) %>%
  summarise(mean = mean(`Weighted F1`),
            sd = sd(`Weighted F1`)) %>%
  ggplot(aes(Dataset, y=mean, fill = Model, groups = Dataset)) + 
  geom_col(position = 'dodge', color = 'white') +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9)) + 
  ylab('Weighted F1') + 
  xlab('') +
  scale_fill_manual(name = 'Model', 
                    values = safepal_v2) +
  coord_cartesian(ylim = c(0.6, 1)) +
  theme_bw(base_size = 14) + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'right',
        axis.text = element_text(color = 'black'),
        legend.text = element_text(color = 'black'))

ggsave(filename = paste0(workingdir,'/plots/',sample_no,"weightedF1_reducedaxis_barplot.pdf"), 
       units = 'in',
       height = 8,
       width = 12,
       dpi = 300)




# Runtime -----------------------------------------------------------------


# Load runtime data
model_eval_time <- read_csv(file = paste0(
  file.path(file.path(workingdir, datadir), 
            '/Model Evaluation - SupervisedModels_Splits_Final.csv')))

# little bit of prepro for plotting
model_eval_time <- model_eval_time %>%
  mutate(Dataset = case_when(Dataset == 'Immune CSF' ~ 'ImmuneCSF',
                             Dataset == 'GSE154567_COVID PBMC' ~ 'GSE154567',
                             Dataset == 'GSE144236_Immune cSCC' ~ 'GSE144236',
                             Dataset == 'Lukassen2020_Lung' ~ 'Lukassen2020',
                             Dataset == 'Mouse HDF (SCP1361)' ~ 'SCP1361'))


# arranging for ploting
model_eval_time$Model <- factor(model_eval_time$Model, 
                           levels = c('scANNA', 
                                      'ACTINN',
                                      'scClassify',
                                      'SingleCellNet',
                                      'CHETAH',
                                      'Random Forest', 
                                      'scPred (SVM)', 
                                      'scPred (NN)'))

# Create a subset of the data for benchmark analysis (time)
model_eval_time <- model_eval_time %>%
  filter(Model %in% c("scANNA", 
                      "ACTINN", 
                      "Random Forest", 
                      "scPred (SVM)", 
                      "scPred (NN)"))

# GPU time
gpu_times <- model_eval_time %>%
  filter(GPU == "GPU")

# CPU time
model_eval_time_cpuonly <- model_eval_time %>%
  filter(!GPU %in% "GPU")

# GPU and CPU time
model_eval_time_recomended <- model_eval_time %>%
  filter(!Model %in% c("scANNA", "ACTINN")) %>% 
  rbind(., gpu_times)


# GPU and CPU - point plot
model_eval_time_recomended %>%
  group_by(Dataset, Model, Interpretability) %>%
  summarise(mean_time = mean(Runtime/60),
            sd_time = sd(Runtime/60),
            mean_train = mean(`Training Support`)) %>%
  ggplot(aes(mean_train, mean_time, color = Model, group = Model, shape = Interpretability)) + 
  geom_line(linetype = 'dashed', size = 1) + 
  geom_point(size = 4) + 
  geom_errorbar(aes(ymin=mean_time-sd_time, ymax=mean_time+sd_time), 
                width=.4,
                position=position_dodge(.9)) +
  xlab('Training Support') + 
  ylab('Average Runtime (Minutes)') +
  scale_y_log10() +
  scale_y_continuous(breaks = seq.int(0,12,2),
                     limits = c(0,12)) +
  scale_color_manual(name = 'Model', 
                    values =  viridis::viridis(5)) +
  theme_bw(base_size = 18) + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'right',
        axis.text = element_text(color = 'black'),
        legend.text = element_text(color = 'black'),
        axis.text.x = element_text(angle = 0, 
                                   vjust = 0))

ggsave(filename = paste0(workingdir,'/plots/',sample_no,"runtime_trainingsuport_recomended.pdf"), 
       units = 'in',
       height = 8,
       width = 12,
       dpi = 300)

# GPU and CPU - ribbon plot
model_eval_time_recomended %>%
  group_by(Dataset, Model, Interpretability) %>%
  summarise(mean_time = mean(Runtime/60),
            sd_time = sd(Runtime/60),
            mean_train = mean(`Training Support`)) %>%
  ggplot(aes(mean_train, 
             mean_time, 
             color = Model, 
             fill = Model,
             group = Model, 
             shape = Interpretability,
             alpha = Model)) + 
  geom_line(linetype = 'dashed', size = 1) + 
  geom_point(size = 3) + 
  geom_ribbon(aes(ymin=mean_time-sd_time,
                  ymax=mean_time+sd_time, 
                  group = Model), 
              alpha=0.15) + 
  xlab('Training Support') + 
  ylab('Average Runtime (Minutes)') +
  scale_y_log10() +
  scale_color_manual(name = 'Model', 
                     values = viridis::viridis(5)) +
  scale_fill_manual(name = 'Model', 
                    values = viridis::viridis(5)) +
  scale_alpha_manual(values=c(1, 
                              0.5, 
                              0.5, 
                              0.5, 
                              0.5, 
                              0.5, 
                              0.5, 
                              0.5)) +
  theme_bw(base_size = 18) + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'right',
        axis.text = element_text(color = 'black'),
        legend.text = element_text(color = 'black'))

ggsave(filename = paste0(workingdir,'/plots/',sample_no,"runtime_trainingsuport_recomended_ribbon.pdf"), 
       units = 'in',
       height = 8,
       width = 12,
       dpi = 300)


# CPU Times only - ribbon plot
model_eval_time_cpuonly %>%
  group_by(Dataset, Model, Interpretability) %>%
  summarise(mean_time = mean(Runtime/60),
            sd_time = sd(Runtime/60),
            mean_train = mean(`Training Support`)) %>%
  ggplot(aes(mean_train, 
             mean_time, 
             color = Model, 
             fill = Model,
             group = Model, 
             shape = Interpretability,
             alpha = Model)) + 
  geom_line(linetype = 'dashed', size = 1) + 
  geom_point(size = 3) + 
  geom_ribbon(aes(ymin=mean_time-sd_time,
                  ymax=mean_time+sd_time, 
                  group = Model), 
              alpha=0.15) + 
  xlab('Training Support') + 
  ylab('Average Runtime (Minutes)') +
  scale_y_log10() +
  scale_color_manual(name = 'Model', 
                     values = viridis::viridis(5)) +
  scale_fill_manual(name = 'Model', 
                     values = viridis::viridis(5)) +
  scale_alpha_manual(values=c(1, 
                              0.5, 
                              0.5, 
                              0.5, 
                              0.5, 
                              0.5, 
                              0.5, 
                              0.5)) +
  theme_bw(base_size = 18) + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'right',
        axis.text = element_text(color = 'black'),
        legend.text = element_text(color = 'black'))


ggsave(filename = paste0(workingdir,'/plots/',sample_no,"runtime_trainingsuport_cputimes_ribbon.pdf"), 
       units = 'in',
       height = 8,
       width = 12,
       dpi = 300)


# CPU Times only - point plot
model_eval_time_cpuonly %>%
  group_by(Dataset, Model, Interpretability) %>%
  summarise(mean_time = mean(Runtime/60),
            sd_time = sd(Runtime/60),
            mean_train = mean(`Training Support`)) %>%
  ggplot(aes(mean_train, mean_time, color = Model, group = Model, shape = Interpretability)) + 
  geom_line(linetype = 'dashed', size = 1) + 
  geom_point(size = 4) + 
  geom_errorbar(aes(ymin=mean_time-sd_time, 
                    ymax=mean_time+sd_time), 
                width=.4,
                position=position_dodge(.9)) +
  xlab('Training Support') + 
  ylab('Average Runtime (Minutes)') +
  scale_y_log10() +
  scale_color_manual(name = 'Model', 
                     values = viridis::viridis(5)) + 
  theme_bw(base_size = 18) + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'right',
        axis.text = element_text(color = 'black'),
        legend.text = element_text(color = 'black'),
        axis.text.x = element_text(angle = 0, 
                                   vjust = 0))

ggsave(filename = paste0(workingdir,'/plots/',sample_no,"runtime_trainingsuport_cputimes.pdf"), 
       units = 'in',
       height = 8,
       width = 12,
       dpi = 300)


