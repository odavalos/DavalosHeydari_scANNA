# This script is used to generate the model evaluation plots for the transfer 
# learning experiments in the scANNA manuscript.

library(tidyverse)
library(RColorBrewer)
library(rcartocolor)
library(patchwork)


# Reading in the data -----------------------------------------------------

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

# read in the data
model_eval <- read_csv(file = paste0(
  file.path(file.path(workingdir, datadir), 
            '/Model Evaluation - TransferLearning_Final.csv')))


# Plotting ----------------------------------------------------------------
base_best <- model_eval %>% 
  filter(Tuning %in% c("Baseline","Best Case"))

model_eval <- model_eval %>% 
  filter(!Tuning %in% c("Baseline","Best Case"))


# arranging for ploting
model_eval$Model <- factor(model_eval$Model, 
                           levels = c('scANNA (Ours)',
                                      "scArches",
                                      'scNym',
                                      'ACTINN-TL'))

# modify rcartocolor's 'Safe' color palette to emphasize scANNA
safepal <- rcartocolor::carto_pal(n=8, 'Safe')
safepal_v2 <- replace(safepal, c(1,5), safepal[c(5,1)])

# weighted f1 plot
model_eval %>%
  ggplot(aes(Tuning, 
             y=`Weighted F1`, 
             fill = Model, 
             groups = Dataset,
             alpha = Model)) + 
  geom_col(width = 0.75,
           position =  position_dodge(0.75), 
           color = "white") +
  ylab('Weighted F1') + 
  xlab('') +
  scale_y_continuous(breaks = seq.int(0,1,0.1)) +
  scale_fill_manual(name = 'Model', 
                    values = safepal_v2) +
  scale_alpha_manual(values=c(1, 
                              0.75, 
                              0.75, 
                              0.75)) +
  facet_grid(~Dataset, 
             scales = "free_x") +
  geom_hline(data = base_best, 
             aes(yintercept = `Weighted F1`, 
                 color = Tuning,
                 linetype = Tuning),
             size = 1) +
  scale_color_manual(name = "", 
                     values = c("orange", 
                                "dodgerblue")) +
  scale_linetype_manual(name = "", 
                        values = c("dotdash",
                                   "dashed")) + 
  theme_bw(base_size = 18) + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'bottom',
        axis.text = element_text(color = 'black'),
        legend.text = element_text(color = 'black'))

ggsave(filename = paste0(workingdir,'/plots/',sample_no,"transferlearning_barplot.pdf"), 
       units = 'in',
       height = 8,
       width = 12,
       dpi = 300)



