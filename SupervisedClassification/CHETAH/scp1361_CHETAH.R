library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(tictoc)
library(CHETAH)
# library(doParallel)


# Reading data ------------------------------------------------------------

sample_no <- 'SCP1361' # change accordingly
tool_name <- 'CHETAH'

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


# Prep scPred -------------------------------------------------------------

# load dataset
scp1361 <- LoadH5Seurat(file = '~/Documents/NACT_Datasets/scp1361_qc_hvg_anno_5k.h5Seurat')

# read in metadata that contains split column
metadata <- read_csv(file = '~/Documents/NACT_Datasets/Metadata_Splits/SCP1361_metadata_splits.csv')

# clean up metadata
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$...1
metadata <- metadata[,-1]

# check to make sure barcodes are ordered properly in data and metadata
table(rownames(scp1361@meta.data) == metadata$barcodes)

# append 5 split columns to the seurat metadata df
scp1361@meta.data <- cbind(scp1361@meta.data, metadata[,c("Split_1","Split_2", "Split_3", "Split_4", "Split_5")])

split_vect <- c("Split_1","Split_2", "Split_3", "Split_4", "Split_5")


# # a function to run through the scClassify process
runSplitCHETAH <- function(seuratdata, split, cluster_col = 'celltypes', dataset = sample_no){
  message(crayon::magenta(paste0("Running CHETAH model on: ", dataset, ' on split: ', split,'\n')))
  
  gc(verbose = FALSE)
  
  # Split the dataset for training the model
  reference <- seuratdata[, seuratdata@meta.data[, split] == "train"]
  query <- seuratdata[, seuratdata@meta.data[, split] == "test"]
  
  ## Make SingleCellExperiments
  reference <- SingleCellExperiment(assays = list(counts = reference@assays$RNA@data),
                                    colData = DataFrame(celltypes = reference@meta.data[[cluster_col]]))
  
  input <- SingleCellExperiment(assays = list(counts = query@assays$RNA@data),
                                reducedDims = SimpleList(UMAP = query@reductions$umap@cell.embeddings))
  
  
  # train the models
  tic("Total Runtime") # start timer
  tic("Training Models") # model training timer
  ## Run CHETAH
  input <- CHETAHclassifier(input = input, ref_cells = reference)
  
  ## Extract celltypes:
  preds <- input$celltype_CHETAH
  
  toc(log = TRUE) # end model training
  toc(log = TRUE) # total run time timer end
  
  # log runtime
  runtime <- unlist(tic.log())
  runtime <- append(paste0('Running CHETAH - ', sample_no, '-', split), runtime)
  runtime <- as.data.frame(runtime)
  colnames(runtime) <- 'time'
  # runtime <- runtime %>%
  #   separate(time, into = c('Part', 'Runtime'), sep = ':') %>%
  #   mutate(data_split = rep(split, nrow(.)))
  
  tic.clearlog()
  
  
  # remove the intermediate predictions and call them unassigned
  preds <- ifelse(test = (!preds %in% unique(query@meta.data[[cluster_col]])), 
                  yes = 'unassigned', 
                  no = preds)
  
  # add the prediction back to query data and convert into a level
  query@meta.data$CHETAH_pred <- as.factor(preds)
  query@meta.data[[cluster_col]] <- factor(query@meta.data[[cluster_col]], 
                                           levels(query$CHETAH_pred))
  
  # model evaluation
  cm <- caret::confusionMatrix(query$CHETAH_pred,
                               reference = query@meta.data[[cluster_col]])
  
  
  model_metrics <- as.data.frame(cm[['byClass']]) %>%
    rownames_to_column(var = 'Class') %>%
    mutate(data_split = rep(split, nrow(.)))
  
  query_metadata <- query@meta.data %>%
    mutate(data_split = rep(split, nrow(.)))
  
  # list to save values
  return_list <- list(Runtime = runtime,
                      Model_Eval = model_metrics,
                      SplitMetadat = query_metadata)
  
  
  return(return_list)
  
}


all_splits <- list()

for(i in split_vect){
  
  
  res_list <- try(runSplitCHETAH(seuratdata = scp1361, split = i))
  
  # repeat {
  #   res_list <- try(runSplitCHETAH(seuratdata = scp1361, split = i))
  #   if (!(inherits(res_list,"try-error")))
  #     break
  # }
  
  all_splits[[i]] <- res_list
  
}


runtime_list <- list()
eval_list <- list()
meta_list <- list()

for(i in split_vect){
  runtime_list[[i]]<- all_splits[[i]][['Runtime']]
  eval_list[[i]]<- all_splits[[i]][['Model_Eval']]
  meta_list[[i]]<- all_splits[[i]][['SplitMetadat']]
}

runtime_df <- as.data.frame(data.table::rbindlist(runtime_list, use.names = TRUE))
eval_df <- as.data.frame(data.table::rbindlist(eval_list, use.names = TRUE))
meta_df <- as.data.frame(data.table::rbindlist(meta_list, use.names = TRUE))


write.csv(runtime_df,
          file = paste0(workingdir,'/',resultsdir, '/',sample_no,'_runtime_',tool_name,'.csv'),
          row.names = TRUE)

write.csv(eval_df,
          file = paste0(workingdir,'/',resultsdir, '/',sample_no,'_modeleval_',tool_name,'.csv'),
          row.names = TRUE)

write.csv(meta_df,
          file = paste0(workingdir,'/',resultsdir, '/',sample_no,'_metadata_',tool_name,'.csv'),
          row.names = TRUE)
