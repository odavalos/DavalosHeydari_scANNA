library(scClassify)
library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(tictoc)
# library(doParallel)


# Reading data ------------------------------------------------------------

sample_no <- 'GSE154567' # change accordingly
tool_name <- 'scClassify'

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
GSE154567 <- LoadH5Seurat(file = '/Volumes/BigDrive/GradSchool/NACT_NewDatasets/GSE154567/GSE154567_qc_hvg_anno_5k_raw.h5Seurat')

# read in metadata that contains split column
metadata <- read_csv(file = '/Volumes/BigDrive/NACT_ComparisonModels/DatasetSplits/Metadata_Splits/GSE154567_metadata_splits.csv')

# clean up metadata
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$index
metadata <- metadata[,-1]
colnames(metadata)[colnames(metadata) == "barcodes"] ="barcode_sub"
colnames(metadata)[colnames(metadata) == "index"] ="barcodes"

# check to make sure barcodes are ordered properly in data and metadata
table(rownames(GSE154567@meta.data) == metadata$barcodes)

GSE154567@meta.data <- metadata

# append 5 split columns to the seurat metadata df
#GSE154567@meta.data <- cbind(GSE154567@meta.data, metadata[,c("Split_1","Split_2", "Split_3", "Split_4", "Split_5")])

split_vect <- c("Split_1","Split_2", "Split_3", "Split_4", "Split_5")


# reference <- scp1361[, scp1361@meta.data[, split_vect[1]] == "train"]
# query <- scp1361[, scp1361@meta.data[, split_vect[1]] == "test"]


GSE154567 <- DietSeurat(GSE154567, 
                        counts = TRUE, 
                        data = TRUE, 
                        scale.data = FALSE,
                        dimreducs = c('umap', 'pca'))


# # a function to run through the scClassify process
runSplitscClassify <- function(seuratdata, split, cluster_col = 'celltypes', dataset = sample_no){
  message(crayon::magenta(paste0("Running scClassify model on: ", dataset, ' on split: ', split,'\n')))
  
  gc(verbose = FALSE)
  
  # Split the dataset for training the model
  reference <- seuratdata[, seuratdata@meta.data[, split] == "train"]
  query <- seuratdata[, seuratdata@meta.data[, split] == "test"]
  row_len <- nrow(query@meta.data)
  data_split <- rep(split, row_len)
  
  
  # train the models
  tic("Total Runtime") # start timer
  tic("Training Models") # model training timer
  scClassify_res <- scClassify(exprsMat_train = reference@assays$RNA@data,
                               cellTypes_train = reference@meta.data$celltypes,
                               exprsMat_test = list(query = query@assays$RNA@data),
                               cellTypes_test = list(query = query@meta.data$celltypes),
                               tree = "HOPACH",
                               algorithm = "WKNN",
                               selectFeatures = c("limma"),
                               similarity = c("pearson"),
                               returnList = FALSE,
                               parallel = TRUE,
                               verbose = FALSE)
  toc(log = TRUE) # end model training
  toc(log = TRUE) # total run time timer end
  
  # log runtime
  runtime <- unlist(tic.log())
  runtime <- append(paste0('Running scClassify - ', sample_no, '-', split), runtime)
  runtime <- as.data.frame(runtime)
  colnames(runtime) <- 'time'
  # runtime <- runtime %>%
  #   separate(time, into = c('Part', 'Runtime'), sep = ':') %>%
  #   mutate(data_split = rep(split, nrow(.)))
  
  tic.clearlog()
  
  # pull put the predictions
  preds <- scClassify_res$testRes$query$pearson_WKNN_limma$predRes
  
  
  # remove the intermediate predictions and call them unassigned
  preds <- ifelse(test = (!preds %in% unique(query@meta.data[[cluster_col]])), 
                  yes = 'unassigned', 
                  no = preds)
  
  # add the prediction back to query data and convert into a level
  query@meta.data$scClassify_pred <- as.factor(preds)
  query@meta.data[[cluster_col]] <- factor(query@meta.data[[cluster_col]], 
                                           levels(query$scClassify_pred))
  
  # model evaluation
  cm <- caret::confusionMatrix(query$scClassify_pred,
                               reference = query@meta.data[[cluster_col]])
  
  
  model_metrics <- as.data.frame(cm[['byClass']]) %>%
    rownames_to_column(var = 'Class') %>%
    mutate(data_split = rep(split, nrow(.)))
  
  # query_metadata <- query@meta.data %>%
  #   mutate(data_split = rep(split, nrow(.)))
  query@meta.data$data_split <- data_split
  query_metadata <- query@meta.data
  
  
  # list to save values
  return_list <- list(Runtime = runtime,
                      Model_Eval = model_metrics,
                      SplitMetadat = query_metadata)
  
  
  return(return_list)
  
}


all_splits <- list()

for(i in split_vect){
  
  res_list <- try(runSplitscClassify(seuratdata = GSE154567, split = i))
  # repeat {
  #   res_list <- try(runSplitscClassify(seuratdata = immunecsf, split = i))
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

