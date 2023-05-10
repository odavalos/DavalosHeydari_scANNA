library(scPred)
library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(tictoc)
library(doParallel)


# Reading data ------------------------------------------------------------

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


# Prep scPred -------------------------------------------------------------

# load dataset
gse144236 <- LoadH5Seurat(file = '~/GSE144236/GSE144236_qc_hvg_anno_5k_raw.h5Seurat')

# read in metadata that contains split column
metadata <- read_csv(file = '~/DatasetSplits/Metadata_Splits/GSE144236_metadata_splits.csv')

# clean up metadata
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$...1
metadata <- metadata[,-1]

# check to make sure barcodes are ordered properly in data and metadata
table(rownames(gse144236@meta.data) == metadata$barcodes)

# append 5 split columns to the seurat metadata df
gse144236@meta.data <- cbind(gse144236@meta.data, metadata[,c("Split_1","Split_2", "Split_3", "Split_4", "Split_5")])

split_vect <- c("Split_1","Split_2", "Split_3", "Split_4", "Split_5")

# a function to run through the scPred process
runSplitscPred <- function(seuratdata, split, cluster_col = 'celltypes', model='svmRadial'){
  message(crayon::magenta(paste0("Running scPred model: ", model, ' on split: ', split,'\n')))
  

  # split data into reference - 'train' and query - 'test'
  reference <- seuratdata[, seuratdata@meta.data[, split] == "train"]
  
  reference@assays$RNA@counts['CASC6', sample(dim(reference)[2],1)] <- 5 # add a count to a random cell to make sure this gene runs through the classifer
  
  query <- seuratdata[, seuratdata@meta.data[, split] == "test"]


  # run in parallel use 4 tasks
  cl <- makePSOCKcluster(4)
  registerDoParallel(cl)
  
  # train the models - SVM
  tic("Total Runtime") # start timer
  tic("Feature space") # Subsection timer
  reference <- getFeatureSpace(reference, "celltypes")
  toc(log = TRUE) # feature space runtime
  tic("Training Models") # model training timer
  reference <- trainModel(reference, 
                          allowParallel = TRUE, 
                          seed = 2022,
                          resampleMethod = 'none',
                          tuneLength = 1) 
  toc(log = TRUE) # end model training
  toc(log = TRUE) # total run time timer end
  
  stopCluster(cl)
  
  runtime <- unlist(tic.log())
  runtime <- append('Parallel run: 4', runtime)
  runtime <- append(paste0('Running scPred - ', model,': ', sample_no, '-', split), runtime)
  runtime <- as.data.frame(runtime) 
  colnames(runtime) <- 'time'
  runtime <- runtime %>%
    separate(time, into = c('Part', 'Runtime'), sep = ':') %>%
    mutate(data_split = rep(split, nrow(.)))
  
  tic.clearlog()
  
  
  # classification of query cells
  query <- NormalizeData(query)
  query <- scPredict(query, reference)
  
  query@meta.data$scpred_prediction <- as.factor(query@meta.data$scpred_prediction)
  query@meta.data[[cluster_col]] <- factor(query@meta.data[[cluster_col]], levels(query$scpred_prediction))
  
  # model evaluation
  cm <- caret::confusionMatrix(query$scpred_prediction, 
                               reference = query$celltypes)
  
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
  
  res_list <- runSplitscPred(seuratdata = gse144236, split = i)
  
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
          file = paste0(workingdir,'/',resultsdir, '/',sample_no,'_runtime_svmRadial.csv'),
          row.names = TRUE)

write.csv(eval_df,
          file = paste0(workingdir,'/',resultsdir, '/',sample_no,'_modeleval_svmRadial.csv'),
          row.names = TRUE)

write.csv(meta_df,
          file = paste0(workingdir,'/',resultsdir, '/',sample_no,'_metadata_model_svmRadial.csv'),
          row.names = TRUE)




