library(scClassify)
library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(tictoc)
# library(doParallel)


# Reading data ------------------------------------------------------------

sample_no <- 'SCP1361' # change accordingly
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
scp1361 <- LoadH5Seurat(file = '/Volumes/BigDrive/GradSchool/NACT_NewDatasets/SCP1361/scp1361_qc_hvg_anno_5k.h5Seurat')

# read in metadata that contains split column
metadata <- read_csv(file = '/Volumes/BigDrive/NACT_ComparisonModels/DatasetSplits/Metadata_Splits/SCP1361_metadata_splits.csv')

# clean up metadata
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$...1
metadata <- metadata[,-1]

# check to make sure barcodes are ordered properly in data and metadata
table(rownames(scp1361@meta.data) == metadata$barcodes)

# append 5 split columns to the seurat metadata df
scp1361@meta.data <- cbind(scp1361@meta.data, metadata[,c("Split_1","Split_2", "Split_3", "Split_4", "Split_5")])

split_vect <- c("Split_1","Split_2", "Split_3", "Split_4", "Split_5")


# reference <- scp1361[, scp1361@meta.data[, split_vect[1]] == "train"]
# query <- scp1361[, scp1361@meta.data[, split_vect[1]] == "test"]


# # a function to run through the scClassify process
runSplitscClassify <- function(seuratdata, split, cluster_col = 'celltypes', dataset = sample_no){
  message(crayon::magenta(paste0("Running scClassify model on: ", dataset, ' on split: ', split,'\n')))
  
  gc(verbose = FALSE)
  
  # Split the dataset for training the model
  reference <- seuratdata[, seuratdata@meta.data[, split] == "train"]
  query <- seuratdata[, seuratdata@meta.data[, split] == "test"]
  
  
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
                               parallel = FALSE,
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

  repeat {
    res_list <- try(runSplitscClassify(seuratdata = scp1361, split = i))
    if (!(inherits(res_list,"try-error")))
      break
  }

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






# scClassify_res <- scClassify(exprsMat_train = reference@assays$RNA@data,
#                              cellTypes_train = reference@meta.data$celltypes,
#                              exprsMat_test = list(query = query@assays$RNA@data),
#                              cellTypes_test = list(query = query@meta.data$celltypes),
#                              tree = "HOPACH",
#                              algorithm = "WKNN",
#                              selectFeatures = c("limma"),
#                              similarity = c("pearson"),
#                              returnList = FALSE,
#                              # parallel = TRUE,
#                              verbose = TRUE)



# plotCellTypeTree(cellTypeTree(scClassify_res$trainRes))


# table(scClassify_res$testRes$query$pearson_WKNN_limma$predRes, query@meta.data$celltypes)

# # pull put the predictions
# preds <- scClassify_res$testRes$query$pearson_WKNN_limma$predRes
# 
# 
# # remove the intermediate predictions and call them unassigned
# preds <- ifelse(test = (!preds %in% unique(query@meta.data[[cluster_col]])), 
#                 yes = 'unassigned', 
#                 no = preds)
# 
# # add the prediction back to query data and convert into a level
# query@meta.data$scClassify_pred <- as.factor(preds)
# query@meta.data[[cluster_col]] <- factor(query@meta.data[[cluster_col]], 
#                                          levels(query$scClassify_pred))
# 
# # model evaluation
# cm <- caret::confusionMatrix(query$scClassify_pred,
#                              reference = query$celltypes)
# 
# 
# model_metrics <- as.data.frame(cm[['byClass']]) %>%
#   rownames_to_column(var = 'Class') %>%
#   mutate(data_split = rep(split, nrow(.)))
# 
# query_metadata <- query@meta.data %>%
#   mutate(data_split = rep(split, nrow(.)))
# 
# 
# # list to save values
# return_list <- list(Runtime = runtime,
#                     Model_Eval = model_metrics,
#                     SplitMetadat = query_metadata)
# 
# 
# return(return_list)


# table(scClassify_res$testRes$query$pearson_WKNN_limma$predRes == scClassify_res$testRes$query$pearson_WKNN_limma$predLabelMat[,4])



# # a function to run through the scClassify process
# runSplitscClassify <- function(seuratdata, split, cluster_col = 'celltypes', mod='nnet'){
#   message(crayon::magenta(paste0("Running scPred model: ", mod, ' on split: ', split,'\n')))
#   
#   rnd_seed <- sample(c(1:4000),1)
#   
#   # split data into reference - 'train' and query - 'test'
#   reference <- seuratdata[, seuratdata@meta.data[, split] == "train"]
#   query <- seuratdata[, seuratdata@meta.data[, split] == "test"]
#   
#   # # Find genes that do not vary - reference
#   # vf_r <- VariableFeatures(FindVariableFeatures(reference, nfeatures = 5000))
#   # nonvargenes_r <- setdiff(rownames(reference), vf_r)
#   # 
#   # if(length(nonvargenes_r)>0){
#   #   # add one pseudo count to allow scPred to keep non var genes
#   #   for(i in nonvargenes_r){
#   #     message(crayon::magenta(paste0(cli::symbol$bullet,' Adding pseudocount to gene: ', i, '\n')))
#   #     reference@assays$RNA@data[i, sample(dim(reference)[2],1)] <- 2
#   #   } 
#   # } else {
#   #   message(crayon::magenta(paste0(cli::symbol$bullet,' All genes are variable!\n')))
#   # }
#   
#   # # # # rerun data processing to get scaled features
#   # reference <- reference %>%
#   #   NormalizeData() %>%
#   #   FindVariableFeatures(nfeatures = 5000) %>%
#   #   ScaleData() %>%
#   #   RunPCA()
#   
#   # # run in parallel use 4 tasks
#   # cl <- makePSOCKcluster(4)
#   # registerDoParallel(cl)
#   
#   # train the models - SVM
#   tic("Total Runtime") # start timer
#   tic("Feature space") # Subsection timer
#   reference <- getFeatureSpace(reference, "celltypes")
#   toc(log = TRUE) # feature space runtime
#   tic("Training Models") # model training timer
#   reference <- trainModel(reference,
#                           model = mod,
#                           allowParallel = TRUE, 
#                           seed = rnd_seed,
#                           resampleMethod = 'none',
#                           tuneLength = 1) 
#   toc(log = TRUE) # end model training
#   toc(log = TRUE) # total run time timer end
#   
#   stopCluster(cl)
#   
#   # log runtime
#   runtime <- unlist(tic.log())
#   runtime <- append(paste0('Random seed: ', rnd_seed), runtime)
#   runtime <- append('Parallel run: 4', runtime)
#   runtime <- append(paste0('Running scClassify - ', mod,': ', sample_no, '-', split), runtime)
#   runtime <- as.data.frame(runtime) 
#   colnames(runtime) <- 'time'
#   runtime <- runtime %>%
#     separate(time, into = c('Part', 'Runtime'), sep = ':') %>%
#     mutate(data_split = rep(split, nrow(.)))
#   
#   tic.clearlog()
#   
#   
#   # classification of query cells
#   query <- NormalizeData(query)
#   query <- scPredict(query, reference, seed = rnd_seed)
#   
#   query@meta.data$scpred_prediction <- as.factor(query@meta.data$scpred_prediction)
#   query@meta.data[[cluster_col]] <- factor(query@meta.data[[cluster_col]], levels(query$scpred_prediction))
#   
#   # model evaluation
#   cm <- caret::confusionMatrix(query$scpred_prediction, 
#                                reference = query$celltypes)
#   
#   model_metrics <- as.data.frame(cm[['byClass']]) %>%
#     rownames_to_column(var = 'Class') %>%
#     mutate(data_split = rep(split, nrow(.)))
#   
#   query_metadata <- query@meta.data %>%
#     mutate(data_split = rep(split, nrow(.)))
#   
#   # list to save values
#   return_list <- list(Runtime = runtime,
#                       Model_Eval = model_metrics,
#                       SplitMetadat = query_metadata)
#   
#   
#   return(return_list)
#   
# }
# 
# 
# 
# 
# all_splits <- list()
# 
# for(i in split_vect){
#   
#   repeat {
#     res_list <- try(runSplitscPred(seuratdata = scp1361, split = i))
#     if (!(inherits(res_list,"try-error"))) 
#       break
#   }
#   
#   all_splits[[i]] <- res_list
#   
# }
# 
# 
# runtime_list <- list()
# eval_list <- list()
# meta_list <- list()
# 
# for(i in split_vect){
#   runtime_list[[i]]<- all_splits[[i]][['Runtime']]
#   eval_list[[i]]<- all_splits[[i]][['Model_Eval']]
#   meta_list[[i]]<- all_splits[[i]][['SplitMetadat']]
# }
# 
# runtime_df <- as.data.frame(data.table::rbindlist(runtime_list, use.names = TRUE))
# eval_df <- as.data.frame(data.table::rbindlist(eval_list, use.names = TRUE))
# meta_df <- as.data.frame(data.table::rbindlist(meta_list, use.names = TRUE))
# 
# 
# write.csv(runtime_df,
#           file = paste0(workingdir,'/',resultsdir, '/',sample_no,'_runtime_NN.csv'),
#           row.names = TRUE)
# 
# write.csv(eval_df,
#           file = paste0(workingdir,'/',resultsdir, '/',sample_no,'_modeleval_NN.csv'),
#           row.names = TRUE)
# 
# write.csv(meta_df,
#           file = paste0(workingdir,'/',resultsdir, '/',sample_no,'_metadata_model_NN.csv'),
#           row.names = TRUE)
# 
