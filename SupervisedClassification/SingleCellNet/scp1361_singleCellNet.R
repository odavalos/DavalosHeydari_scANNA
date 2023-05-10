library(singleCellNet)
library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(tictoc)
# library(doParallel)


# Reading data ------------------------------------------------------------

sample_no <- 'SCP1361' # change accordingly
tool_name <- 'singleCellNet'

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


# # a function to run through the scClassify process
runSplitSCN <- function(seuratdata, split, cluster_col = 'celltypes', dataset = sample_no){
  message(crayon::magenta(paste0("Running singleCellNet model on: ", dataset, ' on split: ', split,'\n')))
  
  gc(verbose = FALSE)
  
  # Split the dataset for training the model
  reference <- seuratdata[, seuratdata@meta.data[, split] == "train"]
  query <- seuratdata[, seuratdata@meta.data[, split] == "test"]
  
  
  #exp_type options can be: counts, normcounts, and logcounts, if they are available in your sce object
  seuratfile_ref = extractSeurat(reference, exp_slot_name = "counts")
  sampTab_ref = seuratfile_ref$sampTab
  expDat_ref = seuratfile_ref$expDat
  # sampTab_ref <- droplevels(sampTab_ref)
  
  seuratfile_query = extractSeurat(query, exp_slot_name = "counts")
  sampTab_query = seuratfile_query$sampTab
  expDat_query = seuratfile_query$expDat
  
  
  # # Split for training and assessment, and transform training data
  set.seed(2023) #can be any random seed number
  stList = splitCommon(sampTab=sampTab_ref, ncells=100, dLevel=cluster_col)
  stTrain = stList[[1]]
  expTrain = expDat_ref[,rownames(stTrain)]
  
  
  
  # train the models
  tic("Total Runtime") # start timer
  tic("Training Models") # model training timer
  
  class_info <- scn_train(stTrain = stTrain, 
                          expTrain = expTrain, 
                          nTopGenes = 10, 
                          nRand = 70, 
                          nTrees = 1000, 
                          nTopGenePairs = 25, 
                          dLevel = cluster_col, 
                          colName_samp = "barcodes")
  
  
  #validate data
  stTestList = splitCommon(sampTab=stList[[2]], ncells=100, dLevel=cluster_col) #normalize validation data so that the assessment is as fair as possible
  stTest = stTestList[[1]]
  expTest = expDat_ref[,rownames(stTest)]
  
  #predict
  classRes_val_all = scn_predict(cnProc=class_info[['cnProc']], expDat=expTest, nrand = 50)

  toc(log = TRUE) # end model training
  toc(log = TRUE) # total run time timer end
  
  # log runtime
  runtime <- unlist(tic.log())
  runtime <- append(paste0('Running singleCellNet - ', sample_no, '-', split), runtime)
  runtime <- as.data.frame(runtime)
  colnames(runtime) <- 'time'
  # runtime <- runtime %>%
  #   separate(time, into = c('Part', 'Runtime'), sep = ':') %>%
  #   mutate(data_split = rep(split, nrow(.)))
  
  tic.clearlog()
  
  preds <- scn_predict(class_info[['cnProc']], expDat_query, nrand=50)
  
  
  # This classifies a cell with  the catgory with the highest classification score or higher than a classification score threshold of your choosing.
  # The annotation result can be found in a column named category in the query sample table.
  
  pred_metadata <- get_cate(classRes = preds, 
                            sampTab = sampTab_query, 
                            dLevel = cluster_col, 
                            sid = "barcodes", 
                            nrand = 50)
  
  # change column names
  colnames(pred_metadata)[colnames(pred_metadata) == "category"] = "scn_preds"
  
  preds <- pred_metadata$scn_preds
  scores <- pred_metadata$scn_score
  
  # # remove the intermediate predictions and call them unassigned
  # preds <- ifelse(test = (!preds %in% unique(query@meta.data[[cluster_col]])), 
  #                 yes = 'unassigned', 
  #                 no = preds)
  
  # add the prediction back to query data and convert into a level
  query@meta.data$singleCellNet_scores <- scores
  query@meta.data$singleCellNet_pred <- as.factor(preds)
  query@meta.data[[cluster_col]] <- factor(query@meta.data[[cluster_col]], 
                                           levels(query$singleCellNet_pred))
  
  # model evaluation
  cm <- caret::confusionMatrix(query$singleCellNet_pred,
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
  
  
  res_list <- try(runSplitSCN(seuratdata = scp1361, split = i))
  # repeat {
  #   res_list <- try(runSplitSCN(seuratdata = scp1361, split = i))
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

