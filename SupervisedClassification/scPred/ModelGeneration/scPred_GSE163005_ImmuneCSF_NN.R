library(scPred)
library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(tictoc)
library(doParallel)


# Reading data ------------------------------------------------------------

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


# Prep scPred -------------------------------------------------------------

# load dataset
immunecsf <- LoadH5Seurat(file = '~/GSE163005_IMMUNE_CSF/GSE_ImmuneCSF_qc_hvg_anno_5k_raw.h5Seurat')

# read in metadata that contains split column
metadata <- read_csv(file = '~/DatasetSplits/Metadata_Splits/GSE_ImmuneCSF_metadata_splits.csv')

# clean up metadata
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$...1
metadata <- metadata[,-1]

# check to make sure barcodes are ordered properly in data and metadata
table(rownames(immunecsf@meta.data) == metadata$barcodes)

# append 5 split columns to the seurat metadata df
immunecsf@meta.data <- cbind(immunecsf@meta.data, metadata[,c("Split_1","Split_2", "Split_3", "Split_4", "Split_5")])

split_vect <- c("Split_1","Split_2", "Split_3", "Split_4", "Split_5")


scPredict2 <- scPred::scPredict
b <- lapply(body(scPredict2), function(a) {
  gsub("(rownames\\(.*\\) <- .*)(;?\\n?)",
       "message('Try to set the row names.'); try(\\1)\\2", deparse(a))
})
body(scPredict2) <- parse(text = c(do.call(c, b), "}"))
environment(scPredict2) <- environment(scPred::scPredict)
assignInNamespace("scPredict", scPredict2, ns = asNamespace("scPred"))




# a function to run through the scPred process
runSplitscPred <- function(seuratdata, split, cluster_col = 'general_celltypes', mod='nnet'){
  message(crayon::magenta(paste0("Running scPred model: ", mod, ' on split: ', split,'\n')))
  
  rnd_seed <- sample(c(1:4000),1)
  
  # split data into reference - 'train' and query - 'test'
  reference <- seuratdata[, seuratdata@meta.data[, split] == "train"]
  query <- seuratdata[, seuratdata@meta.data[, split] == "test"]
  
  # Find genes that do not vary - reference
  vf_r <- VariableFeatures(FindVariableFeatures(reference, nfeatures = 5000))
  nonvargenes_r <- setdiff(rownames(reference), vf_r)
  
  if(length(nonvargenes_r)>0){
    # add one pseudo count to allow scPred to keep non var genes
    for(i in nonvargenes_r){
      message(crayon::magenta(paste0(cli::symbol$bullet,' Adding pseudocount to gene: ', i, '\n')))
      reference@assays$RNA@data[i, sample(dim(reference)[2],1)] <- 2
    } 
  } else {
    message(crayon::magenta(paste0(cli::symbol$bullet,' All genes are variable!\n')))
  }

  
  # run in parallel use 4 tasks
  cl <- makePSOCKcluster(4)
  registerDoParallel(cl)
  
  # train the models - SVM
  tic("Total Runtime") # start timer
  tic("Feature space") # Subsection timer
  reference <- getFeatureSpace(reference, "general_celltypes")
  toc(log = TRUE) # feature space runtime
  tic("Training Models") # model training timer
  reference <- trainModel(reference,
                          model = mod,
                          allowParallel = TRUE, 
                          seed = rnd_seed,
                          resampleMethod = 'none',
                          tuneLength = 1) 
  toc(log = TRUE) # end model training
  toc(log = TRUE) # total run time timer end
  
  stopCluster(cl)
  
  # log runtime
  runtime <- unlist(tic.log())
  runtime <- append(paste0('Random seed: ', rnd_seed), runtime)
  runtime <- append('Parallel run: 4', runtime)
  runtime <- append(paste0('Running scPred - ', mod,': ', sample_no, '-', split), runtime)
  runtime <- as.data.frame(runtime) 
  colnames(runtime) <- 'time'
  runtime <- runtime %>%
    separate(time, into = c('Part', 'Runtime'), sep = ':') %>%
    mutate(data_split = rep(split, nrow(.)))
  
  tic.clearlog()
  
  
  # classification of query cells
  query <- NormalizeData(query)
  query <- scPredict(query, reference, seed = rnd_seed)
  
  query@meta.data$scpred_prediction <- as.factor(query@meta.data$scpred_prediction)
  query@meta.data[[cluster_col]] <- factor(query@meta.data[[cluster_col]], levels(query$scpred_prediction))
  
  # model evaluation
  cm <- caret::confusionMatrix(query$scpred_prediction, 
                               reference = query$general_celltypes)
  
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
    res_list <- try(runSplitscPred(seuratdata = immunecsf, split = i))
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
          file = paste0(workingdir,'/',resultsdir, '/',sample_no,'_runtime_NN.csv'),
          row.names = TRUE)

write.csv(eval_df,
          file = paste0(workingdir,'/',resultsdir, '/',sample_no,'_modeleval_NN.csv'),
          row.names = TRUE)

write.csv(meta_df,
          file = paste0(workingdir,'/',resultsdir, '/',sample_no,'_metadata_model_NN.csv'),
          row.names = TRUE)




