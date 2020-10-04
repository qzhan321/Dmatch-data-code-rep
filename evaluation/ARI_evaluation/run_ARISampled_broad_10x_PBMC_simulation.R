
# Purpose: Primary script to run ARI pipeline

# clear workspace
rm(list=ls())

# source relevant functions
src_dir <- "/home/qizhan/scratch-midway2/others/dmatch_review/things_done_by_me/july17th-july23rd2020/evaluation/ARI_evaluation/"
source(paste0(src_dir, "/ARI_utils/run_ARISampled.R"))
source(paste0(src_dir, "/ARI_utils/ari_calcul_sampled.R"))
source(paste0(src_dir, "/ARI_utils/conclude_ARISampled.R"))

##############################
############ Set arguments

eval_metric <- 'ARI'

# read in files from this_dir
data_dir<-"/project2/mengjiechen/qizhan/qizhan/dmatch/dmatch_response_to_review/others/broad_PBMC_10x_batch_effects_simulation_and_correction/corrected_data/pbmc1_10x_v2_A/"  

# send output to out_dir
out_dir<-"/project2/mengjiechen/qizhan/qizhan/dmatch/dmatch_response_to_review/others/evaluation/broad_PBMC_10x_simulation/"

# create relevant folder within the out_dir
dir.create(paste0(out_dir, eval_metric),showWarnings = FALSE)

##############################
############ Run 'run_ARISampled' function on relevant files
batch_effects <- c("big", "medium", "small")
ratio_names <- c("1_6", "1_3", "1_2")
reps <- 1:10
cpcs <- 1:20
nbiters <- 100
percent_extract <- 0.8
maxiter <- 100
celltypelb='celltype'
batchlb='batch'

for (batch_effect in batch_effects) {
  meta_dir <- paste0("/project2/mengjiechen/qizhan/qizhan/dmatch/dmatch_response_to_review/others/broad_PBMC_10x_batch_effects_simulation_and_correction/simulated_data/pbmc1_10x_v2_A/", batch_effect, "/")
  for (ratio_name in ratio_names) {
    for (rep in reps) {
      methods_use <- 'fastmnn'
      txt_out_dir <- paste0(out_dir, eval_metric, "/", methods_use, "/")
      setwd(paste0(data_dir, methods_use, "/", batch_effect, "/results/"))
      files <- list.files()
      index1 <- grep(ratio_name, files)
      files1 <- files[index1]
      index2 <- grep(paste0("rep",rep, "_"),files1) 
      files2 <- files1[index2]
      index3 <- grep("_pc", files2)
      fn <- files2[index3]
      this_dir <- paste0(data_dir, methods_use, "/", batch_effect, "/")
      #load(paste0(data_dir, batch_effect, "/", filename_prefix, ratio_name, "_rep", rep))
      Rfastmnn<-run_ARISampled(fn,this_dir,out_dir,meta_dir,eval_metric,methods_use, cpcs = cpcs, 
                               batch_effect = batch_effect, ratio_name = ratio_name, rep = rep, 
                               raw = F, celltypelb=celltypelb, batchlb=batchlb, nbiters = nbiters,
                               txt_out_dir = txt_out_dir, maxiter = maxiter)
      
      
      methods_use <- 'harmony'
      txt_out_dir <- paste0(out_dir, eval_metric, "/", methods_use, "/")
      setwd(paste0(data_dir, methods_use, "/", batch_effect, "/results/"))
      files <- list.files()
      index1 <- grep(ratio_name, files)
      files1 <- files[index1]
      index2 <- grep(paste0("rep",rep, "_"),files1) 
      files2 <- files1[index2]
      index3 <- grep("_pc", files2)
      fn <- files2[index3]
      this_dir <- paste0(data_dir, methods_use, "/", batch_effect, "/")
      #load(paste0(data_dir, batch_effect, "/", filename_prefix, ratio_name, "_rep", rep))
      Rharmony<-run_ARISampled(fn,this_dir,out_dir,meta_dir,eval_metric,methods_use, cpcs = cpcs, 
                               batch_effect = batch_effect, ratio_name = ratio_name, rep = rep, 
                               raw = F, celltypelb=celltypelb, batchlb=batchlb, nbiters = nbiters,
                               txt_out_dir = txt_out_dir, maxiter = maxiter)
      
      
      methods_use <- 'dmatch'
      txt_out_dir <- paste0(out_dir, eval_metric, "/", methods_use, "/")
      setwd(paste0(data_dir, methods_use, "/", batch_effect, "/results/"))
      files <- list.files()
      index1 <- grep(ratio_name, files)
      files1 <- files[index1]
      index2 <- grep(paste0("rep",rep, "_"),files1) 
      files2 <- files1[index2]
      index3 <- grep("_pc", files2)
      fn <- files2[index3]
      this_dir <- paste0(data_dir, methods_use, "/", batch_effect, "/")
      #load(paste0(data_dir, batch_effect, "/", filename_prefix, ratio_name, "_rep", rep))
      Rdmatch<-run_ARISampled(fn,this_dir,out_dir,meta_dir,eval_metric,methods_use, cpcs = cpcs, 
                              batch_effect = batch_effect, ratio_name = ratio_name, rep = rep, 
                              raw = F, celltypelb=celltypelb, batchlb=batchlb, nbiters = nbiters,
                              txt_out_dir = txt_out_dir, maxiter = maxiter)
      
      
      methods_use <- 'raw'
      txt_out_dir <- paste0(out_dir, eval_metric, "/", methods_use, "/")
      fn <- paste0("pbmc1_10x_v2_A_", ratio_name, "_", "rep", rep)
      Rraw<-run_ARISampled(fn,this_dir,out_dir,meta_dir,eval_metric,methods_use, cpcs = cpcs, 
                           batch_effect = batch_effect, ratio_name = ratio_name, rep = rep, 
                           raw = T, celltypelb=celltypelb, batchlb=batchlb, nbiters = nbiters,
                           txt_out_dir = txt_out_dir, maxiter = maxiter)
      
      
      methods_use <- 'scmerge'
      txt_out_dir <- paste0(out_dir, eval_metric, "/", methods_use, "/")
      setwd(paste0(data_dir, methods_use, "/", batch_effect, "/results/"))
      files <- list.files()
      index1 <- grep(ratio_name, files)
      files1 <- files[index1]
      index2 <- grep(paste0("rep",rep, "_"),files1) 
      files2 <- files1[index2]
      index3 <- grep("_pc", files2)
      fn <- files2[index3]
      this_dir <- paste0(data_dir, "/", methods_use, "/", batch_effect, "/")
      #load(paste0(data_dir, batch_effect, "/", filename_prefix, ratio_name, "_rep", rep))
      Rscmerge<-run_ARISampled(fn,this_dir,out_dir,meta_dir,eval_metric,methods_use, cpcs = cpcs, 
                               batch_effect = batch_effect, ratio_name = ratio_name, rep = rep, 
                               raw = F, celltypelb=celltypelb, batchlb=batchlb, nbiters = nbiters,
                               txt_out_dir = txt_out_dir, maxiter = maxiter)
      
      
      
      methods_use <- 'seurat3'
      txt_out_dir <- paste0(out_dir, eval_metric, "/", methods_use, "/")
      setwd(paste0(data_dir, methods_use, "/", batch_effect, "/results/"))
      files <- list.files()
      index1 <- grep(ratio_name, files)
      files1 <- files[index1]
      index2 <- grep(paste0("rep",rep, "_"),files1) 
      files2 <- files1[index2]
      index3 <- grep("_pc", files2)
      fn <- files2[index3]
      this_dir <- paste0(data_dir, methods_use, "/", batch_effect, "/")
      #load(paste0(data_dir, batch_effect, "/", filename_prefix, ratio_name, "_rep", rep))
      Rseurat3<-run_ARISampled(fn,this_dir,out_dir,meta_dir,eval_metric,methods_use, cpcs = cpcs, 
                               batch_effect = batch_effect, ratio_name = ratio_name, rep = rep, 
                               raw = F, celltypelb=celltypelb, batchlb=batchlb, nbiters = nbiters,
                               txt_out_dir = txt_out_dir, maxiter = maxiter)
    }
  } 
}

####################################################################

##############################
############ Extracting all data from all methods in dataset 

# Reads files from dir_this 
dir_this<-paste0(out_dir, eval_metric, "/")
#dir_this<-paste0(out_dir, eval_metric, "_OP")
folders <- list.files(dir_this)
# Perform the following function to produce final CSV file
for (i in 1:length(folders)) {
  folder <- folders[i]
  wholedf<-conclude_ARISampled(paste0(dir_this, folder), paste0(folder, "_broad_PBMC_10x_simulation"), nbiters = nbiters) 
  #save.image(paste0(folder, "_broad_PBMC_10x_simulation", "_complete.RData"))
}

####################################################################

##############################
############ Consolidate all raw data into one file

# setwd(dir_this)
# rm(list=ls())
# 
# source("./ARI_utils/ARI_files_consolidate.R")
#  
# ari_consolidate()

### END
