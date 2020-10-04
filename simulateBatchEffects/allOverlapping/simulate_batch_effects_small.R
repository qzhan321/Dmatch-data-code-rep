# Author : Qi Zhan
# Date : 29/06/2020
# Proj : Run dmatch pipeline

########################
#load packages

library(Seurat)  # Seurat 3 version
library(umap)
library(dmatch)
library(ggplot2)
library(cowplot)

rm(list=ls())

########################

src_dir = "/home/qizhan/scratch-midway2/others/dmatch_review/things_done_by_me/July8th-July15th2020/scripts/broad_PBMC/"
working_dir = "/project2/mengjiechen/qizhan/qizhan/dmatch/dmatch_response_to_review/others/broad_PBMC_batch_effects_correction/dmatch/"
read_dir = "/project2/mengjiechen/qizhan/qizhan/dmatch/dmatch_response_to_review/others/potential_dataset/broad_PBMC_datasets/"


########################
# load data 
pbmc_dir <- "pbmc1"
protocol <- "_10x"
outfilename_prefix <- substr(protocol, 2, nchar(protocol))
load(paste0(read_dir,pbmc_dir, "/", pbmc_dir, protocol, "_v2_A"))

sample1 <- pbmc1_10x_v2_A
load(paste0(read_dir, pbmc_dir, "/", pbmc_dir, protocol, "_v2_A", "_celltype_information"))
meta1 <- celltypes

load(paste0(read_dir,pbmc_dir, "/", pbmc_dir, protocol, "_v2_B"))
sample2 <- pbmc1_10x_v2_B
load(paste0(read_dir, pbmc_dir, "/", pbmc_dir, protocol, "_v2_B", "_celltype_information"))
meta2 <- celltypes

load(paste0(read_dir,pbmc_dir, "/", pbmc_dir, protocol, "_v3"))
sample3 <- pbmc1_10x_v3
load(paste0(read_dir, pbmc_dir, "/", pbmc_dir, protocol, "_v3", "_celltype_information"))
meta3 <- celltypes

pbmc_dir <- "pbmc2"
load(paste0(read_dir,pbmc_dir, "/", pbmc_dir, protocol, "_v2"))
sample4 <- pbmc2_10X_v2
load(paste0(read_dir, pbmc_dir, "/", pbmc_dir, protocol, "_v2_celltype_information"))
meta4 <- celltypes





simulate_func <- function(meta, sample, ratio, mean, sd) {
  data1 <- NULL
  data2 <- NULL
  #set.seed(0)
  for (i in 1:length(unique(meta[,1]))) {
    celltype <- unique(meta[,1])[i]
    subset <- sample[, which(meta[,1] == celltype)]
    index <- sample(1:ncol(subset), round(ratio*ncol(subset)), replace = F)
    temp1 <- subset[, index]
    data1 <- cbind(data1, temp1)
    temp2 <- subset[, -index]
    data2 <- cbind(data2, temp2)
  }
  cells1 <- colnames(data1)
  celltypes1 <- meta[cells1,]
  
  cells2 <- colnames(data2)
  celltypes2 <- meta[cells2,]
  
  data2 <- data2 + rnorm(nrow(data2)*ncol(data2), mean = mean, sd = sd)  
  return(list("data1" = data1, "data2" = data2, "meta1" = celltypes1, "meta2" = celltypes2))
}

#### sample1: 1:5, 2:4, 3:3
#
batch_effects <- "small"
batch_mean <- 0.08
batch_sd <- 1

outfilename <- "pbmc1_10x_v2_A"
saveout_dir <- paste0("/project2/mengjiechen/qizhan/qizhan/dmatch/dmatch_response_to_review/others/broad_PBMC_10x_batch_effects_simulation_and_correction/simulated_data/pbmc1_10x_v2_A/", batch_effects, "/") 
meta1 <- as.data.frame(meta1)
rownames(meta1) <- colnames(sample1)
colnames(meta1) <- "celltype"

rep_num <- 10
for (i in 1:rep_num) {
  ratio <- 1/6
  ratio_name <- "1_6"
  sim_list <- simulate_func(meta1, sample1, ratio, mean = batch_mean, sd = batch_sd)
  save(sim_list, file = paste0(saveout_dir, outfilename, "_", ratio_name, "_rep", i))
  
  
  ratio <- 1/3
  ratio_name <- "1_3"
  sim_list <- simulate_func(meta1, sample1, ratio, mean = batch_mean, sd = batch_sd)
  save(sim_list, file = paste0(saveout_dir, outfilename, "_", ratio_name, "_rep", i))
  
  
  ratio <- 1/2
  ratio_name <- "1_2"
  sim_list <- simulate_func(meta1, sample1, ratio, mean = batch_mean, sd = batch_sd)
  save(sim_list, file = paste0(saveout_dir, outfilename, "_", ratio_name, "_rep", i))
}


