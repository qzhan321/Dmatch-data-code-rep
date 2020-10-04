
#### this is for partial overlapping: two cases; 
#### share six out of nine and share four out of nine celltypes; 
#### and two, for the two, share the most common two celltypes

# Author : Qi Zhan
# Date : 29/06/2020
# Proj : Run dmatch pipeline

########################
#load packages

library(Seurat, lib.loc = "/software/R-3.6.1-el7-x86_64/lib64/R/library")  # Seurat 3 version
library(umap)
library(dmatch)
library(ggplot2)
library(cowplot)

rm(list=ls())

########################

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


simulate_func_partial_ovelapping <- function(meta, sample, ratio, mean, sd, shareNum) {
  data1 <- NULL
  data2 <- NULL
  #randomly select shareTypes
  shareType <- sample(unique(meta[,1]), shareNum, replace = F)
  for (i in 1:length(shareType)) {
    celltype <- shareType[i]
    subset <- sample[, which(meta[,1] == celltype)]
    index <- sample(1:ncol(subset), round(ratio*ncol(subset)), replace = F)
    temp1 <- subset[, index]
    data1 <- cbind(data1, temp1)
    temp2 <- subset[, -index]
    data2 <- cbind(data2, temp2)
  }
  
  #for the rest celltype, distribute them to the two datsets randomly. For example
  restType <- setdiff(unique(meta[,1]), shareType)
  restType1 <- sample(restType, ceiling(length(restType)/2), replace = F)
  restType2 <- setdiff(restType, restType1)
  
  for (j in 1:length(restType1)) {
    celltype <- restType1[j]
    subset <- sample[, which(meta[,1] == celltype)]
    data1 <- cbind(data1, subset)
  }
  
  for (k in 1:length(restType2)) {
    celltype <- restType2[k]
    subset <- sample[, which(meta[,1] == celltype)]
    data2 <- cbind(data2, subset)
  }
  
  cells1 <- colnames(data1)
  celltypes1 <- meta[cells1,]
  
  cells2 <- colnames(data2)
  celltypes2 <- meta[cells2,]
  
  data2 <- data2 + rnorm(nrow(data2)*ncol(data2), mean = mean, sd = sd)  
  return(list("data1" = data1, "data2" = data2, "meta1" = celltypes1, "meta2" = celltypes2))
}




simulate_func_partial_ovelapping_two <- function(meta, sample, ratio, mean, sd) {
  data1 <- NULL
  data2 <- NULL
  #select the two most common ones
  shareType <- names(sort(table(meta[,1]))[8:9])
  for (i in 1:length(shareType)) {
    celltype <- shareType[i]
    subset <- sample[, which(meta[,1] == celltype)]
    index <- sample(1:ncol(subset), round(ratio*ncol(subset)), replace = F)
    temp1 <- subset[, index]
    data1 <- cbind(data1, temp1)
    temp2 <- subset[, -index]
    data2 <- cbind(data2, temp2)
  }
  
  #for the rest celltype, distribute them to the two datsets randomly. For example
  restType <- setdiff(unique(meta[,1]), shareType)
  restType1 <- sample(restType, ceiling(length(restType)/2), replace = F)
  restType2 <- setdiff(restType, restType1)
  
  for (j in 1:length(restType1)) {
    celltype <- restType1[j]
    subset <- sample[, which(meta[,1] == celltype)]
    data1 <- cbind(data1, subset)
  }
  
  for (k in 1:length(restType2)) {
    celltype <- restType2[k]
    subset <- sample[, which(meta[,1] == celltype)]
    data2 <- cbind(data2, subset)
  }
  
  cells1 <- colnames(data1)
  celltypes1 <- meta[cells1,]
  
  cells2 <- colnames(data2)
  celltypes2 <- meta[cells2,]
  
  data2 <- data2 + rnorm(nrow(data2)*ncol(data2), mean = mean, sd = sd)  
  return(list("data1" = data1, "data2" = data2, "meta1" = celltypes1, "meta2" = celltypes2))
}




#### sample1: 1:5, 2:4, 3:3


#### sample1: 1:5, 2:4, 3:3
#
batch_effects <- "medium"
batch_mean <- 0.18
batch_sd <- 1


outfilename <- "pbmc1_10x_v2_A"
saveout_dir <- paste0("/project2/mengjiechen/qizhan/qizhan/dmatch/dmatch_response_to_review/others/broad_PBMC_10x_batch_effects_simulation_and_correction/simulated_data_partial_overlapping/pbmc1_10x_v2_A/") 
meta1 <- as.data.frame(meta1)
rownames(meta1) <- colnames(sample1)
colnames(meta1) <- "celltype"

rep_num <- 10
for (i in 1:rep_num) {
  ratio <- 1/6
  ratio_name <- "1_6"
  sim_list <- simulate_func_partial_ovelapping(meta1, sample1, ratio, mean = batch_mean, sd = batch_sd, shareNum = 6)
  save(sim_list, file = paste0(saveout_dir, "/shareSix/", batch_effects, "/", outfilename, "_", ratio_name, "_rep", i))
  sim_list <- simulate_func_partial_ovelapping(meta1, sample1, ratio, mean = batch_mean, sd = batch_sd, shareNum = 4)
  save(sim_list, file = paste0(saveout_dir, "/shareFour/", batch_effects, "/", outfilename, "_", ratio_name, "_rep", i))
  sim_list <- simulate_func_partial_ovelapping_two(meta1, sample1, ratio, mean = batch_mean, sd = batch_sd)
  save(sim_list, file = paste0(saveout_dir, "/shareTwo/", batch_effects, "/", outfilename, "_", ratio_name, "_rep", i))
  
  
  ratio <- 1/3
  ratio_name <- "1_3"
  sim_list <- simulate_func_partial_ovelapping(meta1, sample1, ratio, mean = batch_mean, sd = batch_sd, shareNum = 6)
  save(sim_list, file = paste0(saveout_dir, "/shareSix/", batch_effects, "/", outfilename, "_", ratio_name, "_rep", i))
  sim_list <- simulate_func_partial_ovelapping(meta1, sample1, ratio, mean = batch_mean, sd = batch_sd, shareNum = 4)
  save(sim_list, file = paste0(saveout_dir, "/shareFour/", batch_effects, "/", outfilename, "_", ratio_name, "_rep", i))
  sim_list <- simulate_func_partial_ovelapping_two(meta1, sample1, ratio, mean = batch_mean, sd = batch_sd)
  save(sim_list, file = paste0(saveout_dir, "/shareTwo/", batch_effects, "/", outfilename, "_", ratio_name, "_rep", i))
  
  
  ratio <- 1/2
  ratio_name <- "1_2"
  sim_list <- simulate_func_partial_ovelapping(meta1, sample1, ratio, mean = batch_mean, sd = batch_sd, shareNum = 6)
  save(sim_list, file = paste0(saveout_dir, "/shareSix/", batch_effects, "/", outfilename, "_", ratio_name, "_rep", i))
  sim_list <- simulate_func_partial_ovelapping(meta1, sample1, ratio, mean = batch_mean, sd = batch_sd, shareNum = 4)
  save(sim_list, file = paste0(saveout_dir, "/shareFour/", batch_effects, "/", outfilename, "_", ratio_name, "_rep", i))
  sim_list <- simulate_func_partial_ovelapping_two(meta1, sample1, ratio, mean = batch_mean, sd = batch_sd)
  save(sim_list, file = paste0(saveout_dir, "/shareTwo/", batch_effects, "/", outfilename, "_", ratio_name, "_rep", i))
}


