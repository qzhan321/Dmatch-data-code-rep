

########################
#load packages

#BiocManager::install("scMerge")
library(scMerge)
library(SingleCellExperiment)
library(scater)
library(scran)
library(cowplot)
library(umap)
library(stats)

rm(list=ls())

########################
#settings

filter_genes = F
filter_cells = F
min_cells = 5
min_genes = 300
cosineNorm = FALSE
LogNormalize = FALSE
kmeansK = c(3,3)
replicate_prop = 0.5
npcs=20
visualize = T
outfile_prefix = "broad_pbmc1_10x_v2_A_simulation"
save_obj = F

data("segList", package = "scMerge") 
seg = segList$human$human_scSEG

batch_effects <- "big"
src_dir = "/home/qizhan/scratch-midway2/others/dmatch_review/things_done_by_me/july17th-july23rd2020/broad_PBMC_10x_simulate_batch_effects/scmerge/"
working_dir = paste0("/project2/mengjiechen/qizhan/qizhan/dmatch/dmatch_response_to_review/others/broad_PBMC_10x_batch_effects_simulation_and_correction/corrected_data/pbmc1_10x_v2_A/", batch_effects, "/")
read_dir = paste0("/project2/mengjiechen/qizhan/qizhan/dmatch/dmatch_response_to_review/others/broad_PBMC_10x_batch_effects_simulation_and_correction/simulated_data/pbmc1_10x_v2_A/", batch_effects, "/")


ratio_names <- c("1_6", "1_3", "1_2")
for (i in 1:10) {
  for (j in 1:3) {
    ratio_name <- ratio_names[j]
    file <- paste0("pbmc1_10x_v2_A_", ratio_name, "_rep", i)
    print(file)
    load(paste0(read_dir, file)) 
    ########################
    # load data 
    
    sample1 <- sim_list[[1]]
    sample2 <- sim_list[[2]]
    meta1 <- sim_list[[3]]
    meta2 <- sim_list[[4]]
    
    expr_mat <- cbind(sample1, sample2)
    rownames(expr_mat) <- toupper(rownames(expr_mat))
    batch <- c(rep(1, ncol(sample1)), rep(2, ncol(sample2)))
    celltype <- c(meta1, meta2)
    metadata <- data.frame("batch"=batch, "celltype"=celltype)
    rownames(metadata) <- colnames(expr_mat)
    
    
    
    
    ########################
    # settings for running scmerge
    umapplot_filename = "_scmerge_umap"
    dmatch_out_filename = "_scmerge_out"
    metadata_out_filename = "_scmerge_metadata_out"
    pca_filename = "_scmerge_pca"
    
    batch_label = "batch"
    celltype_label = "celltype"
    npcs = 30
    
    saveout_dir = paste0("/project2/mengjiechen/qizhan/qizhan/dmatch/dmatch_response_to_review/others/broad_PBMC_10x_batch_effects_simulation_and_correction/corrected_data/pbmc1_10x_v2_A/scmerge/", batch_effects, "/results/")
    dmatch_out_filename = "_samples"
    dmatch_out_pcs_filename = "_samples_pcs"
    plotout_dir <- paste0("/project2/mengjiechen/qizhan/qizhan/dmatch/dmatch_response_to_review/others/broad_PBMC_10x_batch_effects_simulation_and_correction/corrected_data/pbmc1_10x_v2_A/scmerge/", batch_effects, "/plots/")
    umapplot_filename <- "_umap_plot"
    heatmap_filename <- "_heatmap"
    pcaplot_filename <- "_pcaplot"
    outfilename_prefix <- file
    visualize = T
    save_obj = F
    
    ########################
    #run
    source(paste0(src_dir,'call_scmerge.R'))
    
    scMerge_obj <- scMerge_preprocess(expr_mat, metadata, 
                                      batch_label = batch_label,
                                      filter_genes = filter_genes, filter_cells = filter_cells,
                                      min_cells = min_cells, min_genes = min_genes, 
                                      cosineNorm = cosineNorm, LogNormalize = LogNormalize)
    
    call_scMerge(scMerge_obj, batch_label, celltype_label, seg = seg, 
                 kmeansK = kmeansK, replicate_prop=replicate_prop, npcs=npcs,
                 plotout_dir = plotout_dir, saveout_dir = saveout_dir, outfilename_prefix = outfilename_prefix, visualize = visualize, save_obj = save_obj)

    
  }
  
}











