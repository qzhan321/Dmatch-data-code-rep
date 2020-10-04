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
# settings for preprocessing
batch_effects <- "medium"
projectName = "dmatch"
filter_genes = F
filter_cells = F
normData = F
Datascaling = T
min_cells = 10 #5
min_genes = 300
norm_method = "LogNormalize"
scale_factor = 10000
visualize = T
outfile_prefix = "broad_pbmc1_10x_v2_A_simulation"
save_obj = F
n_features = 5000
selection_method = "vst"


src_dir = "/home/qizhan/scratch-midway2/others/dmatch_review/things_done_by_me/july17th-july23rd2020/broad_PBMC_10x_simulate_batch_effects/scripts/"
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
    
    # preprocessing
    source(paste0(src_dir,'call_dmatch_seurat_preprocessing.R'))
    
    batches = dmatch_preprocess(
      expr_mat, metadata, projectName,
      filter_genes = filter_genes, filter_cells = filter_cells,
      normData = normData, Datascaling = Datascaling, 
      min_cells = min_cells, min_genes = min_genes, 
      norm_method = norm_method, scale_factor = scale_factor, 
      nfeatures = n_features, selection.method = selection_method)
    
    
    
    ########################
    # settings for running dmath
    umapplot_filename = "_dmatch_umap"
    dmatch_out_filename = "_dmatch_out"
    metadata_out_filename = "_dmatch_metadata_out"
    pca_filename = "_dmatch_pca"
    batch.id <- c(1,2)
    npcs = 30
    Reference_name <- "cell_atlas_ref_panel"
    #Reference_name <- "cell.line.to.ensl.ID"
    hclust_method <- "ward.D"
    TopCellLineNumber = 5
    ShowCellNumber = 10
    cor_threshold <- 0.05
    num_clusters <- 5
    cut_groups_methods <- "ward.D"
    select_clusters_quantile = 0.95
    anchors = c(1,2,3,4,5)
    run_alignment_quantile <- 0.95
    
    saveout_dir = paste0("/project2/mengjiechen/qizhan/qizhan/dmatch/dmatch_response_to_review/others/broad_PBMC_10x_batch_effects_simulation_and_correction/corrected_data/pbmc1_10x_v2_A/dmatch/", batch_effects, "/results/")
    dmatch_out_filename = "_samples"
    dmatch_out_pcs_filename = "_samples_pcs"
    plotout_dir <- paste0("/project2/mengjiechen/qizhan/qizhan/dmatch/dmatch_response_to_review/others/broad_PBMC_10x_batch_effects_simulation_and_correction/corrected_data/pbmc1_10x_v2_A/dmatch/", batch_effects, "/plots/")
    umapplot_filename <- "_umap_plot"
    heatmap_filename <- "_heatmap"
    pcaplot_filename <- "_pcaplot"
    outfilename_prefix <- file
    ########################
    #run
    source(paste0(src_dir,'call_dmatch_run_dmatch.R'))
    call_dmatch(batch_data = batches, npcs = npcs, 
                pcaplot_filename = pcaplot_filename,
                batch.id = batch.id, Reference_name = Reference_name,
                hclust_method = hclust_method,
                TopCellLineNumber = TopCellLineNumber,
                ShowCellNumber = ShowCellNumber,
                cor_threshold = cor_threshold,
                num_clusters = num_clusters,
                cut_groups_methods = cut_groups_methods,
                select_clusters_quantile = select_clusters_quantile,
                anchors = anchors,
                run_alignment_quantile = run_alignment_quantile,
                plotout_dir = plotout_dir, saveout_dir = saveout_dir, 
                outfilename_prefix = outfilename_prefix, 
                dmatch_out_filename = dmatch_out_filename,
                heatmap_filename = heatmap_filename,
                dmatch_out_pcs_filename = dmatch_out_pcs_filename,
                visualize = T, save_obj = T,
                umapplot_filename = umapplot_filename)
  }
  
}
