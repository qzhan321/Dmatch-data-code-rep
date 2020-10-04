
# Author : Kok Siong Ang 
# Date : 18/09/2019
# Proj : fast MNN pipeline

#' fast MNN pipeline
#' 



#' dmatch pipeline
#' 
#' @param expr_mat_dmatch : expression matrix for seurat object, genes by cells
#' @param metadata : dataframe with metadata for seurat object
#' @param projectName : name of the seurat project
#' @param filter_genes : logical for gene filtering (default TRUE)
#' @param filter_cells : logical for cell filtering (default TRUE)
#' @param normData : logical for expression matrix normalization (default TRUE)
#' @param Datascaling : logical for expression matrix scaling (default TRUE)
#' @param min_cells : min cells for gene filtering (default 10)
#' @param min_genes : min genes for cell filtering (default 300)
#' @param scale_factor :  scale factor for normalization (default 10000)


mnn_preprocess <- function(expr_mat_mnn, metadata, projectName,
                              filter_genes = T, filter_cells = T,
                              normData = T, Datascaling = T, regressUMI = F, 
                              nfeatures = 5000, selection.method = "vst",
                              min_cells = 10, min_genes = 300, norm_method = "LogNormalize", 
                              scale_factor = 10000)
{
  
  ##########################################################
  # preprocessing
  
  if(filter_genes == F) {
    min_cells = 0
  }
  if(filter_cells == F) {
    min_genes = 0
  }
  
  b_seurat <- CreateSeuratObject(counts = expr_mat_mnn, meta.data = metadata, project = projectName, 
                                 min.cells = min_cells, min.features = min_genes)
  if (normData) {
    b_seurat <- NormalizeData(object = b_seurat, normalization.method = norm_method, scale.factor = scale_factor, verbose = F)
  } else {
    b_seurat@assays$RNA@data = b_seurat@assays$RNA@counts
  }
  b_seurat <- ScaleData(object = b_seurat, verbose = F)
  b_seurat <- FindVariableFeatures(object = b_seurat, selection.method = selection_method, nfeatures = n_features, verbose = F)
  
  # Get the list of highly variable gens in 2 batches
  hvg_all <- b_seurat@assays$RNA@var.features
  
  #split data into two batches
  meta_data <- b_seurat@meta.data
  
  b1_selected <- row.names(meta_data)[which (meta_data$batch==1)]
  b2_selected <- row.names(meta_data)[which (meta_data$batch==2)]
  b1_meta <- meta_data[which (meta_data$batch==1),]
  b2_meta <- meta_data[which (meta_data$batch==2),]
  
  data_filtered <- b_seurat@assays$RNA@scale.data[hvg_all,]
  dim(data_filtered)
  b1_data <- data_filtered[,b1_selected]
  b2_data <- data_filtered[,b2_selected]
  b1_meta <- b1_meta[b1_selected,]
  b2_meta <- b2_meta[b2_selected,]
  
  return(list( b1_data, b2_data, b1_meta, b2_meta))
}





call_mnn <- function(batch_data, batch_label, celltype_label, npcs = 30, 
                         plotout_dir = "", saveout_dir = "", 
                         outfilename_prefix = "", 
                         visualize = T, save_obj = T)
{

  #plotting
  k_seed = 10
  # tsne_perplex = 30
  # tsneplot_filename = "_fastmnn_tsne"
  mnn_out_filename = "_fastmnn_out"
  metadata_out_filename = "_fastmnn_metadata_out"
  pca_filename = "_fastmnn_pca"

  b1_data <- batch_data[[1]]
  b2_data <- batch_data[[2]]
  b1_meta <- batch_data[[3]]
  b2_meta <- batch_data[[4]]

  metadata <- rbind(b1_meta, b2_meta)

  ##########################################################
  #run

  t1 = Sys.time()

  cosd1 <- cosineNorm(as.matrix(b1_data))
  cosd2 <- cosineNorm(as.matrix(b2_data))

  pcs_total <- batchelor::multiBatchPCA(cosd1, cosd2, d = 30)  
  out_mnn_total <- batchelor::fastMNN(pcs_total[[1]], pcs_total[[2]], pc.input = TRUE)

  t2 = Sys.time()
  print(t2-t1)

  ##########################################################
  #save

  saveRDS(out_mnn_total, file=paste0(saveout_dir, outfilename_prefix, mnn_out_filename, ".RDS"))

  corrected_pca_mnn <- as.data.frame(out_mnn_total$corrected)
  rownames(corrected_pca_mnn) <- colnames(cbind(b1_data,b2_data))
  cells_use <- rownames(corrected_pca_mnn)
  metadata <- metadata[cells_use,]

  corrected_pca_mnn[rownames(metadata), batch_label] <- metadata[ , batch_label]
  #corrected_pca_mnn[rownames(metadata), 'batch'] <- metadata[ , "batch"]
  corrected_pca_mnn[rownames(metadata), celltype_label] <- metadata[ , celltype_label]
  write.table(corrected_pca_mnn, file=paste0(saveout_dir, outfilename_prefix, pca_filename, ".txt"), quote=F, sep='\t', row.names = T, col.names = NA)

  saveRDS(metadata, file = paste0(saveout_dir, outfilename_prefix, metadata_out_filename, ".RDS"))
  #write.table(metadata, file=paste0(saveout_dir, outfilename_prefix, metadata_out_filename, ".txt"), quote=F, sep="\t", row.name=TRUE, col.name=NA)

  if (visualize) {
    
    set.seed(10)
    out_umap <- umap(corrected_pca_mnn[,1:npcs])
    umap_df<- as.data.frame(out_umap$layout)
    
    rownames(umap_df) <- rownames(corrected_pca_mnn)
    dim(umap_df)
    colnames(umap_df) <- c('umap_1', 'umap_2')
    umap_df$batchlb <- as.factor(metadata$batch)     
    umap_df$CellType <- as.factor(metadata$celltype)
    
    p01 <- ggplot(umap_df, aes(x = umap_1, y = umap_2, colour = batchlb)) + geom_point(alpha = 0.6) + theme_bw()
    p01 <- p01 + labs(x='umap_1',y='umap_2',title='fastmnn')
    p01 <- p01 + theme(legend.title = element_text(size=17), 
                       legend.key.size = unit(1.1, "cm"),
                       legend.key.width = unit(0.5,"cm"), 
                       legend.text = element_text(size=14), 
                       plot.title = element_text(color="black", size=20, hjust = 0.5))
    
    p02 <- ggplot(umap_df, aes(x = umap_1, y = umap_2, colour = CellType)) + geom_point(alpha = 0.6) + theme_bw()
    p02 <- p02 + labs(x='umap_1',y='umap_2',title='fastmnn')
    p02 <- p02 + theme(legend.title = element_text(size=17), 
                       legend.key.size = unit(1.1, "cm"),
                       legend.key.width = unit(0.5,"cm"), 
                       legend.text = element_text(size=14), 
                       plot.title = element_text(color="black", size=20, hjust = 0.5))
    
    png(paste0(plotout_dir,outfilename_prefix,umapplot_filename,".png"),width = 2*1000, height = 800, res = 2*72)
    print(plot_grid(p01, p02))
    dev.off()
    
    pdf(paste0(plotout_dir,outfilename_prefix,umapplot_filename,".pdf"),width=15,height=7,paper='special')
    print(plot_grid(p01, p02))
    dev.off()
    
    # set.seed(10)
    # out_tsne <- Rtsne(out_mnn_total$corrected, check_duplicates=FALSE, perplexity = 30)
    # tsne_df<- as.data.frame(out_tsne$Y)
    # 
    # rownames(tsne_df) <- colnames(cbind(b1_data,b2_data))
    # dim(tsne_df)
    # colnames(tsne_df) <- c('tSNE_1', 'tSNE_2')
    # tsne_df$batchlb <- metadata[ , batch_label]
    # #tsne_df$batch <- metadata[ , "batch"]
    # tsne_df$CellType <- metadata[ , celltype_label]
    # 
    # p01 <- ggplot(tsne_df, aes(x = tSNE_1, y = tSNE_2, colour = batchlb)) + geom_point(alpha = 0.6) + theme_bw()
    # p01 <- p01 + labs(x='tSNE_1',y='tSNE_2',title='fast MNN')
    # p01 <- p01 + theme(legend.title = element_text(size=17), 
    #                    legend.key.size = unit(1.1, "cm"),
    #                    legend.key.width = unit(0.5,"cm"), 
    #                    legend.text = element_text(size=14), 
    #                    plot.title = element_text(color="black", size=20, hjust = 0.5))
    # 
    # p02 <- ggplot(tsne_df, aes(x = tSNE_1, y = tSNE_2, colour = CellType)) + geom_point(alpha = 0.6) + theme_bw()
    # p02 <- p02 + labs(x='tSNE_1',y='tSNE_2',title='fast MNN')
    # p02 <- p02 + theme(legend.title = element_text(size=17), 
    #                    legend.key.size = unit(1.1, "cm"),
    #                    legend.key.width = unit(0.5,"cm"), 
    #                    legend.text = element_text(size=14), 
    #                    plot.title = element_text(color="black", size=20, hjust = 0.5))
    # 
    # png(paste0(plotout_dir,outfilename_prefix,tsneplot_filename,".png"),width = 2*1000, height = 800, res = 2*72)
    # print(plot_grid(p01, p02))
    # dev.off()
    # 
    # pdf(paste0(plotout_dir,outfilename_prefix,tsneplot_filename,".pdf"),width=15,height=7,paper='special')
    # print(plot_grid(p01, p02))
    # dev.off()
  }

}







