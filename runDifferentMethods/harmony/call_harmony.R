
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


harmony_preprocess <- function(expr_mat_dmatch, metadata, projectName,
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
  
  b_seurat <- CreateSeuratObject(counts = expr_mat_dmatch, meta.data = metadata, project = projectName, 
                                 min.cells = min_cells, min.features = min_genes)
  if (normData) {
    b_seurat <- NormalizeData(object = b_seurat, normalization.method = norm_method, scale.factor = scale_factor, verbose = F)
  } else {
    b_seurat@assays$RNA@data = b_seurat@assays$RNA@counts
  }
  b_seurat <- ScaleData(object = b_seurat, verbose = F)
  b_seurat <- FindVariableFeatures(object = b_seurat, selection.method = selection_method, nfeatures = n_features, verbose = F)
  
  return(b_seurat)
}






call_harmony_2 <- function(b_seurat, batch_label, celltype_label, npcs = 30, 
                         plotout_dir = "", saveout_dir = "", 
                         outfilename_prefix = "", 
                         visualize = T, save_obj = T)
{

  #Harmony settings
  theta_harmony = 2
  numclust = 50
  max_iter_cluster = 100

  #plotting
  k_seed = 10
  
  umapplot_filename = "_harmony_umap"
  obj_filename = "_harmony_sobj"
  pca_filename = "_harmony_pca"

  ##########################################################
  #run

  t1 = Sys.time()

  b_seurat <- RunPCA(object = b_seurat, npcs = npcs, features = b_seurat@assays$RNA@var.features, verbose = F)

  b_seurat <- RunHarmony(object = b_seurat, batch_label, theta = theta_harmony, plot_convergence = TRUE, 
                          nclust = numclust, max.iter.cluster = max_iter_cluster)

  t2 = Sys.time()
  print(t2-t1)

  ##########################################################
  #save

  harmony_res <- as.data.frame(b_seurat@reductions$harmony@cell.embeddings)
  cells_use <- rownames(harmony_res)

  harmony_res$batch <- b_seurat@meta.data[, batch_label]
  harmony_res$celltype <- b_seurat@meta.data[, celltype_label]
  write.table(harmony_res, file=paste0(saveout_dir,outfilename_prefix,pca_filename,".txt"), quote=F, sep='\t', row.names = T, col.names = NA)

  if(save_obj) {
    saveRDS(b_seurat, file=paste0(saveout_dir,outfilename_prefix,obj_filename,".RDS"))
  }

  if (visualize) {

    ##########################################################
    #preparing plots

    #b_seurat <- RunTSNE(b_seurat, reduction.use = "harmony", dims.use = 1:npcs, k.seed = 10, do.fast = T, check_duplicates = FALSE, perplexity = tsne_perplex)
    b_seurat <- RunUMAP(b_seurat, reduction = "harmony", dims = 1:npcs, k.seed = 10, do.fast = T)

    ##########################################################
    #tSNE plot

    # p11 <- TSNEPlot(b_seurat, do.return = T, pt.size = 0.5, group.by = batch_label)
    # p12 <- TSNEPlot(b_seurat, do.return = T, pt.size = 0.5, group.by = celltype_label)
    # 
    # png(paste0(plotout_dir,outfilename_prefix,tsneplot_filename,".png"),width = 2*1000, height = 800, res = 2*72)
    # print(plot_grid(p11, p12))
    # dev.off()
    # 
    # pdf(paste0(plotout_dir,outfilename_prefix,tsneplot_filename,".pdf"),width=15,height=7,paper='special') 
    # print(plot_grid(p11, p12))
    # dev.off()

    ##########################################################
    #UMAP plot

    p21 <- DimPlot(object = b_seurat, reduction.use = 'umap', group.by = batch_label)
    p22 <- DimPlot(object = b_seurat, reduction.use = 'umap', group.by = celltype_label)   # or orig.ident

    png(paste0(plotout_dir,outfilename_prefix,umapplot_filename,".png"),width = 2*1000, height = 800, res = 2*72)
    print(plot_grid(p21, p22))
    dev.off()

    pdf(paste0(plotout_dir,outfilename_prefix,umapplot_filename,".pdf"),width=15,height=7,paper='special') 
    print(plot_grid(p21, p22))
    dev.off()

  }

  return(b_seurat)
}







