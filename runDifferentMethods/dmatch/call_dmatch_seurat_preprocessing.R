
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


dmatch_preprocess <- function(expr_mat_dmatch, metadata, projectName,
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
