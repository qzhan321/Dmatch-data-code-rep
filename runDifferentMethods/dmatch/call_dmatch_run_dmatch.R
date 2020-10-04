#' dmatch pipeline

#' @param batch_data : the preprocessed data
#' @param npcs : number of principal components (default 30)
#' @param Reference_name: the name of the reference, either "cell.line.to.ensl.ID" for ensembl id, or "cell_atlas_ref_panel" for gene names
#' @param TopCellLineNumber
#' @param hclust_method : 
#' @param ShowCellNumber :  
#' @param cor_threshold
#' @param num_clusters
#' @param select_clusters_quantile
#' @param anchors
#' @param run_alignment_quantile
#' @param plotout_dir : path to saving plot files (default current directory)
#' @param saveout_dir : path to saving output files (default current directory)
#' @param save_obj : logical to save Seurat object (default TRUE)




call_dmatch <- function(batch_data, npcs = 30,
                        pcaplot_filename = "",
                        batch.id, Reference_name,
                        hclust_method = "ward.D",
                        TopCellLineNumber = 5,
                        ShowCellNumber = 20,
                        cor_threshold = 0.1,
                        num_clusters = 5,
                        cut_groups_methods = "ward.D",
                        select_clusters_quantile = 0.98,
                        anchors = c(1),
                        run_alignment_quantile = 0.98,
                        plotout_dir = "", saveout_dir = "", 
                        outfilename_prefix = "", 
                        dmatch_out_filename = "",
                        dmatch_out_pcs_filename = "",
                        visualize = T, save_obj = T,
                        umapplot_filename,
                        heatmap_filename="")
{
  
  b1_data <- as.data.frame(batch_data[[1]])
  b2_data <- as.data.frame(batch_data[[2]])
  b1_meta <- as.data.frame(batch_data[[3]])
  b2_meta <- as.data.frame(batch_data[[4]])
  metadata <- rbind(b1_meta, b2_meta)
  
  ##########################################################
  #run
  t1 = Sys.time()
  PCA<-fastSVD(list(b1_data,b2_data), nPC = npcs)
  dmatch::PCAPlot(PCA = PCA, PCs.to.plot = c(1,2), batchs.to.plot = c(1,2), filename = paste0(plotout_dir,outfilename_prefix,pcaplot_filename, ".png"))
  samples<-CreatedmatchObject(pairwiseSamples.list = list(b1_data,b2_data), batch.id = batch.id, PCA = PCA)
  path<-system.file("extdata", Reference_name, package = "dmatch")
  load(path)
  if (Reference_name == "cell_atlas_ref_panel") {
    rownames(cell.line) <- toupper(rownames(cell.line))
    samples<-projection_to_reference_panel(samples,Reference = cell.line)
    samples<-projection_visualization(samples, hclust.method = hclust_method, 
                                      TopCellLineNumber = TopCellLineNumber, 
                                      ShowCellNumber = ShowCellNumber , 
                                      cor.threshold = cor_threshold, filename = paste0(plotout_dir,outfilename_prefix,heatmap_filename, ".png"))
    samples<-cut_groups(samples,K=num_clusters, method = cut_groups_methods)
    samples<-select_clusters(samples, quantile = select_clusters_quantile)
    
    temp <- samples@select.clusters$cells.num
    print(temp)
    index1 <- which(temp[,1] > 30)
    index2 <- which(temp[,2] > 30)
    anchors <- intersect(index1, index2)
    print(anchors)
    samples<-run_alignment_by_2D(samples, selected = anchors, quantile = run_alignment_quantile)
  } else {
    rownames(cleaned.cell.line) <- toupper(rownames(cleaned.cell.line))
    samples<-projection_to_reference_panel(samples,Reference = cleaned.cell.line)
    samples<-projection_visualization(samples, hclust.method = hclust_method, 
                                      TopCellLineNumber = TopCellLineNumber, 
                                      ShowCellNumber = ShowCellNumber , 
                                      cor.threshold = cor_threshold, filename = paste0(plotout_dir,outfilename_prefix,heatmap_filename, ".pdf"))
    samples<-cut_groups(samples,K=num_clusters, method = cut_groups_methods)
    samples<-select_clusters(samples, quantile = select_clusters_quantile)
    # samples@select.clusters$cells.num
    # samples@select.clusters$shapiro.test.pvalue
    samples<-run_alignment_by_2D(samples, selected = anchors, quantile = run_alignment_quantile)
  }
  t2 = Sys.time()
  print(t2-t1)
  
  ##########################################################
  #save
  saveRDS(samples, file=paste0(saveout_dir, outfilename_prefix, dmatch_out_filename, ".RDS"))
  corrected_pca_dmatch <- as.data.frame(rbind(samples@run_alignment_by_2D.results$Corrected, samples@run_alignment_by_2D.results$Reference))
  saveRDS(corrected_pca_dmatch, file=paste0(saveout_dir, outfilename_prefix, dmatch_out_pcs_filename, ".RDS"))
  
  #umap
  set.seed(10)
  out_umap <- umap(corrected_pca_dmatch)
  umap_df<- as.data.frame(out_umap$layout)
  
  rownames(umap_df) <- rownames(corrected_pca_dmatch)
  dim(umap_df)
  colnames(umap_df) <- c('umap_1', 'umap_2')
  umap_df$batchlb <- as.factor(c(rep(1, nrow(samples@run_alignment_by_2D.results$Corrected)), rep(2, nrow(samples@run_alignment_by_2D.results$Reference))))
  umap_df$CellType <- as.factor(c(b1_meta$celltype, b2_meta$celltype))
  
  p01 <- ggplot(umap_df, aes(x = umap_1, y = umap_2, colour = batchlb)) + geom_point(alpha = 0.6) + theme_bw()
  p01 <- p01 + labs(x='umap_1',y='umap_2',title='dmatch')
  p01 <- p01 + theme(legend.title = element_text(size=17), 
                     legend.key.size = unit(1.1, "cm"),
                     legend.key.width = unit(0.5,"cm"), 
                     legend.text = element_text(size=14), 
                     plot.title = element_text(color="black", size=20, hjust = 0.5))
  
  p02 <- ggplot(umap_df, aes(x = umap_1, y = umap_2, colour = CellType)) + geom_point(alpha = 0.6) + theme_bw()
  p02 <- p02 + labs(x='umap_1',y='umap_2',title='dmatch')
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
} 

  








