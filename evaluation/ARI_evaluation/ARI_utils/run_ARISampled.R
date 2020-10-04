

run_ARISampled <- function(fn, this_dir, out_dir, meta_dir, eval_metric, methods_use, cpcs = 1:20, 
                           celltypelb='celltype', batchlb='batch',batch_effect = "", ratio_name="", 
                           raw = F, rep = NULL, txt_out_dir = "", maxiter = maxiter, nbiters = nbiters){
  setwd(this_dir)
  dataset_no<-paste(rev(strsplit(this_dir, split = "/")[[1]])[2],rev(strsplit(this_dir, split = "/")[[1]])[1], sep = "_")
  
  if (raw) {
    samples <- readRDS(paste0("/project2/mengjiechen/qizhan/qizhan/dmatch/dmatch_response_to_review/others/broad_PBMC_10x_batch_effects_simulation_and_correction/corrected_data/pbmc1_10x_v2_A/dmatch/", batch_effect, "/results/pbmc1_10x_v2_A_", ratio_name, "_rep", rep, "_samples.RDS"))
    thisData <- rbind(samples@run_alignment_by_2D.results$Original, samples@run_alignment_by_2D.results$Reference) 
  } else {
    if (grepl(".txt", fn, fixed = TRUE)) {
      thisData <- read.table(paste0(this_dir,"/results/",fn), head=T, row.names = 1, check.names = FALSE)
    } else if (grepl(".RDS", fn, fixed = TRUE)) {
      thisData <- readRDS(paste0(this_dir,"/results/",fn))
    }
  }
  
  
  temp1 <- strsplit(fn, split = "_")
  temp1 <- temp1[[1]][7]
  temp2 <- substr(fn, 1, 18)
  mfn <- paste(temp2, temp1, sep = "_")
  load(paste0(meta_dir, mfn))
  if (ncol(thisData) == 30) {
    meta_bc <- matrix(NA, nrow = nrow(thisData), ncol = 2)
    thisData <- cbind(thisData, meta_bc)
    meta1 <- data.frame("celltype" = sim_list[[3]], "batch" = rep(1, length(sim_list[[3]])))
    rownames(meta1) <- colnames(sim_list[[1]])
    meta2 <- data.frame("celltype" = sim_list[[4]], "batch" = rep(2, length(sim_list[[4]])))
    rownames(meta2) <- colnames(sim_list[[2]])
    
    meta12 <- rbind(meta1, meta2)

    thisData[,31] <- meta12[rownames(thisData),"celltype"]
    thisData[,32] <- meta12[rownames(thisData),"batch"]
  } else if (ncol(thisData) == 32) {
    meta1 <- data.frame("celltype" = sim_list[[3]], "batch" = rep(1, length(sim_list[[3]])))
    rownames(meta1) <- colnames(sim_list[[1]])
    meta2 <- data.frame("celltype" = sim_list[[4]], "batch" = rep(2, length(sim_list[[4]])))
    rownames(meta2) <- colnames(sim_list[[2]])
    
    meta12 <- rbind(meta1, meta2)
    
    thisData[,31] <- meta12[rownames(thisData),"celltype"]
    thisData[,32] <- meta12[rownames(thisData),"batch"]
  }
  
  colnames(thisData)[31:32] <- c(celltypelb, batchlb)
  setwd(paste0(out_dir, '/', eval_metric))
  rep_name <- paste0("rep_", rep)
  temp<-ari_calcul_sampled(myData=thisData, cpcs=cpcs, isOptimal=FALSE, 
                           method_use = methods_use,  
                           base_name=paste(methods_use, batch_effect, ratio_name, rep_name, eval_metric, sep = "_"), 
                           nbiters = nbiters, percent_extract = percent_extract,
                           maxiter = maxiter, txt_out_dir = txt_out_dir)
  return(temp)
}

