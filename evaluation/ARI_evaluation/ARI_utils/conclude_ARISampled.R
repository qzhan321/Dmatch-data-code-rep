

conclude_ARISampled <- function(dir_this, csv_title, nbiters = 20){
  library(ggplot2)
  
  setwd(dir_this)
  filesneed<-list.files(pattern = "_ARI.txt")
  
  method<-vector()
  medianARIbatch<-vector()
  medianARIcelltype<-vector()
  
  for (x in 1:length(filesneed)){
    temp<-read.table(filesneed[x], header = TRUE, stringsAsFactors = FALSE)
    method[x]<-as.character(temp$use_case[1])
    medianARIbatch[x]<-temp$ari_batch[101]
    medianARIcelltype[x]<-temp$ari_celltype[101]
    rm(temp)
  }
  
  # normalise values to 0 - 1
  min_batch <- min(medianARIbatch)
  max_batch <- max(medianARIbatch)
  min_cell <- min(medianARIcelltype)
  max_cell <- max(medianARIcelltype)
  medianARIbatch_norm <- (medianARIbatch-min_batch)/(max_batch-min_batch)
  medianARIcelltype_norm <- (medianARIcelltype-min_cell)/(max_cell-min_cell)
  
  # produce final fscore ARI, similar to scMerge paper
  medianfscoreARI <- (2 * (1 - medianARIbatch_norm)*(medianARIcelltype_norm))/
                          (1 - medianARIbatch_norm + medianARIcelltype_norm)
  
  sum_xy<-medianARIcelltype_norm+(1-medianARIbatch_norm)
  
  finaldf<-data.frame("ARIMethod" = method, 
                      "ARIbatchMedian" = medianARIbatch, 
                      "ARIbatchMedian_norm" = medianARIbatch_norm, 
                      "ARIcelltypeMedian" = medianARIcelltype, 
                      "ARIcelltypeMedian_norm" = medianARIcelltype_norm, 
                      "ARI_fscore" = medianfscoreARI,
                      "ARI_summedXY" = sum_xy)
  
  write.csv(finaldf, file = paste0("allARI_", csv_title, ".csv"), row.names = FALSE)
  return(finaldf)
}
