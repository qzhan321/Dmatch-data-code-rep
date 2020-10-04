

ari_consolidate<-function(){
  baseName<-rev(unlist(strsplit(getwd(), split = "/", fixed = TRUE)))[1]
  wds<-rev(unlist(strsplit(getwd(), split = "/", fixed = TRUE)))[2]
  filenames<-list.files(pattern = "_ARI.txt")
  bigone<-data.frame(NULL)

  for (x in 1:length(filenames)){
    temp<-read.table(file = filenames[x], header = TRUE, stringsAsFactors = FALSE)
    temp2<-temp[1:(dim(temp)[1]-1),]
    bigone<-rbind(bigone, temp2)
  }
  
  write.table(bigone, file = paste0("AllValues_", baseName, wds, ".txt"), row.names = FALSE, 
              col.names = TRUE, quote = FALSE, sep="\t")
  return(bigone)
}

