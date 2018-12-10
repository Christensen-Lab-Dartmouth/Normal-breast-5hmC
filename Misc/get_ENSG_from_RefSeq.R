get_ENSG_from_RefSeq <- function(refseq_ids, ENSG_model_dataset, ENGS_model){
  ENSG <- c()
  ENSG_temp1 <- c()
  if(length(refseq_ids)==1){
    ENSG_temp1 <- unique(ENSG_model_dataset$ensembl_gene_id[which(ENSG_model_dataset[,c(ENGS_model)] == refseq_ids)])
    if(length(ENSG_temp1)==0){
      ENSG <- NA
      rm(ENSG_temp1)
    } else {
      ENSG <- ENSG_temp1
      rm(ENSG_temp1)
    }
  } 
  else if(length(refseq_ids)>1){
    ENSG_temp2.1 <- list()
    for(j in 1:length(refseq_ids)){
      ENSG_temp2.2 <- unique(ENSG_model_dataset$ensembl_gene_id[which(ENSG_model_dataset[,c(ENGS_model)] == refseq_ids[j])])
      if(length(ENSG_temp2.2)==0){
        ENSG_temp2.1[[j]] <- NA
      } else {
        ENSG_temp2.1[[j]] <- ENSG_temp2.2
      }
      rm(ENSG_temp2.2)
    }
    if(length(unlist(ENSG_temp2.1))==0){
      ENSG <- NA
      rm(ENSG_temp2.1)
    } else {
      ENSG <- unique(unlist(ENSG_temp2.1))
      rm(ENSG_temp2.1)
    }
  } 
  
  else {
    ENSG <- NA
  }
  
  ENSG
}
