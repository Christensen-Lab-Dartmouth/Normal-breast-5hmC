get_ENST_from_RefSeq <- function(refseq_ids, ENST_model_dataset, ENGS_model){
  ENST <- c()
  ENST_temp1 <- c()
  if(length(refseq_ids)==1){
    ENST_temp1 <- unique(ENST_model_dataset$ensembl_transcript_id[which(ENST_model_dataset[,c(ENGS_model)] == refseq_ids)])
    if(length(ENST_temp1)==0){
      ENST <- NA
      rm(ENST_temp1)
    } else {
      ENST <- ENST_temp1
      rm(ENST_temp1)
    }
  } 
  else if(length(refseq_ids)>1){
    ENST_temp2.1 <- list()
    for(j in 1:length(refseq_ids)){
      ENST_temp2.2 <- unique(ENST_model_dataset$ensembl_transcript_id[which(ENST_model_dataset[,c(ENGS_model)] == refseq_ids[j])])
      if(length(ENST_temp2.2)==0){
        ENST_temp2.1[[j]] <- NA
      } else {
        ENST_temp2.1[[j]] <- ENST_temp2.2
      }
      rm(ENST_temp2.2)
    }
    if(length(unlist(ENST_temp2.1))==0){
      ENST <- NA
      rm(ENST_temp2.1)
    } else {
      ENST <- unique(unlist(ENST_temp2.1))
      rm(ENST_temp2.1)
    }
  } 
  
  else {
    ENST <- NA
  }
  
  ENST
}
