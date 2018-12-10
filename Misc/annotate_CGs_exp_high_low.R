annotate_CGs_exp_high_low <- function(ENSG_entires){
  cg_exp_grp <- c()
  cg_exp_grp_temp1 <- c()
  if(length(ENSG_entires)==1){
    cg_exp_grp_temp1 <- genes$exp_grp[which(genes$genes_split == ENSG_entires)]
    if(length(cg_exp_grp_temp1)==0){
      cg_exp_grp <- NA
    } else {
      cg_exp_grp <- cg_exp_grp_temp1
    }
  } else if (length(ENSG_entires)>1){
    cg_exp_grp_temp2.1 <- c()
    for(j in 1:length(ENSG_entires)){
      cg_exp_grp_temp2.2 <- genes$exp_grp[which(genes$genes_split == ENSG_entires[j])]
      if(length(cg_exp_grp_temp2.2)==0){
        cg_exp_grp_temp2.1[j] <- NA
      } else{
        cg_exp_grp_temp2.1[j] <- cg_exp_grp_temp2.2
      }
      rm(cg_exp_grp_temp2.2)
    }
    cg_exp_grp_temp2.1 <- na.omit(cg_exp_grp_temp2.1)
    if(length(cg_exp_grp_temp2.1)>0){
      if(any(cg_exp_grp_temp2.1=="High")){
        cg_exp_grp <- "High"
        rm(cg_exp_grp_temp2.1)
      } else {
        cg_exp_grp <- "Low"
        rm(cg_exp_grp_temp2.1)
      }
    } else {
      rm(cg_exp_grp_temp2.1)
    }
  } else{
    cg_exp_grp <- NA
  }
  cg_exp_grp
}
