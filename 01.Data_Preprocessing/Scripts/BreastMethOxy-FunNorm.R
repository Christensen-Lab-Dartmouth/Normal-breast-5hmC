#######################################################################################################################
#Code from Andy Houseman to calculate 'betas' from FunNorm signals 
#######################################################################################################################

#"Betas" that sum to exactly one (unmethylated, 5mC, 5hmC)
# (Andy's annotation as it is)

# derivative of -log(beta pdf) wrt a parm
diffBeta1 <- function(x,a,b){
  digamma(a)-digamma(a+b)-log(x)
}

# derivative of -log(beta pdf) wrt b parm
diffBeta2 <- function(x,a,b){
  digamma(b)-digamma(a+b)-log(1-x)
}

# beta-likelihood function
likeOxBS <- function(theta, betaBS, betaOxBS, signalBS, signalOxBS){
  theta <- pmin(100,pmax(-100,theta))
  p <- exp(c(0,theta))
  p <- p/sum(p)
  a1  <- signalBS*(p[2]+p[3])
  b1  <- signalBS*(p[1])
  a2 <- signalOxBS*(p[2])
  b2 <- signalOxBS*(p[1]+p[3])  
  -(dbeta(betaBS,a1,b1,log=TRUE)+dbeta(betaOxBS,a2,b2,log=TRUE)) 
}

# Beta score function
scoreOxBS <- function(theta, betaBS, betaOxBS, signalBS, signalOxBS){
  theta <- pmin(100,pmax(-100,theta))
  p <- exp(c(0,theta))
  p <- p/sum(p)
  a1  <- signalBS*(p[2]+p[3])
  b1  <- signalBS*(p[1])
  a2 <- signalOxBS*(p[2])
  b2 <- signalOxBS*(p[1]+p[3])  
  dp <-   matrix(c(
     -p[1]*p[2], p[2]*(p[1]+p[3]), -p[2]*p[3],
     -p[1]*p[3], -p[2]*p[3], p[3]*(p[1]+p[2])
  ),3,2)
  ua1 <- ( diffBeta1(betaBS,a1,b1)*signalBS*(dp[2,1]+dp[3,1]) +
           diffBeta1(betaOxBS,a2,b2)*signalOxBS*(dp[2,1]) )
  ub1 <- ( diffBeta2(betaBS,a1,b1)*signalBS*(dp[1,1]) +
         diffBeta2(betaOxBS,a2,b2)*signalOxBS*(dp[1,1]+dp[3,1]) )

  ua2 <- ( diffBeta1(betaBS,a1,b1)*signalBS*(dp[2,2]+dp[3,2]) +
           diffBeta1(betaOxBS,a2,b2)*signalOxBS*(dp[2,2]) )
  ub2 <- ( diffBeta2(betaBS,a1,b1)*signalBS*(dp[1,2])+
           diffBeta2(betaOxBS,a2,b2)*signalOxBS*(dp[1,2]+dp[3,2]) )
  c(ua1+ub1,ua2+ub2)

}

# Fit one value
fitOneOxBS <- function(betaBS, betaOxBS, signalBS, signalOxBS, eps=1E-5){
  est5mC <- max(min(betaOxBS,1-eps),eps)
  estTotMeth <- max(min(betaBS,1-eps),eps)
  est5hmC <- max(estTotMeth-est5mC,eps)
  
  theta <- log(c(est5mC,est5hmC))-log(1-estTotMeth)
  opt <- try( optim(theta, likeOxBS, method="BFGS", gr=scoreOxBS, 
    betaBS=betaBS, betaOxBS=betaOxBS, 
    signalBS=signalBS, signalOxBS=signalOxBS) )

  if(inherits(opt,"try-error")){
    #print(c(est5mC,est5hmC,estTotMeth))
    out <- c(1-estTotMeth,est5mC,est5hmC)
  }
  else out <- exp(c(0,opt$par))

  #names(out) <- c("C","5mC","5hmC")
  out/sum(out)
}

# Fit multiple values
fitOxBS <- function(betaBS, betaOxBS, signalBS, signalOxBS, eps=1E-5){
  n <- length(betaBS)
  out <- matrix(NA,n,3)
  for(i in 1:n){
    out[i,] <- fitOneOxBS(betaBS[i],betaOxBS[i], signalBS[i], signalOxBS[i],eps=eps)
    if(i %% 5000==0) cat(i,"\n")
  }
  colnames(out) <- c("C","5mC","5hmC")
  rownames(out) <- names(betaBS)
  out
}

##################################
#Processing the manifest file for a second time. If the R workspace was not reset, some of 
#these steps would not be necessary
  setwd("/Users/kevinjohnson/Documents/NormalBreast_Glioma/Preprocessing_Data")
  meta <- read.delim("metaData-breast-150318.txt", sep="\t", head=TRUE, 
  stringsAsFactors=FALSE)
  meta$ArrayID <- with(meta, paste(Beadchip,Position,sep="_"))
  meta$Label.ID.[meta$OX_BS==1] <- gsub("B$","",meta$ID[meta$OX_BS==1])
  meta$Label.ID.[meta$OX_BS==0] <- gsub("A$","",meta$ID[meta$OX_BS==0])
  with(meta, all(Label.ID.==gsub("[A|B]$","",ID))) # Check consistency of IDs
  mrg <- merge(
    data.frame(LabelID=meta$Label.ID.[meta$OX_BS==0],
       ixBS=as.numeric(rownames(meta)[meta$OX_BS==0]),
       stringsAsFactors=FALSE),
    data.frame(LabelID=meta$Label.ID.[meta$OX_BS==1],
       ixOxBS=as.numeric(rownames(meta)[meta$OX_BS==1]),
       stringsAsFactors=FALSE)
)
#Identify the indices of a paired sample (i.e. B10A in row 1, matches B10B in row 2)
  meta$Pair[mrg$ixOxBS] <- mrg$ixBS 
  meta$Pair[mrg$ixBS] <- mrg$ixOxBS 

#'ArrayID' represents the location of the sample in which well and on what slide 
  rownames(meta) <- meta$ArrayID

#Get the Subject ID and location of BS Slide/Well and oxBS Slide/Well for given sample
  tmp <- meta[meta$OX_BS==0,]
  matchedIDs <- data.frame(SubjectID=tmp$Label.ID.,
    BS=tmp$ArrayID, OxBS=meta$ArrayID[tmp$Pair],stringsAsFactors=FALSE)

##################################################
#Load list from "process450Ksignals_kcj.R" file
load("/Users/kevinjohnson/Documents/NormalBreast_Glioma/Preprocessing_Data/OxyBreastSignals-FunNorm.RData")
#The estimation at each probe occurs independently so you don't have to drop poor/SNP probes here
###################################################
#For 30 samples this next step took ~24 hours on the Cluster
#change the 'NCPG' to a small value ~10 or 100 to make sure that it is working
NCPG <- dim(methBS)[1] #NCPG <- 100
NSPEC <- dim(matchedIDs)[1]

fitOxBSfromIDs <- function(idBS, idOxBS, rowsubset=1:NCPG){
  signalBS <- methBS[rowsubset,idBS]+unmethBS[rowsubset,idBS]
  signalOxBS <-methOxBS[rowsubset,idOxBS]+unmethOxBS[rowsubset,idOxBS]
  betaBS <- methBS[rowsubset,idBS]/signalBS
  betaOxBS <- methOxBS[rowsubset,idOxBS]/signalOxBS
  fitOxBS(betaBS,betaOxBS,signalBS,signalOxBS)
}
MethOxy <- array(NA,dim=c(NCPG,NSPEC,3))
dimnames(MethOxy) <- list(rownames(methBS)[1:NCPG],matchedIDs$SubjectID,c("C","5mC","5hmC"))

for(i in 1:NSPEC){
  MethOxy[,i,] <- fitOxBSfromIDs(matchedIDs$BS[i], matchedIDs$OxBS[i], 1:NCPG)
  cat(i,"\n")
}
save(list=c("MethOxy","meta"), file="OxyBreastMethOxy-FunNorm.RData")
