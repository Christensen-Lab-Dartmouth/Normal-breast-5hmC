#######################################################################################################################
# Produce figures of empirical cumulative distribution functions for 5hmC/5mC
#######################################################################################################################
rm(list=ls())
setwd("/Users/Owen 1/Dropbox (Christensen Lab)/NDRI_Breast_5hmC_update/")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load data 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# annotation data 
load("02.Characterization_5hmC_levels/Files/annotation-150214.RData")
# methylation data 
load("02.Characterization_5hmC_levels/Files/OxyBreastMethOxy-FunNorm.RData")
# KEGG key
load("02.Characterization_5hmC_levels/Files/KEGG-Key.RData")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# pre-process data 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Create descriptions for the KEGG pathways
pathDescript <- pathKey$pathway
names(pathDescript) <- pathKey$id 

#Replace the pathKey object with KEGG IDs
pathKey <- names(keggPaths) 
#All the unique CpGs that track to KEGG pathways
uniqueKEGG <- unique(unlist(keggPaths)) 
#Remove the "N_" and "S_"  annotation from the annotation file
geoUCSC <- gsub("^[N|S]_","",longAnnot$RelToIslandUCSC) 

#Labelling all of the 'bad' probes
strat <- rep("bad",dim(longAnnot)[1]) 
#Identifying the 'good' type I probes
strat[lmIndex[["I:0:0:0:0"]]] <- "I" 
#Identifying the 'good' type II probes
strat[lmIndex[["II:0:0:0:0"]]] <- "II" 
#Creating a unique group of CpG Island and Probe-type
strat <- paste(strat, geoUCSC, sep=":") 

#New object name moving forward
strat2 <- strat 

#Extract row information for TSS-associated CpGs
tmp <- grep("TSS",longAnnot$UCSC_RefGene_Group) 
#Adding TSS indication to relevant CpGs
strat2[tmp] <- paste(strat2[tmp],"TSS",sep=":") 
#Matching CpG ID to the stratified names
names(strat2) <- longAnnot$TargetID 

#Drops the 8 elements of 'bad' probes
strat2Index <- split(longAnnot$TargetID, strat2)[-(1:8)] 
#Number of probes per stratification
nPerStrat <- sapply(strat2Index,length) 

#Total number of CpGs
nCpGTot <- dim(MethOxy)[1] 

#(Sel)ecting only the 'good' probes
sel <- unlist(lmIndex[c("I:0:0:0:0","II:0:0:0:0")]) 
names(sel) <- rownames(MethOxy[sel, , 3])

#Attaching names to good probe indices
selIndex <- 1:length(sel)
names(selIndex) <- rownames(MethOxy[sel,,3])

#CG ID and its probe type/location/TSS status
strat2sel <- strat2[sel]
STRAT2SELINDEX <- split(1:length(sel), strat2sel)

#Create an index for just CpG island region based analysis
AnnotSel <- longAnnot[sel, ]
AnnotSel$CpG_Regions <- gsub("^[N|S]_","", AnnotSel$RelToIslandUCSC) 
CPGSTRATSELINDEX <- split(1:length(sel), AnnotSel$CpG_Regions)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Calculate CpG specific statistics
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Calculate the mean, median, and correlations of 5mC and 5hmC
#### mean 
mbeta1 <- apply(MethOxy[sel, , 3], 1, mean, na.rm=TRUE)
mbeta2 <- apply(MethOxy[sel, , 2], 1, mean, na.rm=TRUE)
#### median
mdbeta1 <- apply(MethOxy[sel, , 3], 1, median, na.rm=TRUE)
mdbeta2 <- apply(MethOxy[sel,,2],1,median,na.rm=TRUE)
#### correlations
corbeta12 <- apply(MethOxy[sel,,2:3], 1, cor, use="pair")[2,] # pearson
corsbeta12 <- apply(MethOxy[sel,,2:3], 1, cor, use="pair",method="spearman")[2,] # spearman

# calculate empirical cumulatiev distribution function (ECDF)
#### using mean betas 
ecdfMean1 <- ecdf(mbeta1)
ecdfMean2 <- ecdf(mbeta2)
#### using median betas 
ecdfMed1 <- ecdf(mdbeta1)
ecdfMed2 <- ecdf(mdbeta2)
# for CpG specific corrlations 
ecdfCorrP <- ecdf(corbeta12)
ecdfCorrS <- ecdf(corsbeta12)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# generate ECDF plots
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Produce ECDF median figure for 5hmC and 5mC in Breast
ppi = 300
png("02.Characterization_5hmC_levels/Figures/ECDF_median_Breast.png",  width=6.5*ppi, height=5*ppi, res=ppi)
plot(ecdfMed2,main="Cumulative Density 5(h)mC", xlab="Median average beta",ylab="Cumulative proportion", col="blue", cex.lab = 1.3, cex.axis = 1.15, las = 1)
plot(ecdfMed1 ,add=TRUE,col="red")
legend(0.7,0.3,c("5mC","5hmC"),lwd=3,col=c("blue","red"))
dev.off()

#Produce ECDF correlations figure for 5hmC and 5mC in Breast
ppi = 300
png("02.Characterization_5hmC_levels/Figures/ECDF_corr_Breast.png", width=6.5*ppi, height=5*ppi, res=ppi)
plot(ecdfCorrP, main="Cumulative Density Correlation", xlab="Correlation coefficient", ylab="Cumulative proportion", col="black", cex.lab = 1.3, cex.axis = 1.15, las = 1)
plot(ecdfCorrS, add=TRUE,col="purple")
legend(0.35, 0.3, c("Pearson","Spearman"), lwd=3, col=c("black","purple"))
abline(v=0, lty=3)
dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# generate 5hmC/5mC methylation percentile plots 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Set arbitrary thresholds to plot 5hmC and 5hmC by cut-offs (as opposed to ECDF curve)
bThresh1 <- c(0.01,0.025,0.05,0.1,0.2,0.1) 
bThresh2 <- c(0,0.25,0.5,0.75,1.0)
qtile <- c(0.1,0.2,0.3,0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)

#For each quantile
qtileStratMedian1 <- sapply(STRAT2SELINDEX, function(u)quantile(mdbeta1[u],prob=qtile))
qtileStratMedian2 <- sapply(STRAT2SELINDEX, function(u)quantile(mdbeta2[u],prob=qtile))
qtile_CpG_StratMedian1 <- sapply(CPGSTRATSELINDEX, function(u)quantile(mdbeta1[u],prob=qtile))
qtile_CpG_StratMedian2 <- sapply(CPGSTRATSELINDEX, function(u)quantile(mdbeta2[u],prob=qtile))

#Colors for plots
STRATPALETTE <- c("Island"="wheat","OpenSea"="blue",
  "Shelf"="cyan", "Shore"="green")

# function to generate plots 
plotQtile <- function(x,legendcex=0.9){
  nq <- dim(qtile_CpG_StratMedian1)[1]
  ns <- dim(qtile_CpG_StratMedian1)[2]
  labs <- colnames(x)
  col <- rep("",ns)
  sym <- rep(0,ns)
  cx <- lt <- rep(1,ns)
  tss <- type2 <- rep(FALSE,ns)

  tss[grep("TSS$",labs)] <- TRUE
  type2[grep("^II[:]",labs)] <- TRUE 

  for(j in 1:4) {
    col[grep(names(STRATPALETTE)[j],labs)] <- STRATPALETTE[j]
  }
  sym[tss & !type2] <- 15
  sym[tss & type2] <- 22
  sym[!tss & !type2] <- 19
  sym[!tss & type2] <- 21
  cx[tss] <- 1.3
  lt[!tss] <- 2

  plot(c(1,nq),range(x),xlab="Quantile",ylab="Median beta-value",type="n",xaxt="n", las = 1, cex.lab = 1.3, cex.axis = 1.15)
  axis(1,1:nq,rownames(x))
  for(i in 1:ns) {
    points(1:nq, x[,i], col=col[i], pch=sym[i], cex=cx[i])
    lines(1:nq, x[,i], col=col[i], lwd=2, lty=lt[i])
  }
  legend(1.1,max(x),labs,col=col,pch=sym,lty=lt,cex=legendcex)
}

# apply function to output plots 
#### 5hmC
ppi = 300
png("02.Characterization_5hmC_levels/Figures/Median_5hmC_Quantiles_CpGStrat2.png", width=6.5*ppi, height=5*ppi, res=ppi)
plotQtile(qtile_CpG_StratMedian1) ; title("Median 5hmC")
dev.off()
#### 5mC 
png("02.Characterization_5hmC_levels/Figures/Median_5mC_Quantiles_CpGStrat2.png", width=6.5*ppi, height=5*ppi, res=ppi)
plotQtile(qtile_CpG_StratMedian2) ; title("Median 5mC")
dev.off()

