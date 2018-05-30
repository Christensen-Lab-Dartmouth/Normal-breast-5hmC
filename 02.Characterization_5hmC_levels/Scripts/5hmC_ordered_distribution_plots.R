#######################################################################################################################
# Generate ordered distribution plot for 5hmC abundance (averaged accross samples) 
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
# pre-process data (remove low quality probes and add KEGG pathway info)
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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Calculate CpG specific statistics + generate ordered distribution plots
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Calculate the median 5hmC values
mdbeta1 <- apply(MethOxy[sel, , 3], 1, median, na.rm=TRUE)

# order CpGs basd on beta values 
median_5hmc_ordered <- mdbeta1[order(mdbeta1, decreasing=FALSE)]

# produce ordered distribution plot w/ lines to indicate beta value signals of 5, 10 & 15%
ppi =300
png("02.Characterization_5hmC_levels/Figures/ordered_5hmC_distribution.png", width=5*ppi, height=5*ppi, res=ppi)
plot(median_5hmc_ordered, pch=16, xaxt="n", bty='l', las = 1, ylab = "Median 5hmC beta-value", xlab = "Ordered CpG loci", cex.lab = 1.3, cex.axis = 1.15)
abline(a = 0.05, b = 0, col = "purple", lwd = 2, lty = "dotdash")
abline(a = 0.1, b = 0, col = "blue", lwd = 2, lty = "dotdash")
abline(a = 0.15, b = 0, col = "red", lwd = 2, lty = "dotdash")
dev.off()

# save CpG number with minimum beta value among those w/ the top 1% 5hmC beta values 
CpG_begin = 387617 - 3876 
# save CpG number of CpG w/ highest 5hmC value 
CpG_end = 387617 
# produce ordered distribution plot, highlighting the high 5hmC CpGs (top 1% median)
png("02.Characterization_5hmC_levels/Figures/ordered_5hmC_distribution_high_5hmC.png", width=5*ppi, height=5*ppi, res=ppi)
plot(median_5hmc_ordered, pch=16, xaxt="n", bty='l', las = 1, ylab = "Median 5hmC beta-value", xlab = "Ordered CpG loci", cex.lab = 1.3, cex.axis = 1.15)
points(x=CpG_begin:CpG_end, median_5hmc_ordered[CpG_begin:CpG_end], col='red', pch=16)
abline(a = 0.141132, b = 0, col = "green", lwd = 2, lty = "dotdash")
dev.off()

