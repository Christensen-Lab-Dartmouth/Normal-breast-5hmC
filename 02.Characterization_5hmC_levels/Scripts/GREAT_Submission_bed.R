#######################################################################################################################
# produce BED file for high 5hmC CpGs 
#######################################################################################################################
rm(list=ls())
library(GenomicRanges)
setwd("/Users/Owen/Dropbox (Christensen Lab)/NDRI_Breast_5hmC copy/")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load Data 
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
STRAT2SELINDEX <- split(1:length(sel), strat2sel)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# generate GRanges object and use this to produce BED file
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
MethOxy_breast <- MethOxy
#Calculate the average (mean) of 5hmC across 18 breasts to be plotted
mdbeta1_5hmC_breast <- as.data.frame(apply(MethOxy_breast[sel, , 3], 1, median, na.rm=TRUE))
colnames(mdbeta1_5hmC_breast) = "Median_5hmC"
#Check the number of CpGs that meet the criteria for high 5hmC
sum(mdbeta1_5hmC_breast$Median_5hmC>= 0.1) #14,733 CpG sites have a median value at least of 0.1
#Create a subset of "high mean 5hmC" probes
High_Breast_5hmC = subset(mdbeta1_5hmC_breast, Median_5hmC>= 0.1)
Hydro_CpGs = rownames(High_Breast_5hmC)
rm(High_Breast_5hmC)

# high 5hmC probes 
high_5hmC_CpGs_2 <- readRDS("02.Characterization_5hmC_levels/Files/high_5hmc_top1%_5hmC.rds")

# Load annotation file that contains MAPINFO for CpG coordinates
load("02.Characterization_5hmC_levels/Files/HumanMethylation450_15017482_v1-2.Rdata")
# Restrict the annotation file to only those CpGs used in the Glioma analysis
annot_breast <- annot[annot$Name %in% high_5hmC_CpGs_2$id, ] 
# Create a new variable with MAPINFO for CpG start
annot_breast$hg19_start = annot_breast$MAPINFO
# Create another new variable with MAPINFO for CpG end
annot_breast$hg19_end = annot_breast$MAPINFO + 1
# Create a new column that contans chromosome as 'chr' variable
annot_breast$chr = paste("chr", annot_breast$CHR, sep='')
annot_breast_sub <- annot_breast[, c("Name", 'UCSC_RefGene_Name','UCSC_RefGene_Group', 'Relation_to_UCSC_CpG_Island', 'hg19_start', 'hg19_end', 'chr')]

# Use GRanges package to transform data frame
high_5hmc_gr <- makeGRangesFromDataFrame(annot_breast_sub, keep.extra.columns=TRUE, ignore.strand=TRUE, seqnames.field = "chr", 
                                              start.field="hg19_start", end.field="hg19_end")
#Convert GenomicRange objects to BED files
Breast_High_5hmC = data.frame(cbind(chrom=as.vector(seqnames(high_5hmc_gr)),start=start(high_5hmc_gr),end=end(high_5hmc_gr)))
write.table(Breast_High_5hmC, file="02.Characterization_5hmC_levels/Files/Bed_files_GREAT/Breast_genes_high_5hmC.bed", quote=F, sep="\t", row.names=F, col.names=F)
save(high_5hmc_gr, file = "02.Characterization_5hmC_levels/Files/high_5hmc_gr.Rdata")

# load in names of CGs w/ differential meth. sig. value <0.05 in TCGA data set
high_5hmC_sub <- readRDS("06.Cancer_Comparison/Files/high_5hmC_CGs_below_bonf_05_in_TCGA.RData")
# generate BED file to be used in enrichment analysis for these CGs to 
high_5hmc_gr_sub <- high_5hmc_gr[match(high_5hmC_sub, high_5hmc_gr@elementMetadata$Name)]
#Convert GenomicRange objects to BED files
Breast_High_5hmC_sub = data.frame(cbind(chrom=as.vector(seqnames(high_5hmc_gr_sub)),start=start(high_5hmc_gr_sub),end=end(high_5hmc_gr_sub)))
write.table(Breast_High_5hmC_sub, file="02.Characterization_5hmC_levels/Files/Breast_genes_high_5hmC_low_P_TCGA.bed", quote=F, sep="\t", row.names=F, col.names=F)



