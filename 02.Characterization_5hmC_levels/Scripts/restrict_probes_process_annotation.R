#######################################################################################################################
# restrict 450K data to high quality probes 
# clean up probe annotation from annotation file to assist w/ analysis 
#######################################################################################################################
library(GenomicRanges)
rm(list = ls())
setwd("/Users/Owen 1/Dropbox (Christensen Lab)/NDRI_Breast_5hmC copy/")

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
# pre-process data to remove low quality probes and add KEGG pathway info
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
# add strat column to annotation 
longAnnot$probe_type <- strat

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

# are anty Chen and SNP probes in remaining annotation data 
table(AnnotSel$excludeChen)
table(AnnotSel$SNPinProbe)
# no. so none need be removed 

# add MAPINFO to annotation file using Illumina 450K file 
load("05.Nano_String/Files/annot.RData")
annot2 <- annot[match(AnnotSel$TargetID, annot$Name),]
all(AnnotSel$TargetID == annot2$Name)
AnnotSel$MAPINFO <- annot2$MAPINFO

# index methylation data for good probes selected for above 
MethOxy_2 <- MethOxy[sel, , ]

# save annotation file 
save(AnnotSel, file = "02.Characterization_5hmC_levels/Files/annotation_good_probes.Rdata")
# save methylation data indexed for good probes 
save(MethOxy_2, file = "02.Characterization_5hmC_levels/Files/MethOxy_FunNorm_good_probes.Rdata")

