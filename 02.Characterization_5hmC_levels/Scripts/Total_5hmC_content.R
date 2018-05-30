#######################################################################################################################
# explore relation of total 5hmC content with: 
# - total 5mC content, 
# - proportions of putative cell types (identified by RefFreeEWAS), and 
# - epigenetic enzyme expression
#######################################################################################################################
rm(list = ls())
library(ggplot2)
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
# estimated cell proportionsfrom RefFreeEWAS
cp <- read.csv("02.Characterization_5hmC_levels/Files/NDRI_Cell_Proportions.csv")

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
# visualize total 5hmC/5mC content across samples 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Determine the average level of 5hmC/5mC per CpG across all subjects 
#### hmC
Hydroxy_sum <- apply(MethOxy[sel, , 3], 2, sum, na.rm=TRUE)
Avg_5hmC <- (Hydroxy_sum/387617)
#saveRDS(Hydroxy_sum, file = "02.Characterization_5hmC_levels/Files/subject_total_5hmC_content.rds")
#### mC
methly_sum <- apply(MethOxy[sel, , 2], 2, sum, na.rm=TRUE)
Avg_5mC <- (methly_sum/387617)
#### C
free_sum <- apply(MethOxy[sel, , 1], 2, sum, na.rm=TRUE)

# visualize 5hmC across subjects 
ppi = 300
png("02.Characterization_5hmC_levels/Figures/Total_5hmC_levels.png", height=5*ppi, width=4.5*ppi, res=ppi)
barplot(Avg_5hmC, ylim=c(0,0.06), las=2, col = "indianred1", ylab = "Average 5hmC beta-value", xlab = "Samples", cex.lab = 1.2)
dev.off()
png("02.Characterization_5hmC_levels/Figures/Total_5mC_levels.png", height=5*ppi, width=4.5*ppi, res=ppi)
barplot(Avg_5mC, ylim=c(0,0.06), las=2, col = "indianred1", ylab = "Average 5mC beta-value", xlab = "Samples", cex.lab = 1.2)
dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# examine relationship between total 5hmC and 5mC content 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# calculate % 5hmC of total cytosine content for eacvh sample 
percent_5hmC <- rep(NA, length(Hydroxy_sum))
percent_totals <- matrix(NA, nrow = length(Hydroxy_sum), ncol = 3)
for(i in 1:nrow(percent_totals)){
  percent_totals[i, 1] <- (Hydroxy_sum[i])/(Hydroxy_sum[i]+methly_sum[i]+free_sum[i])
  percent_totals[i, 2] <- (methly_sum[i])/(Hydroxy_sum[i]+methly_sum[i]+free_sum[i])
  percent_totals[i, 3] <- (free_sum[i])/(Hydroxy_sum[i]+methly_sum[i]+free_sum[i])
  rownames(percent_totals) <- names(Hydroxy_sum)
  colnames(percent_totals) <- c("hmC", "mC", "C")
  percent_totals <- as.data.frame(percent_totals)
}

# test the correlation and assocition between 5hmC and 5mC 
cor.test(percent_totals[,1], percent_totals[,2])
lm1 <- summary(lm(percent_totals[,1] ~ percent_totals[,2]))
lm1$coefficients[2,1]
lm1$coefficients[2,4]
sd(percent_totals[,2])

# barplot of total 5hmC accross al subjects vs 5mC in all subjects 
png("02.Characterization_5hmC_levels/Figures/total_5hmc_5mC_barplot.png", height=5*ppi, width=4.5*ppi, res=ppi)
barplot(apply(percent_totals[,c(1:2)], 2, mean), las = 1, col = c("dodgerblue", "indianred1"), ylab = "Total 5(h)mC across all samples")
dev.off()
# histogram of total 5hmC distribution across all subjects 
png("02.Characterization_5hmC_levels/Figures/percent_5hmC_hist.png", height=5*ppi, width=4.5*ppi, res=ppi)
hist(Avg_5hmC, col = "dodgerblue", las = 1, xlab = "Total 5hmC per CpG", main = NULL, cex.lab = 1.2)
dev.off()
# histogram of total 5mC distribution across all subjects 
png("02.Characterization_5hmC_levels/Figures/percent_5mC_hist.png", height=5*ppi, width=4.5*ppi, res=ppi)
hist(Avg_5mC, col = "indianred1", las = 1, xlab = "Totla 5mC per CpG", main = NULL, cex.lab = 1.2)
dev.off()

# write short function to establish desired theme that can be re-used for ggplot plots 
my_theme <- function() {
  p <- theme(legend.key.size = unit(1, "cm"), 
             panel.grid.major = element_blank(), 
             panel.border = element_rect(fill = NA, colour = "black", size = 0.6, linetype = "solid"), 
             panel.background = element_blank(), 
             axis.text.x=element_text(colour="black", size = 11, hjust = 1), 
             axis.text.y=element_text(colour="black", size = 11),         
             axis.title.x=element_text(colour="black", size = 11),         
             axis.title.y=element_text(colour="black", size = 11),         
             legend.key = element_blank()
  )
  return (p)
}

# plot % of total C for 5mC vs 5hmC 
png("02.Characterization_5hmC_levels/Figures/Total_5hmC_vs_total_5mC.png", height=5*ppi, width=5*ppi, res=ppi)
p <- ggplot(percent_totals, aes(hmC, mC)) + 
  geom_smooth(method = "lm") + 
  geom_point() + 
  xlab("Total 5hmC") + 
  ylab("Total 5mC")+ 
  #annotate("text", x = 0.04, y = 0.445, label = "Beta = -0.28") +
  my_theme()
dev.off()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# examine relationship between total 5hmC proportions of putative cell types (identified by RefFreeEWAS)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# linear regression of total 5hmC per sample and cell type proportions 
summary(lm(Avg_5hmC ~ cp$X1))
summary(lm(Avg_5hmC ~ cp$X2))

# plot cell type 1 across subjects 
png("03.Cellular_proportions/Figures/NDRI_breast_cell_type_props_1_barplot.png", height=5*ppi, width=5*ppi, res=ppi)
barplot(cp$X1, ylim = c(0,1), col = "indianred", main = "Cell type 1", xlab = "Sample", ylab = "Cellular proportion", las = 1)
dev.off()

png("03.Cellular_proportions/Figures/NDRI_breast_cell_type_props_2_barplot.png", height=5*ppi, width=5*ppi, res=ppi)
barplot(cp$X2, ylim = c(0,1), col = "indianred", main = "Cell type 2", xlab = "Sample", ylab = "Cellular proportion", las = 1)
dev.off()

# plot % total 5hmC against cellular proportions for all subjects 
#### putative cell type 1
png("02.Characterization_5hmC_levels/Total_5hmC_vs_cell_type_1.png", height=5*ppi, width=5*ppi, res=ppi)
p <- ggplot(percent_totals, aes(cp$X1, Avg_5hmC)) + 
  geom_smooth(method = "lm") + 
  geom_point() + 
  ylab("Total 5hmC") + 
  xlab("Cellular proportion - Putative cell type 1") + 
  #xlim(0,1) +
  ylim(0, 0.045) +
  #annotate("text", x = 0.04, y = 0.445, label = "Beta = -0.28") +
  my_theme
dev.off()

#### putative cell type 2
png("02.Characterization_5hmC_levels/Total_5hmC_vs_cell_type_2.png", height=5*ppi, width=5*ppi, res=ppi)
p <- ggplot(percent_totals, aes(cp$X2, Avg_5hmC)) + 
  geom_smooth(method = "lm") + 
  geom_point() + 
  ylab("Total 5hmC") + 
  xlab("Cellular proportion - Putative cell type 2")+ 
  #xlim(0,1) +
  ylim(0, 0.045) +
  #annotate("text", x = 0.04, y = 0.445, label = "Beta = -0.28") +
  my_theme()
dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# examine relationship between total 5hmC proportions and epigenetic enzyme expression
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# read in normalized expression data from Nanostring panel 
Breast_Nanostring_Norm = read.csv("04.Nano_String/Files/NormalBreastTissue_NormalizedData.csv", header=T, row.names=1)
cols = c(2:19)    
Breast_Nanostring_Norm[ , cols] = apply(Breast_Nanostring_Norm[ ,cols], 2, function(x) as.numeric(x))

# Load in covariate file for NDRI samples to examine association between Age and BMI
NDRI_Covars <- read.csv("01.Data_Preprocessing/Files/NDRI_Normal_Samples.csv", header=T, row.names=1)

#Order the vector and data fram to be in the same order
Covars_ordered <- NDRI_Covars[order(rownames(NDRI_Covars)), ]
Breast_Nanostring_Norm <- Breast_Nanostring_Norm[, c(1, order(colnames(Breast_Nanostring_Norm)))]
all(rownames(Covars_ordered)==colnames(Breast_Nanostring_Norm)[2:19])

# index expression data for expression of epigenetic enzymes 
exp <- Breast_Nanostring_Norm[Breast_Nanostring_Norm$Gene_ID == "TET1" |
                                Breast_Nanostring_Norm$Gene_ID == "TET2" |
                                Breast_Nanostring_Norm$Gene_ID == "TET3" |
                                Breast_Nanostring_Norm$Gene_ID == "DNMT1" |
                                Breast_Nanostring_Norm$Gene_ID == "DNMT3A" |
                                Breast_Nanostring_Norm$Gene_ID == "DNMT3B",]

# remove geneID columns at start and end 
exp <- exp[,2:19]

# average expression across DNMT3A transcripts 
exp[dim(exp)[1]+1,] <- apply(exp[c("DNMT3A_v2","DNMT3A_v3","DNMT3A_v4"),], 2, median)
rownames(exp)[dim(exp)[1]] <- "DNMT3A"
exp

# quick visualizations of relationships 
plot(percent_totals[,"hmC"], as.numeric(exp["TET1",]), main = "TET1")
plot(percent_totals[,"hmC"], as.numeric(exp["TET2",]), main = "TET2")
plot(percent_totals[,"hmC"], as.numeric(exp["TET3",]), main = "TET3")
plot(percent_totals[,"hmC"], as.numeric(exp["DNMT3A",]), main = "DNMT3A")
plot(percent_totals[,"hmC"], as.numeric(exp["DNMT3B",]), main = "DNMT3B")
plot(percent_totals[,"hmC"], as.numeric(exp["DNMT1",]), main = "DNMT1")

# test correlations between enzymes expression and total 5hmC levels 
cor.test(percent_totals[,"hmC"], as.numeric(exp["TET1",]), method = "spearman")
cor.test(percent_totals[,"hmC"], as.numeric(exp["TET2",]), method = "spearman")
cor.test(percent_totals[,"hmC"], as.numeric(exp["TET3",]), method = "spearman")
cor.test(percent_totals[,"hmC"], as.numeric(exp["DNMT3A",]), method = "spearman")
cor.test(percent_totals[,"hmC"], as.numeric(exp["DNMT3B",]), method = "spearman")
cor.test(percent_totals[,"hmC"], as.numeric(exp["DNMT1",]), method = "spearman")

# fit linear models to test association between expression and 5hmC levels 
summary(lm(percent_totals[,"hmC"] ~ as.numeric(exp["TET1",])))
summary(lm(percent_totals[,"hmC"] ~ as.numeric(exp["TET2",])))
summary(lm(percent_totals[,"hmC"] ~ as.numeric(exp["TET3",])))
summary(lm(percent_totals[,"hmC"] ~ as.numeric(exp["DNMT3A",])))
summary(lm(percent_totals[,"hmC"] ~ as.numeric(exp["DNMT3B",])))
summary(lm(percent_totals[,"hmC"] ~ as.numeric(exp["DNMT1",])))
