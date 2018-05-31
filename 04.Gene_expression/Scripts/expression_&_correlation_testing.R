#######################################################################################################################
# Analyze correlations between 5hmC/5mC abundance and gene expression in normal breast tissue 
# (using available gene expression data collected on NanoString nCounter panel)
#######################################################################################################################
rm(list = ls())
setwd("/Users/Owen 1/Dropbox (Christensen Lab)/NDRI_Breast_5hmC_update/")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load data 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# gene expression data 
Breast_Nanostring_Norm = read.csv("04.Gene_expression/Files/NormalBreastTissue_NormalizedData.csv", header=T, row.names=1)
cols = c(2:19)    
# convert columns to numeric
Breast_Nanostring_Norm[ , cols] = apply(Breast_Nanostring_Norm[ ,cols], 2, function(x) as.numeric(x))

# covariate file
NDRI_Covars <- read.csv("01.Data_Preprocessing/Files/NDRI_Normal_Samples.csv", header=T, row.names=1)

# Order the vector and data fram to be in the same order
Covars_ordered <- NDRI_Covars[order(rownames(NDRI_Covars)), ]
Breast_Nanostring_Norm <- Breast_Nanostring_Norm[, c(1, order(colnames(Breast_Nanostring_Norm)))]
all(rownames(Covars_ordered)==colnames(Breast_Nanostring_Norm)[2:19])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Selection for CpGs of Interest
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Maximum Liklihood Estimation data and annotation files
load("02.Characterization_5hmC_levels/Files/annotation-150214.RData")
load("02.Characterization_5hmC_levels/Files/OxyBS-Breast-FunNorm.RData")
load("02.Characterization_5hmC_levels/Files/KEGG-Key.RData") 

# Create descriptions for the KEGG pathways
pathDescript <- pathKey$pathway
names(pathDescript) <- pathKey$id 

# Replace the pathKey object with KEGG IDs
pathKey <- names(keggPaths) 
# All the unique CpGs that track to KEGG pathways
uniqueKEGG <- unique(unlist(keggPaths)) 
# Remove the "N_" and "S_"  annotation from the annotation file
geoUCSC <- gsub("^[N|S]_","",longAnnot$RelToIslandUCSC) 

# Labelling all of the 'bad' probes
strat <- rep("bad",dim(longAnnot)[1]) 
# Identifying the 'good' type I probes
strat[lmIndex[["I:0:0:0:0"]]] <- "I" 
# Identifying the 'good' type II probes
strat[lmIndex[["II:0:0:0:0"]]] <- "II" 
# Creating a unique group of CpG Island and Probe-type
strat <- paste(strat, geoUCSC, sep=":") 
# New object name moving forward
strat2 <- strat 

# Extract row information for TSS-associated CpGs
tmp <- grep("TSS",longAnnot$UCSC_RefGene_Group) 
# Adding TSS indication to relevant CpGs
strat2[tmp] <- paste(strat2[tmp],"TSS",sep=":") 
# Matching CpG ID to the stratified names
names(strat2) <- longAnnot$TargetID 

# Drops the 8 elements of 'bad' probes
strat2Index <- split(longAnnot$TargetID, strat2)[-(1:8)] 
# Number of probes per stratification
nPerStrat <- sapply(strat2Index,length) 

# Total number of CpGs
nCpGTot <- dim(MethOxy)[1] 

# (Sel)ecting only the 'good' probes
sel <- unlist(lmIndex[c("I:0:0:0:0","II:0:0:0:0")]) 
names(sel) <- rownames(MethOxy[sel, , 3])

# Attaching names to good probe indices
selIndex <- 1:length(sel)
names(selIndex) <- rownames(MethOxy[sel,,3])

# CG ID and its probe type/location/TSS status
strat2sel <- strat2[sel]
STRAT2SELINDEX <- split(1:length(sel), strat2sel)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# test correlation between average 5hmC content and epigenetic enzyme expression 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Determine the average level of 5hmC/5mC for ALL CpGs
Hydroxy_sum <- apply(MethOxy[sel, , 3], 2, sum, na.rm=TRUE)
Avg_5hmC <- (Hydroxy_sum/387617)

# visualize total levels accross al subjects 
barplot(Avg_5hmC, ylim=c(0,0.1), las=2)

# Summary statistics to accompany graph of 5hmC
summary(Avg_5hmC)

# Test associations between Epigenetic enzyme expression levels and 'total' 5hmC
cor.test(Avg_5hmC, as.numeric(Breast_Nanostring_Norm["TET1", 2:19])) 
cor.test(Avg_5hmC, as.numeric(Breast_Nanostring_Norm["TET2", 2:19])) 
cor.test(Avg_5hmC, as.numeric(Breast_Nanostring_Norm["TET3", 2:19])) 
cor.test(Avg_5hmC, as.numeric(Breast_Nanostring_Norm["DNMT1", 2:19])) 
cor.test(Avg_5hmC, as.numeric(Breast_Nanostring_Norm["IDH1", 2:19])) 
cor.test(Avg_5hmC, as.numeric(Breast_Nanostring_Norm["IDH2", 2:19])) 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# identify CpGs tracking to Nanostring genes and index 5hmC data accordiungly 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load annotation file that contains MAPINFO for CpG coordinates
load("04.Gene_expression/Files/annot.RData")

# Subset the annotation file to the CpGs used in Johnson et al
annot_full <- annot[rownames(MethOxy[sel, , 3]), ]

# get 5hmC data 
All_5hmC = t(MethOxy[, , 3])
All_5hmC_nano = All_5hmC[colnames(Breast_Nanostring_Norm)[2:19], ]
cpg = (t(All_5hmC_nano))
CpGgenes <- annot_full
CpGgenes$UCSC_RefGene_Name = as.character(CpGgenes$UCSC_RefGene_Name)

# Define function to select CpGs of interest based on gene name in 450k annotation file
grepM = function(pattern, varname, data, exact = T, identifier = ";") {
  if(exact) {
    x = data[,varname]
    tmp = grep(pattern, x)
    x.sub = data[tmp,]
    s = apply(as.matrix(x.sub[,varname]), 1, function(x) strsplit(x, identifier))
    ind = lapply(s, function(x) sum(pattern == x[[1]])>=1)
    ind.1 = unlist(ind)
    tmp.1 = tmp[ind.1]
    tmp.1
    
  }
  else { 
    tmp
  }
}
nano_string_genes <- as.character(unique(Breast_Nanostring_Norm$Gene_ID))
CpGs_Nanostring <- unlist(lapply(nano_string_genes, grepM, varname = "UCSC_RefGene_Name", data = annot_full, 
                                 exact = T, identifier = ";"))
CpGs_of_interest <- rownames(annot_full[CpGs_Nanostring, ])
CpGgenes = annot_full[CpGs_Nanostring, ]
CpGgenes$UCSC_RefGene_Name = as.character(CpGgenes$UCSC_RefGene_Name)

# calculate median expression over multiple transcripts for DNMT3A + RASSF1 
# DNMT3A
Breast_Nanostring_Norm[dim(Breast_Nanostring_Norm)[1]+1,] <- apply(Breast_Nanostring_Norm[c("DNMT3A_v2","DNMT3A_v3","DNMT3A_v4"),], 2, median)
rownames(Breast_Nanostring_Norm)[dim(Breast_Nanostring_Norm)[1]] <- "DNMT3A"
# RASSF1
Breast_Nanostring_Norm[dim(Breast_Nanostring_Norm)[1]+1,] <- apply(Breast_Nanostring_Norm[c("RASSF1_vB","RASSF1_vC","RASSF1_vH"),], 2, median)
rownames(Breast_Nanostring_Norm)[dim(Breast_Nanostring_Norm)[1]] <- "RASSF1"

# index 5hmC data for CpGs mapping to genes with expression data 
All_5hmC = t(MethOxy[CpGs_of_interest, , 3])
All_5hmC_nano = All_5hmC[colnames(Breast_Nanostring_Norm)[2:19], ]
write.csv(All_5hmC_nano, file="04.Gene_expression/Files/Normal_Breast_5hmC_CpGs_at_Nanostring_genes.csv")
# do the same for 5mC
All_5mC = t(MethOxy[CpGs_of_interest, , 2])
All_5mC_nano = All_5mC[colnames(Breast_Nanostring_Norm)[2:19], ]
write.csv(All_5mC_nano, file="04.Gene_expression/Files/Normal_Breast_5mC_CpGs_at_Nanostring_genes.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# test correlations between CpG specific 5hmC/5mC and gene expression 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Listing the output file where R will dump the results
myoutf1 = "04.Gene_expression/Files/corr_5hmC_nano_normalized_AllCpGs_breast.txt"

conOut = file(myoutf1, "w")
for(h in 1:nrow(cpg)){
  mycpg = row.names(cpg)[h]
  xx = CpGgenes[mycpg, "UCSC_RefGene_Name"]
  xx = unlist(strsplit(xx, ";"))[1]
  if(is.na(xx)) next    
  if(xx=="") next  
  se = which(Breast_Nanostring_Norm$Gene_ID==xx)
  if(length(se)==0) {
    next
  }
  else {
    for(i in 1:length(se))
    {
      corr1 = cor.test(as.numeric(Breast_Nanostring_Norm[se[i], 2:19]), as.numeric(cpg[h,]), method="s", exact = F)  
      curLine = paste(mycpg, rownames(Breast_Nanostring_Norm)[se[i]], CpGgenes[mycpg, 21], CpGgenes[mycpg, 11], CpGgenes[mycpg, 12], CpGgenes[mycpg, 23], 
                      median(as.numeric(cpg[h,])),round(corr1$estimate, 6), round(corr1$p.value, 6), sep="\t")
      writeLines(curLine, conOut)  
    }
  }
}
close(conOut)  
#writeLines(curLine, conOut)  

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# now re-run for 5mC 
All_5mC = t(MethOxy[CpGs_of_interest, , 2])
All_5mC_nano = All_5mC[colnames(Breast_Nanostring_Norm)[2:19], ]
cpg = (t(All_5mC_nano))

# Listing the output file where R will dump the results
myoutf1 = "04.Gene_expression/Files/corr_5mC_nano_normalized_AllCpGs_breast.txt"

#
conOut = file(myoutf1, "w")
for(h in 1:nrow(cpg)){
  mycpg = row.names(cpg)[h]
  xx = CpGgenes[mycpg, "UCSC_RefGene_Name"]
  xx = unlist(strsplit(xx, ";"))[1]
  if(is.na(xx)) next    
  if(xx=="") next  
  se = which(Breast_Nanostring_Norm$Gene_ID==xx)
  if(length(se)==0) {
    next
  }
  else {
    for(i in 1:length(se))
    {
      corr1 = cor.test(as.numeric(Breast_Nanostring_Norm[se[i], 2:19]), as.numeric(cpg[h,]), method="s", exact = F)  
      curLine = paste(mycpg, rownames(Breast_Nanostring_Norm)[se[i]], CpGgenes[mycpg, 21], CpGgenes[mycpg, 11], CpGgenes[mycpg, 12], CpGgenes[mycpg, 23],
                      median(as.numeric(cpg[h,])),round(corr1$estimate, 6), round(corr1$p.value, 6), sep="\t")
      writeLines(curLine, conOut)  
    }
  }
}
close(conOut)  

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# restrict to high 5hmC CpGs and breast specific genes 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
hmc <- read.table("04.Gene_expression/Files/corr_5hmC_nano_normalized_AllCpGs_breast.txt", sep = "\t", header = F, quote="", stringsAsFactors = FALSE)
mc <- read.table("04.Gene_expression/Files/corr_5mC_nano_normalized_AllCpGs_breast.txt", sep = "\t", header = F, quote="", stringsAsFactors = FALSE)
colnames(hmc) <- c("Illumina_ID", "Transcript_Variant",	"UCSC_Gene_Name",	"CHR",	"MAPINFO", "UCSC_Gene_Group",	"Median_5hmC",	"Spearman_Cor",	"Spearman_Pval")
colnames(mc) <- c("Illumina_ID", "Transcript_Variant", "UCSC_Gene_Name",	"CHR",	"MAPINFO", "UCSC_Gene_Group",	"Median_5hm",	"Spearman_Cor",	"Spearman_Pval")

# index hmc for those with high (>0.1 beta) 5hmC
top1_perc <- readRDS("02.Characterization_5hmC_levels/Files/high_5hmc_top1%_5hmC.rds")
hmc_sub <- hmc[hmc$Illumina_ID %in% top1_perc$id,]
hmc_sub2 <- hmc_sub[hmc_sub$Spearman_Pval < 0.05,]

# restrict to only breast specific genes and epigenetic enzymes 
table(hmc_sub$Transcript_Variant)
hmc_sub_2 <- hmc_sub[hmc_sub$Transcript_Variant=="RAB32" |
                       hmc_sub$Transcript_Variant=="RASSF1" |
                       hmc_sub$Transcript_Variant=="RASSF1_vB" |
                       hmc_sub$Transcript_Variant=="RASSF1_vC" |
                       hmc_sub$Transcript_Variant=="RASSF1_vH" |
                       hmc_sub$Transcript_Variant=="TWIST1" | 
                       hmc_sub$Transcript_Variant=="DNMT3A" | 
                       hmc_sub$Transcript_Variant=="DNMT3A_v2" | 
                       hmc_sub$Transcript_Variant=="DNMT3A_v3" | 
                       hmc_sub$Transcript_Variant=="DNMT3A_v4",]

# index 5mc cpgs for those that have high 5hmC
mc_sub <- mc[mc$Illumina_ID %in% top1_perc$id,]
mc_sub_2 <- mc_sub[mc_sub$Transcript_Variant=="RAB32" |
                     mc_sub$Transcript_Variant=="RASSF1" |
                     mc_sub$Transcript_Variant=="RASSF1_vB" |
                     mc_sub$Transcript_Variant=="RASSF1_vC" |
                     mc_sub$Transcript_Variant=="RASSF1_vH" |
                     mc_sub$Transcript_Variant=="TWIST1" | 
                     mc_sub$Transcript_Variant=="DNMT3A" | 
                     mc_sub$Transcript_Variant=="DNMT3A_v2" | 
                     mc_sub$Transcript_Variant=="DNMT3A_v3" |
                     mc_sub$Transcript_Variant=="DNMT3A_v4",]

# out put these as lists 
write.csv(hmc_sub_2, file = "04.Gene_expression/Files/high_5hmC_CpG_expression_correlations_5hmC.csv")
write.csv(mc_sub_2, file = "04.Gene_expression/Files/high_5mC_CpG_expression_correlations_5mC.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# make expression vs methylation plots for all high 5hmC CpGs located in genes with expression data 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# write function to generate plots given gene and cg 
exp_vs_meth_plot <- function(gene, CG, meth_type_df, meth_type_name){
  exp <- as.data.frame(t(Breast_Nanostring_Norm[rownames(Breast_Nanostring_Norm) == gene, 2:19]))
  meth <- as.data.frame(t(as.data.frame(t(meth_type_df[, colnames(meth_type_df) == CG]))))
  z <- cbind(exp, meth)
  colnames(z) <- c("expression", "methylation")
  z$expression <- as.numeric(z$expression)
  z$methylation <- as.numeric(z$methylation)
  # make plot
  ggplot(z, aes(methylation, expression)) + 
    geom_smooth(method = "lm")  + 
    geom_point() + xlab(paste0(meth_type_name, " (", CG,")")) + 
    ylab(paste0("-Log2 ", "(", gene, ")", " expression")) + 
    theme(legend.key.size = unit(1, "cm"), 
        panel.grid.major = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "black", size = 0.6, linetype = "solid"), 
        panel.background = element_blank(),
        axis.text.x=element_text(colour="black", size = 12, hjust = 1), 
        axis.text.y=element_text(colour="black", size = 12), 
        axis.title.x=element_text(colour="black", size = 12), 
        axis.title.y=element_text(colour="black", size = 12), 
        legend.key = element_blank())
}

# apply this to make plots for all CGs of interest 
library(ggpubr)
ppi = 300

# save vectors of all gene names and GpGs of interest 
genes <- c(rep("RAB32", 6), 
           rep(c("RASSF1", "RASSF1_vB", "RASSF1_vC", "RASSF1_vH"), 2),
           rep("TWIST1", 4),
           rep(c("DNMT3A", "DNMT3A_v2", "DNMT3A_v3", "DNMT3A_v4"), 6))
CpGs <- c("cg02664328", "cg01892997", "cg24744430", "cg18987220", "cg23267550", "cg01915609", 
          "cg19854901", "cg19854901", "cg19854901", "cg19854901", "cg24049629", "cg24049629", "cg24049629", "cg24049629", 
          "cg26279021", "cg14391419", "cg10126205", "cg27334919", 
          "cg23009818", "cg23009818", "cg23009818", "cg23009818", 
          "cg17207266", "cg17207266", "cg17207266", "cg17207266", 
          "cg17742416", "cg17742416", "cg17742416", "cg17742416", 
          "cg00050692", "cg00050692", "cg00050692", "cg00050692", 
          "cg20702417", "cg20702417", "cg20702417", "cg20702417", 
          "cg05896193", "cg05896193", "cg05896193", "cg05896193")
genes[19]
CpGs[19]
# loop over CpGs/genes they are located witin, and plot 5hmC + 5mC vs expression for each 
for(i in 1:length(genes)){
  p1 <- exp_vs_meth_plot(genes[i], CpGs[i], All_5hmC_nano, "5hmC")
  p2 <- exp_vs_meth_plot(genes[i], CpGs[i], All_5mC_nano, "5mC")
  png(paste0("04.Gene_expression/Figures/", genes[i], " vs ", CpGs[i], ".png"), width=7*ppi, height=3*ppi, res=ppi)
  #print(ggarrange(p1, p2, labels = c("A", "B"), ncol = 2, nrow = 1))
  print(ggarrange(p1, p2, ncol = 2, nrow = 1))
  dev.off()
}
