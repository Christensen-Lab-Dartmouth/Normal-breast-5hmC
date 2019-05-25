#######################################################################################################################
# Test high 5hmC CpGs for statsitical enrichment among genomic features (such as enahncer regions - H3K4me1) as 
# identified from Roadmap Epigeniomics project ChIP-seq experiments, performed in various breast/breast derived 
# tissues (myoepithelial cells, human mammary epithelial cells (HMECs), variants HMECs (vHMECs))
#######################################################################################################################
library(GenomicRanges)
library(rtracklayer)
rm(list = ls())

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load Data 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# annotation file
setwd("/Users/Owen 1/Dropbox (Christensen Lab)/NDRI_Breast_5hmC_update/")
load("02.Characterization_5hmC_levels/Files/annotation_good_probes.Rdata")
annot <- AnnotSel
# load methylation data
load("02.Characterization_5hmC_levels/Files/MethOxy_FunNorm_good_probes.Rdata")
# high 5hmC probes 
top1_perc <- readRDS("02.Characterization_5hmC_levels/Files/high_5hmc_top1%_5hmC.rds")
high_5hmC <- top1_perc

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# make Grange object for good probes from 450K array  
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Generate new variables for MAPINFO. A probe lists one location, to get genomic ranges we need two locations
annot$hg19_start = annot$MAPINFO
annot$hg19_end = annot$MAPINFO + 1
#Make the 'CHR' integer variable a character
annot$chr = as.character(annot$chr)
names(annot) #To identify columns of interest so that they do not conflict (multiple genomic locations)
annot_sub = annot[, c(1, 3, 27, 28)]

#Create a 'GRanges' object from the Illumina annotation file 
Illumina_gr <- makeGRangesFromDataFrame(annot_sub, keep.extra.columns=TRUE, ignore.strand=TRUE, seqnames.field = "chr", start.field="hg19_start", end.field="hg19_end")
#Print the GRange object to see what it looks like
names(Illumina_gr) <- annot_sub$TargetID
Illumina_gr

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# make Grange object for high 5hmC probes
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# add high 5hmC indicator 
AnnotSel$top_5hmC <- NA
AnnotSel$top_5hmC[AnnotSel$TargetID %in% high_5hmC$id] <- 1
AnnotSel$top_5hmC[!AnnotSel$TargetID %in% high_5hmC$id] <- 0
table(AnnotSel$top_5hmC)
table(AnnotSel$top_5hmC, AnnotSel$probe_type)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Perform enrichment test on each Roadmap Epigenomics data set
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# prepare a function that will use M-H test to test 5hmC enirhcment with each genomic feature
test_enrichment <- function(mark_object, mark_name){
  # index for sites w/ -log10P > 2 as roadmap advise this provides good signal/noice separation 
  mark_bs <- mark_object[score(mark_object)>=2, ]
  # drop XY & MT coords 
  non_auto_coords <- which(seqnames(mark_bs)=="chrY" | seqnames(mark_bs)=="chrX" | seqnames(mark_bs)=="chrM")  
  mark_bs <- mark_bs[-non_auto_coords,]
  # find the overlapping regions between the high confidence peaks and the CpGs on the 450k array 
  overlaps <- findOverlaps(mark_bs, Illumina_gr)
  # get the indicies of the overlapping sites in the 450K annotation file 
  indicies <- subjectHits(overlaps)
  # add dummy variable to annotation file for CpGs overlapping with this mark 
  AnnotSel[,paste0(mark_name)] <- 0
  AnnotSel[,paste0(mark_name)][indicies] <- 1
  # make 3D table of 5hmC sites vs mark vs CpG regions 
  MH_table <- table(AnnotSel$top_5hmC, AnnotSel[,paste0(mark_name)], AnnotSel$CpG_Regions)
  # perform mantel-hansel test to test enrichment of 5hmC within mark 
  mantelhaen.test(MH_table, exact=TRUE)
}

########################  E027 - Breast myoepithelial cells ########################  

#Note on consolidated ChIP-seq data:
#Signal processing engine of MACSv2.0.10 peak caller to generate genome-wide signal coverage tracks
#Negative log10 of the Poisson p-value of ChIP-seq or DNase counts relative to expected background counts
#These signal confidence scores provides a measure of statistical significan of the observed enrichment.


# get ChIP-seq data for available chromatin marks: save paths to data online
#E027_paths <- list("http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E027-H3K4me1.pval.signal.bigwig", 
                   #"http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E027-H3K4me3.pval.signal.bigwig",
                   #"http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E027-H3K9ac.pval.signal.bigwig",
                   #"http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E027-H3K9me3.pval.signal.bigwig",
                   #"http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E027-H3K27me3.pval.signal.bigwig",
                   #"http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E027-H3K36me3.pval.signal.bigwig")

# import BIGWIG files from consolidated epigenomes from ROADMAP website 
#E027_data <- lapply(E027_paths, import.bw)

# save data so we don't have to download again 
#save(E027_data, file = "/Users/Owen/Dropbox (Christensen Lab)/NDRI_Breast_5hmC_update/05.Enrichment_Analyses/roadmap_breast_myoepith_H3_ChIP_seq.Rdata")
load("04.Enrichment_Analyses/Roadmap_Epigenomics/Files/roadmap_breast_myoepith_ChIP_seq.Rdata")

# check annotation and IlluminGR object are in same order so that code in function below will work 
all(annot_sub$TargetID == names(Illumina_gr))

# apply function for each mark to test the enrichment 
#MH_h3k4me1 <- test_enrichment(E027_data[[1]], "h3k4me1")
#Error in uniroot(function(t) mn2x2xk(1/t) - x, c(.Machine$double.eps,  : 
   #f() values at end points not of opposite sign
   #Called from: uniroot(function(t) mn2x2xk(1/t) - x, c(.Machine$double.eps, 1))
MH_h3k4me3 <- test_enrichment(E027_data[[2]], "h3k4me3")
MH_h3k9ac <- test_enrichment(E027_data[[3]], "h3k9ac")
MH_h3k9me3 <- test_enrichment(E027_data[[4]], "h3k9me3")
MH_h3k27me3 <- test_enrichment(E027_data[[5]], "h3k27me3")
MH_h3k36me3 <- test_enrichment(E027_data[[6]], "h3k36me3")

# generate and output results table of 5hmC enrichment with each genomic feature 
tab <- matrix(NA, ncol = 5, nrow = 5)
colnames(tab) <- c("P-value", "Odds ratio (95% CI)", "OR", "LB", "UB")
rownames(tab) <- c("H3K4me3", "H3K36me3", "H3K9ac", "H3K9me3", "H3K27me3")
results <- list(MH_h3k4me3, MH_h3k36me3, MH_h3k9ac, MH_h3k9me3, MH_h3k27me3)
for(i in 1:length(results)){
  tab[i,1] <- paste0(results[[i]]$p.value)
  tab[i,2] <- paste0(format(round(results[[i]]$estimate, digits = 2), nsmall = 2), " (", 
                     format(round(results[[i]]$conf.int[1], digits = 2), nsmall = 2), "-", 
                     format(round(results[[i]]$conf.int[2], digits = 2), nsmall = 2), ")")
  tab[i,3] <- results[[i]]$estimate
  tab[i,4] <- results[[i]]$conf.int[1]
  tab[i,5] <- results[[i]]$conf.int[2]
}
write.csv(tab, file = "04.Enrichment_Analyses/Roadmap_Epigenomics/Files/roadmap_breast_myoepith_5hmC_enrichment_results.csv")

# clean workspace
rm(E027_paths, E027_data, tab, results, MH_h3k4me3, MH_h3k9ac, MH_h3k9me3, MH_h3k27me3, MH_h3k36me3)

########################  E028 - Variant Human Mammary Epithelial cells (vHMEC) ########################  

# E028 = numeric epigenome identifier for consolidated data from variant HMEC
#E028_paths <- list("http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E028-DNase.pval.signal.bigwig",
                   #"http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E028-H3K4me1.pval.signal.bigwig",
                   #"http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E028-H3K4me3.pval.signal.bigwig",
                   #"http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E028-H3K9me3.pval.signal.bigwig",
                   #"http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E028-H3K27me3.pval.signal.bigwig",
                   #"http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E028-H3K36me3.pval.signal.bigwig")
#E028_data <- lapply(E028_paths, import.bw)
#save(E028_data, file = "03.Enrichment_Analyses/roadmap_vHMEC_ChIP_seq.Rdata")
load("04.Enrichment_Analyses/Roadmap_Epigenomics/Files/roadmap_vHMEC_ChIP_seq.Rdata")

# apply function for each mark to test the enrichment 
MH_DNase <- test_enrichment(E028_data[[1]], "DNase")
MH_h3k4me1 <- test_enrichment(E028_data[[2]], "h3k4me1")
MH_h3k4me3 <- test_enrichment(E028_data[[3]], "h3k4me3")
MH_h3k9me3 <- test_enrichment(E028_data[[4]], "h3k9me3")
MH_h3k27me3 <- test_enrichment(E028_data[[5]], "h3k27me3")
MH_h3k36me3 <- test_enrichment(E028_data[[6]], "h3k36me3")

# generate and output results table of 5hmC enrichment with each genomic feature 
tab <- matrix(NA, ncol = 5, nrow = 6)
colnames(tab) <- c("P-value", "Odds ratio (95% CI)", "OR", "LB", "UB")
rownames(tab) <- c("DNase", "h3k4me1", "H3K4me3", "H3K36me3", "H3K9me3", "H3K27me3")
results <- list(MH_DNase, MH_h3k4me1, MH_h3k4me3, MH_h3k36me3, MH_h3k9me3, MH_h3k27me3)
for(i in 1:length(results)){
  tab[i,1] <- paste0(results[[i]]$p.value)
  tab[i,2] <- paste0(format(round(results[[i]]$estimate, digits = 2), nsmall = 2), " (", 
                     format(round(results[[i]]$conf.int[1], digits = 2), nsmall = 2), "-", 
                     format(round(results[[i]]$conf.int[2], digits = 2), nsmall = 2), ")")
  tab[i,3] <- results[[i]]$estimate
  tab[i,4] <- results[[i]]$conf.int[1]
  tab[i,5] <- results[[i]]$conf.int[2]
}

write.csv(tab, file = "04.Enrichment_Analyses/Roadmap_Epigenomics/Files/roadmap_vHMEC_5hmC_enrichment_results.csv")

# clean workspace
rm(i, tab, results, MH_DNase, MH_h3k4me1, MH_h3k4me3, MH_h3k9me3, MH_h3k27me3, MH_h3k36me3)

######################## E119 - Human Mammary Epithelial cells (HMEC) ########################  

# E119 = numeric epigenome identifier for consolidated data from HMEC
#E119_paths <- list("http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E119-DNase.pval.signal.bigwig",
                   #"http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E119-H2A.Z.pval.signal.bigwig",
                   #"http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E119-H3K4me1.pval.signal.bigwig",
                   #"http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E119-H3K4me2.pval.signal.bigwig",
                   #"http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E119-H3K4me3.pval.signal.bigwig",
                   #"http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E119-H3K9ac.pval.signal.bigwig",
                   #"http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E119-H3K9me3.pval.signal.bigwig",
                   #"http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E119-H3K27ac.pval.signal.bigwig",
                   #"http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E119-H3K27me3.pval.signal.bigwig",
                   #"http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E119-H3K36me3.pval.signal.bigwig", 
                   #"http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E119-H3K79me2.pval.signal.bigwig",
                   #"http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E119-H4K20me1.pval.signal.bigwig")

#E119_data <- lapply(E119_paths, import.bw)
#save(E119_data, file = "03.Enrichment_Analyses/roadmap_HMEC_ChIP_seq.Rdata")
load("04.Enrichment_Analyses/Roadmap_Epigenomics/Files/roadmap_HMEC_ChIP_seq.Rdata")

# apply function for each mark to test the enrichment 
MH_DNase <- test_enrichment(E119_data[[1]], "DNase")
MH_H2A.Z <- test_enrichment(E119_data[[2]], "H2A.Z")
MH_h3k4me1 <- test_enrichment(E119_data[[3]], "h3k4me1")
MH_h3k4me2 <- test_enrichment(E119_data[[4]], "h3k4me2")
MH_h3k4me3 <- test_enrichment(E119_data[[5]], "h3k4me3")
MH_h3k9ac <- test_enrichment(E119_data[[6]], "h3k9ac")
MH_h3k9me3 <- test_enrichment(E119_data[[7]], "h3k9me3")
MH_h3k27ac <- test_enrichment(E119_data[[8]], "h3k27ac")
MH_h3k27me3 <- test_enrichment(E119_data[[9]], "h3k27me3")
MH_h3k36me3 <- test_enrichment(E119_data[[10]], "h3k36me3")
MH_h3k79me2 <- test_enrichment(E119_data[[11]], "h3k79me2")
MH_h4k20me1 <- test_enrichment(E119_data[[12]], "h3k20me1")

# generate and output results table of 5hmC enrichment with each genomic feature 
tab <- matrix(NA, ncol = 5, nrow = 12)
colnames(tab) <- c("P-value", "Odds ratio (95% CI)", "OR", "LB", "UB")
rownames(tab) <- c("DNase", "H3K4me1", "H3K4me2", "H3K4me3", "H3K9ac", "H3K27ac", "H3K36me3", "H3K79me2", "H3K20me1", 
                   "H3K9me3",  "H3K27me3", 
                   "H2A.Z")
results <- list(MH_DNase, MH_h3k4me1, MH_h3k4me2, MH_h3k4me3, MH_h3k9ac, MH_h3k27ac, MH_h3k36me3, MH_h3k79me2, MH_h4k20me1,
                MH_h3k9me3, MH_h3k27me3, 
                MH_H2A.Z)
for(i in 1:length(results)){
  tab[i,1] <- paste0(results[[i]]$p.value)
  tab[i,2] <- paste0(format(round(results[[i]]$estimate, digits = 2), nsmall = 2), " (", 
                     format(round(results[[i]]$conf.int[1], digits = 2), nsmall = 2), "-", 
                     format(round(results[[i]]$conf.int[2], digits = 2), nsmall = 2), ")")
  tab[i,3] <- results[[i]]$estimate
  tab[i,4] <- results[[i]]$conf.int[1]
  tab[i,5] <- results[[i]]$conf.int[2]
}
write.csv(tab, file = "04.Enrichment_Analyses/Roadmap_Epigenomics/Files/roadmap_HMEC_5hmC_enrichment_results.csv")

