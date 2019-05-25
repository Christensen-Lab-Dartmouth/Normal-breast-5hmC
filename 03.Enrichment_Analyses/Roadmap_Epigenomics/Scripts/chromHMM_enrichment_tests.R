#######################################################################################################################
# 5hmC enrichment among 15-state model ChromHMM chromatin states
#######################################################################################################################

library(LOLA)
library(GenomicRanges)
library(rtracklayer)
rm(list = ls())

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load Data 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("/Users/Owen 1/Dropbox (Christensen Lab)/NDRI_Breast_5hmC_update/")
dir1 <- "02.Characterization_5hmC_levels/Files/"
dir2 <- "03.Enrichment_Analyses/Roadmap_Epigenomics/Files/"

# annotation file
load(paste0(dir1, "annotation_good_probes.Rdata"))
annot <- AnnotSel
# load methylation data
load(paste0(dir1, "MethOxy_FunNorm_good_probes.Rdata"))
# high 5hmC probes 
top1_perc <- readRDS(paste0(dir1, "high_5hmc_top1%_5hmC.rds"))
high_5hmC <- top1_perc

E027_bed <- LOLA::readBed(paste0(dir2, "E027_15_coreMarks_mnemonics.bed"))
E119_bed <- LOLA::readBed(paste0(dir2, "E119_15_coreMarks_mnemonics.bed"))

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
AnnotSel$top_5hmC[AnnotSel$TargetID %in% high_5hmC$ID] <- 1
AnnotSel$top_5hmC[!AnnotSel$TargetID %in% high_5hmC$ID] <- 0
table(AnnotSel$top_5hmC)
table(AnnotSel$top_5hmC, AnnotSel$probe_type)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define chromatin states and run enrichment for each one 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# drop XY coords 
#### E027
E027_non_auto_coords <- which(seqnames(E027_bed)=="chrY" | 
                                seqnames(E027_bed)=="chrX" | 
                                seqnames(E027_bed)=="chrM")
E027_bed <- E027_bed[-E027_non_auto_coords,]
#### E119
E119_non_auto_coords <- which(seqnames(E119_bed)=="chrY" | 
                                seqnames(E119_bed)=="chrX" | 
                                seqnames(E119_bed)=="chrM")
E119_bed <- E119_bed[-E119_non_auto_coords,]

# define chromHMM states 
states <- c("1_TssA", "2_TssAFlnk", "3_TxFlnk", "4_Tx", "5_TxWk", "6_EnhG", "7_Enh", 
            "8_ZNF/Rpts", "9_Het", "10_TssBiv", "11_BivFlnk", "12_EnhBiv", "13_ReprPC", 
            "14_ReprPCWk", "15_Quies")

# define function to run enrichment test using MH test
test_enrichment_MH <- function(state, state_coords){
  state_coords_sub <- state_coords[which(names(state_coords)==state)]
  overlaps <- findOverlaps(state_coords_sub, Illumina_gr)
  indicies <- subjectHits(overlaps)
  AnnotSel[,paste0(state)] <- 0
  AnnotSel[,paste0(state)][indicies] <- 1
  MH_table <- table(factor(AnnotSel$top_5hmC, levels = c("1","0")), 
                    factor(AnnotSel[,paste0(state)], levels = c("1","0")), 
                    AnnotSel$CpG_Regions)
  mantelhaen.test(MH_table, exact=TRUE)
}
# apply to calculate5hmC enrichment for each state 
res_E027 <- lapply(states, test_enrichment_MH, E027_bed)
res_E119 <- lapply(states, test_enrichment_MH, E119_bed)

# generate and output results table of 5hmC enrichment with each genomic feature 
results_summary <- function(res_list){
  tab <- matrix(NA, ncol = 5, nrow = 15)
  colnames(tab) <- c("P-value", "Odds ratio (95% CI)", "OR", "LB", "UB")
  rownames(tab) <- c("1_TssA", "2_TssAFlnk", "3_TxFlnk", "4_Tx", "5_TxWk", "6_EnhG", "7_Enh", 
                     "8_ZNF/Rpts", "9_Het", "10_TssBiv", "11_BivFlnk", "12_EnhBiv", "13_ReprPC", 
                     "14_ReprPCWk", "15_Quies")
  # look over each state (list elements) and add results to table 
  for(i in 1:length(res_list)){
    tab[i,1] <- paste0(res_list[[i]]$p.value)
    tab[i,2] <- paste0(format(round(res_list[[i]]$estimate, digits = 2), nsmall = 2), " (", 
                       format(round(res_list[[i]]$conf.int[1], digits = 2), nsmall = 2), "-", 
                       format(round(res_list[[i]]$conf.int[2], digits = 2), nsmall = 2), ")")
    tab[i,3] <- res_list[[i]]$estimate
    tab[i,4] <- res_list[[i]]$conf.int[1]
    tab[i,5] <- res_list[[i]]$conf.int[2]
  }
  tab
}
# apply to each results set of interest 
res_sum_E027 <- results_summary(res_E027)
res_sum_E119 <- results_summary(res_E119)

# write to csvs 
write.csv(res_sum_E027, file = paste0(dir2, "E027_ChromHMM_5hmC_enrichment_results.csv"))
write.csv(res_sum_E119, file = paste0(dir2, "E119_ChromHMM_5hmC_enrichment_results.csv"))


