#######################################################################################################################
# Test 5hmC enrichment among 15-state model ChromHMM chromatin states, specifically for chromatin states 
# non-overlappiung between HMECs & v-HMECs 
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

E119_bed <- LOLA::readBed(paste0(dir2, "E119_15_coreMarks_mnemonics.bed"))
E028_bed <- LOLA::readBed(paste0(dir2, "E028_15_coreMarks_mnemonics.bed"))

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
# determine non-overlapping HMEC and vHMEC regions for each chromatin state 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# drop XY coords 
#### E028
E028_non_auto_coords <- which(seqnames(E028_bed)=="chrY" | 
                                seqnames(E028_bed)=="chrX" | 
                                seqnames(E028_bed)=="chrM")
E028_bed <- E028_bed[-E028_non_auto_coords,]
#### E119
E119_non_auto_coords <- which(seqnames(E119_bed)=="chrY" | 
                                seqnames(E119_bed)=="chrX" | 
                                seqnames(E119_bed)=="chrM")
E119_bed <- E119_bed[-E119_non_auto_coords,]

# define chromHMM states 
states <- c("1_TssA", "2_TssAFlnk", "3_TxFlnk", "4_Tx", "5_TxWk", "6_EnhG", "7_Enh", 
            "8_ZNF/Rpts", "9_Het", "10_TssBiv", "11_BivFlnk", "12_EnhBiv", "13_ReprPC", 
            "14_ReprPCWk", "15_Quies")

# determine regions of chromatin states that are non-overlapping 
get_overlap <- function(state, samp1_coords, samp2_coords, type){
  #state <- states[1]
  #samp1_coords <- E119_bed
  #samp2_coords <- E028_bed
  samp1_coords_sub <- samp1_coords[which(names(samp1_coords)==state)]
  samp2_coords_sub <- samp2_coords[which(names(samp2_coords)==state)]
  # get samp1 ranges not found in samp2 (chromatin regions lost by vHMEC relative to HMECs)
  if(type=="lost") {
    o <- subsetByOverlaps(samp1_coords_sub, samp2_coords_sub, invert = TRUE)
  }
  # get samp2 ranges not found in samp1 (chromatin regions gained by vHMEC relative to HMECs)
  if(type=="gained") {
    o <- subsetByOverlaps(samp2_coords_sub, samp1_coords_sub, invert = TRUE)
  }  
  o
}
non_overlapping_lost <- lapply(states, get_overlap, E119_bed, E028_bed, "lost")
non_overlapping_gained <- lapply(states, get_overlap, E119_bed, E028_bed, "gained")
non_overlapping_joint <- mapply(function(x,y) c(x,y), non_overlapping_lost, non_overlapping_gained)

# calculate length of elements in each state 
lost_n <- unlist(lapply(non_overlapping_lost, length))
gained_n <- unlist(lapply(non_overlapping_gained, length))
joint_n <- unlist(lapply(non_overlapping_joint, length))

# add names to lists 
names(non_overlapping_lost) <- states
names(non_overlapping_gained) <- states
names(non_overlapping_joint) <- states

# calculate % of non-overlapping sites 
#### lost 
mapply(function(x,y) ((x)/y)*100, lost_n, joint_n)
#### gained
mapply(function(x,y) ((x)/y)*100, gained_n, joint_n)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define chromatin states and run enrichment for each one 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# define function to run enrichment test using MH test
test_enrichment_MH <- function(state_coords, state){
  overlaps <- findOverlaps(state_coords, Illumina_gr)
  indicies <- subjectHits(overlaps)
  temp <- rep(0, length(AnnotSel$top_5hmC))
  temp[indicies] <- 1
  MH_table <- table(factor(AnnotSel$top_5hmC, levels = c("1","0")), 
                    factor(temp, levels = c("1","0")), 
                    AnnotSel$CpG_Regions)
  mh <- mantelhaen.test(MH_table, exact=TRUE)
}
# apply to calculate 5hmC enrichment for each state 
res_lost <- lapply(non_overlapping_lost, test_enrichment_MH)
res_gained <- lapply(non_overlapping_gained, test_enrichment_MH)
res_joint <- lapply(non_overlapping_joint, test_enrichment_MH)




overlaps <- findOverlaps(non_overlapping_gained[[8]], Illumina_gr)
indicies <- subjectHits(overlaps)
temp <- rep(0, length(AnnotSel$top_5hmC))
temp[indicies] <- 1
MH_table <- table(factor(AnnotSel$top_5hmC, levels = c("1","0")), 
                  factor(temp, levels = c("1","0")), 
                  AnnotSel$CpG_Regions)
mh <- mantelhaen.test(MH_table, exact=TRUE)



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
res_lost_sum <- results_summary(res_lost)
res_gained_sum <- results_summary(res_gained)
res_joint_sum <- results_summary(res_joint)

# save results as csv files 
write.csv(res_lost_sum, file = paste0(dir2, "HMEC_vs_vHMEC_lost_states_ChromHMM_5hmC_enrichment_results.csv"))
write.csv(res_gained_sum, file = paste0(dir2, "HMEC_vs_vHMEC_gained_states_ChromHMM_5hmC_enrichment_results.csv"))
write.csv(res_joint_sum, file = paste0(dir2, "HMEC_vs_vHMEC_all_differential_states_ChromHMM_5hmC_enrichment_results.csv"))

