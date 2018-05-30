#######################################################################################################################
# Determine relation of high 5hmC CpGs to canonical transcription start sites (TSSs)
#######################################################################################################################
library(GenomicRanges) ; library(rtracklayer) ; library(genomation)
rm(list = ls())
setwd("/Users/Owen 1/Dropbox (Christensen Lab)/NDRI_Breast_5hmC_update/")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load Data 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# annotation file
load("02.Characterization_5hmC_levels/Files/annotation_good_probes.Rdata")
annot <- AnnotSel
# load methylation data
load("02.Characterization_5hmC_levels/Files/MethOxy_FunNorm_good_probes.Rdata")
# high 5hmC probes 
high_5hmC <- readRDS("02.Characterization_5hmC_levels/Files/high_5hmc_top1%_5hmC.rds")

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
annot_sub$top_5hmC <- NA
annot_sub$top_5hmC[annot_sub$TargetID %in% high_5hmC$id] <- 1
annot_sub$top_5hmC[!annot_sub$TargetID %in% high_5hmC$id] <- 0
table(annot_sub$top_5hmC)

annot_5hmC <- annot_sub[annot_sub$top_5hmC==1,]
annot_5hmC_gr <- makeGRangesFromDataFrame(annot_5hmC, keep.extra.columns=TRUE, ignore.strand=TRUE, seqnames.field = "chr", start.field="hg19_start", end.field="hg19_end")
annot_5hmC_gr

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# calculate distance and orientation of 5hmC relative to nearest canonical TSS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# read in BED file to get coordinates for unique promoters,TSSs, introns, exons, etc. 
transcFeat = readTranscriptFeatures("02.Characterization_5hmC_levels/Files/UCSC_hg19_refGene.bed", up.flank = 2000, down.flank = 2000, unique.prom = TRUE)
names(transcFeat)

# get indicies of high 5hmC CpGs in annotation file that precede transcription start sites (TSSs)
tss_precede <- precede(annot_5hmC_gr, transcFeat$TSSes)
# do the same for those that follow TSSs
tss_follow <- follow(annot_5hmC_gr, transcFeat$TSSes)

# identify the nearest TSS for each of the high 5hmC CpGs 
nearest_tss <- nearest(annot_5hmC_gr, transcFeat$TSSes)
# index for these TSSs only
transcFeat_tss_sub <- transcFeat$TSSes[nearest_tss]

# calculate distance to nearest TSS for each high 5hmC CpG
dist_to_TSS <- distance(annot_5hmC_gr, transcFeat_tss_sub)
# divide by 1000 to convert to kilobases
dist_to_TSS <- dist_to_TSS/1000

# identify which of CpGs are up- or down-stream of the nearest 
upstream_indicies <- which(nearest_tss==tss_precede)
downstream_indicies <- which(nearest_tss!=tss_precede)
dist_to_TSS_upstream <- dist_to_TSS[upstream_indicies]
dist_to_TSS_downstream <- dist_to_TSS[downstream_indicies]

# visualize distributions of each
hist(dist_to_TSS_upstream, 20)
hist(dist_to_TSS_downstream, 20)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# calculate and visualize proportions of high 5hmC in various distance ranges of nearest TSS 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# write function that will calculate proportions of high 5hmC CpGs within specified 'bins' for genomic range
# (assuming distnace to nearest TSS, up and downstream, have been determined)
calc_prop <- function(dist_up, dist_down, bins){
  p_up <- c()
  p_down <- c()
  for(i in 1:(length(bins)-1)){
    t1 <- table(dist_up>bins[i] & dist_up<=bins[i+1])
    if(length(t1)>1){
      p_up[i] <- t1[[2]]/length(annot_5hmC_gr)*100
    } else {
      p_up[i] <- 0
    }
  }
  for(i in 1:(length(bins)-1)){
    t2 <- table(dist_down>bins[i] & dist_down<=bins[i+1])
    if(length(t2)>1){
      p_down[i] <- t2[[2]]/length(annot_5hmC_gr)*100
    } else {
      p_down[i] <- 0
    }
  }
  t_up_end <- table(dist_up>bins[length(bins)])
  if(length(t_up_end)>1){p_up_end <- t_up_end[[2]]/length(annot_5hmC_gr)*100}
  else{p_up_end <- 0}
  
  t_down_end <- table(dist_down>bins[length(bins)])
  if(length(t_down_end)>1){p_down_end <- t_down_end[[2]]/length(annot_5hmC_gr)*100}
  else{p_down_end <- 0}
  
  p_combo <- c(p_up_end, rev(p_up), p_down, p_down_end)
  p_combo
}

# caclulate & plot proportions of high 5hmC CpGs within 0-5, 5-50, 50-500, or >500kb of the nearest TSS 
props <- calc_prop(dist_to_TSS_upstream, dist_to_TSS_downstream, c(0, 5, 50, 500))
names(props) <- c("< -500", "-500 to 50", "-50 to -5", "-5 to 0", "0 to 5", "5 to 50", "50 to 500", "> 500")
ppi=300
png("02.Characterization_5hmC_levels/Figures/high_5hmC_CpGs_TSS_relation_-5-+5kb_barplot.png", width=6*ppi, height=5*ppi, res=ppi)
mp <- barplot(props, col = "indianred", las = 1, ylim = c(0, 30),
        ylab = "Proportion of high 5hmC CpGs (%)", xlab = "", 
        cex.lab = 1.3, cex.axis = 1.15, xaxt = "n")
labels <- c("< -500", "-500 to 50", "-50 to -5", "-5 to 0", "0 to 5", "5 to 50", "50 to 500", "> 500")
text(mp, par("usr")[3], labels = labels, srt = 45, adj = c(1.1,1.1), xpd = T, cex=1.3)
mtext("Distance to TSS (kb)", side=1, line=4, cex = 1.15)
dev.off()

# caclulate & plot proportions of high 5hmC CpGs +/-5kb of nearest TSS, giving proportions for every 1kb window 
props <- calc_prop(dist_to_TSS_upstream, dist_to_TSS_downstream, seq(0, 5, 1))
# drop 1st & last element as these relate to the proportion of CpGs outside of this region 
props <- props[-c(1, length(props))]
names(props) <- c("-5 to -4", "-4 to -3", "-3 to -2", "-2 to -1", "-1 to 0", "0 to 1", "1 to 2", "2 to 3", "3 to 4", "4 to 5")
png("02.Characterization_5hmC_levels/Figures/high_5hmC_CpGs_TSS_relation_-5-+5kb_lineplot_10bp.png", width=6*ppi, height=5*ppi, res=ppi)
mp <- barplot(props, col = "indianred", las = 1, ylim = c(0, 16),
        ylab = "Proportion of high 5hmC CpGs (%)", xlab = "", 
        cex.lab = 1.3, cex.axis = 1.15, xaxt = "n")
text(mp, par("usr")[3], labels = names(props), srt = 45, adj = c(1.1,1.1), xpd = T, cex=1.3)
mtext("Distance to TSS (kb)", side=1, line=4, cex = 1.15)
dev.off()

