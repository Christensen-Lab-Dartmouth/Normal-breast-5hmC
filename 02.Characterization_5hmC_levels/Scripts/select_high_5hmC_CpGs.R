#######################################################################################################################
# select 'high5hmC CpGs' to be used for downstream analysis & generate summary statistics for these CpGs 
#######################################################################################################################
rm(list = ls())
setwd("/Users/Owen 1/Dropbox (Christensen Lab)/NDRI_Breast_5hmC_update/")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load Data 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# high5hmC data 
load("02.Characterization_5hmC_levels/Files/MethOxy_FunNorm_good_probes.Rdata")
MethOxy <- MethOxy_2
rm(MethOxy_2)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# select high 5hmC CpGs
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# calculate median 5hmC values for CpGs across subjects 
median1 <- apply(MethOxy[, , 3], 1, median, na.rm=TRUE)
mean1 <- apply(MethOxy[, , 3], 1, mean, na.rm=TRUE)
sd1 <- apply(MethOxy[, , 3], 1, sd, na.rm=TRUE)

# group into df 
dat1 <- data.frame(I(names(median1)), median1, mean1, sd1)
rownames(dat1) <- NULL
colnames(dat1) = c("ID", "Median_5hmC", "Mean_5hmC", "SD")
head(dat1)

#Order CpGs the based on median beta value 
dat2 <- dat1[order(dat1$Median_5hmC, decreasing = FALSE),]

# select those CpGs amon the top 1% median across all subjects 
top_1 <- round(dim(dat2)[1]/100, digits = 0)
site_start <- dim(dat2)[1]-top_1+1
site_end <- dim(dat2)[1]
dat3 <- dat2[site_start:site_end,]

# save ordered CpG IDs w/ their 5hmC beta values
saveRDS(dat3, file = "02.Characterization_5hmC_levels/Files/high_5hmc_top1%_5hmC.rds")

# how many CpGs w/ 5hmC signals above thresholds of interest 
table(dat2$Median_5hmC>0.05)
table(dat2$Median_5hmC>0.10)
table(dat2$Median_5hmC>0.15)
table(dat2$Median_5hmC>0.20)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# identify most/least variable CpGs across subjects and calculate proportions of CpGs w/ 'recurrent' 5hmC 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Order based on sd values
sd_2 <- sd1[order(sd1, decreasing = FALSE)]

# get top 1% sd CpGs (most variable)
top_1_sd <- round(length(sd_2)/100, digits = 0)
site_start <- length(sd_2)-top_1_sd+1
site_end <- length(sd_2)
sd_3 <- sd_2[site_start:site_end]

# get bottom 1% sd CpGs (least variable)
site_start <- 1
site_end <- top_1_sd
sd_4 <- sd_2[site_start:site_end]

# overlap between most variable (highest sd) CGs and high 5hmC CGs
table(dat3$ID %in% names(sd_3))
500/3376 # 14.8% 
# overlap between least variable (highest sd) CGs and high 5hmC CGs
table(dat3$ID %in% names(sd_4))
# 0%

# how high 5hmC CGs were recurrent (consistently elevated 5hmC across subjects)
# ie (many CGs had at least half the samples displaying the minimal median value for high 5hmC CGs - 0.1411354)
Meth_high_5hmC <- MethOxy[dat3$ID, , 3]
recurrent_q25 <- c()
recurrent_q50 <- c()
recurrent_q75 <- c()
quantiles_medians <- quantile(dat3$Median_5hmC)
quantiles_medians
for(i in 1:nrow(Meth_high_5hmC)){
  recurrent_q25[i] <- sum(Meth_high_5hmC[i,] > quantiles_medians[[2]])
  recurrent_q50[i] <- sum(Meth_high_5hmC[i,] > quantiles_medians[[3]])
  recurrent_q75[i] <- sum(Meth_high_5hmC[i,] > quantiles_medians[[4]])
}
table(recurrent_q25>=9) # 3281/3876 = 84.6%
table(recurrent_q50>=9) # 2229/3876 = 57.5%
table(recurrent_q75>=9) # 1137/3876 = 29.3%

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# generate summary table of all high 5hmC CGs (Supp. Data 1 for manuscript)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# import annotation data 
load("02.Characterization_5hmC_levels/Files/annotation_good_probes.Rdata")
annot <- AnnotSel

# subset for high 5hmC CGs 
annot_2 <- annot[match(dat3$ID, annot$TargetID),]

# group annotation data and summary statistics into one data frame 
dat4 <- data.frame(I(dat3$ID), I(annot_2$chr), I(annot_2$MAPINFO), 
                   I(annot_2$UCSC_RefGene_Name), I(annot_2$UCSC_RefGene_Group), 
                   round(dat3$Median_5hmC, digits = 3), round(dat3$Mean_5hmC, digits = 3), round(dat3$SD, digits = 3))
str(dat4)

# remove "chr" prefix from chromosome no.
chromo <- as.character(dat4$annot_2.chr)
chromo2 <- sapply(chromo, function (x) strsplit(x, "r")[[1]][2])
dat4$annot_2.chr <- chromo2

# update colnames 
colnames(dat4) <- c("Illumina CpG ID", "Chromosome (hg19)", 
                    "Genomic Coordinate (hg19)", "UCSC RefGene Name", "UCSC RefGene Group", 
                    "Median 5hmC", "Mean 5hmC", "Standard deviation 5hmC")

# write out to csv 
write.csv(dat4, file = "02.Characterization_5hmC_levels/Files/high_5hmC_CpG_summary.csv")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# generate scatter plot of all high 5hmC CGs plotted against 5mC value in each subject 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# have a look at the 14733 'high 5hmC' CpGs plotted against their 5mC levels 
MethOxy_2 <- MethOxy[dat3$ID, , ]
dpi=300
png("02.Characterization_5hmC_levels/Figures/high_5hmC_vs_5mC.png", width=5*dpi, height=5*dpi, res=dpi)
smoothScatter(MethOxy_2[, , 3]~MethOxy_2[, , 2], colramp = colorRampPalette(c("white", "blue", "orange"), space = "Lab"), nrpoints = 0, postPlotHook= NULL,
              xlab = "5mC beta-value", ylab = "5hmC beta-value", bty = "n", las = 1, cex.lab = 1.3, cex.axis = 1.15, ylim = c(0, 0.6), xaxs = "i", xlim = c(0,1), yaxs = "i")
dev.off()

