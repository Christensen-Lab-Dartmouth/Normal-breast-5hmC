#######################################################################################################################
# Pre-process and calculate expression percentiles for breast tissue samples in GTEx RNA-seq dataset 
# NOTE: Only transcripts mapping to CpGs on 450K array were used to calculate expression percentiles 
#######################################################################################################################

# load libraries 
library(hexbin)
library(readr)
library(vsn)
library(doParallel)
library(data.table)
library(ggplot2)

# set directories
dir_1 <- "~/Dropbox (Christensen Lab)/NDRI_Breast_5hmC_update/04.Gene_expression/GTeX/Files/"
dir_2 <- "~/Dropbox (Christensen Lab)/NDRI_Breast_5hmC_update/04.Gene_expression/GTeX/Figures/"

########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load data 
########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# expression from RNA-seq analyses (GTEx)
exp <- fread(paste0(dir_1, "GTEx_Analysis_2016-01-15_v7_RSEMv1.2.22_transcript_tpm.txt"), header=TRUE, stringsAsFactors = F)

# available covariate/sample data 
covs <- read_delim(paste0(dir_1, "GTEx_v7_Annotations_SampleAttributesDS.txt"), delim = "\t", col_names = TRUE)

# ENST IDs mapping to any 450K CG
ENST_450K <- readRDS(paste0(dir_1, "ENST_accessions_mapping_to_450K.rds"))

########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# subset to breast tissue 
########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# create separate data frame for transcript to gene mapping 
anno <- data.frame(exp$transcript_id, exp$gene_id)

# subset to breast tissues 
covs_breast <- covs[which(covs$SMTS=="Breast"),]
ind1 <- colnames(exp) %in% covs_breast$SAMPID
exp_breast <- data.frame(exp[, ..ind1])

# add transcript IDs as rownames 
rownames(exp_breast) <- exp$transcript_id

# get covariate data for the 16 samples that are missing from the expression data set
covs_missing <- covs_breast[!covs_breast$SAMPID %in% colnames(exp), ]
# NOTE: all seem to be SMAFRZW status==EXCLUDE, indicating they were not identified as 
# the samples best suited for use in analysis, with a specific focus on eQTL analysis 

# clean environment 
rm(covs, covs_missing)

# subset covs for subjects w/ SMAFRZW status==EXCLUDE 
covs_breast_2 <- covs_breast[covs_breast$SAMPID %in% colnames(exp_breast),]

########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# pre-process and visualize expression count data 
########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# log2 transform
exp_breast_log <- log2(exp_breast+1)

# evaluate heteroscadisicty and effect of transformation 
#### non-log transformed 
ppi=300
png(paste0(dir_2, "GTEx_breast_mean_sd_plot_untransformed.png"), width=6*ppi, height=5*ppi, res=ppi)
mds <- meanSdPlot(as.matrix(exp_breast))
mds$gg + ggtitle("GTEx breast un-transformed, all transcripts")
dev.off()
#### log-transformed 
png(paste0(dir_2, "GTEx_breast_mean_sd_plot_log_transformed.png"), width=6*ppi, height=5*ppi, res=ppi)
mds <- meanSdPlot(as.matrix(exp_breast_log))
mds$gg + ggtitle("GTEx breast log transformed, all transcripts")
dev.off()

# check distribution for random genes 

####gene 1
png(paste0(dir_2, "GTEx_breast_density_plots_random_genes.png"), width=7*ppi, height=12*ppi, res=ppi)
par(mfrow=c(3,2))
plot(density(as.numeric(exp[155, 3:dim(exp)[2]])), main="counts - random gene 1")
plot(density(as.numeric(exp_breast_log[155,])), main="log2 - counts - random gene 1")
#### gene 2
plot(density(as.numeric(exp[7656, 3:dim(exp)[2]])), main="counts - random gene 2")
plot(density(as.numeric(exp_breast_log[7656,])), main="log2 - random gene 2")
#### gene 3
plot(density(as.numeric(exp[2855, 3:dim(exp)[2]])), main="counts - random gene 3")
plot(density(as.numeric(exp_breast_log[2855,])), main="log2 - random gene 3")
dev.off()

# split rownames into stable ENST IDs 
rownames(exp_breast_log) <- as.character(sapply(as.character(rownames(exp_breast_log)), function(x) strsplit(x, "[.]")[[1]][1]))

# restrict to only transcripts mapping to CpGs on the 450K array
# and passing QC for the breast tissue
exp_breast_log_450K <- exp_breast_log[rownames(exp_breast_log) %in% ENST_450K,]

# check the meansdplot for these 
png(paste0(dir_2, "GTEx_breast_mean_sd_plot_log_transformed_450K_transcripts_only.png"), width=6*ppi, height=5*ppi, res=ppi)
mds <- meanSdPlot(as.matrix(exp_breast_log_450K))
mds$gg + ggtitle("GTEx breast 450K transcripts only")
dev.off()

# calculate median expression for each of these transcripts  
medians <- apply(exp_breast_log_450K, 1, median)

# check distribution of medians
png(paste0(dir_2, "GTEx_breast_median_expression_density.png"), width=5*ppi, height=5*ppi, res=ppi)
hist(medians, main = "GTEx breast tissues median ENST expression")
dev.off()

########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# calculate expression percentiles for transcripts 
########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# calculate percentile distribution   
perc.rank <- ecdf(medians)
percentile <- sapply(medians, perc.rank)

# order by percentile and add ranks to df 
percentile_2 <- percentile[order(percentile, decreasing = TRUE)]
percentile_df <- data.frame(percentile_2, 1:length(percentile_2))
rownames(percentile_df) <- NULL
percentile_df$transcript_id <- names(medians)

# rename columns 
colnames(percentile_df) <- c("percentile", "rank", "transcript")

# reorder columns 
percentile_df <- percentile_df[,c("transcript", "percentile", "rank")]

# write out the genes df 
saveRDS(percentile_df, paste0(dir_1, "GTEx_breast_expression_percentiles.rds"))


