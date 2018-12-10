
library(hexbin)
library(readr)
library(vsn)
library(doParallel)
library(data.table)

# directories
dir_1 <- "~/Dropbox (Christensen Lab)/NDRI_Breast_5hmC_update/GTeX/Files/"

# expression from RNA-seq analyses 
exp <- fread(paste0(dir_1, "GTEx_Analysis_2016-01-15_v7_RSEMv1.2.22_transcript_tpm.txt"), header=TRUE, stringsAsFactors = F)

# expression from RNA-seq analyses 
#exp <- read_delim(paste0(dir_1, "GTEx_Analysis_2016-01-15_v7_RSEMv1.2.22_transcript_tpm.txt"), delim = "\t", col_names = TRUE)
# available covariate/sample data 
covs <- read_delim(paste0(dir_1, "GTEx_v7_Annotations_SampleAttributesDS.txt"), delim = "\t", col_names = TRUE)



# look at the expression data 
exp[1:10, 1:10]
exp2 <- as.data.frame(exp)

# create separate data frame for transcript to gene mapping 
anno <- data.frame(exp$transcript_id, exp$gene_id)

# look at those covs 
covs[1:10, 1:10]

# check tissue types included 
table(covs$SMTS)

# check sample overlap between expression file 
table(covs_breast$SAMPID %in% colnames(exp))

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






# log2 transform
exp_breast_log <- log2(exp_breast+1)

# evaluate heteroscadisicty and effect of transformation 
meanSdPlot(as.matrix(exp_breast))
meanSdPlot(as.matrix(exp_breast_log))

####gene 1
plot(density(as.numeric(exp[155, 3:dim(exp)])), main="counts - random gene 1")
plot(density(as.numeric(exp_breast_log[155,])), main="log2 - counts - random gene 1")
#### gene 2
plot(density(as.numeric(exp[7656, 3:dim(exp)[2]])), main="counts - random gene 2")
plot(density(as.numeric(exp_breast_log[7656,])), main="log2 - random gene 2")
#### gene 3
plot(density(as.numeric(exp[2855, 3:dim(exp)[2]])), main="counts - random gene 3")
plot(density(as.numeric(exp_breast_log[2855,])), main="log2 - random gene 3")





# find duplicate gene models (due to transcript variants)
ind2 <- which(duplicated(exp$gene_id))
# drop duplicates 
percentile_df_2 <- percentile_df[-ind2,]



combine_transcripts <- function(x){
  
  which(exp$gene_id == x)
  
  if()
  
}




# calculate median expression for each gene 
meds <- apply(exp_breast, 1, median)
meds_log <- apply(exp_breast_log, 1, median)

# check distribution of medians
par(mfrow=c(1,2))
hist(meds)
hist(meds_log)





# check quantile distribution 
quantile(meds_log, probs = c(0, 0.25, 0.5, 0.75, 1)) 

# calculate percentiles  
perc.rank <- ecdf(meds_log)
percentile <- sapply(meds_log, perc.rank)

# order by percentile and add ranks to df 
percentile_2 <- percentile[order(percentile, decreasing = TRUE)]
percentile_df <- data.frame(percentile_2, 1:length(percentile_2))
rownames(percentile_df) <- NULL
percentile_df$transcript_id <- exp$transcript_id
percentile_df$gene_id <- exp$gene_id

# rename columns 
colnames(percentile_df) <- c("percentile", "rank", "transcript", "gene")

# reorder columns 
percentile_df <- percentile_df[,c("gene", "transcript", "percentile", "rank")]

# look at it 
head(percentile_df)

# find duplicate gene models (due to transcript variants)
ind2 <- which(duplicated(percentile_df$gene))
# drop duplicates 
percentile_df_2 <- percentile_df[-ind2,]

# look at it 
head(percentile_df)

# write out the genes df 
saveRDS(percentile_df, paste0(dir_1, "GTEx_breast_expression_percentiles.rds"))



