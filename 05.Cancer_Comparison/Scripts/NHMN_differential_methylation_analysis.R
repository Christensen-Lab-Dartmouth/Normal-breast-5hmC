#######################################################################################################################
# Test enrichment of high 5hmC CpGs among ductal carcinoma in-situ (DCIS) vs adjacent-normal differentially 
# methylated loci in New Hampshire Mammography Network (NHMN) data set 
#######################################################################################################################
rm(list=ls())
library(RefFreeEWAS)
library(data.table)
library(limma)
library(ggplot2)
library(ggthemes)
library(RnBeads)
library(RnBeads.hg19)
library(dplyr)
library(doParallel)
setwd("/Users/Owen 1/Dropbox (Christensen Lab)/NDRI_Breast_5hmC_update/")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load Data & pre-process before analysis 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load processed and QC'd (by Lucas) TCGA breast betas 
#save.dir <- "/Users/Owen 1/Dropbox (Christensen Lab)/NDRI_Breast_5hmC_update/06.TCGA_Cancer_Comparison/Files/"
#rns <- load.rnb.set(path=paste0(save.dir, "NHMN_rnb_set_filtered.zip"))
# get betas 
#betas <- meth(rns)
#saveRDS(betas, file = "NHMN_rnb_betas_filtered.Rdata")
betas <- readRDS("05.Cancer_Comparison/Files/NHMN_rnb_betas_filtered.Rdata")

# index betas for high 5hmC CpGs 
#annot <- annotation(rns)
#saveRDS(annot, file = "NHMN_rnb_betas_filtered_annotation.Rdata")
annot <- readRDS("05.Cancer_Comparison/Files/NHMN_rnb_betas_filtered_annotation.Rdata")
rownames(betas) <- rownames(annot)
betas[1:10,1:10]

# load high 5hmC CpG list 
high_5hmC <- readRDS("02.Characterization_5hmC_levels/Files/high_5hmc_top1%_5hmC.rds")

# load covariate data 
covariates <- read.csv("05.Cancer_Comparison/Files/GSE66313/NHMN_samples.csv", header = T, sep = ",", stringsAsFactors = F)

# remove subjects from covariate file that were removed during filtering 
covs_2 <- covariates[covariates$Sample_ID %in% colnames(betas),]
all(covs_2$Sample_ID==colnames(betas))
covs_2$Sample_ID[7]
colnames(betas)[7]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# run epigenome-wide association study (EWAS) approach to test CpG loci for differential methylation status
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# add sample ID as rownames to covs
rownames(covs_2) <- covs_2$Sample_ID

# Create a design matrix for covariates of interest
table(covs_2$tissue.ch1)
covs_2$tissue.ch1 <- as.factor(covs_2$tissue.ch1)
levels(covs_2$tissue.ch1)
table(covs_2$tissue.ch1)
XX <- model.matrix(~tissue.ch1+subject.age.ch1, data = covs_2)

# Convert beta-values to M-values for gaussian consideration
betas_2 <- ifelse(betas>=1,1-1E-6,ifelse(betas<=0,1E-6,betas))
Betas_NHMN <- log(betas_2)-log(1-betas_2)
all(colnames(Betas_NHMN)==rownames(XX))

# Apply limma
lf_Null <- eBayes(lmFit(Betas_NHMN, XX))
#save(lf_Null, file="06.TCGA_Cancer_Comparison/Files/NHMN_breast_limma_models.RData")
#load("05.Cancer_Comparison/Files/NHMN_breast_limma_models.RData")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# visualize results 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
str(lf_Null)

# index results for high 5hmC CpGs available in the TCGA data set 
overlap_CpGs <- high_5hmC$id[which(high_5hmC$id%in%rownames(lf_Null$p.value))]

# explore epigenome-wide results
x<-lf_Null$p.value
hist(lf_Null$p.value[ , 'tissue.ch1 DCIS']) # P-value distribution 
table(lf_Null$p.value[ , 'tissue.ch1 DCIS'] < (0.05/length(overlap_CpGs))) # how many above Bonf. threshold
table(is.na(lf_Null$p.value[ , 'tissue.ch1 DCIS'] < (0.05/length(overlap_CpGs)))) # + w/o NAs

# explore high 5hmC CpG specific results
table(lf_Null$p.value[overlap_CpGs, 'tissue.ch1 DCIS'] < (0.05/length(overlap_CpGs))) # how many above Bonf. threshold
coefs <- lf_Null$coefficients[overlap_CpGs, ]
ind1 <- which(lf_Null$p.value[overlap_CpGs, 'tissue.ch1 DCIS'] < (0.05/length(overlap_CpGs)))
coefs_sig <- coefs[ind1,]
sum(coefs_sig[, 'tissue.ch1 DCIS'] < 0)
sum(coefs_sig[, 'tissue.ch1 DCIS'] > 0)

# restrict annotation to high 5hmC CpGs  
annot_overlap <- annot[overlap_CpGs, ] 
annot_results <- annot[rownames(lf_Null), ] 

# how many high 5hmC CpGs are havce P<0.05 
sum(lf_Null$p.value[overlap_CpGs, 'tissue.ch1 DCIS'] < (0.05)) 

# check dist of CpG island relation among high 5hmC CpG and all results 
table(annot_overlap$`CGI Relation`)
table(is.na(annot_overlap$`CGI Relation`))
table(annot_results$`CGI Relation`)
table(is.na(annot_results$`CGI Relation`))
# Collapse "North" and "South" nomenclature for CpG islands
annot_overlap$`CGI Relation` <- gsub("^[N|S].....", "", annot_overlap$`CGI Relation`)
annot_results$`CGI Relation` <- gsub("^[N|S].....", "", annot_results$`CGI Relation`)
# get the number of CpGs in each context from high 5hmC set 
tab1 <- table(annot_overlap$`CGI Relation`)
# check no of CpGs in overall set in each context 
table(annot_results$`CGI Relation`)
# create subset of annot_results that doesn't contain high 5hmC CpGs
annot_results_sub <- annot_results[-match(rownames(annot_overlap),rownames(annot_results)),]
table(is.na(match(rownames(annot_overlap),rownames(annot_results))))

# Select the results that meet the threshold for significance
Limma_Results <- as.data.frame(cbind((lf_Null$coef[overlap_CpGs, 'tissue.ch1 DCIS']), lf_Null$p.value[overlap_CpGs, 'tissue.ch1 DCIS']))
results = mutate(Limma_Results, sig=ifelse(Limma_Results$V2<0.05/(length(overlap_CpGs)), "P-value < 0.05", "Not Sig"))
table(results[,3])
colnames(results) = c("Coefficient", "P_value", "limma_model") 

# extract coef and P-values from results for 5hmC set 
coef_hmc_sites <- Limma_Results$V1
log_pval_hmc_sites <- -log10(Limma_Results$V2)

# randomly sample 1000 CpG sets with the same CpG island proportion as high 5hmC set 
set.seed(174)
mat1 <- matrix(NA, nrow = length(overlap_CpGs), ncol = 1000)
for(i in 1:1000){
  Island = sample(rownames(annot_results_sub[annot_results_sub$`CGI Relation`=="Island", ]), tab1[["Island"]])
  OpenSea = sample(rownames(annot_results_sub[annot_results_sub$`CGI Relation`=="Open Sea", ]), tab1[["Open Sea"]])
  Shelf = sample(rownames(annot_results_sub[annot_results_sub$`CGI Relation`=="Shelf", ]), tab1[["Shelf"]])
  Shore = sample(rownames(annot_results_sub[annot_results_sub$`CGI Relation`=="Shore", ]), tab1[["Shore"]])
  NHMN_random_CpGs = c(Island, OpenSea, Shelf, Shore)
  plist <- lf_Null$p.value[NHMN_random_CpGs, 'tissue.ch1 DCIS']
  mat1[,i] <- plist[order(plist)]
}
#saveRDS(mat1, file = "NHMN_random_CpG_Pvalue_matrix_seed_465.rds")
#mat1 <- readRDS("05.Cancer_Comparison/Files/NHMN_random_CpG_Pvalue_matrix_seed_465.rds")

# calculate the average p-value across each row 
log_pval_means_random_sites <- -log10(apply(mat1, 1, mean))
log_pval_medians_random_sites <- -log10(apply(mat1, 1, median))

# visualize to check
par(mfrow=c(1,2))
plot(log_pval_means_random_sites) 
plot(log_pval_medians_random_sites) 
plot(log_pval_hmc_sites) 
dev.off()

# run KS-test 
ks.test(log_pval_hmc_sites, log_pval_means_random_sites)
ks.test(log_pval_hmc_sites, log_pval_medians_random_sites)

# Plot the cdf of the age-related and random P-values
cdf_hmc_pval = ecdf(log_pval_hmc_sites)
cdf_random_pval = ecdf(log_pval_means_random_sites)
#quants <- quantile(-log10(pval_hmc_sites), probs = seq(0, 1, 1/3572))
#tab2 <- as.data.frame(cbind(quants, -log10(pval_hmc_sites)))

# Figure for locus-specific TCGA vs normal
png("05.Cancer_Comparison/Figures/NHMN_Volcano.png", height=8*250, width=8*300, pointsize = 16, res=300)
p = ggplot(results, aes(Coefficient, -log10(P_value))) +
  geom_point(aes(col=limma_model))  + xlim(-2.5, 2.5) +
  xlab("coefficient") + ylab("-log10(P-value)") +
  scale_color_manual(values=c("black", "red")) + 
  #geom_hline(yintercept = -log10(.05/dim(lf_Null$p.value)[1]), color = "red", linetype = "dashed", size = 1.3) + 
  geom_hline(yintercept = -log10(.05/length(overlap_CpGs)), color = "red", linetype = "dashed", size = 1.3) + 
  labs(color="limma model \n significance") + 
  theme(legend.key = element_blank()) + 
  theme(legend.key.size = unit(1, "cm"), 
        #panel.grid.major = element_line(colour = "#d3d3d3"), 
        panel.border = element_rect(fill = NA, colour = "black", size = 0.6, linetype = "solid"), 
        panel.background = element_blank(),
        axis.text.x=element_text(colour="black", size = 20, hjust = 1), 
        axis.text.y=element_text(colour="black", size = 20), 
        axis.title.x=element_text(colour="black", size = 20), 
        axis.title.y=element_text(colour="black", size = 20), 
        legend.key = element_blank(), 
        legend.text = element_text(colour="black", size = 20), 
        legend.title = element_text(colour="black", size = 20))
print(p)
dev.off()

# plot P-value cumulative distribution plots
png("05.Cancer_Comparison/Figures/NHMN_Pvalue_ecdf.png", height=7*300, width=7*275, pointsize = 12, res=300)
plot(cdf_hmc_pval, main="DCIS-Adjacent Normal Differential \n DNA methylation (NHMN)", 
     xlab="-log10(P-value)", ylab="Cumulative proportion - P-value", col="red", 
     cex.axis = 1.45, cex.lab = 1.45, las = 1)
points(log_pval_hmc_sites[order(log_pval_hmc_sites)], col="red", cdf_hmc_pval(log_pval_hmc_sites)[order(cdf_hmc_pval(log_pval_hmc_sites))], cex = 0.3)
plot(cdf_random_pval, add=TRUE, col="black") 
points(log_pval_means_random_sites[order(log_pval_means_random_sites)], col="black", cdf_random_pval(log_pval_means_random_sites)[order(cdf_random_pval(log_pval_means_random_sites))], cex = 0.3)
legend(2, 0.5, c("High 5hmC CpGs", "Random CpGs"), lwd=3, bty="n", col=c("red","black"), cex = 1.3)
text(3, 0.3, "Kolmogorov-Smirnov test P = 0.017", cex = 1.3) 
dev.off()

