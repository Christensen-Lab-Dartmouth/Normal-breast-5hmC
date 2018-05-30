#######################################################################################################################
# Test enrichment of high 5hmC CpGs among tumor vs normal differentially methylated loci in TCGA data set 
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
# load Data & ore-process before analysis 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load processed and QC'd (by Lucas) TCGA breast betas 
#save.dir <- "/Users/Owen/Dropbox (Christensen Lab)/BRCA TCGA/BRCAfiltered"
#rns <- load.rnb.set(path=paste0(save.dir, ".zip"))
# get betas 
#betas <- meth(rns)
#saveRDS(betas, file = "Files/rnb_betas_filtered.Rdata")
betas <- readRDS("05.Cancer_Comparison/Files/rnb_betas_filtered.Rdata")

# index betas for high 5hmC CpGs 
#annot <- annotation(rns)
#saveRDS(annot, file = "Files/rnb_betas_filtered_annotation.Rdata")
annot <- readRDS("05.Cancer_Comparison/Files/rnb_betas_filtered_annotation.Rdata")
rownames(betas) <- rownames(annot)
betas[1:10,1:10]

# load high 5hmC CpG list 
high_5hmC <- readRDS("02.Characterization_5hmC_levels/Files/high_5hmc_top1%_5hmC.rds")

# load covariate data 
covariates <- read.csv("05.Cancer_Comparison/Files/BRCA_samples filtered3.csv", header = T, sep = ",", stringsAsFactors = F)

# check for missing subjects between betas and covariate data 
table(covariates$Sample_ID %in% colnames(betas))
# seems to be 4 missing, likely due to removal during array QA/QC

# organize rownames 
rownames(covariates) <- covariates$X
covariates$X <- NULL
# index covariate file 
table(covariates$casecontrol, is.na(covariates$age))
# drop subjects w/o age data 
covariates <- covariates[!is.na(covariates$age), ]
# drop subvjects with mets 
Mets_indices <- which(covariates$ajcc_metastasis_pathologic_pm=="M1")
covariates <- covariates[-Mets_indices, ]
table(covariates$casecontrol)
dim(betas) ; dim(covariates)
table(is.na(match(covariates$Sample_ID, colnames(betas))))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# run epigenome-wide association study (EWAS) approach to test CpG loci for differential methylation status
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# index beta for the subjects reminaing in covariate data 
betas_2 <- betas[, na.omit(match(rownames(covariates), colnames(betas)))]
covariates_2 <- covariates[-which(is.na((match(rownames(covariates), colnames(betas))))),]
dim(betas_2) ; dim(covariates_2)
# Are the samples the same in both the TCGA betas and the covariates
all(colnames(betas_2)==rownames(covariates_2))

# Create a design matrix for covariates of interest
table(covariates_2$casecontrol)
covariates_2$casecontrol[covariates_2$casecontrol==1] <- "Primary Tumor"
covariates_2$casecontrol[covariates_2$casecontrol==11] <- "Solid Tissue Normal"
covariates_2$casecontrol <- as.factor(covariates_2$casecontrol)
levels(covariates_2$casecontrol)
table(covariates_2$casecontrol)
covariates_2$casecontrol <- factor(covariates_2$casecontrol, levels = c("Solid Tissue Normal", "Primary Tumor"))
levels(covariates_2$casecontrol)
table(covariates_2$casecontrol)
XX <- model.matrix(~casecontrol+age, data = covariates_2)

# Convert beta-values to M-values for gaussian consideration
betas_3 <- ifelse(betas_2>=1,1-1E-6,ifelse(betas_2<=0,1E-6,betas_2))
Betas_TCGAM <- log(betas_3)-log(1-betas_3)
all(colnames(Betas_TCGAM)==rownames(XX))

# Apply limma
lf_Null <-  eBayes(lmFit(Betas_TCGAM, XX))
#save(lf_Null, file="05.Cancer_Comparison/Files/TCGA_breast_limma_models.RData")
load("05.Cancer_Comparison/Files/TCGA_breast_limma_models.RData")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# visualize results 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# explore epigenome-wide results 
hist(lf_Null$p.value[ , 'casecontrolPrimary Tumor']) # P-value distribution 
table(lf_Null$p.value[ , 'casecontrolPrimary Tumor'] < (0.05/dim(lf_Null$p.value)[1])) # how many above Bonf. threshold
table(is.na(lf_Null$p.value[ , 'casecontrolPrimary Tumor'] < (0.05/dim(lf_Null$p.value)[1]))) # + w/o NAs

# index results for high 5hmC CpGs available in the TCGA data set 
overlap_CpGs <- high_5hmC$id[which(high_5hmC$id%in%rownames(lf_Null$p.value))]

# restrict annotation to high 5hmC CpGs  
annot_overlap <- annot[overlap_CpGs, ] 
annot_results <- annot[rownames(lf_Null), ] 

# Select the results that meet the threshold for significance
Limma_Results <- as.data.frame(cbind((lf_Null$coef[overlap_CpGs, 'casecontrolPrimary Tumor']), lf_Null$p.value[overlap_CpGs, 'casecontrolPrimary Tumor']))
results = mutate(Limma_Results, sig=ifelse(Limma_Results$V2<0.05/(dim(Limma_Results)[1]), "P-value < 0.05 (1,712 CpGs)", "Not Sig (1,860 CpGs)"))
table(results[,3])
colnames(results) = c("Coefficient", "P_value", "limma_model") 

# extract coef and P-values from results for 5hmC set 
coef_hmc_sites <- Limma_Results$V1
log_pval_hmc_sites <- -log10(Limma_Results$V2)

# how many high 5hmC CpGs are havce P<0.05 
sum(lf_Null$p.value[overlap_CpGs, 'casecontrolPrimary Tumor'] < (0.05/dim(results)[1])) 
high_5hmC_low_P <- names(which(lf_Null$p.value[overlap_CpGs, 'casecontrolPrimary Tumor'] < (0.05/dim(results)[1])))
#saveRDS(high_5hmC_low_P, file = "05.Cancer_Comparison/Files/high_5hmC_CGs_below_bonf_05_in_TCGA.RData")

# check how many sig CGs have +ve or -ve coefficients 
coefs <- lf_Null$coefficients[overlap_CpGs, ]
ind1 <- which(lf_Null$p.value[overlap_CpGs, 'casecontrolPrimary Tumor'] < (0.05/dim(results)[1]))
coefs_sig <- coefs[ind1,]
sum(coefs_sig[, 'casecontrolPrimary Tumor'] < 0)
sum(coefs_sig[, 'casecontrolPrimary Tumor'] > 0)
plot(coefs[coefs[ind1, 'casecontrolPrimary Tumor'] < 0, 2])





# determine coefficients for CpGs located within breast specific TSGs of interest 
TSG_annot <- readRDS("02.Characterization_5hmC_levels/Files/TSG_annot.rds")

coefs[match(TSG_annot$MBNL1$TargetID, rownames(coefs)),]
coefs[match(TSG_annot$ARID1B$TargetID, rownames(coefs)),]
coefs[match(TSG_annot$DNMT3A$TargetID, rownames(coefs)),]
coefs[match(TSG_annot$FOXO3$TargetID, rownames(coefs)),]

plot(coefs[match(TSG_annot$MBNL1$TargetID, rownames(coefs)),][,2]) ; abline(h=0, col = "red")
plot(coefs[match(TSG_annot$ARID1B$TargetID, rownames(coefs)),]) ; abline(h=0, col = "red")
plot(coefs[match(TSG_annot$DNMT3A$TargetID, rownames(coefs)),]) ; abline(h=0, col = "red")
plot(coefs[match(TSG_annot$FOXO3$TargetID, rownames(coefs)),]) ; abline(h=0, col = "red")

coefs_sig[match(TSG_annot$MBNL1$TargetID, rownames(coefs_sig)),]
coefs_sig[match(TSG_annot$ARID1B$TargetID, rownames(coefs_sig)),]
coefs_sig[match(TSG_annot$DNMT3A$TargetID, rownames(coefs_sig)),]
coefs_sig[match(TSG_annot$FOXO3$TargetID, rownames(coefs_sig)),]

plot(coefs_sig[match(TSG_annot$MBNL1$TargetID, rownames(coefs_sig)),][,2]) ; abline(h=0, col = "red")
plot(coefs_sig[match(TSG_annot$ARID1B$TargetID, rownames(coefs_sig)),]) ; abline(h=0, col = "red")
plot(coefs_sig[match(TSG_annot$DNMT3A$TargetID, rownames(coefs_sig)),]) ; abline(h=0, col = "red")
plot(coefs_sig[match(TSG_annot$FOXO3$TargetID, rownames(coefs_sig)),]) ; abline(h=0, col = "red")






# do the same thing with the set of high 5hmC CGs that overlaps with the solo-WCGWs
annot_5hmC_PMD_soloCGs_IDs <- readRDS("Misc/Solo-WCGW-comparison/Data/high_5hmC_CGs_breast_solo_WCGWs_overlap.rds")
# check among all tests 
coefs_soloCGs_5hmC <- coefs[na.omit(match(annot_5hmC_PMD_soloCGs_IDs, rownames(coefs))),]
dim(coefs_soloCGs_5hmC)
plot(coefs[match(annot_5hmC_PMD_soloCGs_IDs, rownames(coefs)),][,2]) ; abline(h=0, col = "red")
table(coefs_soloCGs_5hmC[,2]<0)[[2]][1]/dim(coefs_soloCGs_5hmC)[1]*100

# check among those tests that were significant after bonferonni correction 
sig_coefs_soloCGs_5hmC <- coefs_sig[na.omit(match(annot_5hmC_PMD_soloCGs_IDs, rownames(coefs_sig))),]
dim(sig_coefs_soloCGs_5hmC)
plot(coefs_sig[match(annot_5hmC_PMD_soloCGs_IDs, rownames(coefs_sig)),][,2]) ; abline(h=0, col = "red")
table(sig_coefs_soloCGs_5hmC[,2]<0)
# all of the solo_WCGWs that overlap with high 5hmC CpGs in breast have negative coefficients 





# what about all of the soloWCGWs (regardless of overlap with high 5hmC)
hm450k_PMD_soloCGs <- read.table("Misc/Solo-WCGW-comparison/Data/hm450.comPMD.probes.tsv", header = TRUE, stringsAsFactors = FALSE, sep = "\t")$x

coefs_full <- lf_Null$coefficients
ind2 <- which(lf_Null$p.value[, 'casecontrolPrimary Tumor'] < (0.05/dim(lf_Null)[1]))
coefs_full_sig <- coefs_full[ind2,]

coefs_soloCGs <- coefs_full[na.omit(match(hm450k_PMD_soloCGs, rownames(coefs_full))),]
plot(coefs_soloCGs[,2]) ; abline(h=0, col = "red")
table(coefs_soloCGs[,2]<0)[[2]][1]/dim(coefs_soloCGs)[1]*100

sig_coefs_soloCGs <- coefs_full_sig[na.omit(match(hm450k_PMD_soloCGs, rownames(coefs_full_sig))),]
plot(sig_coefs_soloCGs[,2]) ; abline(h=0, col = "red")
table(sig_coefs_soloCGs[,2]<0)[[2]][1]/dim(sig_coefs_soloCGs)[1]*100







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


# randomly sample 1000 CpG sets with the same CpG island proportion as high 5hmC set 
set.seed(100)
mat1 <- matrix(NA, nrow = length(overlap_CpGs), ncol = 1000)
for(i in 1:1000){
  Island = sample(rownames(annot_results_sub[annot_results_sub$`CGI Relation`=="Island", ]), tab1[["Island"]])
  OpenSea = sample(rownames(annot_results_sub[annot_results_sub$`CGI Relation`=="Open Sea", ]), tab1[["Open Sea"]])
  Shelf = sample(rownames(annot_results_sub[annot_results_sub$`CGI Relation`=="Shelf", ]), tab1[["Shelf"]])
  Shore = sample(rownames(annot_results_sub[annot_results_sub$`CGI Relation`=="Shore", ]), tab1[["Shore"]])
  TCGA_random_CpGs = c(Island, OpenSea, Shelf, Shore)
  plist <- lf_Null$p.value[TCGA_random_CpGs, 'casecontrolPrimary Tumor']
  mat1[,i] <- plist[order(plist)]
}

# calculate the average p-value across each row 
pval_means_random_sites <- apply(mat1, 1, mean)
log_pval_means_random_sites <- -log10(apply(mat1, 1, mean))
log_pval_medians_random_sites <- -log10(apply(mat1, 1, median))

# visualize to check
par(mfrow=c(1,2))
plot(log_pval_means_random_sites) 
plot(log_pval_medians_random_sites) 
plot(log_pval_hmc_sites) 

# run KS-test 
ks.test(log_pval_hmc_sites, log_pval_means_random_sites)
#data:  log_pval_hmc_sites and log_pval_means_random_sites
#D = 0.033035, p-value = 0.04056
#alternative hypothesis: two-sided
ks.test(log_pval_hmc_sites, log_pval_medians_random_sites)

# Plot the cdf of the age-related and random P-values
cdf_hmc_pval = ecdf(log_pval_hmc_sites)
cdf_random_pval = ecdf(log_pval_means_random_sites)
#quants <- quantile(-log10(pval_hmc_sites), probs = seq(0, 1, 1/3572))
#tab2 <- as.data.frame(cbind(quants, -log10(pval_hmc_sites)))

# Figure for locus-specific TCGA vs normal
png("05.Cancer_Comparison/Figures/TCGA_Normal_Volcano2.png", height=8*250, width=10*300, pointsize = 16, res=300)
p = ggplot(results, aes(Coefficient, -log10(P_value))) +
  geom_point(aes(col=limma_model))  + xlim(-2.5, 2.5) +
  xlab("coefficient") + ylab("-log10(P-value)") +
  scale_color_manual(values=c("black", "red")) + 
  geom_hline(yintercept = -log10(.05/dim(results)[1]), color = "red", linetype = "dashed", size = 1.3) + 
  labs(color="limma model \n significance") + 
  theme(legend.key = element_blank()) + 
  #ggtitle("Tumor(n=753) -Normal(n=95) Differential \n DNA methylation (TCGA)") + 
  geom_vline(xintercept = 0, colour = "black", linetype="dotted") + 
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
1428-887

# plot P-valye cumulative distribution plots
png("05.Cancer_Comparison/Figures/TCGA_Normal_Pvalue_ecdf.png", height=7*300, width=7*275, pointsize = 12, res=300)
plot(cdf_hmc_pval, main="Tumor-Normal Differential \n DNA methylation (TCGA)", 
     xlab="-log10(P-value)", ylab="Cumulative proportion - P-value", col="red", 
     cex.axis = 1.45, cex.lab = 1.45, las = 1)
points(log_pval_hmc_sites[order(log_pval_hmc_sites)], col="red", cdf_hmc_pval(log_pval_hmc_sites)[order(cdf_hmc_pval(log_pval_hmc_sites))], cex = 0.3)
plot(cdf_random_pval, add=TRUE, col="black") 
points(log_pval_means_random_sites[order(log_pval_means_random_sites)], col="black", cdf_random_pval(log_pval_means_random_sites)[order(cdf_random_pval(log_pval_means_random_sites))], cex = 0.3)
legend(50, 0.5, c("High 5hmC CpGs", "Random CpGs"), lwd=3, bty="n", col=c("red","black"), cex = 1.3)
text(70, 0.3, "Kolmogorov-Smirnov test P = 0.046", cex = 1.3) 
dev.off()

png("05.Cancer_Comparison/Figures/x.png", height=7*300, width=5*275, pointsize = 12, res=300)
boxplot(log_pval_hmc_sites, log_pval_means_random_sites)
dev.off()

