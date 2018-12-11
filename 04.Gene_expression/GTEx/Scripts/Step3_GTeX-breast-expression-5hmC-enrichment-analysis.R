#######################################################################################################################
# Test enrichment of 5hmC among highly expressed genes in normal breast tissue 
#######################################################################################################################
rm(list=ls())
setwd("/Users/Owen 1/Dropbox (Christensen Lab)/NDRI_Breast_5hmC_update/")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load requird data-sets 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# set directories
base_dir <- "~/Dropbox (Christensen Lab)/NDRI_Breast_5hmC_update/"
dir_1 <- paste0(base_dir, "04.Gene_expression/GTEx/Files/")
dir_2 <- paste0(base_dir, "02.Characterization_5hmC_levels/Files/")

# high 5hmC probes 
high_5hmC <- readRDS(paste0(dir_2, "high_5hmc_top1%_5hmC.rds"))

# GTEx breast expression percentiles 
genes <- readRDS(paste0(dir_1, "GTEx_breast_expression_percentiles.rds"))

# ENST IDs mapping to any 450K CG
ENST_450K <- readRDS(paste0(dir_1, "ENST_accessions_mapping_to_450K.rds"))

# ENST gene IDs mapping to high 5hmC CGs
ENST_high_hmC <- readRDS(paste0(dir_1, "ENST_accession_with_high_5hmC_CGs.rds"))

########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# pre-process percentile dataset
# subset to only 450K associated transcripts
# add indicator for transcripts mapping to high 5hmC CpGs 
########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# split into stable ENST IDs 
genes$transcript_split <- sapply(as.character(genes$transcript), function(x) strsplit(x, "[.]")[[1]][1])
genes$transcript_split <- as.character(genes$transcript_split)

# how many of the 450K genes are present in the dataset
table(genes$transcript_split %in% ENST_450K)

# subset 'genes' to only genes represented on 450K array 
ind1 <- genes$transcript_split %in% ENST_450K
genes_450K <- genes[ind1,]

# create indicator variable in 'genes' for if each gene maps to a high 5hmC CG 
genes_450K$hmC_CGs <- NA
genes_450K$hmC_CGs[genes_450K$transcript_split %in% ENST_high_hmC] <- "1"
genes_450K$hmC_CGs[!genes_450K$transcript_split %in% ENST_high_hmC] <- "0"
genes_450K$hmC_CGs <- factor(genes_450K$hmC_CGs, levels = c("1", "0"))

########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# label genes into desired groups based on expression percentile 
########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# function to annotate genes df based on specified expression range 
annotate_exp_percentiles <- function(x, p1){
  temp <- c()
  if(x >= p1){
    temp <- "High"
  } else{
    temp <- "Low"
  }
  temp <- factor(temp, levels=c("High", "Low"))
}
# apply for splits of interest 
genes_450K$exp_75_25 <- sapply(genes_450K$percentile, annotate_exp_percentiles, 0.75)
genes_450K$exp_50_50 <- sapply(genes_450K$percentile, annotate_exp_percentiles, 0.50)

########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# run enrichment tests 
########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# make and look at contingency tables 
tab1 <- table(genes_450K$exp_75_25, genes_450K$hmC_CGs)
tab2 <- table(genes_450K$exp_50_50, genes_450K$hmC_CGs)

# run enrichment tests
res1 <- fisher.test(tab1)
res2 <- fisher.test(tab2)

# make results table 
results <- matrix("", nrow = 6, ncol = 6)
colnames(results) <- c("High 5hmC associated transcripts", 
                       "Other 450K transcripts", 
                       "OR", "LB.95.CI", "UB.95.CI", "P.value")
rownames(results) <- c("Q1+2+3 vs Q4", "High.1", "Low.1", 
                       "Q1+2 vs Q3+4", "High.2", "Low.2")
results[2, c(1:2)] <- c(tab1[1,1], tab1[1,2])
results[3, c(1:6)] <- c(tab1[1,1], tab1[2,2], 
                        round(res1$estimate[1], digits=2), round(res1$conf.int[1], digits=2), 
                        round(res1$conf.int[2], digits=2), res1$p.value)
results[5, c(1:2)] <- c(tab2[1,1], tab2[1,2])
results[6, c(1:6)] <- c(tab2[1,1], tab2[2,2], 
                        round(res2$estimate[1], digits=2), round(res2$conf.int[1], digits=2), 
                        round(res2$conf.int[2], digits=2), res2$p.value)

# write out results 
write.csv(results, file = paste0(dir_1, "GTEx_expression_5hmC_enrichment_results.csv"))
