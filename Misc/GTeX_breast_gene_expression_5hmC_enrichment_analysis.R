



# directories
base_dir <- "~/Dropbox (Christensen Lab)/NDRI_Breast_5hmC_update/"
dir_1 <- paste0(base_dir, "GTEx/Files/")
dir_2 <- paste0(base_dir, "02.Characterization_5hmC_levels/Files/")

# complete 450K annotation for all CGs w/ GTEx breast high vs low annotated expression data 
annot <- readRDS(paste0(dir_1, "450_annotation_with_GTEx_Breast_expression_annotation.rds"))

# probes that passed QC is normal breast tissue for 5hmC study 
load(paste0(dir_2, "annotation_good_probes.Rdata"))
annot_good <- AnnotSel

# high 5hmC probes 
high_5hmC <- readRDS(paste0(dir_2, "high_5hmc_top1%_5hmC.rds"))
# GTEx breast expression percentiles 
genes <- readRDS(paste0(dir_1, "GTEx_breast_expression_percentiles.rds"))

# ENSG IDs mapping to any 450K CG
ENSG_450K <- readRDS(paste0(dir_1, "ENSG_accessions_mapping_to_450K.rds"))

# ENSG gene IDs mapping to high 5hmC CGs
ENSG_high_hmC <- readRDS(paste0(dir_1, "ENSG_accession_with_high_5hmC_CGs.rds"))

########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# split into stable ENSG IDs 
genes$genes_split <- sapply(as.character(genes$gene), function(x) strsplit(x, "[.]")[[1]][1])
genes$genes_split <- as.character(genes$genes_split)

# how many of the 450K genes are present in the dataset
table(genes$genes_split %in% ENSG_450K)

# subset 'genes' to only genes represented on 450K array 
ind1 <- genes$gene[genes$genes_split %in% ENSG_450K]
genes_450K <- genes[ind1,]

# create indicator variable in 'genes' for if each gene maps to a high 5hmC CG 
genes_450K$hmC_CGs <- NA
genes_450K$hmC_CGs[genes_450K$genes_split %in% ENSG_high_hmC] <- "1"
genes_450K$hmC_CGs[!genes_450K$genes_split %in% ENSG_high_hmC] <- "0"
genes_450K$hmC_CGs <- factor(genes_450K$hmC_CGs, levels = c("1", "0"))

# check the counts 
table(genes_450K$hmC_CGs)

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

genes$exp_50_50 <- sapply(genes$percentile, annotate_exp_percentiles, 0.50)
genes$exp_60_40 <- sapply(genes$percentile, annotate_exp_percentiles, 0.60)
genes$exp_70_30 <- sapply(genes$percentile, annotate_exp_percentiles, 0.70)
genes$exp_80_20 <- sapply(genes$percentile, annotate_exp_percentiles, 0.80)
genes$exp_90_10 <- sapply(genes$percentile, annotate_exp_percentiles, 0.90)


########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# run enrichment tests 
########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# run MH test 
res1 <- fisher.test(genes_450K$exp_50_50, genes_450K$hmC_CGs)
res1

res1 <- fisher.test(genes_450K$exp_50_50, genes_450K$hmC_CGs)
res1

res1 <- fisher.test(genes_450K$exp_50_50, genes_450K$hmC_CGs)
res1

res1 <- fisher.test(genes_450K$exp_50_50, genes_450K$hmC_CGs)
res1
# print summary of results 
res_summary <- function(res_in, split){
  P_val <- round(res_in$p.value, digits = 4)
  OR <- round(res_in$estimate, digits = 2)
  CI_95 <- paste0("(", round(res_in$conf.int[1], digits = 2), "-", round(res_in$conf.int[2], digits = 2), ")")
  paste0("Gene Expression split High-Low ", split, ": ", "P-value = ", P_val, ", ", OR, " ", CI_95)
}
res_summary(res1, "50-50")





# re-make genes file with expression split into tertiles 
##### if you think this is the best way... 
##### high vs low based on various splits of genes may make more sense... or quartiles and do vs rest 
##### how do you get an OR if you make 
# subset genes file for genes on 450K array (ENSG)
# create indicator variable in genes for high 5hmC CG or no high 5hmC CG
# use genes data frame to create contingency tables for enrichment tests 


# need to consider if you are doing this analysis right... 
# you are testing here, for enrichment of high 5hmC CGs amoung CpGs that track to highly expressed genes 
# maybe it makes more sense to test for enirhment of high 5hmC CpGs among highly expressed genes
# the latter probs requires either chi-sq or fisher test... check chi-sq OR 
# would be nice to rpeort OR rather than %s


# also need a string of all gene possible on 450K as this needs to be your background!!! VERY IMPORTANT 

# should also check that ORs for your other work (ENCODE enrichment analyses etc.) are being calculated right...
# fir fisher test, i think it matters what you set as which label maybe??? 
# you can test this by changing labels!!!! and seeing if you get same results 
# also check their directionality 





