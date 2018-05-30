###########################
# Genomic Regions Enrichment of Annotations Tool (GREAT) analysis of high 5hmC loci, using rGREAT package 
# Author: Owen Wilkins, Kevin Johnson
###########################
rm(list = ls())
setwd("/Users/Owen/Dropbox (Christensen Lab)/NDRI_Breast_5hmC_update/")
#devtools::install_github("jokergoo/rGREAT")
library(rGREAT)
library(LOLA)
library(GenomicRanges)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# rGREAT analysis of high 5hmC loci
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# load BED file containing high 5hmC CpGs 
bed <- readBed('02.Characterization_5hmC_levels/Files/Breast_high_5hmC.bed') 
background <- readBed('02.Characterization_5hmC_levels/Files/Breast_450K_Array_5hmC.bed') 

# submit GREAT job (specifying background set as 450K probes that passed QC)
job = submitGreatJob(gr = bed, bg = background, species = "hg19")
job

# get tables for each GO set ordered by variable of interest 
get_GO_tables <- function(results, GO_table_to_get, order_by){
  tb = getEnrichmentTables(results)
  # restrict to which GP processes are desired 
  tb <- tb[[GO_table_to_get]]
  # calculate FDR from raw P-values 
  tb$Hyper_FDR_Q_val <- p.adjust(tb$Hyper_Raw_PValue, method = "BH")
  # restrict to GO terms w/ FDR<0.05
  tb <- tb[tb$Hyper_FDR_Q_val<0.05,]
  # order results by FDR or enrichment value 
  if(order_by == "qval") 
    tb <- tb[order(tb$Hyper_FDR_Q_val, decreasing = FALSE),]
  if(order_by == "enrichment") 
    tb <- tb[order(tb$Hyper_Fold_Enrichment, decreasing = TRUE),]
  # return processed table 
  tb
}
# GO Biological Process
tb_bp_qval <- get_GO_tables(job, "GO Biological Process", "qval")
tb_bp_enrichment <- get_GO_tables(job, "GO Biological Process", "enrichment")
# GO Molecular Function 
tb_mf_qval <- get_GO_tables(job, "GO Molecular Function", "qval")
tb_mf_enrichment <- get_GO_tables(job, "GO Molecular Function", "enrichment")
# GO Cellular Component
tb_cc_qval <- get_GO_tables(job, "GO Cellular Component", "qval")
tb_cc_enrichment <- get_GO_tables(job, "GO Cellular Component", "enrichment")

# check avaiable ontologies & categories 
availableOntologies(job)
availableCategories(job)

# plot association graphs for genomic regions and genes 
par(mfrow = c(1, 3))
res = plotRegionGeneAssociationGraphs(job)

# write results files to csv 
write.csv(tb_bp_qval, file = "04.Enrichment_Analyses/rGREAT/Files/rGREAT_high_5hmC_GO_biological_process.csv")
write.csv(tb_mf_qval, file = "04.Enrichment_Analyses/rGREAT/Files/rGREAT_high_5hmC_GO_molecular_function.csv")
write.csv(tb_cc_qval, file = "04.Enrichment_Analyses/rGREAT/Files/rGREAT_high_5hmC_GO_cellular_component.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# perform additional analysis on only the high 5hmC CGs w/ sig. differential meth. status between 
# tumor and normal in TCGA data set 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# load in BED file for high 5hmC CGs w/ differential meth status below P-threshold in TCGA data set  
bed_sub <- readBed('02.Characterization_5hmC_levels/Files/Breast_genes_high_5hmC_low_P_TCGA.bed') 

# submit new job 
job = submitGreatJob(gr = bed_sub, bg = background, species = "hg19")
job

# GO Biological Process
tb_bp_qval_sub <- get_GO_tables(job, "GO Biological Process", "qval")
tb_bp_enrichment_sub <- get_GO_tables(job, "GO Biological Process", "enrichment")
# GO Molecular Function 
tb_mf_qval_sub <- get_GO_tables(job, "GO Molecular Function", "qval")
tb_mf_enrichment_sub <- get_GO_tables(job, "GO Molecular Function", "enrichment")
# GO Cellular Component
tb_cc_qval_sub <- get_GO_tables(job, "GO Cellular Component", "qval")
tb_cc_enrichment_sub <- get_GO_tables(job, "GO Cellular Component", "enrichment")

# check avaiable ontologies & categories 
availableOntologies(job)
availableCategories(job)

# plot association graphs for genomic regions and genes 
par(mfrow = c(1, 3))
res = plotRegionGeneAssociationGraphs(job)

# write results files to csv 
write.csv(tb_bp_qval_sub, file = "04.Enrichment_Analyses/rGREAT/Files/rGREAT_high_5hmC_GO_biological_process_low_P_TCGA.csv")
write.csv(tb_mf_qval_sub, file = "04.Enrichment_Analyses/rGREAT/Files/rGREAT_high_5hmC_GO_molecular_function_low_P_TCGA.csv")
write.csv(tb_cc_qval_sub, file = "04.Enrichment_Analyses/rGREAT/Files/rGREAT_high_5hmC_GO_cellular_component_low_P_TCGA.csv")

