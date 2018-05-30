#######################################################################################################################
# Locus Overlap Analysis (LOLA) to identify genomic features enriched for high 5hmC CpGs
# Analysis script 
#######################################################################################################################
setwd("/Users/Owen 1/Dropbox (Christensen Lab)/NDRI_Breast_5hmC copy/")
rm(list = ls())
library(LOLA)
library(qvalue)
library(GenomicRanges) 
#devtools::install_github("databio/simpleCache") # install simplecache 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Locus Overlap Analysis of high 5hmC loci
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Define the universe & the gene set 
array <- readBed('02.Characterization_5hmC_levels/Files/Breast_450K_Array_5hmC.bed') 
hydroxymethyl <- readBed('02.Characterization_5hmC_levels/Files/Breast_high_5hmC.bed') 
universe <- GRanges(array)
geneset <- GRanges(hydroxymethyl)
str(hydroxymethyl)

# Compute gene enrichment set (see document "Using LOLA Core")
# Load the LOLA core (cached version; hg19/38, ENCODE TFBS, UCSC CGIs, Citrome epigenome)
lolaDB <- loadRegionDB("04.Enrichment_Analyses/LOLA/Files/LOLACore/hg19/")
locResults <- runLOLA(geneset, universe, lolaDB)

# check out results 
dim(locResults)
locResults[1:10,] #list the top 10 results
locResults[order(meanRnk, decreasing=FALSE),][1:20,]
locResults[order(maxRnk, decreasing=FALSE),][1:20,]

# look through results by collection 
locResults[collection=="sheffield_dnase", ][order(meanRnk, decreasing=FALSE),]
x <- locResults[collection=="encode_tfbs", ][order(meanRnk, decreasing=FALSE),]
locResults[collection=="codex", ][order(meanRnk, decreasing=FALSE),]
locResults[collection=="cistrome_epigenome", ][order(meanRnk, decreasing=FALSE),]
locResults[collection=="ucsc_features", ][order(meanRnk, decreasing=FALSE),]

# write results to tsv files
#writeCombinedEnrichment(locResults, outFolder = "04.Enrichment_Analyses/LOLA/Files/", includeSplits=TRUE)

# save LOLA results as Rdata
#save(locResults, file = "04.Enrichment_Analyses/LOLA/Files/lola_results.Rdata")
