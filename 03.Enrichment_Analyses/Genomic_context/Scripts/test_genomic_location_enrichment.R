#######################################################################################################################
# Determine genomic location (promoter, intron, exon, intergenic) of high 5hmC CpGs 
# Test enrichment of 5hmC within these regions against the 450K background set 
#######################################################################################################################
library(genomation)
library(GenomicRanges)
rm(list = ls())
setwd("/Users/Owen 1//Dropbox (Christensen Lab)/NDRI_Breast_5hmC_update//")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load Data 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# annotation data 
load("02.Characterization_5hmC_levels/Files/annotation_good_probes.Rdata")
annot <- AnnotSel
# high 5hmC probes 
top1_perc <- readRDS("02.Characterization_5hmC_levels/Files/high_5hmc_top1%_5hmC.rds")
high_5hmC_CpGs_2 <- top1_perc

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# determine genomic regions of high 5hmC CpGs  
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# add high 5hmC indicator to annotation data 
annot$top_5hmC <- NA
annot$top_5hmC[annot$TargetID %in% high_5hmC_CpGs_2$id] <- 1
annot$top_5hmC[!annot$TargetID %in% high_5hmC_CpGs_2$id] <- 0
table(annot$top_5hmC)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# make Grange object for good probes from 450K array  
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Generate new variables for MAPINFO. A probe lists one location, to get genomic ranges we need two locations
annot$hg19_start = annot$MAPINFO
annot$hg19_end = annot$MAPINFO + 1
#Make the 'CHR' integer variable a character
annot$chr = as.character(annot$chr)
names(annot) #To identify columns of interest so that they do not conflict (multiple genomic locations)
annot_sub = annot[, c(1, 3, 28, 29)]

#Create a 'GRanges' object from the Illumina annotation file 
Illumina_gr <- makeGRangesFromDataFrame(annot_sub, keep.extra.columns=TRUE, ignore.strand=TRUE, seqnames.field = "chr", start.field="hg19_start", end.field="hg19_end")
#Print the GRange object to see what it looks like
names(Illumina_gr) <- annot_sub$TargetID
Illumina_gr

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# make Grange object for good probes from 450K array  
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# read in UCSC ReF Gene file that contains assigments for promoters, introns, exons for hg19 
transcFeat = readTranscriptFeatures("02.Characterization_5hmC_levels/Files/UCSC_hg19_refGene.bed", up.flank = 2000, down.flank = 2000)
transcFeat

# The 'genomation' function annotates GRanges object as overlapping with promoter, exon, intron, or intergenic regions
ann_450k <- annotateWithGeneParts(Illumina_gr, transcFeat) 
head(ann_450k@members)
membership <- ann_450k@members
# return %s for features 
annotStats1 = getTargetAnnotationStats(ann_450k, percentage=TRUE, precedence=TRUE)
annotStats1

# annotate CGs using followin precedence: promoter>exon>intron (as some are associated w/ more than 1)
dim(membership)
membership2 <- membership
for(i in 1:nrow(membership2)){
  if(membership[i,1] == 1){membership2[i,] <- c(1,0,0)}
  if(membership[i,1] == 0 & membership[i,2] == 1){membership2[i,] <- c(0,1,0)}
  if(membership[i,1] == 0 & membership[i,2] == 0 & membership[i,3] == 1){membership2[1,] <- c(0,0,1)}
}
# check they all eqal total no of CGs
apply(membership2, 2, table)

# add counts for tyranscription features to annotation 
annot$promoters <- membership2[,"prom"]
annot$exon <- membership2[,"exon"]
annot$intron <- membership2[,"intron"]
annot$intergenic <- 0
annot$intergenic[which(annot$promoters==0 & annot$exon==0 & annot$intron==0)] <- 1
# check new variable counts 
table(annot$promoters)
table(annot$exon)
table(annot$intron)
table(annot$intergenic)

# make sure TRUE totals to all probes 
table(annot$promoters)[[2]] + table(annot$exon)[[2]] + table(annot$intron)[[2]] + table(annot$intergenic)[[2]]

# calculate % of CpGs in regions of interest 
pro_exo_int_high_hmC <- table(annot$top_5hmC, annot$intron)[4]+table(annot$top_5hmC, annot$promoters)[4]+table(annot$top_5hmC, annot$exon)[4]
intergenic_high_hmC <- table(annot$top_5hmC, annot$intergenic)[4]
pro_exo_int_high_hmC/3876*100

# save annotation containing gene part annotation 
saveRDS(annot, file = "02.Characterization_5hmC_levels/Files/annotation_good_probes_w_gene_parts.rds")

# test enrichmentr of high 5hmC CpGs among these genomic features against 450K background set 
transc_features <- c("promoters", "exon", "intron", "intergenic")
test_enrichment <- function(feature){
  MH_table <- table(annot$top_5hmC, annot[,paste0(feature)], annot$CpG_Regions)
  mantelhaen.test(MH_table, exact=TRUE)
}
res <- lapply(transc_features, test_enrichment)
names(res) <- c("promoters", "exon", "intron", "intergenic")
res

table(annot$top_5hmC, annot[,"exon"])
fisher.test(table(annot$top_5hmC, annot[,"exon"]))

table(factor(annot$top_5hmC, levels = c("1", "0")), factor(annot[,"exon"], levels = c("1", "0")))
fisher.test(table(factor(annot$top_5hmC, levels = c("1", "0")), factor(annot[,"exon"], levels = c("1", "0"))))


table(annot$top_5hmC, annot[,"intergenic"], annot$CpG_Regions)
mantelhaen.test(table(annot$top_5hmC, annot[,"intergenic"], annot$CpG_Regions), exact=TRUE)

# generate and output results table of 5hmC enrichment with each genomic feature 
tab <- matrix(NA, ncol = 5, nrow = 4)
colnames(tab) <- c("P-value", "Odds ratio (95% CI)", "OR", "LB", "UB")
rownames(tab) <- c("Promoters", "Exons", "Introns", "Intergenic")
results <- list(res$promoters, res$exon, res$intron, res$intergenic)
for(i in 1:length(results)){
  tab[i,1] <- paste0(results[[i]]$p.value)
  tab[i,2] <- paste0(format(round(results[[i]]$estimate, digits = 2), nsmall = 2), " (", 
                     format(round(results[[i]]$conf.int[1], digits = 2), nsmall = 2), "-", 
                     format(round(results[[i]]$conf.int[2], digits = 2), nsmall = 2), ")")
  tab[i,3] <- results[[i]]$estimate
  tab[i,4] <- results[[i]]$conf.int[1]
  tab[i,5] <- results[[i]]$conf.int[2]
}
write.csv(tab, file = "03.Enrichment_Analyses/Genomic_context/FIles/transcription_feature_hg19_enrichment.csv")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Test enrichment at CpG Island context 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# create one-hot encoding variables for CpGI context
### island
annot$Island <- 0
annot$Island[annot$CpG_Regions=="Island"] <- 1
### shore
annot$Shore <- 0
annot$Shore[annot$CpG_Regions=="Shore"] <- 1
### shelf
annot$Shelf <- 0
annot$Shelf[annot$CpG_Regions=="Shelf"] <- 1
### sea
annot$OpenSea <- 0
annot$OpenSea[annot$CpG_Regions=="OpenSea"] <- 1


# test enrichmentr of high 5hmC CpGs among these genomic features against 450K background set 
cgi_features <- c("Island", "Shore", "Shelf", "OpenSea")
test_enrichment_2 <- function(feature){
  cont_tab <- table(annot$top_5hmC, annot[,paste0(feature)])
  fisher.test(cont_tab)
}
res <- lapply(cgi_features, test_enrichment_2)
names(res) <- c("Island", "Shore", "Shelf", "OpenSea")
res


# generate and output results table of 5hmC enrichment with each genomic feature 
tab <- matrix(NA, ncol = 5, nrow = 4)
colnames(tab) <- c("P-value", "Odds ratio (95% CI)", "OR", "LB", "UB")
rownames(tab) <- c("Island", "Shore", "Shelf", "OpenSea")
results <- list(res$Island, res$Shore, res$Shelf, res$OpenSea)
for(i in 1:length(results)){
  tab[i,1] <- paste0(results[[i]]$p.value)
  tab[i,2] <- paste0(format(round(results[[i]]$estimate, digits = 2), nsmall = 2), " (", 
                     format(round(results[[i]]$conf.int[1], digits = 2), nsmall = 2), "-", 
                     format(round(results[[i]]$conf.int[2], digits = 2), nsmall = 2), ")")
  tab[i,3] <- results[[i]]$estimate
  tab[i,4] <- results[[i]]$conf.int[1]
  tab[i,5] <- results[[i]]$conf.int[2]
}
write.csv(tab, file = "03.Enrichment_Analyses/Genomic_context/FIles/CpGI_feature_hg19_enrichment.csv")

# test 5hmC enrichment among CGs located at promoter associated CpG islands 
annot$promoter_Island <- 0
annot$promoter_Island[annot$promoters==1 & annot$CpG_Regions=="Island"] <- 1
fisher.test(annot$top_5hmC, annot[,"promoter_Island"])
