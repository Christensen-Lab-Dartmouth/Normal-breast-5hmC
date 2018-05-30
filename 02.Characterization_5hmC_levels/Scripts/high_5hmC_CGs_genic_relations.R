#######################################################################################################################
# Identify genes w/ the most 'high 5hmC CpGs' associated with them, and summarize these data 
#######################################################################################################################
rm(list = ls())
setwd("/Users/Owen 1/Dropbox (Christensen Lab)/NDRI_Breast_5hmC_update/")
library(genomation)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load Data 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# annotation data 
annot <- readRDS("02.Characterization_5hmC_levels/Files/annotation_good_probes_w_gene_parts.rds")

# high 5hmC probes 
top1_perc <- readRDS("02.Characterization_5hmC_levels/Files/high_5hmc_top1%_5hmC.rds")
high_5hmC_CpGs_2 <- top1_perc

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# generate GRanges object for high 5hmC CpGs 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# add high 5hmC indicator to annotation data 
annot$top_5hmC <- NA
annot$top_5hmC[annot$TargetID %in% high_5hmC_CpGs_2$id] <- 1
annot$top_5hmC[!annot$TargetID %in% high_5hmC_CpGs_2$id] <- 0
table(annot$top_5hmC)

#Generate new variables for MAPINFO. A probe lists one location, to get genomic ranges we need two locations
annot$hg19_start = annot$MAPINFO
annot$hg19_end = annot$MAPINFO + 1
#Make the 'CHR' integer variable a character
annot$chr = as.character(annot$chr)
names(annot) #To identify columns of interest so that they do not conflict (multiple genomic locations)
annot_sub = annot[, c("TargetID", "chr", "hg19_start", "hg19_end", "top_5hmC")]

# convert to GRanges object 
annot_5hmC <- annot_sub[annot_sub$top_5hmC==1,]
annot_5hmC_gr <- makeGRangesFromDataFrame(annot_5hmC, keep.extra.columns=TRUE, ignore.strand=TRUE, seqnames.field = "chr", start.field="hg19_start", end.field="hg19_end")
annot_5hmC_gr

# index for high 5hmC CGs
annot_2 <- annot[annot$top_5hmC==1,]

# replace those w/o refgene names w/ NA
annot_2$UCSC_RefGene_Name[annot_2$UCSC_RefGene_Name==""] <- NA

# count these to see how many are in intergenic regions 
table(is.na(annot_2$UCSC_RefGene_Name))

# find out which cgs have multiple genes associated with them 
pb = txtProgressBar(min = 0, max = dim(annot)[1], initial = 0) 
mul_genes <- list()
for(i in 1:3876){
  setTxtProgressBar(pb,i)  
  split <- strsplit(annot_2$UCSC_RefGene_Name, ";")[[i]]
  if(length(unique(split))==1){ mul_genes[[i]]<- unique(split) }
  else{ mul_genes[[i]] <- split }
}
close(pb)
count1 <- lapply(mul_genes, table)
count2 <- unlist(count1)
tab1 <- table(names(count2))
gene_counts <- data.frame(tab1)
colnames(gene_counts) <- c("Gene", "Count")
gene_counts <- gene_counts[order(gene_counts$Count, decreasing = TRUE),]
gene_counts$Gene <- as.character(gene_counts$Gene)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# look at CGs for  specific genes 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# create subsets to view genes with at least X many high 5hmC CGs
g3 <- gene_counts[gene_counts$Count >= 3,]
g4 <- gene_counts[gene_counts$Count >= 4,]
g5 <- gene_counts[gene_counts$Count >= 5,]
g6 <- gene_counts[gene_counts$Count >= 6,]
g7 <- gene_counts[gene_counts$Count >= 7,]

# find out how what the max number of genes associated with one CG is 
x1 <- c()
for(i in 1:length(count1)){
  x1[i] <- length(names(count1[[i]]))
}
table(x1)
# it is 5 genes, so need to make 5 new columns in the annotation data for each CG to 
# include the associated genes for each CG 
annot_2$gene_1 <- 0
annot_2$gene_2 <- 0
annot_2$gene_3 <- 0
annot_2$gene_4 <- 0
annot_2$gene_5 <- 0
annot_2$gene_1[is.na(annot_2$UCSC_RefGene_Name)] <- NA
annot_2$gene_2[is.na(annot_2$UCSC_RefGene_Name)] <- NA
annot_2$gene_3[is.na(annot_2$UCSC_RefGene_Name)] <- NA
annot_2$gene_4[is.na(annot_2$UCSC_RefGene_Name)] <- NA
annot_2$gene_5[is.na(annot_2$UCSC_RefGene_Name)] <- NA

annot3<- annot_2
# fill the new gene columns with the relevant gene names 
for(i in 1:length(count1)){
  if(length(names(count1[[i]]))==0){
    annot_2[i,34] <- NA
    } else if(length(names(count1[[i]]))>1){
    for(j in 1:length(count1[[i]])){
      annot_2[i,34+j] <- names(count1[[i]])[j]
    }
  } else {
    annot_2[i,34] <- names(count1[[i]])
  }
}

# check the counts in each new column 
gene_tab_count <- apply(annot_2[,34:38], 2, table)

# check each gene variable for gene of interest 
which(names(gene_tab_count$gene_1)==g7$Gene[2])
which(names(gene_tab_count$gene_2)==g7$Gene[2])
which(names(gene_tab_count$gene_3)==g7$Gene[2])
which(names(gene_tab_count$gene_4)==g7$Gene[2])

# index for the CGs in those genes
annot3 <- annot_2[which(annot_2$gene_1==g7$Gene[1]),]
annot4 <- annot_2[which(annot_2$gene_1==g7$Gene[2]),]
annot5 <- annot_2[which(annot_2$gene_1==g7$Gene[3] | annot_2$gene_2==g7$Gene[3]),]
annot6 <- annot_2[which(annot_2$gene_1==g7$Gene[4]),]
annot7 <- annot_2[which(annot_2$gene_1==g7$Gene[5]),]

# index for CGs from genes of interest 
mb1 <- annot_2[which(annot_2$gene_1==g7$Gene[2]),]
mb2 <- annot_2[which(annot_2$gene_2==g7$Gene[2]),]
mb3 <- annot_2[which(annot_2$gene_3==g7$Gene[2]),]
mb_com <- rbind(mb1, mb2, mb3)
arid1b_1 <- annot_2[which(annot_2$gene_1==g6$Gene[2]),]
dnmt3a_1 <- annot_2[which(annot_2$gene_1==g6$Gene[5]),]
foxo3_1 <- annot_2[which(annot_2$gene_1==g5$Gene[13]),]
TSG_annot <- list(mb_com, arid1b_1, dnmt3a_1, foxo3_1)
names(TSG_annot) <- c("MBNL1", "ARID1B", "DNMT3A", "FOXO3")
saveRDS(TSG_annot, file = "02.Characterization_5hmC_levels/Files/TSG_annot.rds")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# generate file w/ annotation info for CGs located in genes w/ at least 5 high 5hmC CGs
# for use as a supplemental figure 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
annot_w_context <- readRDS("02.Characterization_5hmC_levels/Files/annotation_good_probes_w_gene_parts.rds")
annot_w_context_2 <- annot_w_context[match(annot_2$TargetID, annot_w_context$TargetID),]
head(annot_2$TargetID)
head(annot_w_context_2$TargetID)
all(annot_2$TargetID==annot_w_context_2$TargetID)

# add genomic context vector to annotation file 
annot_2$genomic_context <- NA
for(i in 1:dim(annot_2)[1]){
  if(annot_2$promoters[i] ==1) {annot_2$genomic_context[i] <- "promoter"}
  if(annot_2$intron[i] ==1) {annot_2$genomic_context[i] <- "intron"}
  if(annot_2$exon[i] ==1) {annot_2$genomic_context[i] <- "exon"}
  if(annot_2$intergenic[i] ==1) {annot_2$genomic_context[i] <- "intergenic"}
}
table(annot_2$genomic_context)

# get annotation data for CGs associated with all genes with 5 or more 'high 5hmC CGs' 
annot_5_CGs <- list()
for(i in 1:nrow(g5)){
  annot_5_CGs[[i]] <- annot_2[which(annot_2$gene_1==g5$Gene[i] 
                               | annot_2$gene_2==g5$Gene[i] 
                               | annot_2$gene_3==g5$Gene[i]
                               | annot_2$gene_4==g5$Gene[i]
                               | annot_2$gene_5==g5$Gene[i]),]
}
table(annot_5_CGs[[1]]$genomic_context)

# bind all of the list elements together into one data frame 
annot_5_CGs_df <- do.call(rbind, annot_5_CGs)

# ouput 
#write.csv(annot_5_CGs_df, "02.Characterization_5hmC_levels/Files/annotation_info_5_high-5hmC-CGs_or_higher.csv")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# generate genomic context count table for all genes with 5 or more high 5hmC CpGs 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# split all gene/transcript/context data and oput into 1 data frame
annot_5_CGs_split <- data.frame(I(unlist(strsplit(annot_5_CGs_df$UCSC_RefGene_Name, ";"))), 
                            I(unlist(strsplit(annot_5_CGs_df$UCSC_RefGene_Accession, ";"))), 
                            unlist(strsplit(annot_5_CGs_df$UCSC_RefGene_Group, ";")))
colnames(annot_5_CGs_split) <- c("Gene", "Accession", "Context")

# generate matrix to fill with with gene/transcript/count data 
context_counts <- matrix(NA, nrow = length(unique(annot_5_CGs_split$Gene)), ncol = 7) 
colnames(context_counts) <- c("Gene", "No. high 5hmC CpGs", "I450K_CGs_in_gene", "Promoter", "Exon", "Intron", "Intergenic")
context_counts[, "Gene"] <- unique(annot_5_CGs_split$Gene)
#context_counts[, "No. high 5hmC CpGs"] <- unique(annot_5_CGs_split$Gene)

# remove the genes that dont have at least 5 high 5hmC CpGs but got included as a CG associated w/ a gene that does 
# have at least 5 high 5hmC CGs associated w/ it also has a relation to a transcript of that gene 
ind_to_drop <- which(context_counts[,"Gene"] %in% unique(annot_5_CGs_split$Gene)[!unique(annot_5_CGs_split$Gene)%in%g5$Gene])
context_counts <- context_counts[-ind_to_drop,]

# identify the genomic contexts present that need to be accounted for 
table(annot_5_CGs_split$Context)
genomic_context_counts <- lapply(annot_5_CGs, function(x) table(factor(x$genomic_context, levels=c("intron", "promoter", "exon", "intergenic"))))
table(names(unlist(genomic_context_counts)))

# iteratively fill in table for each gene 
for(i in 1:nrow(context_counts)){
  context_counts[i, 2] <- g5$Count[which(g5$Gene == context_counts[i, "Gene"])]
  if(genomic_context_counts[[i]]["promoter"]!=0) {context_counts[i, "Promoter"] <- genomic_context_counts[[i]]["promoter"][[1]]}
  if(genomic_context_counts[[i]]["exon"]!=0) {context_counts[i, "Exon"] <- genomic_context_counts[[i]]["exon"][[1]]}
  if(genomic_context_counts[[i]]["intron"]!=0) {context_counts[i, "Intron"] <- genomic_context_counts[[i]]["intron"][[1]]}
  if(genomic_context_counts[[i]]["intergenic"]!=0) {context_counts[i, "Intergenic"] <- genomic_context_counts[[i]]["intergenic"][[1]]}
  if(genomic_context_counts[[i]]["promoter"]==0) {context_counts[i, "Promoter"] <- 0}
  if(genomic_context_counts[[i]]["exon"]==0) {context_counts[i, "Exon"] <- 0}
  if(genomic_context_counts[[i]]["intron"]==0) {context_counts[i, "Intron"] <- 0}
  if(genomic_context_counts[[i]]["intergenic"]==0) {context_counts[i, "Intergenic"] <- 0}
}
head(context_counts)

# now count CGs associated w/ each gene accross array 
pb = txtProgressBar(min = 0, max = dim(annot)[1], initial = 0) 
mul_genes_genome_wide <- list()
for(i in 1:dim(annot)[1]){
  setTxtProgressBar(pb,i)  
  split <- strsplit(annot$UCSC_RefGene_Name, ";")[[i]]
  if(length(unique(split))==1){ mul_genes_genome_wide[[i]]<- unique(split) }
  else{ mul_genes_genome_wide[[i]] <- split }
}
close(pb)

# extract CpG count data for each gene 
count1_genome_wide <- lapply(mul_genes_genome_wide, table)
count2_genome_wide <- unlist(count1_genome_wide)
tab1_genome_wide <- table(names(count2_genome_wide))
# index for genes with more than 5 high 5hmC CGs associated with them 
tab1_genome_wide_2 <- tab1_genome_wide[match(as.character(context_counts$Gene), names(tab1_genome_wide))]

# add these counts to context_counts 
context_counts$I450K_CGs_in_gene <- tab1_genome_wide_2

# ouput file as csv 
write.csv(context_counts, "02.Characterization_5hmC_levels/Files/genomic_context_counts_info_5_high-5hmC-CGs_or_higher.csv")

