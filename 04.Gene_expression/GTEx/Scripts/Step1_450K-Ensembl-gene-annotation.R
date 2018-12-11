#######################################################################################################################
# Annotate 450K CpGs with associated transcripts using Ensembl gene identifiers
#######################################################################################################################

rm(list=ls())
# packages 
library(genomation)
library(biomaRt)
library(ensembldb)
library(FDb.InfiniumMethylation.hg19)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(doParallel)
library(data.table)
# set directories
base_dir <- "~/Dropbox (Christensen Lab)/NDRI_Breast_5hmC_update/"
dir_1 <- paste0(base_dir, "04.Gene_expression/GTEx/Files/")
dir_2 <- paste0(base_dir, "02.Characterization_5hmC_levels/Files/")
dir_3 <- paste0(base_dir, "04.Gene_expression/GTEx/Scripts/")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load requird data-sets 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# expression from RNA-seq analyses (GTEx)
exp <- fread(paste0(dir_1, "GTEx_Analysis_2016-01-15_v7_RSEMv1.2.22_transcript_tpm.txt"), header=TRUE, stringsAsFactors = F)

# get 450K annotation data 
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# probes that passed QC is normal breast tissue for 5hmC study 
load(paste0(dir_2, "annotation_good_probes.Rdata"))
annot_good <- AnnotSel

# high 5hmC probes 
high_5hmC <- readRDS(paste0(dir_2, "high_5hmc_top1%_5hmC.rds"))

########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# subset annotation for good probes & add indicators for high 5hmC 
########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# subset complete annotation for good probes 
ann <- ann[rownames(ann) %in% annot_good$TargetID,]
rm(annot_good)

# add high 5hmC indicator to annotation data 
ann$top_5hmC <- NA
ann$top_5hmC[rownames(ann) %in% high_5hmC$ID] <- "1"
ann$top_5hmC[!rownames(ann) %in% high_5hmC$ID] <- "0"
ann$top_5hmC <- factor(ann$top_5hmC, levels = c("1", "0"))

########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# use BioMart to map ENST ids to RefSeq ids
########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# split into stable ENST IDs 
transcript_split <- as.character(sapply(as.character(exp$transcript), function(x) strsplit(x, "[.]")[[1]][1]))

# select Mart to use from ENSEMBL (GRCH37/hg19)
ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL",
                  host="grch37.ensembl.org",
                  path="/biomart/martservice",
                  dataset="hsapiens_gene_ensembl")
# select human dataset  
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

# store attributes of this Mart for viewing in df 
ensembl_genes_atts <- listAttributes(mart)

# extract Refseq mapping for ENSEBL genes from this Mart 
ensembl_transcripts <- transcript_split
# do first for RefSeq mRNAs (due to unknown error when trying to pull to many attributes at once)
ensb_refseq_mrna <- getBM(filters = "ensembl_transcript_id", 
                     attributes = c("ensembl_transcript_id", "refseq_mrna"),
                     values = ensembl_transcripts,
                     mart = mart)
ensb_refseq_mrna$refseq_mrna[ensb_refseq_mrna$refseq_mrna==""] <- NA

# and for RefSeq ncRNAs
ensb_refseq_ncRNA <- getBM(filters = "ensembl_transcript_id", 
                     attributes = c("ensembl_transcript_id", "refseq_ncrna"),
                     values = ensembl_transcripts,
                     mart = mart)
ensb_refseq_ncRNA$refseq_ncrna[ensb_refseq_ncRNA$refseq_ncrna==""] <- NA

########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# link UCSC refGene accession from 450K annotation to ENST ids 
########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# function to split entries from UCSC_RefGene_Accession in 450K annotation and return as list
split_funct1 <- function(x) {
  out <- list()
  if(x==""){
    out <- NA
  } else {
    out <- strsplit(x, ";")[[1]]
  }
  out
}
# apply to 450K annotation file 
split.1 <- lapply(ann$UCSC_RefGene_Accession, split_funct1)

# source function to link UCSC_RefGene_Accessions to ENST ids using ENSEMBL db data from BioMart
source(paste0(dir_3, "get_ENST_from_RefSeq.R"))
# initiate cores 
num_cores <- detectCores()-1
# setup cluster 
cl <- makeCluster( num_cores)
# export items to cluster 
clusterExport(cl, c("split.1", "get_ENST_from_RefSeq", "ensb_refseq_mrna", "ensb_refseq_ncRNA"))
# apply the function to get ENSTs
### mRNA
ENST_mrna <- parLapply(cl, split.1, get_ENST_from_RefSeq, ensb_refseq_mrna, "refseq_mrna")
### ncRNA
ENST_ncrna <- parLapply(cl, split.1, get_ENST_from_RefSeq, ensb_refseq_ncRNA, "refseq_ncrna")
# close cluster 
stopCluster(cl)

########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# add ENST ids to annotation file 
########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cat_custom <- function(x, collapse){
  out <- list()
  t1 <- c()
  if(length(x)==1){
    out <- x
  } else if(length(x)>1){
    for(j in 1:length(x)){
      t1[j] <- x[j]
    }
    if(collapse == TRUE){
      out <- paste(t1, collapse = ";")
    } else {
      out <- t1
    }
  } else{
    out <- NA
  }
}
# apply to ENST mrna 
ENST_mrna_2.1 <- lapply(ENST_mrna, cat_custom, FALSE)
ENST_mrna_2.2 <- lapply(ENST_mrna, cat_custom, TRUE)
# add names for CGs to ENST_mrna_2.1
names(ENST_mrna_2.1) <- rownames(ann)
# add ENST_mrna_2.2 to annotation file 
ann$ENST_gene_accession_mRNA <- unlist(ENST_mrna_2.2)

# apply to ENST ncRNA
ENST_ncrna_2.1 <- lapply(ENST_ncrna, cat_custom, FALSE)
ENST_ncrna_2.2 <- lapply(ENST_ncrna, cat_custom, TRUE)
# add names for CGs to ENST_mrna_2.1
names(ENST_ncrna_2.1) <- rownames(ann)
# add ENST_mrna_2.2 to annotation file 
ann$ENST_gene_accession_ncRNA <- unlist(ENST_ncrna_2.2)

# save annotation with this info
saveRDS(ann, file = paste0(dir_1, "450K_annotation_for_high_quality_probes_with_ENST_accessions.rds"))

########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# identify ENST associated with 450K CpGs and also high 5hmC CGs 
########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# get ENSTs MRNA associated with CG on 450K
ENST_mrna_2.1_all <- unlist(ENST_mrna_2.1)
# remove NAs, to return vector of ENSTs with high 5hmC CGs
ENST_mrna_2.1_all <- ENST_mrna_2.1_all[-which(is.na(ENST_mrna_2.1_all))]
names(ENST_mrna_2.1_all) <- NULL
ENST_mrna_2.1_all <- unique(ENST_mrna_2.1_all)


# get ENSTs NCRNA associated with CG on 450K
ENST_ncrna_2.1_all <- unlist(ENST_ncrna_2.1)
# remove NAs, to return vector of ENSTs with high 5hmC CGs
ENST_ncrna_2.1_all <- ENST_ncrna_2.1_all[-which(is.na(ENST_ncrna_2.1_all))]
names(ENST_ncrna_2.1_all) <- NULL
ENST_ncrna_2.1_all <- unique(ENST_ncrna_2.1_all)

ENST_450K <- c(ENST_mrna_2.1_all, ENST_ncrna_2.1_all)

# save this list for use in enrichment testing 
saveRDS(ENST_450K, file = paste0(dir_1, "ENST_accessions_mapping_to_450K.rds"))

#### 
# subset ENST mrna for CGs associated with high 5hmC CGs 
ENST_mrna_2.1_hmC <- ENST_mrna_2.1[names(ENST_mrna_2.1) %in% high_5hmC$ID]
# drop names and unlist
names(ENST_mrna_2.1_hmC) <- NULL
ENST_mrna_2.1_hmC_2 <- unique(unlist(ENST_mrna_2.1_hmC))
# remove NAs, to return vector of ENSTs with high 5hmC CGs
ENST_mrna_2.1_hmC_2 <- ENST_mrna_2.1_hmC_2[-which(is.na(ENST_mrna_2.1_hmC_2))]

# subset ENST ncrna for CGs associated with high 5hmC CGs 
ENST_ncrna_2.1_hmC <- ENST_ncrna_2.1[names(ENST_ncrna_2.1) %in% high_5hmC$ID]
# drop names and unlist
names(ENST_ncrna_2.1_hmC) <- NULL
ENST_ncrna_2.1_hmC_2 <- unique(unlist(ENST_ncrna_2.1_hmC))
# remove NAs, to return vector of ENSTs with high 5hmC CGs
ENST_ncrna_2.1_hmC_2 <- ENST_ncrna_2.1_hmC_2[-which(is.na(ENST_ncrna_2.1_hmC_2))]

ENST_hmC <- c(ENST_mrna_2.1_hmC_2, ENST_ncrna_2.1_hmC_2)

# save this list for use in enrichment testing 
saveRDS(ENST_hmC, file = paste0(dir_1, "ENST_accession_with_high_5hmC_CGs.rds"))

