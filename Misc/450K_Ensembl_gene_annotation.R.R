
library(genomation)
library(biomaRt)
library(ensembldb)
library(FDb.InfiniumMethylation.hg19)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(doParallel)

ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# directories
base_dir <- "~/Dropbox (Christensen Lab)/NDRI_Breast_5hmC_update/"
dir_1 <- paste0(base_dir, "GTEx/Files/")
dir_2 <- paste0(base_dir, "02.Characterization_5hmC_levels/Files/")
dir_3 <- paste0(base_dir, "GTEx/Scripts/")

# GTEx breast high vs low annotated expression data 
genes <- readRDS(paste0(dir_1, "GTEx_breast_expression_percentiles.rds"))

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
table(ann$top_5hmC)

########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# use BioMart to map ENSG ids to RefSeq ids
########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# split into stable ENSG IDs 
#genes$genes_split <- sapply(as.character(genes$gene), function(x) strsplit(x, "[.]")[[1]][1])
#genes$genes_split <- as.character(genes$genes_split)
genes$transcript_split <- sapply(as.character(genes$transcript), function(x) strsplit(x, "[.]")[[1]][1])
genes$transcript_split <- as.character(genes$transcript_split)

# retrun exp_grp as.character()
#genes$exp_grp <- as.character(genes$exp_grp)

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
#ensembl_genes <- genes$genes_split
ensembl_transcripts <- genes$transcript_split
# do first for RefSeq mRNAs (due to unknown error when trying to pull to many attributes at once)
ensb_refseq_mrna <- getBM(filters = "ensembl_transcript_id", 
                     attributes = c("ensembl_transcript_id", "refseq_mrna"),
                     values = ensembl_transcripts,
                     mart = mart)
ensb_refseq_mrna$refseq_mrna[ensb_refseq_mrna$refseq_mrna==""] <- NA
table(is.na(ensb_refseq_mrna$refseq_mrna))

# and for RefSeq ncRNAs
ensb_refseq_ncRNA <- getBM(filters = "ensembl_transcript_id", 
                     attributes = c("ensembl_transcript_id", "refseq_ncrna"),
                     values = ensembl_transcripts,
                     mart = mart)
ensb_refseq_ncRNA$refseq_ncrna[ensb_refseq_ncRNA$refseq_ncrna==""] <- NA
table(is.na(ensb_refseq_ncRNA$refseq_ncrna))

########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# link UCSC refGene accession from 450K annotation to ENSG ids 
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

# source function to link UCSC_RefGene_Accessions to ENSG ids using ENSEMBL db data from BioMart
source(paste0(dir_3, "get_ENSG_from_RefSeq.R"))
# initiate cores 
num_cores <- detectCores()-1
# setup cluster 
cl <- makeCluster(num_cores)
# export items to cluster 
clusterExport(cl, c("split.1", "get_ENSG_from_RefSeq", "ensb_refseq_mrna", "ensb_refseq_ncRNA"))
# apply the function and add results to genes df 
### mRNA
ENSG_mrna <- parLapply(cl, split.1, get_ENSG_from_RefSeq, ensb_refseq_mrna, "refseq_mrna")
### ncRNA
ENSG_ncrna <- parLapply(cl, split.1, get_ENSG_from_RefSeq, ensb_refseq_ncRNA, "refseq_ncrna")
# close cluster 
stopCluster(cl)
#save(ENSG_mrna, ENSG_ncrna, file = paste0(dir_1, "ENSG_temp.rds"))

########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# add ENSG ids to annotation file 
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
# apply to ENSG mrna 
ENSG_mrna_2.1 <- lapply(ENSG_mrna, cat_custom, FALSE)
ENSG_mrna_2.2 <- lapply(ENSG_mrna, cat_custom, TRUE)
# add names for CGs to ENSG_mrna_2.1
names(ENSG_mrna_2.1) <- rownames(ann)
# add ENSG_mrna_2.2 to annotation file 
ann$ENSG_gene_accession_mRNA <- unlist(ENSG_mrna_2.2)

# apply to ENSG ncRNA
ENSG_ncrna_2.1 <- lapply(ENSG_ncrna, cat_custom, FALSE)
ENSG_ncrna_2.2 <- lapply(ENSG_ncrna, cat_custom, TRUE)
# add names for CGs to ENSG_mrna_2.1
names(ENSG_ncrna_2.1) <- rownames(ann)
# add ENSG_mrna_2.2 to annotation file 
ann$ENSG_gene_accession_ncRNA <- unlist(ENSG_ncrna_2.2)

# save annotation with this info
saveRDS(ann, file = paste0(dir_1, "450K_annotation_for_high_quality_probes_with_ENSG_accessions.rds"))

########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# identify ENSG associated with 450K CpGs and also high 5hmC CGs 
########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# get ENSGs MRNA associated with CG on 450K
ENSG_mrna_2.1_all <- unlist(ENSG_mrna_2.1)
# remove NAs, to return vector of ENSGs with high 5hmC CGs
ENSG_mrna_2.1_all <- ENSG_mrna_2.1_all[-which(is.na(ENSG_mrna_2.1_all))]
names(ENSG_mrna_2.1_all) <- NULL
ENSG_mrna_2.1_all <- unique(ENSG_mrna_2.1_all)


# get ENSGs NCRNA associated with CG on 450K
ENSG_ncrna_2.1_all <- unlist(ENSG_ncrna_2.1)
# remove NAs, to return vector of ENSGs with high 5hmC CGs
ENSG_ncrna_2.1_all <- ENSG_ncrna_2.1_all[-which(is.na(ENSG_ncrna_2.1_all))]
names(ENSG_ncrna_2.1_all) <- NULL
ENSG_ncrna_2.1_all <- unique(ENSG_ncrna_2.1_all)

ENSG_450K <- c(ENSG_mrna_2.1_all, ENSG_ncrna_2.1_all)

# save this list for use in enrichment testing 
saveRDS(ENSG_450K, file = paste0(dir_1, "ENSG_accessions_mapping_to_450K.rds"))

#### 
# subset ENSG mrna for CGs associated with high 5hmC CGs 
ENSG_mrna_2.1_hmC <- ENSG_mrna_2.1[names(ENSG_mrna_2.1) %in% high_5hmC$ID]
# drop names and unlist
names(ENSG_mrna_2.1_hmC) <- NULL
ENSG_mrna_2.1_hmC_2 <- unique(unlist(ENSG_mrna_2.1_hmC))
# remove NAs, to return vector of ENSGs with high 5hmC CGs
ENSG_mrna_2.1_hmC_2 <- ENSG_mrna_2.1_hmC_2[-which(is.na(ENSG_mrna_2.1_hmC_2))]

# subset ENSG ncrna for CGs associated with high 5hmC CGs 
ENSG_ncrna_2.1_hmC <- ENSG_ncrna_2.1[names(ENSG_ncrna_2.1) %in% high_5hmC$ID]
# drop names and unlist
names(ENSG_ncrna_2.1_hmC) <- NULL
ENSG_ncrna_2.1_hmC_2 <- unique(unlist(ENSG_ncrna_2.1_hmC))
# remove NAs, to return vector of ENSGs with high 5hmC CGs
ENSG_ncrna_2.1_hmC_2 <- ENSG_ncrna_2.1_hmC_2[-which(is.na(ENSG_ncrna_2.1_hmC_2))]

ENSG_hmC <- c(ENSG_mrna_2.1_hmC_2, ENSG_ncrna_2.1_hmC_2)

# save this list for use in enrichment testing 
saveRDS(ENSG_hmC, file = paste0(dir_1, "ENSG_accession_with_high_5hmC_CGs.rds"))

