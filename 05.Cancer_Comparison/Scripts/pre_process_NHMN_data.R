#######################################################################################################################
# Read in and pre-process methylation data (from iDATS) for normal and DCIS breast tissue samples from NH mammography network
# GSE66313 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE66313)
#######################################################################################################################
rm(list = ls())
library(minfi)
library(GEOquery)
library(ENmix)
library(RnBeads)
library(RnBeads.hg19)
library(doParallel)
library(IlluminaHumanMethylation450kmanifest) # required for funnorm
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
setwd("/Users/Owen 1/Dropbox (Christensen Lab)/NDRI_Breast_5hmC_update/")
#setwd("/Users/Owen 1/Dropbox (Christensen Lab)/NDRI_Breast_5hmC_update/05.Cancer_Comparison/Files/")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# prepare idats and sample sheet for import 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# get GSE files 
#getGEOSuppFiles("GSE66313")
#untar("GSE66313/GSE66313_RAW.tar", exdir = "05.Cancer_Comparison/Files/GSE66313/idat")
idat_names <- list.files("05.Cancer_Comparison/Files/GSE66313/idat", pattern = "idat")
head(idat_names)

# unzip the iDATs
#idatFiles <- list.files("05.Cancer_Comparison/Files/GSE66313/idat", pattern = "idat.gz$", full = TRUE)
#sapply(idatFiles, gunzip, overwrite = TRUE)

# getGEO("GSE66313")
geoMat <- readRDS("05.Cancer_Comparison/Files/GSE66313.rds")
sample_annot <- pData(geoMat[[1]])
head(sample_annot)

# extract GSM ids 
GSM_id <- sapply(strsplit(idat_names, "_"), function(x) x[1])
# check they are the same (order and characters) as rownames for pheno data 
all(rownames(sample_annot)==unique(GSM_id))

# index idat_names to green idats only and extract the sentrix ID + POS from these, then add to sample annot
idat_names_Grn <- list.files("05.Cancer_Comparison/Files/GSE66313/idat", pattern = "_Grn.idat")
idat_names_Red <- list.files("05.Cancer_Comparison/Files/GSE66313/idat", pattern = "_Red.idat")
sample_annot$Sentrix_ID <- sapply(strsplit(idat_names_Grn, "_"), function(x) paste0(x[1], "_", x[2]))
sample_annot$Sentrix_Position <- sapply(strsplit(idat_names_Grn, "_"), function(x) x[3])
sample_annot$Sample_ID <- sapply(1:nrow(sample_annot), function(x) paste0(sample_annot$Sentrix_ID[x], "_", sample_annot$Sentrix_Position[x]))
# add complete iDAT file name too sample annot too 
sample_annot$idat_Grn <- idat_names_Grn
sample_annot$idat_Red <- idat_names_Red

# reorder columns to desired order 
sample_annot_2 <- sample_annot[,c("Sample_ID", "tissue:ch1", "Sentrix_ID", "Sentrix_Position", "subject age:ch1", "idat_Grn", "idat_Red", colnames(sample_annot)[1:37])]
# save file as csv 
#write.csv(sample_annot_2, "GSE66313/NHMN_samples.csv")
rm(sample_annot, sample_annot_2)

###################################################################################################
# DNA METHYLATION DATA PREPROCESSING AND QC USING minfi or RnBeads
# Author: Lucas A Salas lucas.a.salas.diaz@dartmouth.edu
#
# Notes:
# 1) Before start prepare your idat files on a single directory, if idats are located in several files use a recursive call
# 2) Sample sheet information (*.csv) differs between minfi (include 7 lines of Genome Studio header) and RnBeads (no header), use the appropiate wrapper
# 3) The main difference between the algorithms is the filtering approach, when comparing datasets between packages use a common filtering approach
###################################################################################################
#APPROACH 1: minfi (Hansen & Aryee, 2013) 
#source("https://bioconductor.org/biocLite.R")
#biocLite("minfi")
getwd()
#Set your idat directory
setwd("05.Cancer_Comparison/Files/GSE66313/idat")
#Your sample sheet might be in the working directory if it is not change it accordingly
sheet<-read.metharray.sheet(getwd())
#Warning: There should not be only this csv in your target directory and any other subfiles. If not the process will stop.
RGset <- read.metharray.exp(targets = sheet, extended = TRUE) #S4 object 27K, 450K or EPIC arrays are imported automatically 
names <- pData(RGset)$Sample_Name #or equivalent
groups <- pData(RGset)$Sample_Group#or equivalent

#plot controls
plotCtrl(RGset,IDorder=NULL)
qcscore1<-QCinfo(RGset, detPthre=0.000001, nbthre=3, CpGthre=0.05, samplethre=0.05,outlier=TRUE, distplot=T)
qcscore2<-QCinfo(RGset, detPthre=0.000001, nbthre=3, CpGthre=0.05, samplethre=0.10,outlier=TRUE, distplot=T)
qcscore1$badsample # 17 low quality samples 
qcscore2$badsample # 6 low quality samples 

#Original threshold
#mdat=preprocessNoob(epicData)
#rgSet=QCfilter(mdat,qcinfo=qcscore2,outlier=TRUE)
#beta=bmiq.mc(rgSet,nCores=20)
#pheno<-as.data.frame(pData(rgSet))
#densityPlot(beta, sampGroups = groups, pal = c16)
#mdsPlot(beta, sampGroups = pheno$Slide, pal = c16, legendPos = "bottomright", legendNCol= sunglasses,
          #cov<-data.frame(Subject=factor(pheno$Sample_Name),
            #CellType=factor(pheno$CellType),
            #slide=factor(pheno$Slide)),
        #array=factor(pheno$Array))
        #pcrplot(beta, cov, npc=8)

# drop these 6 samples that have low bisulfite conversion intensities & high % of low quality CpGs
RGset2 <- RGset[, !colnames(RGset) %in% qcscore1$badsample]
sheet2 <- sheet[!sheet$Sample_Name %in% qcscore1$badsample, ]
table(sheet2$status)

#Density plots. First QC
dpi = 300
png("raw_density_bean_plot.png", width = dpi*15, height = dpi*10, res = dpi)
par(mfrow=c(1,2))
densityBeanPlot(RGset2, sampGroups=sheet2$status)
densityPlot(RGset2, sampGroups=sheet2$status)
dev.off()
#Plotting individual samples
#plotBetasByType(minfi::getBeta(RGset2)[,1])
#plotBetasByType(minfi::getBeta(RGset)[,1])
#Generate the simple QC report as a pdf. Sometimes the Beanplots and density plots are not properly printed use the functions above
qcReport(RGset2)
#Sample mix ups
png("snp_beta_heatmap.png", width = dpi*10, height = dpi*10, res = dpi)
heatmap(getSnpBeta(RGset2))
dev.off()

# check BS conversion controls 
controlStripPlot(RGset2, controls="BISULFITE CONVERSION II")
#OOB (Out of band plots)
#convert to a Methylset
Mset<-preprocessRaw(RGset2)
#Fix outliers 
#You must have the proper annotation
#Mset@annotation
Mset2<-minfiQC(Mset, fixOutliers = T, verbose = T)
#Warning message:
  #In .getSex(CN = CN, xIndex = xIndex, yIndex = yIndex, cutoff = cutoff) :
  #An inconsistency was encountered while determining sex. One possibility is that only one sex is present. We recommend further checks, for example with the plotSex function.

#Check extreme outliers not everything should be discarded
png("QC_plot.png", width = dpi*10, height = dpi*10, res = dpi)
plotQC(Mset2$qc)
dev.off()
#Check predicted sex 
png("sex_QC_plot.png", width = dpi*10, height = dpi*10, res = dpi)
plotSex(Mset2$qc)
dev.off()
#Check for batch
sheet2$slide_2 <- sapply(1:nrow(sheet2), function(x) strsplit(sheet2$Slide[x], "_")[[1]][2])
png("bacth_MDS_QC_plot.png", width = dpi*15, height = dpi*10, res = dpi)
par(mfrow=c(1,2))
mdsPlot(Mset2$object, numPositions = 1000, sampGroups=sheet2$Array, pch = 16)
mdsPlot(Mset2$object, numPositions = 1000, sampGroups=sheet2$slide_2, pch = 16)
#mdsPlot(Mset2$object, numPositions = 1000, sampGroups=sheet2$status, pch = 16)
dev.off()
addQC(Mset, qc= Mset2$qc)
#Repeat the density plots. Second QC. The object is no longer an RGset2 use getBeta
densityPlot(minfi::getBeta(Mset),  sampGroups=sheet2$status)
#Convert into a GRset Map to Genome
GRset<-mapToGenome(Mset)
#Normalization (Funnorm), background correction (methylumi.noob), dye correction. Change according to your dataset
GrsetFUNnorm<-preprocessFunnorm(RGset2, nPCs=2, sex = NULL, bgCorr = TRUE, dyeCorr = TRUE, verbose = TRUE)

#Repeat the density plots. Second QC. The object is no longer an RGset2 use getBeta
png("normalized_density_bean_plot_5cutoff.png", width = dpi*15, height = dpi*10, res = dpi)
par(mfrow=c(1,2))
densityBeanPlot(minfi::getBeta(GrsetFUNnorm), sampGroups=sheet2$status)
densityPlot(minfi::getBeta(GrsetFUNnorm), sampGroups=sheet2$status)
dev.off()
#Plotting individual samples
densityPlot(as.matrix(minfi::getBeta(GrsetFUNnorm)[,3]))
#Extract your phenotype
pData(GrsetFUNnorm)
# extract betas 
betas_minfi <- minfi::getBeta(GrsetFUNnorm)
head(betas_minfi)
colnames(betas_minfi)
head(rownames(betas_minfi))

# save betas
saveRDS(betas_minfi, file = "betas_bgrcorr_funnorm_minfi.rds")
saveRDS(GrsetFUNnorm, file = "GrsetFUNnorm_bgrcorr_funnorm_minfi.rds")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# import idats + save as rnbeads object 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# set data directory 
data.dir <- getwd()

# set idat directory 
idat.dir <- file.path(data.dir, "GSE66313/idat/")

# save directory of sample sheet2 and initailze directory for report + analysis 
# NOTE: Case-control should be the second column to run correctly the differential methylation analyses
sample.annotation <- file.path(data.dir, "GSE66313/NHMN_samples.csv")
analysis.dir <- file.path(data.dir, "analysis")
report.dir <- file.path(analysis.dir, "reportsNHMN")
rnb.initialize.reports(report.dir)
logger.start(fname=NA)

# set data source directory for import function 
data.source <- c(idat.dir, sample.annotation)

# Import idats 
result <- rnb.run.import(data.source=data.source, data.type="infinium.idat.dir", dir.reports=report.dir)
rnb_NHMN <- result$rnb.set

# save the set (NOTE: this did not work writing directly to dropbox)
save.rnb.set(rnb_NHMN, "/Users/Owen 1/Thesis/NHMN_rnb_set", archive = TRUE)
rm(result, rnb_NHMN, GSM_id, idat_Red, idat_Grn, idat_names, idat_names_Red, idat_names_Grn, idat.dir, sample.annotation)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# index rnb set for subjects filtered in minfi QC
# index rnbset for SNP betas 
# replace rnbset betas w/ backgnd corrected funnormized betas from minfi 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("/Users/Owen 1/Dropbox (Christensen Lab)/NDRI_Breast_5hmC_update/")

# read in rnb set data 
rnb.set.unfiltered <- load.rnb.set("05.Cancer_Comparison/Files/NHMN_rnb_set.zip")
betas_rnb <- meth(rnb.set.unfiltered, row.names=T)

# read in minfi betas 
betas_minfi <- readRDS("05.Cancer_Comparison/Files/GSE66313/idat/betas_bgrcorr_funnorm_minfi.rds")
snp_betas_minfi <- readRDS("snp_betas_from_mset2_minfi.rds")

# check dimensions of both betas 
dim(betas_rnb); dim(betas_minfi)

# why are there some that are missing from minfi betas
ind <- which(!rownames(betas_rnb) %in% rownames(betas_minfi))
betas_rnb_ind <- betas_rnb[ind,]
# seems that they are SNPs
# so SNPs in rnbset need to be replaced w/ minfi snp betas or removed from rnbset 

# filter out SNP betas from rnb set 
rnb.set.f1 <- remove.sites(rnb.set.unfiltered, ind)
betas.rnb.set.f1 <- meth(rnb.set.f1, row.names=T)
dim(betas.rnb.set.f1)

# drop samples removed during minfi QC
s2r <- samples(rnb.set.f1)[!samples(rnb.set.f1) %in% colnames(betas_minfi)]
rnb.set.f2 <- remove.samples(rnb.set.f1, s2r)
betas.rnb.set.f2 <- meth(rnb.set.f2, row.names=T)
#saveRDS(betas.rnb.set.f2, file = "betas.rnb.set.f2")

# make sure minfi betas and filtered rnbset are same dimensions and CG order 
dim(betas_minfi)
dim(betas.rnb.set.f2)
all(rownames(betas_minfi)==rownames(betas.rnb.set.f2))
all(samples(rnb.set.f2)==colnames(betas_minfi))

# replace rnbeads betas w/ betas prodiced by minfi QC and funnorm  
rnb.set.f2.2 <- rnb.set.f2
rnb.set.f2.2@meth.sites <- betas_minfi

# check they have been replaced properly 
head(meth(rnb.set.f2, row.names=T))
head(meth(rnb.set.f2.2, row.names=T))
meth(rnb.set.f2, row.names=T)[1:10,1:5]
meth(rnb.set.f2.2, row.names=T)[1:10,1:5]
#saveRDS(rnb.set.f2.2, file = "betas.rnb.set.f2.2.rds")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# QA/QC probe filtering 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

logger.start(fname=NA)
# set data directory 
data.dir <- "/Users/Owen 1/Thesis/rnbeads_temp/"
analysis.dir <- file.path(data.dir, "analysis")
report.dir <- file.path(analysis.dir, "reportsNHMN")
rnb.initialize.reports(report.dir)

# set up parallel computing 
parallel.isEnabled()
num.cores <- 10
parallel.setup(num.cores)
parallel.isEnabled()
options(fftempdir="/Users/Owen 1/fftempdir")

# set options for filtering of probe sets 
rnb.options(analysis.name ="NHMN_5hmC_breast_project", 
            gz.large.files = TRUE,
            normalization.method="none",
            normalization.background.method = "none",
            filtering.snp="yes",
            filtering.cross.reactive = TRUE,
            filtering.greedycut = TRUE,
            filtering.greedycut.pvalue.threshold = 0.05,
            filtering.sex.chromosomes.removal=TRUE, 
            identifiers.column="Sample_ID",
            covariate.adjustment.columns = NULL,
            differential = FALSE,
            export.to.csv = TRUE,
            export.to.ewasher = TRUE,
            disk.dump.big.matrices = TRUE,
            enforce.memory.management = TRUE)
report.dir <- file.path(analysis.dir, "reports_details")
rnb.initialize.reports(report.dir)
logger.start(fname=NA)

# perform filtering 
result <- rnb.run.preprocessing(rnb.set.f2.2, dir.reports = report.dir)
rnb.set.f3 <- result$rnb.set
gc()

# save the filtered rnB set 
save.rnb.set(rnb.set.f3, "/Users/Owen 1/Thesis/rnbeads_temp/NHMN_rnb_set_filtered")
gc()

# run some exploratory analysis
rnb.run.exploratory(rnb.set.f3, report.dir)
rm(list=ls())
gc()

