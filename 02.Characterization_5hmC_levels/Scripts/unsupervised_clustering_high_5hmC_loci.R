#######################################################################################################################
# Unsupervised analyses of high 5hmC CpGs based on age, BMI, cell type, & total 5hmC content 
#######################################################################################################################
library(data.table) ; library(gplots) ; library(plyr)
rm(list = ls())
setwd("/Users/Owen 1/Dropbox (Christensen Lab)/NDRI_Breast_5hmC_update/")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load Data 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# annotation 
load("02.Characterization_5hmC_levels/Files/annotation_good_probes.Rdata")
annot <- AnnotSel
# betas
load("02.Characterization_5hmC_levels/Files/MethOxy_FunNorm_good_probes.Rdata")
# high 5hmC probes 
high_5hmC <- readRDS("02.Characterization_5hmC_levels/Files/high_5hmc_top1%_5hmC.rds")
# cellular proportions from RefFreeEWAS
cell_props <- read.csv("03.Cellular_proportions/Files/NDRI_Cell_Proportions.csv", stringsAsFactors = F)
# NDRI covariate data 
NDRI_Covars <- read.csv("01.Data_Preprocessing/Files/NDRI_Normal_Samples.csv", row.names=1)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# prep. data sets 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# order subjects identically in betas and covars files 
all(colnames(MethOxy_2[, , 3])==rownames(NDRI_Covars))
all(colnames(MethOxy_2[, , 3][order(colnames(MethOxy_2[, , 3]))])==rownames(NDRI_Covars)[order(rownames(NDRI_Covars))])
betas <- MethOxy_2[, , 3][,order(colnames(MethOxy_2[, , 3]))]
covars <- NDRI_Covars[order(rownames(NDRI_Covars)),]
all(colnames(betas)==rownames(covars))

# add total 5hmC content to covars
total_hmC <- readRDS("02.Characterization_5hmC_levels/Files/subject_total_5hmC_content.rds")
total_hmC_2 <- total_hmC[order(names(total_hmC))]
all(names(total_hmC_2)==rownames(covars))
covars$total_5hmC <- total_hmC_2

# index betas to high 5hmC CGs 
betas_hmC <- betas[match(high_5hmC$id, rownames(betas)),]

# load heatmap function 
source("02.Characterization_5hmC_levels/Scripts/heatmap3.R")

# index annot for high 5hmC CGs
annot2 <- annot[match(high_5hmC$id, annot$TargetID),]

# modify levels
levels(annot2$RelToIslandUCSC)[levels(annot2$RelToIslandUCSC)==""] <- "Open Sea"
# Collapse "North" and "South" nomenclature for CpG islands
CpGs_Annotated <- gsub("^[N|S]_","",annot2$RelToIslandUCSC) 

# generate hsv vectors of colors for each variable based on min, median and max values
## Age 
colorAge <- hsv(0.01, covars$Age/max(covars$Age), 1)
## BMI 
colorBMI <- hsv(0.1, covars$BMI/max(covars$BMI), 1)
## cell type 
cell_props_2 <- cell_props[order(cell_props$Sample_Name),]
covars$cell_type_1 <- cell_props_2$X1
covars$cell_type_2 <- cell_props_2$X2
colorcelltype1 <- hsv(0.7, covars$cell_type_1/max(covars$cell_type_1), 1)
colorcelltype2 <- hsv(0.72, covars$cell_type_2/max(covars$cell_type_2), 1)
## total 5hmC content 
color5hmC <- hsv(0.8, covars$total_5hmC/max(covars$total_5hmC), 1)

# bind colors together into one table 
c.lab <- cbind(color5hmC, colorcelltype2, colorcelltype1, colorBMI, colorAge)
colnames(c.lab) <- c("Total 5hmC", "Cell type 2", "Cell type 1", "BMI", "Age")

# generate vector of colors for CpG island regions 
## Column colors 
colregions <- c("black", "blue", "ghostwhite", "darkgray")
## Get CpG Island Colors
Regions <- c("1","2","3","4")
colorRegionsVector <- c(as.character(as.numeric(as.factor(CpGs_Annotated))))
for(i in 1:length(colregions))
{
  colorRegionsVector[colorRegionsVector == Regions[i]] <- colregions[i]
}
r.lab <- matrix(NA, nrow = 1, ncol = dim(betas_hmC)[1])
r.lab[1,] <- colorRegionsVector
colnames(r.lab) <- rownames(betas_hmC)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# implement clustering analysis 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Define the color palette
heatcols <- colorRampPalette(c("yellow","black",'blue'))(n = 1000)
png("02.Characterization_5hmC_levels/Figures/high_5hmC_Related_CpGs_Heatmap_2.png", height=8*300, width=6*300, res=300)
heatmap.3(betas_hmC, col=heatcols, dendrogram='none', 
          RowSideColors = r.lab, ColSideColors = c.lab, ColSideColorsSize=3.5, RowSideColorsSize=1.2, 
          main = "", margins=c(5,8), 
          xlab="", labRow=FALSE, labCol=FALSE, trace = "none", KeyValueName="5hmC")
par(lend = 1) 
#legend(0.8, 1, legend = c("CpG Island", "Shore", "Shelf", "Open Sea"),
       #col = c("black", "darkgray", "ghostwhite",  "blue"), 
       #fill = c("black", "darkgray", "ghostwhite",  "blue"), 
       #bty = "n")
legend(1e-10, 0.65, legend = c("Max", "Median", "Min"),
       col = hsv(0.01, c(max(covars$Age), median(covars$Age), min(covars$Age))/max(covars$Age), 1), 
       fill = hsv(0.01, c(max(covars$Age), median(covars$Age), min(covars$Age))/max(covars$Age), 1), 
       bty = "n", title = "Age", cex = 0.7)
legend(1e-10, 0.52, legend = c("Max", "Median", "Min"),
       col = hsv(0.1, c(max(covars$Age), median(covars$Age), min(covars$Age))/max(covars$Age), 1), 
       fill = hsv(0.1, c(max(covars$Age), median(covars$Age), min(covars$Age))/max(covars$Age), 1), 
       bty = "n", title = "BMI", cex = 0.7)
legend(1e-10, 0.39, legend = c("Max", "Median", "Min"),
       col = hsv(0.7, c(max(covars$cell_type_1), median(covars$cell_type_1), min(covars$cell_type_1))/max(covars$cell_type_1), 1), 
       fill = hsv(0.7, c(max(covars$cell_type_1), median(covars$cell_type_1), min(covars$cell_type_1))/max(covars$cell_type_1), 1), 
       bty = "n", title = "Cell type 1", cex = 0.7)
legend(1e-10, 0.26, legend = c("Max", "Median", "Min"),
       col = hsv(0.72, c(max(covars$cell_type_2), median(covars$cell_type_2), min(covars$cell_type_2))/max(covars$cell_type_2), 1), 
       fill = hsv(0.72, c(max(covars$cell_type_2), median(covars$cell_type_2), min(covars$cell_type_2))/max(covars$cell_type_2), 1), 
       bty = "n", title = "Cell type 2", cex = 0.7)
legend(1e-10, 0.13, legend = c("Max", "Median", "Min"),
       col = hsv(0.8, c(max(covars$Age), median(covars$Age), min(covars$Age))/max(covars$Age), 1), 
       fill = hsv(0.8, c(max(covars$Age), median(covars$Age), min(covars$Age))/max(covars$Age), 1), 
       bty = "n", title = "Total 5hmC", cex = 0.7)
legend("topright",      
       legend = c("CpG Island", "Shore", "Shelf", "Open Sea"),
       col = c("black", "darkgray", "ghostwhite",  "blue"), 
       bty ='n',
       lwd = 5, cex = 0.7)
dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# implement PCA analysis 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# group subjects by quartiles of variables of interest 
## age
quartiles_age <- quantile(covars$Age, probs = c(0.25, 0.5, 0.75, 1))
covars$age_bin <- NA
covars$age_bin[which(covars$Age<=quartiles_age[[1]])] <- 1
covars$age_bin[which(covars$Age>quartiles_age[[1]] & covars$Age<=quartiles_age[[2]])] <- 2
covars$age_bin[which(covars$Age>quartiles_age[[2]] & covars$Age<=quartiles_age[[3]])] <- 3
covars$age_bin[which(covars$Age>quartiles_age[[3]] & covars$Age<=quartiles_age[[4]])] <- 4
covars$age_bin
##BMI
quartiles_bmi <- quantile(covars$BMI, probs = c(0.25, 0.5, 0.75, 1))
covars$bmi_bin <- NA
covars$bmi_bin[which(covars$BMI<=quartiles_bmi[[1]])] <- 1
covars$bmi_bin[which(covars$BMI>quartiles_bmi[[1]] & covars$BMI<=quartiles_bmi[[2]])] <- 2
covars$bmi_bin[which(covars$BMI>quartiles_bmi[[2]] & covars$BMI<=quartiles_bmi[[3]])] <- 3
covars$bmi_bin[which(covars$BMI>quartiles_bmi[[3]] & covars$BMI<=quartiles_bmi[[4]])] <- 4
covars$bmi_bin
##cell type 1
quartiles_ct1 <- quantile(covars$cell_type_1, probs = c(0.25, 0.5, 0.75, 1))
covars$ct1_bin <- NA
covars$ct1_bin[which(covars$cell_type_1<=quartiles_ct1[[1]])] <- 1
covars$ct1_bin[which(covars$cell_type_1>quartiles_ct1[[1]] & covars$cell_type_1<=quartiles_ct1[[2]])] <- 2
covars$ct1_bin[which(covars$cell_type_1>quartiles_ct1[[2]] & covars$cell_type_1<=quartiles_ct1[[3]])] <- 3
covars$ct1_bin[which(covars$cell_type_1>quartiles_ct1[[3]] & covars$cell_type_1<=quartiles_ct1[[4]])] <- 4
covars$ct1_bin
##cell type 2
quartiles_ct2 <- quantile(covars$cell_type_2, probs = c(0.25, 0.5, 0.75, 1))
covars$ct2_bin <- NA
covars$ct2_bin[which(covars$cell_type_2<=quartiles_ct2[[1]])] <- 1
covars$ct2_bin[which(covars$cell_type_2>quartiles_ct2[[1]] & covars$cell_type_2<=quartiles_ct2[[2]])] <- 2
covars$ct2_bin[which(covars$cell_type_2>quartiles_ct2[[2]] & covars$cell_type_2<=quartiles_ct2[[3]])] <- 3
covars$ct2_bin[which(covars$cell_type_2>quartiles_ct2[[3]] & covars$cell_type_2<=quartiles_ct2[[4]])] <- 4
covars$ct2_bin
##total 5hmC content
quartiles_hmC <- quantile(covars$total_5hmC, probs = c(0.25, 0.5, 0.75, 1))
covars$hmC_bin <- NA
covars$hmC_bin[which(covars$total_5hmC<=quartiles_hmC[[1]])] <- 1
covars$hmC_bin[which(covars$total_5hmC>quartiles_hmC[[1]] & covars$total_5hmC<=quartiles_hmC[[2]])] <- 2
covars$hmC_bin[which(covars$total_5hmC>quartiles_hmC[[2]] & covars$total_5hmC<=quartiles_hmC[[3]])] <- 3
covars$hmC_bin[which(covars$total_5hmC>quartiles_hmC[[3]] & covars$total_5hmC<=quartiles_hmC[[4]])] <- 4
covars$hmC_bin

# perform PCA on high 5hmC CG betas 
PCA <- prcomp(scale(betas_hmC))

# plot PC1 vs PC2 w/ color coding for variables, and visualize proportion of variance explained by PCs 
png("02.Characterization_5hmC_levels/Figures/PCA_high_5hmC_loci_PC1-vs-PC2.png", height=8*300, width=12.5*300, res=300)
par(mfrow=c(2,3))
plot(PCA$x[,1:2], col = covars$age_bin, main = "Age", las = 1) ; abline(h = 0, lty = 2, lwd = 1.5) ; abline(v = 0, lty = 2, lwd = 1.5)
legend("topright", legend = c("<= 38 yrs", ">38 yrs - <= 54.5 yrs", ">54.50 yrs - <= 65.75 yrs", ">65.75 yrs - <= 80 yrs"), col = c(1, 2, 3, 4), bty ='n', pch = 1, cex =1)
plot(PCA$x[,1:2], col = covars$bmi_bin, main = "BMI", las = 1); abline(h = 0, lty = 2, lwd = 1.5) ; abline(v = 0, lty = 2, lwd = 1.5)
legend("topright", legend = c("<= 23.01", ">23.01 - <= 24.86", ">24.86 - <= 31.64", ">31.64 - <= 62.73"), col = c(1, 2, 3, 4), bty ='n', pch = 1, cex =1)
plot(PCA$x[,1:2], col = covars$ct1_bin, main = "Cell type 1", las = 1); abline(h = 0, lty = 2, lwd = 1.5) ; abline(v = 0, lty = 2, lwd = 1.5)
legend("topright", legend = c("<= 0.09", ">0.09 - <= 0.25", ">0.25 - <= 0.39", ">0.39 - <= 0.75"), col = c(1, 2, 3, 4), bty ='n', pch = 1, cex =1)
plot(PCA$x[,1:2], col = covars$ct2_bin, main = "Cell type 2", las = 1); abline(h = 0, lty = 2, lwd = 1.5) ; abline(v = 0, lty = 2, lwd = 1.5)
legend("topright", legend = c("<= 0.57", ">0.57 - <= 0.73", ">0.73 - <= 0.91", ">0.91 - <= 0.98"), col = c(1, 2, 3, 4), bty ='n', pch = 1, cex =1)
plot(PCA$x[,1:2], col = covars$hmC_bin, main = "Total 5hmC", las = 1); abline(h = 0, lty = 2, lwd = 1.5) ; abline(v = 0, lty = 2, lwd = 1.5)
legend("topright", legend = c("Q1", "Q2", "Q3", "Q4"), col = c(1, 2, 3, 4), bty ='n', pch = 1, cex =1)
barplot(summary(PCA)$importance[2,], las = 1, ylab = "Proportion of Variance explained")
dev.off()

