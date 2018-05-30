#######################################################################################################################
# Explore cell type proprtion estimates accross subjects 
#######################################################################################################################
rm(list = ls())
setwd("/Users/Owen 1/Dropbox (Christensen Lab)/NDRI_Breast_5hmC_update/")

# read in putative cellular proportions 
x1 <- read.csv("02.Characterization_5hmC_levels/Files/NDRI_Cell_Proportions.csv")

# generate exploratory barplots 
par(mfrow=c(2,1))
barplot(x1$X1)
barplot(x1$X2)
dev.off()

# output individual barplot to png 
ppi = 300
png("02.Characterization_5hmC_levels/Files/NDRI_breast_cell_type_props.png", height=7*ppi, width=5*ppi, res=ppi)
par(mfrow=c(2,1))
barplot(x1$X1, ylim = c(0,1), col = "indianred", main = "Cell type 1", xlab = "Sample", ylab = "proportion", las = 1)
barplot(x1$X2, ylim = c(0,1), col = "indianred", main = "Cell type 2", xlab = "Sample", ylab = "proportion", las = 1)
dev.off()
