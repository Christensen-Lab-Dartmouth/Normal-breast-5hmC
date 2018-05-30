#######################################################################################################################
# Locus Overlap Analysis (LOLA) to identify genomic features enriched for high 5hmC CpGs
# Results visualization
#######################################################################################################################
rm(list = ls())
library(LOLA)
library(ggplot2)
setwd("/Users/Owen 1/Dropbox (Christensen Lab)/NDRI_Breast_5hmC_update/")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load pre-computed results 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load("04.Enrichment_Analyses/LOLA/Files/lola_results.Rdata")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# pre-process & explore results 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# look through results by collection 
cistrome_epi <- locResults[locResults$collection=="cistrome_epigenome", ]
encode <- locResults[locResults$collection=="encode_tfbs", ]
cistrome <- locResults[locResults$collection=="cistrome_cistrome", ]
shef <- locResults[locResults$collection=="sheffield_dnase", ]
codex <- locResults[locResults$collection=="codex", ]
ucsc <- locResults[locResults$collection=="ucsc_features", ]
rm(locResults)

# order by mean Rnk
cistrome_epi <- cistrome_epi[order(cistrome_epi$meanRnk, decreasing=FALSE),]
encode <- encode[order(encode$meanRnk, decreasing=FALSE),]
cistrome <- cistrome[order(cistrome$meanRnk, decreasing=FALSE),]
shef <- shef[order(shef$meanRnk, decreasing=FALSE),]

# subset to those w/ qvalue <0.05
cistrome_epi <- cistrome_epi[cistrome_epi$qValue < 0.05,]
encode <- encode[encode$qValue < 0.05,]
cistrome <- cistrome[cistrome$qValue < 0.05,]
shef <- shef[shef$qValue < 0.05,]

##### index for only tissues derived from breast  
# cistrome_epi
table(cistrome_epi$cellType)
cistrome_epi <- cistrome_epi[cistrome_epi$cellType == "MCF-7" | cistrome_epi$cellType == "T47D",]
# encode
table(encode$cellType)
encode <- encode[encode$cellType == "MCF-7" | encode$cellType == "MCF10A-Er-Src",]
# cistrome
table(cistrome$cellType)
cistrome <- cistrome[cistrome$cellType == "MCF-7" | cistrome$cellType == "MDA-MB-231 Cells" | cistrome$cellType == "T47D",]
cistrome$cellType[which(cistrome$cellType=="MDA-MB-231 Cells")] <- "MDA-MB-231"

# quick exploration of summary statistics 
par(mfrow=c(2,2))
plot(-log10(cistrome_epi$qValue), bty = "n")
plot(-log10(encode$qValue), bty = "n")
plot(-log10(cistrome$qValue), bty = "n")
dev.off()

par(mfrow=c(2,2))
plot(cistrome_epi$logOddsRatio, bty = "n")
plot(encode$logOddsRatio, bty = "n")
plot(cistrome$logOddsRatio, bty = "n")
dev.off()

rm(shef, codex, ucsc)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot & write out results from each collection of interest
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# plot results for TFBS enrichment using data from cistrome database 
ppi=300
png("04.Enrichment_Analyses/LOLA/Figures/lola_cistrome_TFBS.png", width=9*ppi, height=8*ppi, res=ppi)
fill = c("chartreuse3", "lightskyblue1", "sienna1")
p1 <- ggplot(cistrome, aes(x = antibody, y = -log10(qValue), size = logOddsRatio, fill = cellType))+
  geom_point(alpha = 0.5, shape = 21) +
  guides(size=guide_legend(override.aes = list(fill="black", alpha=1)))+
  labs(x = "Transcription factor binding sites (cistrome database)", y = expression("-log"[10]*"(Q-value)"),
       size = "logOdds ratio", fill = "Cell type") +
  scale_size(range = c(1, 12)) +
  scale_fill_manual(values = fill) + 
  theme(legend.key.size = unit(1, "cm"),
        panel.grid.major = element_line(colour = "#d3d3d3"), 
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "black", size = 0.6, linetype = "solid"), 
        panel.background = element_blank(),
        axis.text.x=element_text(angle = 45, colour="black", size = 20, hjust = 1), 
        axis.text.y=element_text(colour="black", size = 20), 
        axis.title.x=element_text(colour="black", size = 20), 
        axis.title.y=element_text(colour="black", size = 20), 
        legend.key = element_blank(), 
        legend.text = element_text(colour="black", size = 20), 
        legend.title = element_text(colour="black", size = 20))
p1
dev.off()

# write csv files for results tables from each collection 
write.csv(cistrome, file = "04.Enrichment_Analyses/LOLA/Files/lola_results_high_5hmC_cistrome.csv")
write.csv(encode, file = "04.Enrichment_Analyses/LOLA/Files/lola_results_high_5hmC_encode.csv")
write.csv(cistrome_epi, file = "04.Enrichment_Analyses/LOLA/Files/lola_results_high_5hmC_cistrome_epi.csv")



