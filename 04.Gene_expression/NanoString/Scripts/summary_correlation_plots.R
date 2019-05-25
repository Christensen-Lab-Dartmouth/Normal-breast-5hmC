#######################################################################################################################
# Generate plots of correlations between 5hmC/5mC and breast specific gene expression 
#######################################################################################################################
rm(list = ls())
setwd("/Users/Owen 1/Dropbox (Christensen Lab)/NDRI_Breast_5hmC_update/")
library(ggplot2)

# read in data 
hmC <- read.csv("04.Nano_String/Files/high_5hmC_CpG_expression_correlations_5hmC.csv", stringsAsFactors = F)
mC <- read.csv("04.Nano_String/Files/high_5mC_CpG_expression_correlations_5mC.csv", stringsAsFactors = F)
str(hmC)
str(mC)

# drop the median calculated expression entries for RASSF1 and DNMT3A to leave the values for each transcript, 
hmC <- hmC[-which(hmC$Transcript_Variant=="RASSF1" | hmC$Transcript_Variant=="DNMT3A"),]
mC <- mC[-which(mC$Transcript_Variant=="RASSF1" | mC$Transcript_Variant=="DNMT3A"),]
head(hmC)
# convert transcripts variable to factors so you can specify order they appear in plot 
hmC$Transcript_Variant <- factor(hmC$Transcript_Variant, levels = c("RAB32", "TWIST1", "RASSF1_vB", "RASSF1_vC", "RASSF1_vH", "DNMT3A_v3", "DNMT3A_v4"))
head(hmC)
mC$Transcript_Variant <- factor(mC$Transcript_Variant, levels = c("RAB32", "TWIST1", "RASSF1_vB", "RASSF1_vC", "RASSF1_vH", "DNMT3A_v3", "DNMT3A_v4"))
head(hmC)

# plot correlations between high 5hmC CpGs and gene expression 
dpi = 300
png("04.Nano_String/Figures/high_5hmC_CpG_expression_correlations_5hmC.png", width = dpi*10, height = dpi*10, res = dpi)
plot_5hmc <- ggplot(hmC, aes(x = Transcript_Variant, y = Spearman_Cor, size = -log10(Spearman_Pval), fill = Transcript_Variant)) +  
  geom_point(alpha = 0.5, shape = 21) + 
  expand_limits(y=c(-0.8,0.8)) + 
  scale_size(range = c(1, 20)) +
  geom_hline(yintercept = 0) + 
  xlab("") + 
  ylab("Correlation coefficient (5hmC & expression)") + 
  labs(size = expression("0.001")) + 
  #labs(size = expression("1e-03")) +
  guides(size=guide_legend(override.aes = list(fill="black", alpha=1)))+
  theme(legend.key.size = unit(1, "cm"), 
       # panel.grid.major = element_line(colour = "#d3d3d3"), 
        panel.border = element_rect(fill = NA, colour = "black", size = 0.6, linetype = "solid"), 
        panel.background = element_blank(),
        axis.text.x=element_text(angle = 45, colour="black", size = 25, hjust = 1), 
        axis.text.y=element_text(colour="black", size = 25), 
        axis.title.x=element_text(colour="black", size = 25), 
        axis.title.y=element_text(colour="black", size = 25), 
        legend.key = element_blank(), 
        legend.text = element_text(colour="black", size = 20), 
        legend.title = element_text(colour="black", size = 20))
plot_5hmc 
dev.off()

# same for 5mC at same CpGs 
png("04.Nano_String/Figures/high_5hmC_CpG_expression_correlations_5mC.png", width = dpi*10, height = dpi*10, res = dpi)
plot_5mc <- ggplot(mC, aes(x = Transcript_Variant, y = Spearman_Cor, size = -log10(Spearman_Pval), fill = Transcript_Variant))+ 
  geom_point(alpha = 0.5, shape = 21) + 
  expand_limits(y=c(-0.8,0.8)) + 
  scale_size(range = c(1, 20)) +
  geom_hline(yintercept = 0) + 
  xlab("") + 
  ylab("Correlation coefficient (5mC & expression)") + 
  labs(size = expression("P-value")) + 
  guides(size=guide_legend(override.aes = list(fill="black", alpha=1)))+
  theme(legend.key.size = unit(1, "cm"), 
        #panel.grid.major = element_line(colour = "#d3d3d3"), 
        panel.border = element_rect(fill = NA, colour = "black", size = 0.6, linetype = "solid"), 
        panel.background = element_blank(),
        axis.text.x=element_text(angle = 45, colour="black", size = 25, hjust = 1), 
        axis.text.y=element_text(colour="black", size = 25), 
        axis.title.x=element_text(colour="black", size = 25), 
        axis.title.y=element_text(colour="black", size = 25), 
        legend.key = element_blank(), 
        legend.text = element_text(colour="black", size = 20), 
        legend.title = element_text(colour="black", size = 20))
plot_5mc
dev.off()

