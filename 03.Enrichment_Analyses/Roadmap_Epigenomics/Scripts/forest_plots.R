#######################################################################################################################
# Generate forest plots w/ tables of 5hmC enrichment tests against ChIP-seq data from Roadmap Epigenomics project 
#######################################################################################################################
rm(list = ls())
setwd("/Users/Owen 1/Dropbox (Christensen Lab)/NDRI_Breast_5hmC_update/")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load Data 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

HMEC <- read.csv("04.Enrichment_Analyses/Roadmap_Epigenomics/Files/roadmap_HMEC_5hmC_enrichment_results.csv", stringsAsFactors = F)
myo <- read.csv("04.Enrichment_Analyses/Roadmap_Epigenomics/Files/roadmap_breast_myoepith_5hmC_enrichment_results.csv", stringsAsFactors = F)
gen <- read.csv("04.Enrichment_Analyses/Genomic_context//Files/transcription_feature_hg19_enrichment.csv")
vHMEC <- read.csv("04.Enrichment_Analyses/Roadmap_Epigenomics/Files/roadmap_vHMEC_5hmC_enrichment_results.csv", stringsAsFactors = F)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# pre-process data so rows appear in plots in the desired order
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#### basic genomic features (promoters, introns, etc.)
# rename columns 
colnames(gen) <- c("Genomic feature", "P-value", "Odds ratio (95% CI)", "OR", "LB", "UB")
# reverse order of features then refactor so that they appear in desired order 
gen$`Genomic feature` <- rev(gen$`Genomic feature`)
gen$`Genomic feature` <- factor(gen$`Genomic feature`, levels = gen$`Genomic feature`)
gen$`P-value` <- rev(gen$`P-value`)
gen$`Odds ratio (95% CI)` <- rev(gen$`Odds ratio (95% CI)`)
gen$OR <- rev(gen$OR)
gen$LB <- rev(gen$LB)
gen$UB <- rev(gen$UB)

#### breast myoepithelial cells 
# rename columns 
colnames(myo) <- c("Genomic feature", "P-value", "Odds ratio (95% CI)", "OR", "LB", "UB")
# reverse order of features then refactor so that they appear in desired order 
myo$`Genomic feature` <- rev(myo$`Genomic feature`)
myo$`Genomic feature` <- factor(myo$`Genomic feature`, levels = myo$`Genomic feature`)
myo$`P-value` <- rev(myo$`P-value`)
myo$`Odds ratio (95% CI)` <- rev(myo$`Odds ratio (95% CI)`)
myo$OR <- rev(myo$OR)
myo$LB <- rev(myo$LB)
myo$UB <- rev(myo$UB)

#### HMECs
# rename columns 
colnames(HMEC) <- c("Genomic feature", "P-value", "Odds ratio (95% CI)", "OR", "LB", "UB")
# reverse order of features then refactor so that they appear in desired order 
HMEC$`Genomic feature` <- rev(HMEC$`Genomic feature`)
HMEC$`Genomic feature` <- factor(HMEC$`Genomic feature`, levels = HMEC$`Genomic feature`)
HMEC$`P-value` <- rev(HMEC$`P-value`)
HMEC$`Odds ratio (95% CI)` <- rev(HMEC$`Odds ratio (95% CI)`)
HMEC$OR <- rev(HMEC$OR)
HMEC$LB <- rev(HMEC$LB)
HMEC$UB <- rev(HMEC$UB)

# bind these data sets together so they can be plotted together 
comb <- rbind(HMEC, myo, gen)

#### vHMECs
# rename columns 
colnames(vHMEC) <- c("Genomic feature", "P-value", "Odds ratio (95% CI)", "OR", "LB", "UB")
vHMEC[2,1] <- "H3k4me1"
# reverse order of features then refactor so that they appear in desired order 
vHMEC$`Genomic feature` <- rev(vHMEC$`Genomic feature`)
vHMEC$`Genomic feature` <- factor(vHMEC$`Genomic feature`, levels = vHMEC$`Genomic feature`)
vHMEC$`P-value` <- rev(vHMEC$`P-value`)
vHMEC$`Odds ratio (95% CI)` <- rev(vHMEC$`Odds ratio (95% CI)`)
vHMEC$OR <- rev(vHMEC$OR)
vHMEC$LB <- rev(vHMEC$LB)
vHMEC$UB <- rev(vHMEC$UB)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# generate plots
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# generate forest plot for genomic context (promoters, introns, etc.), myoepthelial cells, HMECs
ppi=300
png("04.Enrichment_Analyses/Roadmap_Epigenomics/Figures/normal_cells_5hmC_enrichment_forest_plot.png", width=5*ppi, height=5*ppi, res=ppi)
p <- ggplot(comb, aes(y = OR, x = `Genomic feature`)) +
  geom_point(shape = 16, size = 2) +
  geom_errorbar(aes(ymin = LB, ymax = UB), width=.1) +
  geom_hline(yintercept = 1, linetype = 2, size = 0.4) +
  coord_flip() +
  labs(title = "") +
  ylab("Odds Ratio") +
  xlab("") +
  theme_classic() +
  theme(strip.background = element_blank(), legend.position="none", 
        axis.text.y = element_text(color = "black", size = 10), 
        axis.text.x = element_text(color = "black", size = 10), 
        panel.background = element_rect(colour = "black", size=0.7), 
        line = element_blank())
p
dev.off()

# generate forest plot for vHMECs only 
ppi=300
png("04.Enrichment_Analyses/Roadmap_Epigenomics/Figures/vHMEC_5hmC_enrichment_forest_plot.png", width=5*ppi, height=2.3*ppi, res=ppi)
p <- ggplot(vHMEC, aes(y = OR, x = `Genomic feature`)) +
  geom_point(shape = 16, size = 2) +
  geom_errorbar(aes(ymin = LB, ymax = UB), width=.1) +
  geom_hline(yintercept = 1, linetype = 2, size = 0.4) +
  coord_flip() +
  labs(title = "") +
  ylab("Odds Ratio") +
  xlab("") +
  theme_classic() +
  theme(strip.background = element_blank(), legend.position="none", 
        axis.text.y = element_text(color = "black", size = 14), 
        axis.text.x = element_text(color = "black", size = 14), 
        panel.background = element_rect(colour = "black", size=0.7), 
        line = element_blank())
p
dev.off()
