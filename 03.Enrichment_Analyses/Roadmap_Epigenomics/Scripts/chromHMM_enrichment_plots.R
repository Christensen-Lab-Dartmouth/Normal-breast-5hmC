#######################################################################################################################
# Generate forest plots w/ tables of 5hmC enrichment tests against ChIP-seq data from Roadmap Epigenomics project 
#######################################################################################################################
rm(list = ls())
library(ggplot2)
setwd("/Users/Owen 1/Dropbox (Christensen Lab)/NDRI_Breast_5hmC_update/")
dir2 <- "03.Enrichment_Analyses/Roadmap_Epigenomics/Figures/"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load Data 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

E027 <- read.csv("03.Enrichment_Analyses/Roadmap_Epigenomics/Files/E027_ChromHMM_5hmC_enrichment_results.csv", stringsAsFactors = F)
E119 <- read.csv("03.Enrichment_Analyses/Roadmap_Epigenomics/Files/E119_ChromHMM_5hmC_enrichment_results.csv", stringsAsFactors = F)
E119_vs_E028_lost <- read.csv("03.Enrichment_Analyses/Roadmap_Epigenomics/Files/HMEC_vs_vHMEC_lost_states_ChromHMM_5hmC_enrichment_results.csv", stringsAsFactors = F)
E119_vs_E028_gained <- read.csv("03.Enrichment_Analyses/Roadmap_Epigenomics/Files/HMEC_vs_vHMEC_gained_states_ChromHMM_5hmC_enrichment_results.csv", stringsAsFactors = F)
E119_vs_E028_joint <- read.csv("03.Enrichment_Analyses/Roadmap_Epigenomics/Files/HMEC_vs_vHMEC_all_differential_states_ChromHMM_5hmC_enrichment_results.csv", stringsAsFactors = F)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# pre-process data so rows appear in plots in the desired order
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pre_process_1 <- function(x){
  #### basic genomic features (promoters, introns, etc.)
  # rename columns 
  colnames(x) <- c("state", "P-value", "Odds ratio (95% CI)", "OR", "LB", "UB")
  # reverse order of features then refactor so that they appear in desired order 
  x$state <- rev(x$state)
  x$state <- factor(x$state, levels = x$state)
  x$`P-value` <- rev(x$`P-value`)
  x$`Odds ratio (95% CI)` <- rev(x$`Odds ratio (95% CI)`)
  x$OR <- rev(x$OR)
  x$LB <- rev(x$LB)
  x$UB <- rev(x$UB)
  x
}

E027 <- pre_process_1(E027)
E119 <- pre_process_1(E119)
E119_vs_E028_lost <- pre_process_1(E119_vs_E028_lost)
E119_vs_E028_gained <- pre_process_1(E119_vs_E028_gained)
E119_vs_E028_joint <- pre_process_1(E119_vs_E028_joint)

process_names <- function(x){
  x$state <- factor(sapply(as.character(x$state), function(x1) strsplit(x1, "_")[[1]][2]), 
                    levels = rev(c("TssA", "TssAFlnk", "TxFlnk", "Tx", "TxWk", "EnhG", "Enh", 
                               "ZNF/Rpts", "Het", "TssBiv", "BivFlnk", "EnhBiv", "ReprPC", 
                               "ReprPCWk", "Quies")))
  x
}
E027 <- process_names(E027)
E119 <- process_names(E119)
E119_vs_E028_lost <- process_names(E119_vs_E028_lost)
E119_vs_E028_gained <- process_names(E119_vs_E028_gained)
E119_vs_E028_joint <- process_names(E119_vs_E028_joint)


# convert odds ratios below 1 to -ve values 
#process_OR <- function(x){
#  for(i in 1:nrow(x)){
#    if(x$OR[i]<1){
#      x$OR[i] <- -(1/x$OR[i])
#      x$LB[i] <- -(1/x$LB[i])
#      x$UB[i] <- -(1/x$UB[i])
#    }
#  }
#  x
#}
#E027 <- process_OR(E027)
#E119 <- process_OR(E119)
#E119_vs_E028_lost <- process_OR(E119_vs_E028_lost)
#E119_vs_E028_gained <- process_OR(E119_vs_E028_gained)
#E119_vs_E028_joint <- process_OR(E119_vs_E028_joint)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# generate plots
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

f_plot <- function(dat, title, ylim_min, ylim_max, out_dir, file_name){
  ppi=300
  png(paste0(out_dir, file_name), width=4.2*ppi, height=4.4*ppi, res=ppi)
  p <- ggplot(dat, aes(y = OR, x = state)) +
    geom_point(shape = 16, size = 2) +
    ylim(ylim_min, ylim_max) +
    geom_errorbar(aes(ymin = LB, ymax = UB), width=.1) +
    geom_hline(yintercept = 1, linetype = 2, size = 0.4) +
    coord_flip() +
    labs(title = paste0(title)) +
    ylab("Odds Ratio") +
    xlab("") +
    theme_classic() +
    theme(strip.background = element_blank(), legend.position="none", 
          axis.text.y = element_text(color = "black", size = 11.5), 
          axis.text.x = element_text(color = "black", size = 11.5), 
          panel.background = element_rect(colour = "black", size=0.7), 
          line = element_blank())
  print(p)
  dev.off()
}
# apply to each data set 
f_plot(E027, "Breast myoepithelial cells", 0, 5, dir2, "E027_5hmC_enrichment_forest_plot.png")
f_plot(E119, "HMEC", 0, 5, dir2, "E119_5hmC_enrichment_forest_plot.png")
f_plot(E119_vs_E028_lost, "HMEC vs vHMEC - lost states", 0, 11.5, dir2, "HMEC_vs_vHMEC_5hmC_lost_chromatin_states_enrichment_forest_plot.png")
f_plot(E119_vs_E028_gained, "HMEC vs vHMEC - gained states", 0, 11.5, dir2, "HMEC_vs_vHMEC_5hmC_gained_chromatin_states_enrichment_forest_plot.png")
f_plot(E119_vs_E028_joint, "HMEC vs vHMEC - all states", 0, 5, dir2, "HMEC_vs_vHMEC_5hmC_all_chromatin_states_enrichment_forest_plot.png")



