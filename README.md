# 5-hydroxymethylcytosine (5hmC) in normal breast tissue

## *Genome-wide abundance of 5hmC in breast tissue reveals unique function in dynamic gene regulation & carcinogenesis*

Owen M. Wilkins, Kevin C. Johnson, E. Andres Houseman, Jessica E. King, Carmen J. Marsit, Brock C. Christensen <br />
May 2018

BioRxiv manuscript: https://www.biorxiv.org/content/early/2018/06/07/339069

## Project description

Recent studies have demonstrated 5hmC to be a stable epigenetic mark that contributes to a range of cellular processes, ranging from embryonic development to somatic cell reprogramming. However, genome-scale maps of 5hmC in human tissues are lacking. We utilized a paired bisulfite & oxidative-bisulfite DNA treatment approach followed by DNA methylation profiling (450K array) to study the distribution and abundance of 5hmC in normal breast tissue.

Our work reveals a striking enrichment of 5hmC among breast-specific DNA regulatory elements positively associated with transcriptional activity, while being simultaneously depleted among regulatory regions involved in transcriptional repression. In addition, we demonstrate that 5hmC is positively correlated with gene expression in normal breast tissue. Finally, we find 5hmC in enriched among CpG loci differentially methylated between premalignant/malignant breast tissue and adjacent-normal, suggested 5hmC may contribute to carcinogenesis.

This repository contains all code used to process and generate our results. An overview of the analytical steps and data availability are provided below. This work is currently under submission for publication and this repository will be updated when the associated manuscript is published. 

## Analysis overview

* **01.Data_Preprocessing**
    1. Processing of raw DNA methylation intensity data (IDATs)
    2. Estimation of 5hmC/5mC proportions in normal breast tissue

* **02.Characterization_5hmC_levels**
    1. Genomic characterization of 5hmC distribution & abundance
    2. Selection of *high 5hmC CpGs* (CpG loci among top 1% in 5hmC abundance)

* **03.Enrichment_Analyses**
    1. Gene set enrichment of *high 5hmC CpGs*
    2. Tests of *high 5hmC CpG* enrichment among genomic features (promoters, exons, etc.)
    3. Tests of *high 5hmC CpG* enrichment among DNA regulatory regions (e.g. histone modifications)

* **04.Gene_expression**
    1. Gene-level enrichment analysis of 5hmc among highly expressed genes in breast tissues (GTEx)
    2. CpG-specific correlation analysis between 5hmC abundance & gene expression

* **05.Cancer_Comparison**
    1. Tests of *high 5hmC CpG* enrichment among differentially methylated CpGs between DCIS & adjacent-normal
    2. Tests of *high 5hmC CpG* enrichment among differentially methylated CpGs between tumor & adjacent-normal

## Data availability

* Normal breast (*n*=18) tissue DNA methylation data (IDATs, 450K, bisulfite & oxidative-bisulfite treated DNA):
    - NDRI, GEO accession number: GSE73895
* DCIS (*n*=40) & adjacent-normal (*n*=15) tissue DNA methylation data (IDATs, 450K, bisulfite-treated DNA): GEO,
    - NHMN, GEO accession number: GSE66313
* Breast tumor & adjacent-normal DNA methylation data:
    - TCGA project, NCI Genomic Data Commons (https://portal.gdc.cancer.gov/)
* Breast-specific DNase hypersensitivity site & histone modification coordinates (from ChIP-seq experiments):
    - NIH Roadmap Epigenomics Project (https://www.ncbi.nlm.nih.gov/geo/roadmap/epigenomics/)
* Gene expression (RNA-seq) in normal breast tissues:
    - GTEx Project (https://gtexportal.org/home/)

> **Abbreviations**:  <br />
> NDRI - National Disease Research Interchange <br />
> NHMN - New Hampshire Mammography Network <br />
> TCGA - The Cancer Genome Atlas <br />
> NCI - National Cancer Institute <br />
> GTEx - Genotype-Tissue Expression project <br />
