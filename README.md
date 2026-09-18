# C4-EAS-Brain-Network | East Asian Brain C4 Regulatory Atlas

> Prenatal brain gene co-expression networks link complement signaling in psychiatric disorders
> Central South University
> Manuscript in submission · 2026

![Pipeline Status](https://img.shields.io/badge/Pipeline-Stable-success)
![License](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey)
![Platform](https://img.shields.io/badge/platform-SLURM_HPC-blue)
![R](https://img.shields.io/badge/R-%E2%89%A54.0.0-blue)

### What this is

This repository contains the analytical pipeline used to elucidate the spatiotemporal regulation of the *C4* locus (*C4A* and *C4B*) across human brain development, with a specific focus on East Asian (EAS) populations. 

The workflow integrates genetic structural variation (CNV) imputation, transcriptomic profiling (bulk and deconvoluted cell-type expression), stage-specific seed network analysis, and disease heritability enrichment to demonstrate how early-established *C4A* genetic risk converges on postnatal synaptic networks in psychiatric disorders like schizophrenia (SCZ) and bipolar disorder (BIP).

### Pipeline overview

The complete analysis is structured into modular directories. The pipeline consists of the following core stages:

| # | Stage | Directory | Tools / Methods |
| :--- | :--- | :--- | :--- |
| 01 | **C4 RNA Expression Profiling** | `01_C4_RNA_Expression` | R, Differential Expression, Sex & Cell-type Specificity Analysis |
| 02 | **C4 Structural Variation & Imputation** | `02_C4_region_CNV` | Osprey, PLINK, bcftools, R (Haplotype Refinement) |
| 03 | **C4 CNV-Expression QTL Analysis** | `03_C4_CNV-exp_QTL` | Bulk eQTL, Cell-type Deconvolution (e.g., Dtangle/Bmind), Cell-type eQTL |
| 04 | **C4A Seed Co-expression Networks** | `04_C4A_network` | Seed-based Co-expression, Pathway Enrichment |
| 05 | **Disease Heritability Enrichment** | `04_C4A_network` | MAGMA, GWAS Summary Statistics Integration |

### Directory Structure Details

* **`01_C4_RNA_Expression/`**: Scripts evaluating the baseline expression of *C4* paralogs across development. Includes analyses for sex-biased expression (`01a`) and cell-type specific enrichment patterns (`01b`).
* **`02_C4_region_CNV/`**: 
    * `genotype_SNP_calling/`: Upstream processing for joint variant calling and pre-imputation formatting.
    * `CNV_imputation/`: Core pipeline utilizing the `Osprey` algorithm and an extended 1000 Genomes (1KG) EAS reference panel to accurately infer complex *C4A* and *C4B* copy numbers and haplotypes.
* **`03_C4_CNV-exp_QTL/`**: Scripts linking the imputed *C4* structural variations to transcriptomic readouts to establish eQTL effects. This includes both bulk tissue analysis (`02`) and high-resolution deconvoluted cell-type specific eQTL mapping (`01`, `03`).
* **`04_C4A_network/`**: Generates stage-specific (prenatal vs. postnatal) *C4A* positive and negative seed co-expression networks. Subsequent scripts map genetic liability from psychiatric GWAS (e.g., SCZ, BIP) onto these developmental trajectories using MAGMA to identify critical windows of vulnerability.

### Prerequisites & Dependencies

* **Command-Line Tools**: `Osprey` (for *C4* complex CNV imputation), `bcftools` (v1.9+), `PLINK 2.0`, `MAGMA`.
* **R Packages**: `DESeq2` / `edgeR`, `clusterProfiler`, `ggplot2`, cell-type deconvolution libraries.

