# Bioinformatics Omics Analysis (R)

A portfolio of bioinformatics analyses demonstrating applied **proteomics and RNA-seq workflows** using R and Bioconductor.  
This repository is intended to showcase skills relevant to **industry and academic bioinformatics and data scientist roles**.

---

## Project Overview

This repository contains reproducible R scripts for analyzing multi-omics datasets, including:
- Differential expression analysis
- Statistical modeling
- Data normalization and QC
- Biological interpretation and visualization

These workflows reflect real research pipelines used in translational and discovery biology settings.

---

## Tools & Technologies

- **Languages:** R
- **Libraries:** tidyverse, DESeq2, limma, edgeR, ggplot2
- **Data Types:** RNA-seq, total proteomics
- **Skills Demonstrated:**
  - Experimental design awareness
  - Statistical testing & multiple-hypothesis correction
  - Data visualization for biological insight
  - Modular and reusable R code

---

## Repository Structure

scripts/
  proteomics/
  rnaseq/
functions/
README.md

| Script | Description |
|------|------------|
| `J_R_OGT KO MICE BRAIN TOTAL PROTEOME ANALYSIS.R` | Differential proteomics analysis in OGT knockout mouse brain tissue |
| `J_R_TOTAL PROTEOME ANALYSIS SY5Y HUMAN.R` | Differential total proteomics analysis in OGT knockdown Human neuroblastoma SY5Y cell lines, following serum re-activation for kinase activation  |
| `J_STAT GENOMICS FINAL RNASEQ.R` | RNA-seq differential expression workflow |
| `functions2_FROM_DR_THOMPSON_GITHUB.R` | Utility functions supporting analysis pipelines |

---

## Example Workflow

Typical analysis steps include:
1. Data import and preprocessing  
2. Quality control and normalization  
3. Differential expression analysis  
4. Statistical validation  
5. Visualization (volcano plots, PCA, heatmaps)  
6. Biological interpretation  

---

## Notes for Recruiters

This repository showcases **end-to-end bioinformatics workflows** for proteomics and RNA-seq data.  
All scripts were developed as part of my **Ph.D. research**, demonstrating skills relevant to **both industry and academic bioinformatics roles**, including:

- Data preprocessing and normalization  
- Statistical analysis  
- Visualization and data interpretation  

Scripts are structured for **reproducibility** and can be applied to new datasets, highlighting the ability to deliver **reliable and interpretable results** in real-world research or industry projects.

---

## Reproducibility

- All scripts are designed to be run in **RStudio**.
- Required R packages can be installed using:

```r
install.packages("tidyverse")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("DESeq2", "limma", "edgeR"))
```
Scripts include clear step-by-step comments to guide users through each workflow.

---

## Author

**Sophiya J Hanigan, PhD**  
Bioinformatics Scientist | Omics Data Analysis  
Work in this repository was completed during Ph.D. research  
GitHub: [https://github.com/sjohnephrame](https://github.com/sjohnephrame)
