# Bioinformatics Omics Analysis (R)

A portfolio of bioinformatics analyses demonstrating applied **proteomics and RNA-seq workflows** using R and Bioconductor. Work in this repository was completed during my Ph.D. research. This repository is intended to showcase skills relevant to **industry and academic bioinformatics and data scientist roles**.

---

## Project Overview

This repository contains reproducible R scripts for analyzing proteomics and transcriptomics datasets, including:
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

| Script | Description |
|------|------------|
| `J_Proteomics/J_complete_proteomics_pipeline_human_sy5y_cells_ctrl_vs_ogt_kd.R` | Differential total proteomics analysis in O-GlcNAc transferase gene (OGT) knockdown Human neuroblastoma SY5Y cell lines, following serum re-activation for kinase activation  |
| `J_Proteomics/J_complete_proteomics_pipeline_mouse_brain_ctrl_vs_ogt_ko.R` | Differential proteomics analysis in OGT knockout mouse brain tissue |
| `J_rnaseq/J_complete_rnaseq_pipeline_mouse_liver_saline_vs_tmg.R` | RNA-seq differential expression workflow for OGT knockout mouse liver tissue |
| `J_utils/J_analysis_utils.R` | Utility functions supporting analysis pipelines |

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
Work in this repository was completed during my Ph.D. research  
GitHub: [https://github.com/sjohnephrame](https://github.com/sjohnephrame)
