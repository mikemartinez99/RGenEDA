# RGenEDA 
**Clean, unified, streamlined, and reproducible frameworks for genomic exploratory data analyses**  


 <!-- badges: start -->
 ![status](https://img.shields.io/badge/status-in--development-orange) 
 ![Version](https://img.shields.io/badge/version-2.0.0-blue)
 ![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)
  [![R-CMD-check](https://github.com/mikemartinez99/RGenEDA/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mikemartinez99/RGenEDA/actions/workflows/R-CMD-check.yaml)
  [![Codecov test coverage](https://codecov.io/gh/mikemartinez99/RGenEDA/graph/badge.svg)](https://app.codecov.io/gh/mikemartinez99/RGenEDA)
  <!-- badges: end -->


### Development <!-- omit in toc -->

> [!IMPORTANT]
> **RenEDA is currently under active development. If there is a feature you would like implemented, please submit a feature request on the Issues page of submit a pull request. For information on new additions see the** **[Change Log](#change-log)**


> [!WARNING]
> **RGenEDA has updated to v2.0.0 which introduced some breaking changes. Please ensure you are using the most up-to-date version**

--- 


![Alt text](img/examples.png)

## Installation
Install the latest version from Github

```r
install.packages("devtools") # if not already instaled

library(devtools)
devtools::install_github("mikemartinez99/RGenEDA")

library(RGenEDA)

```

## Usage
For a full demo, see [Snail1 KO Dataset Demo](https://github.com/mikemartinez99/RGenEDA/blob/main/vignettes/Snail1_Vignette.md)  
[Link to data](https://www.ebi.ac.uk/gxa/experiments/E-MTAB-5244/Downloads); data described in [Matsuri et al., 2018](https://pubmed.ncbi.nlm.nih.gov/29729076/)


## Change-log

**![Version](https://img.shields.io/badge/version-0.0.0.900-blue)** **![Release](https://img.shields.io/badge/Initial-grey)**  
*April 2025: "The poor-man's EDA"*
  - Implemented core, standalone functions
    - Variance calculation, PCA calculation, and Eucliden distances, and Eigenvector correlations

**![Version](https://img.shields.io/badge/version-1.0.0-blue)** **![Release](https://img.shields.io/badge/Major-red)**  
*October 16th, 2025: "The GenEDA Object Update*
  - Introduced S4 **GenEDA** object for streamlined data handling
  - Added **methods** for HVG calculation and PCA storage within S4 slots.
  - Added **visualization functions**
    - PCA plots
    - Variance plots for HVGs
    - Eigenvector (gene-loading) plots
  - Revamped eigencorrelation plots to only use top number of HVGs identified in HVG selection/PCA
  - Enhanced plot output flexibility (return **ggplot2, pheatmap**, or **ComplexHeatmap** objects for full customization
  - Added **unit tests**
  - Added **Pasilla dataset vignette**

**![Version](https://img.shields.io/badge/version-1.0.1-blue)** **![Release](https://img.shields.io/badge/Patch-green)**    
*November 3rd, 2025*  
  - Rely soley on **Pheatmap** for all heatmap visualizations
  - Fixed top-gene-loading bug for selected PCs in PlotEigenHeatmap function
  - Added optional **DEGs slot** to GenEDA object to interface with **DESeq2**
  - Added **MA plot** functionality

**![Version](https://img.shields.io/badge/version-1.0.2-blue)** **![Release](https://img.shields.io/badge/Patch-green)**      
*November 10th, 2025*  
  - Fixed example blocks in roxygen headers
  - Addressed all R cmd build errors, warnings, and notes

**![Version](https://img.shields.io/badge/version-2.0.0-blue)** **![Release](https://img.shields.io/badge/Major-red)**      
*November 15th, 2025: "The DEG update"*  
  - Added **assay** argument to **SetDEGs** function
  - Added **assay** and **saveAssay** argument to **FilterDEGs** function to specify which assay to filter and what to save filtered assay as
  - Added **assay** argument to **PlotMA** to specify which DEG slot to plot
  - Added **PlotVolcano** function for direct volcano plotting functionality
  - Added **FindHVDEGs** to intersect DEGs with HVGs calculated during EDA functions
  - Added **GenSave** function to save save pheatmap objects in a manner similar to ggplpot2::ggsave()
  - Renamed eigencorr, ordcorr, and distanceHeatmap functions to **PlotEigenCorr**, **PlotOrdCorr**, and **PlotDistances** respectively to match other plotting function name conventions
  - Replaced Pasilla dataset demo with **Snail1 KO dataset** from **Matsuri et al., 2018**

## Citation  
If you use RGenEDA in your work, please cite: **DOI: dx.doi.org/10.17504/protocols.io.bp2l6z6rdgqe/v1**

## 📬 Contact
Mike Martinez M.S. - Dartmouth Genomic Data Science Core

f007qps@dartmouth.edu








