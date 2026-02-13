# RGenEDA 
**Clean, unified, streamlined, and reproducible frameworks for genomic exploratory data analyses**  


 <!-- badges: start -->
 ![status](https://img.shields.io/badge/status-stable-green) 
 ![Version](https://img.shields.io/badge/version-2.0.2-blue)
  [![pkgdown documentation](https://img.shields.io/badge/docs-pkgdown-blue.svg)](https://mikemartinez99.github.io/RGenEDA/)
 ![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)
  [![R-CMD-check](https://github.com/mikemartinez99/RGenEDA/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mikemartinez99/RGenEDA/actions/workflows/R-CMD-check.yaml)
  [![Codecov test coverage](https://codecov.io/gh/mikemartinez99/RGenEDA/graph/badge.svg)](https://app.codecov.io/gh/mikemartinez99/RGenEDA)
 
  <!-- badges: end -->


> [!WARNING]
> **RGenEDA has updated to v2.0.0 which introduced some breaking changes. Please ensure you are using the most up-to-date version**

--- 


![Alt text](https://raw.githubusercontent.com/mikemartinez99/RGenEDA/main/img/examples.png)


## Installation
Install the latest version from Github

```r
install.packages("devtools") # if not already instaled

library(devtools)
devtools::install_github("mikemartinez99/RGenEDA")

library(RGenEDA)

```

## Usage
For a full demo, see [Snail1 KO Dataset Demo](https://mikemartinez99.github.io/RGenEDA/articles/Snail1_Vignette.html)  

## Change-log

**![Version](https://img.shields.io/badge/version-0.0.0.900-blue)** **![Release](https://img.shields.io/badge/Initial-grey)**
*April 2025: "The poor-man's EDA"*<br>  
  - Implemented core, standalone functions<br>  
    - Variance calculation, PCA calculation, and Eucliden distances, and Eigenvector correlations<br>

**![Version](https://img.shields.io/badge/version-1.0.0-blue)** **![Release](https://img.shields.io/badge/Major-red)**
*October 16th, 2025: "The GenEDA Object Update*<br>  
  - Introduced S4 **GenEDA** object for streamlined data handling<br>    
  - Added **methods** for HVG calculation and PCA storage within S4 slots<br>    
  - Added **visualization functions**<br>    
    - PCA plots<br>    
    - Variance plots for HVGs<br>    
    - Eigenvector (gene-loading) plots<br>    
  - Revamped eigencorrelation plots to only use top number of HVGs identified in HVG selection/PCA<br>    
  - Enhanced plot output flexibility (return **ggplot2, pheatmap**, or **ComplexHeatmap** objects for full customization<br>    
  - Added **unit tests**<br>    
  - Added **Pasilla dataset vignette**<br>  

**![Version](https://img.shields.io/badge/version-1.0.1-blue)** **![Release](https://img.shields.io/badge/Patch-green)**   
*November 3rd, 2025*<br> 
  - Rely soley on **Pheatmap** for all heatmap visualizations<br>  
  - Fixed top-gene-loading bug for selected PCs in PlotEigenHeatmap function<br>    
  - Added optional **DEGs slot** to GenEDA object to interface with **DESeq2**<br>   
  - Added **MA plot** functionality<br>  

**![Version](https://img.shields.io/badge/version-1.0.2-blue)** **![Release](https://img.shields.io/badge/Patch-green)**      
*November 10th, 2025*<br>  
  - Fixed example blocks in roxygen headers<br>    
  - Addressed all R cmd build errors, warnings, and notes<br>  

**![Version](https://img.shields.io/badge/version-2.0.0-blue)** **![Release](https://img.shields.io/badge/Major-red)**  
*November 15th, 2025: "The DEG update"*<br>  
  - Added **assay** argument to **SetDEGs** function<br>  
  - Added **assay** and **saveAssay** argument to **FilterDEGs** function to specify which assay to filter and what to save filtered assay as<br>    
  - Added **assay** argument to **PlotMA** to specify which DEG slot to plot<br>    
  - Added **PlotVolcano** function for direct volcano plotting functionality<br>    
  - Added **FindHVDEGs** to intersect DEGs with HVGs calculated during EDA functions<br>    
  - Added **GenSave** function to save save pheatmap objects in a manner similar to ggplpot2::ggsave()<br>    
  - Renamed eigencorr, ordcorr, and distanceHeatmap functions to **PlotEigenCorr**, **PlotOrdCorr**, and **PlotDistances** respectively to match other plotting function name conventions<br>    
  - Replaced Pasilla dataset demo with **Snail1 KO dataset** from **Matsuri et al., 2018**<br>

**![Version](https://img.shields.io/badge/version-2.0.1-blue)** **![Release](https://img.shields.io/badge/Patch-green)**      
*January 21st, 2026*<br>  
  - Added in **PlotScree** function<br>
  - Updated PCA utilities to dynamically generate PCs up until PC10
  - Renamed the **DimReduction$Loadings** slot to **DimReduction$Scores** to more accurately reflect the information contained in this slot
  - Added Scree Plot to vignette<br>

**![Version](https://img.shields.io/badge/version-2.0.2-blue)** **![Release](https://img.shields.io/badge/Patch-green)**      
*February 13th, 2026*<br>  
  - Fixed bug in RunPCA where PCs beyond 10 are named as NA <br>
  - Fixed bug in PlotScree where PCs beyond 10 are shown as NA <br>

## Citation  
If you use RGenEDA in your work, please cite: **dx.doi.org/10.17504/protocols.io.bp2l6z6rdgqe/v1**

## 📬 Contact
Mike Martinez M.S. - Dartmouth Genomic Data Science Core

f007qps@dartmouth.edu








