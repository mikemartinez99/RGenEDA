# RGenEDA (v1.0.0)

![status](https://img.shields.io/badge/status-in--development-orange)
![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)
![R Version](https://img.shields.io/badge/R-4.4.3-blue)

**Unified, reproducible frameworks for genomic exploratory data analyses.**


## 🚧 Development

RGenEDA is currently under development. If there is a feature you would like implemented, please submit a feature request on the Issues page or submit a pull request!


![Alt text](img/examples.png)

## Installation
To intall the RGenEDA package through R, run the following command:

```r
library(devtools)
devtools::install_github("mikemartinez99/RGenEDA")

library(RGenEDA)

```

## Usage
For a full demo, see [Pasilla Dataset Demo](https://github.com/mikemartinez99/RGenEDA/blob/main/vignettes/introduction.md)


## Change-log

**v0.0.0.900** *April 2025*
  - Initial implementation (stand-alone functions)
  - Variance calculation, PCA calculation, eigen-vector correlation, Euclidean distance calculation

**v1.0.0:** *October 16th, 2025* 

  - Created S4 `GenEDA` object to streamline usage.
  - Added `GenEDA` methods for HVG calculation and PCA stored as S4 slots.
  - Added function for count distribution and eigen vector analysis
  - Added PCA plotting function, eigen vector plotting function
  - Updated HVG plotting function
  - Updated eigencorrelation plotting function
  - Return all plots as either ggplot2 objects or heatmap objects (ComplexHeatmap, Pheatmap) to allow robust customization
  - Added unit tests
  - Added demo to vignettes

**v1.0.1** *November 3rd, 2025*

  - Modified `PlotEigenHeatmap` to use pheatmap instead of `ComplexHeatmap`
  - Fixed bug to correctly index top genes from selected PC

## Contact
Mike Martinez M.S. - Dartmouth Genomic Data Science Core

f007qps@dartmouth.edu








