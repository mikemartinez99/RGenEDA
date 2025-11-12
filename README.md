# RGenEDA (v1.0.2)

![status](https://img.shields.io/badge/status-in--development-orange)
![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)
![R Version](https://img.shields.io/badge/R-4.4.3-blue)  <!-- badges: start -->
  [![R-CMD-check](https://github.com/mikemartinez99/RGenEDA/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mikemartinez99/RGenEDA/actions/workflows/R-CMD-check.yaml)
  [![Codecov test coverage](https://codecov.io/gh/mikemartinez99/RGenEDA/graph/badge.svg)](https://app.codecov.io/gh/mikemartinez99/RGenEDA)
  <!-- badges: end -->

**Unified, reproducible frameworks for genomic exploratory data analyses.**

# Development <!-- omit in toc -->

> [!IMPORTANT]
> RenEDA is currently under active development. If there is a feature you would like implemented, please submit a feature request on the Issues page of submit a pull request. For information on new additions see the Change Log!

--- 


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

**v1.0.0:** *October 16th, 2025* (Major)

  - Created S4 `GenEDA` object to streamline usage.
  - Added `GenEDA` methods for HVG calculation and PCA stored as S4 slots.
  - Added function for count distribution and eigen vector analysis
  - Added PCA plotting function, eigen vector plotting function
  - Updated HVG plotting function
  - Updated eigencorrelation plotting function
  - Return all plots as either ggplot2 objects or heatmap objects (ComplexHeatmap, Pheatmap) to allow robust customization
  - Added unit tests
  - Added demo to vignettes

**v1.0.1** *November 3rd, 2025* (Minor)

  - Modified `PlotEigenHeatmap` to use pheatmap instead of `ComplexHeatmap`
  - Fixed bug to correctly index top genes from selected PC
  - Added DEGs slot to GenEDA object and support for MA plots

**v1.0.2** *November 10th, 2025* (Patch)

  - Fixed example blocks in roxygen headers
  - Addressed R cmd build errors, warnings, and notes

## Contact
Mike Martinez M.S. - Dartmouth Genomic Data Science Core

f007qps@dartmouth.edu








