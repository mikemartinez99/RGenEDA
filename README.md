# Genomic Explotatory Data Analysis (EDA) 
This package contains functions to help explore genomic data of any type. All you require is a counts matrix containing samples and feature information (i.e., genes, microbes, functions, etc...)

![Alt text](/img/RGenEDA_hex.png)

# Table of Contents
- [Installation](#installation)
- [Examples](#examples)
- [Contact](#contact)
- [License](#license)

## Installation
To intall the RGenEDA package through R, run the following command:

```r
library(devtools)
devtools::install_github("mmartinez99/RGenEDA/")

library(RGenEDA)

```

RGenEDA requires a few dependencies. To install these dependencies, run the following code: 
To install these packages
```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")
install.packages(c("RColorBrewer", "grid, "magick", "scales", "dendextend", "pheatmap"))
```
## Examples
To demonstrate RGenEDA, we will utilize the pasilla dataset from R package `pasilla`. We will follow the standard [DESeq2 Workflow](https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html), but genomic data generated from other software is also acceptable. See [Introduction vignette](https://github.com/mikemartinez99/RGenEDA/blob/main/vignettes/introduction.Rmd) for more information. These are very basic examples and can be scaled accordingly to any number of metadata features.

### Sample Distance heatmap
The sample distance heatmap is helpful for visualizing sample-to-sample Euclidean distances. By overlaying metadata features of interest, you can guage how your samples are clustering to identify features driving sample differences, or even identify batch effects. Here, we visualize two metadata categories, treatment and library-type. 
![Alt text](img/Sample_Distance_HM.tiff)

### Variable Feature Selection
When performing principal component analysis, the first step is to identify highly variable features. The `plotVariance` function ranks genes by decreasing variance to help you identify a justifiable cutoff for the number of features to include in principal component analysis. This function also saves the ranked variances as a vector to be easily passed to the `generatePCs` function. This function will output results for principal component analysis and allows users to save results to a dataframe for highly customizable plotting with ggplot2.
![Alt text](img/Variable_Features.tiff)

### Eigenvector correlations
Eigenvectors (PCs) represent individual axes of variation in data. We can correlate metadata features to these eigenvectors to identify which features are highly correlated with major axes of variation. 
![Alt text](img/EigenCorrelations.tiff)

## Contact
Michael Martinez M.S.

f007qps@dartmouth.edu

## License
This package is licensed under the **GNU General Public License v3.0 (GPL-3.0)**

You are free to:
- Share - Copy and redistribute the material in any medium or format.
- Adapt - Modify, remix, and build upon the material for any purpose, even commercially.

**Under the following terms**
- Copyleft - If you distribute versions of this package, you must also license them under GPL-3.0
- Attribution - You must give appropriate credit to the original authors.
- No additional restrictions - You may not impose legal terms that prevent others from using the freedoms granted by this license.

