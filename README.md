# Genomic Explotatory Data Analysis (EDA) 
This package contains functions to help explore genomic data of any type. All you require is a counts matrix containing samples and feature information (i.e., genes, microbes, functions, etc...)

![Alt text](/img/RGenEDA_hex.png)

# Table of Contents
- [Installation](#installation)
- [Examples] (#examples)
- [Contact](#contact)

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
To demonstrate RGenEDA, we will utilize an example dataset of mouse mammary RNASeq. Here, we will follow the standard [DESeq2 Workflow](https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html).

## Contact
Michael Martinez M.S.

f007qps@dartmouth.edu

