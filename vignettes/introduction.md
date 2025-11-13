---
title: "Introduction to RGenEDA"
output:
  rmarkdown::html_document:
vignette: >
  %\VignetteIndexEntry{Introduction to RGenEDA}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



# Overview

`RGenEDA` is designed to provide a unified, reproducible framework for exploratory data analysis across multiple omics data types.

This vignette introduces the `RGenEDA` package using bulk RNA-seq data from the pasilla dataset.
`RGenEDA` facilitates Exploratory Data Analysis (EDA) for any omics dataset, provided a counts matrix and metadata are available.

We’ll follow a simplified version of the [deseq2 framework](https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) for normalization and use RGenEDA’s functionality to assess variance structure, dimensionality reduction, and sample relationships.

## Table of contents
- [Load and inspect the data](#load-and-inspect-the-data)
- [Define metadata](#define-metadata)
- [Processing and normalization](#processing-and-normalization)
- [Create a GenEDA object](#create-a-geneda-object)
- [Count distrubutions across samples](#count-distributions-across-samples)
- [Sample Eucliden distances with hierarchical clustering](#sample-euclidean-distances-with-hierarchical-clustering)
- [Identify highly variable genes](#identify-highly-variable-genes)
- [Principal component analysis](#principal-component-analysis)
- [Extract and visualize PCA results](#extract-and-visualize-pca-results)
- [Explore Eigen vectors of individual PCs](#explore-eigen-vectors-of-individual-pcs)
- [Correlate PCs with metadata](#correlate-pcs-with-metadata)
- [Explore DEGs](#explore-degs)

## Load and inspect the data

We begin by loading the `pasilla` dataset, which contains gene-level RNA-seq counts for *Drosophila melanogaster*.
For your own data, you can use `read_csv()` or `fread()` to import counts from a tabular file.


``` r
datafile <-  system.file( "extdata/pasilla_gene_counts.tsv", package="pasilla" )
count.table <-  read.table( datafile, header=TRUE, row.names=1, quote="", comment.char="" )
head(count.table)
#>             untreated1 untreated2 untreated3 untreated4 treated1 treated2 treated3
#> FBgn0000003          0          0          0          0        0        0        1
#> FBgn0000008         92        161         76         70      140       88       70
#> FBgn0000014          5          1          0          0        4        0        0
#> FBgn0000015          0          2          1          2        1        0        0
#> FBgn0000017       4664       8714       3564       3150     6205     3072     3334
#> FBgn0000018        583        761        245        310      722      299      308
```

## Define metadata

Next, we define the associated sample metadata.
In `pasilla`, samples differ by treatment (condition) and library type (single-end or paired-end).

``` r

# Set metadata
cond.type <- c("untreated", "untreated", "untreated", "untreated",
               "treated", "treated", "treated")
lib.type  <- c("single-end", "single-end", "paired-end", "paired-end",
               "single-end", "paired-end", "paired-end")

# Create metadata dataframe
meta <- data.frame(
  condition = cond.type,
  library = lib.type
)
rownames(meta) <- colnames(count.table)

# Color palettes for plotting
colorList <- list(
  condition = c("untreated" = "#E41A1C", "treated" = "#377EB8"),
  library   = c("single-end" = "#FF7F00", "paired-end" = "#4DAF4A")
)

head(meta)
#>            condition    library
#> untreated1 untreated single-end
#> untreated2 untreated single-end
#> untreated3 untreated paired-end
#> untreated4 untreated paired-end
#> treated1     treated single-end
#> treated2     treated paired-end
```

## Processing and normalization

We next build a `DESeqDataSet` and perform variance-stabilizing transformation (VST) for downstream EDA.
Lowly expressed genes are filtered out before running `DESeq2.`


``` r

dds <- DESeqDataSetFromMatrix(
  countData = count.table,
  colData = meta,
  design = ~ condition + library + condition:library
)
#> Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in design formula are characters, converting to factors

# Set reference levels
dds$condition <- relevel(dds$condition, ref = "untreated")
dds$library <- relevel(dds$library, ref = "single-end")

# Prefilter: keep genes with at least 10 counts in ≥3 samples
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep,]

# Run DESeq2
dds <- DESeq(dds)

# Variance-stabilizing transform
vsd <- vst(dds)
mat <- assay(vsd)

```

## Create a GenEDA object

With the normalized counts and metadata prepared, we can create a `GenEDA` object.
This object will store all components of your analysis, from normalized data and metadata (bare minimum requirements) to PCA and HVG results (downstream.) Raw counts passed to `DESeq2` can optionally be stored.



``` r

obj <- GenEDA(
  normalized = mat,
  metadata = meta)

# View object summary
obj
#> geneda object
#>   features: 8148
#>   samples:  7
#>   HVGs: 0
#>   counts: NULL
#>   DEGs: NULL
```

## Count distributions across samples

To visualize normalized count distributions across samples, the `PlotCountDist` function can be used. This is a helpful way to visualize effectiveness of normalization, as the overall distributions should be similar across samples. Samples with very low or very high overall counts compared to others might indicate problematic samples, technical artifacts, or batch effect. 


``` r

PlotCountDist(obj, split_by = "condition")
```

<div class="figure" style="text-align: center">
<img src="figure/count-dist-1.png" alt="plot of chunk count-dist"  />
<p class="caption">plot of chunk count-dist</p>
</div>


## Sample Eucliden distances with hierarchical clustering

To visualize replicate similarity, we can plot Euclidean distances between samples using the `distanceHeatmap` function.
Darker colors indicate higher similarity, while lighter colors represent dissimilar samples. This provides a quick assessment of replicate quality and metadata features that drive clustering. 


``` r

hm <- distanceHeatmap(
  obj,
  meta_cols = c("condition", "library"),
  palettes = colorList,
  return = "plot"
)
hm$heatmap

```

<div class="figure" style="text-align: center">
<img src="figure/clustering-1.png" alt="plot of chunk clustering" width="60%" />
<p class="caption">plot of chunk clustering</p>
</div>

## Identify highly variable genes

Next, we assess variance across all genes to identify those most variable across samples as these drive much of the biological signal. You can visualize the full variance profile for all genes profiled with `plotHVGVariance`.


``` r

#----- Plot variance profile
PlotHVGVariance(obj)

```

<div class="figure" style="text-align: center">
<img src="figure/hvgs-plot-1.png" alt="plot of chunk hvgs-plot" width="60%" />
<p class="caption">plot of chunk hvgs-plot</p>
</div>

Based on the plot, we’ll retain the top 2,000 most variable genes using `FindVariableFeatures` which stores the gene names in the HVGs slot of the `GenEDA` object


``` r

#----- Add HVGs to object
obj <- FindVariableFeatures(obj, 2000)
head(HVGs(obj))
#> [1] "FBgn0039155" "FBgn0029856" "FBgn0003360" "FBgn0053909" "FBgn0085787" "FBgn0025111"
```

## Principal component analysis

Using the identified HVGs, we perform PCA with `RunPCA`.
This function stores PCA results in the DimReduction slot, including:

•	`$Loadings` (sample scores)
	
•	`$Eigenvectors` (gene contributions)
	
•	`$percent_var` (Percent variance explained per component, up to PC5)

If `FindVariableFeatures` was not ran beforehand, `RunPCA` will calculate HVGs by default with 2000 features. This argument can be overriden using the `nfeatures` argument. 


``` r

obj <- RunPCA(obj)
#> Calculating principal components from top 2000 HVGs
#> Percent variations:
#>       PC1       PC2       PC3       PC4       PC5 
#> "45.35 %" "29.18 %"  "15.5 %"  "5.17 %"     "3 %"

# Inspect PCA outputs
head(obj@DimReduction$Loadings)
#>                  PC1       PC2        PC3        PC4        PC5          PC6          PC7
#> untreated1 -5.244188  3.267803 -7.3173995  0.8956240  0.6560211  0.263014919 4.443147e-14
#> untreated2 -5.741324  4.054147  4.4148996 -1.6852478  2.2255095 -0.157432838 4.412048e-14
#> untreated3 -5.300851 -4.045810  3.1725971  3.6410897 -0.7680523  0.555854065 4.431683e-14
#> untreated4 -4.404109 -4.693799 -0.6643761 -3.0558063 -2.1096526 -0.686050207 4.444540e-14
#> treated1    5.941176  8.466444  1.4317275  0.3751014 -2.0255952 -0.006515056 4.370513e-14
#> treated2    7.325773 -3.284794 -0.6180029  0.9908693  1.2177454 -2.128977372 4.237626e-14
head(obj@DimReduction$Eigenvectors)
#>                       PC1          PC2           PC3          PC4          PC5
#> FBgn0011764 -0.0001939999  0.011784784 -0.0002716843  0.026109530 -0.036764223
#> FBgn0002441 -0.0050539094 -0.005417428 -0.0025994524  0.005128475  0.017970971
#> FBgn0001276 -0.0037799827  0.012863950  0.0022284691 -0.053304421  0.009856952
#> FBgn0025864  0.0069592766 -0.034519739 -0.0074099283  0.033151524  0.009045152
#> FBgn0000063  0.0140917988  0.012765500 -0.0057142017 -0.007164571  0.024494828
#> FBgn0023507 -0.0248617605 -0.051313236 -0.0154158330  0.016058148 -0.023410164
head(obj@DimReduction$percent_var)
#>       PC1       PC2       PC3       PC4       PC5 
#> "45.35 %" "29.18 %"  "15.5 %"  "5.17 %"     "3 %"
```

## Extract and visualize PCA results

You can easily extract PCA results merged with metadata using `ExtractPCA`.
This enables flexible downstream visualization with ggplot2 or other frameworks.

To quickly plot PCA results, the `PlotPCA` function can be used. 


``` r

pcaDF <- ExtractPCA(obj)
head(pcaDF)
#>                  PC1       PC2        PC3        PC4        PC5          PC6          PC7 condition    library
#> untreated1 -5.244188  3.267803 -7.3173995  0.8956240  0.6560211  0.263014919 4.443147e-14 untreated single-end
#> untreated2 -5.741324  4.054147  4.4148996 -1.6852478  2.2255095 -0.157432838 4.412048e-14 untreated single-end
#> untreated3 -5.300851 -4.045810  3.1725971  3.6410897 -0.7680523  0.555854065 4.431683e-14 untreated paired-end
#> untreated4 -4.404109 -4.693799 -0.6643761 -3.0558063 -2.1096526 -0.686050207 4.444540e-14 untreated paired-end
#> treated1    5.941176  8.466444  1.4317275  0.3751014 -2.0255952 -0.006515056 4.370513e-14   treated single-end
#> treated2    7.325773 -3.284794 -0.6180029  0.9908693  1.2177454 -2.128977372 4.237626e-14   treated paired-end

# Plot PCA
PlotPCA(object = obj,
        x = 1,
        y = 2,
        color_by = "condition",
        colors = c("untreated" = "red", "treated" = "blue"),
        shape_by = "library")
#> Warning: `aes_string()` was deprecated in ggplot2 3.0.0.
#> ℹ Please use tidy evaluation idioms with `aes()`.
#> ℹ See also `vignette("ggplot2-in-packages")` for more information.
#> ℹ The deprecated feature was likely used in the RGenEDA package.
#>   Please report the issue to the authors.
#> This warning is displayed once every 8 hours.
#> Call `lifecycle::last_lifecycle_warnings()` to see where this warning was generated.
```

<div class="figure" style="text-align: center">
<img src="figure/extract-pca-1.png" alt="plot of chunk extract-pca" width="70%" />
<p class="caption">plot of chunk extract-pca</p>
</div>
## Explore Eigen vectors of individual PCs

We can explore the individual Eigen vectors that comprise a particular component of interest (usually PC1 and PC2) using `extractEigen`.

Similarly, to visually explore the top Eigen vectors, the `PlotEigenHeatmap` function can be ran directly. Heatmap values are the normalized expression values scaled and Z-scored.


``` r

pc1_eigen <- extractEigen(object = obj,
                          component = "PC1")
head(pc1_eigen)
#>          Gene   EigenVector       PctVar
#> 1 FBgn0011764 -0.0001939999 3.763595e-06
#> 2 FBgn0002441 -0.0050539094 2.554200e-03
#> 3 FBgn0001276 -0.0037799827 1.428827e-03
#> 4 FBgn0025864  0.0069592766 4.843153e-03
#> 5 FBgn0000063  0.0140917988 1.985788e-02
#> 6 FBgn0023507 -0.0248617605 6.181071e-02


PlotEigenHeatmap(obj,
                 pc = "PC1",
                 n = 25,
                 annotate_by = "condition",
                 annotate_colors = list(condition = c("untreated" = "red", 
                                                            "treated" = "blue")))
```

<div class="figure" style="text-align: center">
<img src="figure/eigenvecs-1.png" alt="plot of chunk eigenvecs"  />
<p class="caption">plot of chunk eigenvecs</p>
</div>

```
#> $topGenes
#>                    Gene EigenVector    PctVar
#> FBgn0003360 FBgn0003360 -0.24780515 6.1407392
#> FBgn0025111 FBgn0025111  0.21706116 4.7115549
#> FBgn0026562 FBgn0026562 -0.18774337 3.5247574
#> FBgn0000071 FBgn0000071  0.16430773 2.6997029
#> FBgn0001226 FBgn0001226  0.13402389 1.7962404
#> FBgn0003501 FBgn0003501  0.12895193 1.6628601
#> FBgn0023479 FBgn0023479 -0.11983389 1.4360162
#> FBgn0024288 FBgn0024288 -0.11828805 1.3992062
#> FBgn0011260 FBgn0011260  0.11292362 1.2751745
#> FBgn0000406 FBgn0000406 -0.11034756 1.2176584
#> FBgn0001224 FBgn0001224  0.10569481 1.1171393
#> FBgn0027279 FBgn0027279 -0.09480187 0.8987395
#> FBgn0016715 FBgn0016715 -0.09433320 0.8898753
#> FBgn0000116 FBgn0000116  0.09212341 0.8486722
#> FBgn0024315 FBgn0024315 -0.09159774 0.8390146
#> FBgn0001225 FBgn0001225  0.09134209 0.8343378
#> FBgn0002868 FBgn0002868 -0.08920522 0.7957571
#> FBgn0001137 FBgn0001137 -0.08531569 0.7278766
#> FBgn0023549 FBgn0023549 -0.08466929 0.7168889
#> FBgn0004396 FBgn0004396  0.08292638 0.6876784
#> FBgn0003748 FBgn0003748  0.08059367 0.6495339
#> FBgn0026376 FBgn0026376 -0.07749462 0.6005415
#> FBgn0000567 FBgn0000567 -0.07406692 0.5485909
#> FBgn0024984 FBgn0024984  0.07399523 0.5475294
#> FBgn0003137 FBgn0003137  0.07392526 0.5464943
#> 
#> $expression
#>             untreated1 untreated2 untreated3 untreated4  treated1  treated2  treated3
#> FBgn0000071   7.298675   7.616443   7.673259   7.732185  9.441005  9.591838  9.669924
#> FBgn0000116   8.981496   8.858281   9.142795   8.908643  9.643627 10.377361 10.178400
#> FBgn0000406   8.228731   9.001807   8.991604   9.316908  7.441680  7.419333  7.699209
#> FBgn0000567   8.656902   8.752507   8.996029   8.745890  7.992790  7.947280  7.742959
#> FBgn0001137   9.777190   9.415891   9.449740   9.591886  8.708348  8.399583  8.472375
#> FBgn0001224   8.728100   8.692523   8.618221   8.991123  9.999053  9.934651 10.158735
#> FBgn0001225   8.414694   8.299953   8.412998   8.632434  9.415161  9.502773  9.680591
#> FBgn0001226   9.154509   9.529798   9.613386   9.938920 11.001805 11.189031 11.301815
#> FBgn0002868   8.198075   8.021239   7.576221   7.779286  6.922341  6.760905  6.779488
#> FBgn0003137  11.719768  11.702291  11.819422  11.635386 12.540373 12.738203 12.575383
#> FBgn0003360  13.040072  12.908337  12.525933  12.662331  9.791321  9.811007  9.720418
#> FBgn0003501   7.401548   7.773660   7.795655   7.611522  9.189572  9.326955  9.142274
#> FBgn0003748   9.825406   9.764532   9.641737   9.633788 10.754013 10.710866 10.652714
#> FBgn0004396  10.453501  11.139993  10.889284  10.470948 12.077549 11.634709 11.655947
#> FBgn0011260   7.298675   7.860262   7.355831   7.245390  9.212307  8.698174  8.663440
#> FBgn0016715   8.106460   8.353993   8.330111   8.396878  7.031744  7.221243  7.160907
#> FBgn0023479  11.767226  12.268267  12.247328  12.075607 10.673837 10.684594 10.555510
#> FBgn0023549   8.193639   8.723803   8.760993   8.794971  7.469244  7.631016  7.625931
#> FBgn0024288   7.385010   7.557842   7.661529   7.600928  5.993462  6.192237  6.119027
#> FBgn0024315   7.857037   7.792984   7.917194   7.833633  6.726166  6.714770  6.759058
#> FBgn0024984   7.280619   7.219320   7.228245   7.314733  7.996406  8.215415  8.232886
#> FBgn0025111   8.940574   8.849470   9.048080   8.979631 11.382618 11.657189 11.704143
#> FBgn0026376  10.498716  11.217161  11.020944  10.684737 10.155819  9.806767  9.863900
#> FBgn0026562  15.644058  15.643428  16.245294  16.440073 13.165698 13.853773 13.921332
#> FBgn0027279  11.910754  12.039229  12.011273  11.894614 10.891415 10.799751 10.765910
#> 
#> $heatmap
```

## Correlate PCs with metadata

To interpret principal components, we can correlate them with sample metadata using `eigencorr`.
This function computes Pearson correlations and displays them as a heatmap, helping to reveal which metadata features are most associated with major axes of variation. This function returns a list of 4 elements:

•	`$cor_matrix` (Pearson correlation values)
	
•	`$pval_matrix` (Associated correlation p-values)
	
•	`$stars` (asterisk representations of p-values)
	
•	`$plot` (Eigencorr plot, as a `ggplot2` object)

**Note:** `ordcorr` can be used for microbiome data as it correlates metadata features with NMDS beta values rather than PCs.



``` r

ec <- eigencorr(obj, num_pcs = 5)
ec$plot

```

<div class="figure" style="text-align: center">
<img src="figure/eigencorr-1.png" alt="plot of chunk eigencorr" width="70%" />
<p class="caption">plot of chunk eigencorr</p>
</div>

## Explore DEGs

We can also explore the differentially expressed genes by appending these to our object. We can also quickly filter out data and save it as a new DEG slot with defined assay name. 


``` r

res <- results(dds) |> 
  as.data.frame()

obj <- SetDEGs(obj, res)
obj <- FilterDEGs(obj,
                  padj_thresh = 0.05,
                  log2FC_thresh = 1,
                  assayName = "padj05_lfc1")

x <- PlotMA(obj,
            alpha = 0.05,
            fc = 1)
#> Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
#> ℹ Please use `linewidth` instead.
#> ℹ The deprecated feature was likely used in the RGenEDA package.
#>   Please report the issue to the authors.
#> This warning is displayed once every 8 hours.
#> Call `lifecycle::last_lifecycle_warnings()` to see where this warning was generated.
#> Scale for size is already present.
#> Adding another scale for size, which will replace the existing scale.
x
```

<div class="figure" style="text-align: center">
<img src="figure/DEGs-1.png" alt="plot of chunk DEGs" width="70%" />
<p class="caption">plot of chunk DEGs</p>
</div>
