# RGenEDA 0.99.1

## Breaking changes

* `extractEigen()` and `PlotEigenHeatmap()`: the `PctVar` column and row
  annotation are renamed to `Contribution`. The value is each gene's percent
  share of the PC's squared loadings (`100 * v^2 / sum(v^2)`), not the
  percent of variance it explains. The values are unchanged.
* `DimReduction$percent_var` is now a named numeric vector of percents
  (e.g. `85.43`) instead of strings such as `"85.43 %"`. `PlotPCA()` and
  `PlotScree()` use it directly, and the validity check rejects
  non-numeric values. Re-run `RunPCA()` on objects saved with older
  versions.
* Imports use `importFrom()` only. `grDevices` was added to `Imports`, and
  `scales` was removed. `vegan` moved to `Suggests`: `PlotOrdCorr()` needs
  it installed.

## Behaviour changes

* `PlotEigenCorr()` adjusts p-values for multiple testing with
  `p.adjust()` (new argument `p_adjust_method`, default `"BH"`). It returns
  `padj_matrix` and uses it for the stars and the plot. Raw p-values are
  still in `pval_matrix`.
* `PlotOrdCorr()` uses the same multiple-testing correction
  (`p_adjust_method`, default `"BH"`) and returns `padj_matrix`.
* `p_adjust_method` in `PlotEigenCorr()` and `PlotOrdCorr()` no longer
  accepts `"none"`. Raw p-values and padj are always both returned, and the
  plots always use padj.

## New features

* Added `as(object, "SummarizedExperiment")` to convert a `geneda` object
  (with its HVGs, PCA and DEG tables) to a `SummarizedExperiment`.
  `as.geneda()` reads that layout back, so the round trip is lossless.
* All main functions accept a `SummarizedExperiment` (including
  `DESeqDataSet` and `DESeqTransform`) directly. `FindVariableFeatures()`,
  `RunPCA()`, `SetDEGs()` and `FilterDEGs()` return the input object with
  the results added.
* `as.geneda()` uses an assay named `"normalized"` when there is one.
* Added `inst/scripts/make-Snail1KO.R`, which rebuilds the `Snail1KO`
  dataset from Expression Atlas (E-MTAB-5244). `?Snail1KO` now gives the
  dimensions, organism and cell line, gene ID type, processing and licence
  (CC BY 4.0).
* Added a package man page (`?RGenEDA`) and `inst/CITATION`
  (`citation("RGenEDA")`). Cleaned up the accessor titles and removed
  references to functions that do not exist from the `generatePCs()` and
  `FindVariableFeatures()` docs.

## Bug fixes and documentation

* Consistent PCA terms: "scores" for samples (`DimReduction$Scores`,
  `ExtractPCA()`) and "loadings" (rotation) for genes
  (`DimReduction$Eigenvectors`). Fixed stale references to a `Loadings`
  element in the validity check, the `PlotEigenHeatmap()` error and the
  `DimReduction()` / `generatePCs()` docs.
* `PlotMA()`: `padj` is now a required column. `NA` values in `padj` are
  treated as not significant (smallest point, never labelled). The duplicate
  `scale_size()` was removed, and the legend now shows the inclusive
  thresholds.
* `PlotPCA()` no longer uses the deprecated `aes_string()`.
* `GenSave()` checks its input before it opens a graphics device, and always
  closes the device (even on error). It also accepts a plain `pheatmap`
  object.
* `SetDEGs()` now rejects DEG tables with automatic row names (`1`, `2`,
  ...). The old `is.null(rownames())` check could never be true.
* Gene variances are computed with `MatrixGenerics::rowVars()`
  (`MatrixGenerics` added to Imports).
* `PlotDistances()`: the heatmap now clusters on the sample distances, so
  the plotted order matches the returned `order`. `reorder = FALSE` now
  turns off clustering in the plot. The docs are corrected: it is a
  pheatmap, and `meta_cols` are drawn as annotation strips.
* Documentation now states the scope as bulk RNA-seq (DESeq2 workflow)
  instead of "multiple omics data types". Outdated vignette text and typos
  fixed.

# RGenEDA 0.99.0

* Renumbered to 0.99.0 for Bioconductor submission (continues from 2.1.0).
* Added `as.geneda()` to build a `geneda` object from a `SummarizedExperiment`,
  `DESeqDataSet`, or `DESeqTransform`, carrying over counts, normalized
  values, `colData`, and optional DEG tables. Normalized values can come from
  an assay, a transformed object, a matrix, or be computed with DESeq2 `vst`
  or log2 CPM.
* Added `biocViews`, `URL`, and `BugReports` fields to `DESCRIPTION`.
* Removed `LazyData`; load the example data with `data(Snail1KO)`.
* Added runnable examples for all exported functions using the bundled
  `Snail1KO` data.
* Documented return values for `FindHVDEGs()` and the `show()` method.
* Replaced `sapply()` with `vapply()` and `1:n` with `seq_len()`.
* Vignette now uses accessors instead of `@` slot access and reports
  `sessionInfo()`.
* Restyled code to Bioconductor style (4-space indents, 80-character lines).

# RGenEDA 2.1.0

* Fixed cosmetics of printing statements.
* Minor bug fixes in PCA functions and accessors.

# RGenEDA 2.0.2

* Fixed bug in `RunPCA()` where PCs beyond 10 were named `NA`.
* Fixed bug in `PlotScree()` where PCs beyond 10 were shown as `NA`.

# RGenEDA 2.0.1

* Added `PlotScree()`.
* PCA utilities now generate PCs dynamically up to PC10.
* Renamed `DimReduction$Loadings` to `DimReduction$Scores`.
* Added scree plot to the vignette.

# RGenEDA 2.0.0

* Added `assay` argument to `SetDEGs()`, `PlotMA()`, and `FilterDEGs()`
  (with `saveAssay`).
* Added `PlotVolcano()`, `FindHVDEGs()`, and `GenSave()`.
* Renamed `eigencorr`, `ordcorr`, and `distanceHeatmap` to `PlotEigenCorr()`,
  `PlotOrdCorr()`, and `PlotDistances()`.
* Replaced the pasilla demo with the Snail1 KO dataset (Maturi et al., 2018).

# RGenEDA 1.0.2

* Fixed example blocks in roxygen headers.
* Addressed all `R CMD check` errors, warnings, and notes.

# RGenEDA 1.0.1

* All heatmaps now use pheatmap.
* Fixed top-gene-loading bug in `PlotEigenHeatmap()`.
* Added optional `DEGs` slot to interface with DESeq2.
* Added MA plot functionality.

# RGenEDA 1.0.0

* Introduced the S4 `geneda` object.
* Added methods for HVG calculation and PCA storage.
* Added PCA, HVG variance, and eigenvector plots.
* Eigencorrelation plots now use the top HVGs from HVG selection/PCA.
* Added unit tests and a pasilla vignette.

# RGenEDA 0.0.0.900

* Initial standalone functions: variance, PCA, Euclidean distances, and
  eigenvector correlations.
