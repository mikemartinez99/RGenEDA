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
