# Bioconductor submission fixes

A running log of the changes made to prepare RGenEDA for Bioconductor
submission. These fixes are released as version 0.99.1 (see `NEWS.md`). Add
new entries at the top of the section they belong to: "Statistical and naming
fixes", "Bug fixes and deprecations" or "Packaging and style fixes".

This file is listed in `.Rbuildignore` and is not part of the package build.

## Statistical and naming fixes

### `percent_var` stored as numbers — **BREAKING CHANGE**

> **Breaking change.** `DimReduction$percent_var` is now numeric, not a
> string. Code that parses it with `gsub()` will break, and objects saved
> with an older version fail `validObject()` until you re-run `RunPCA()`.

**Problem.** `generatePCs()` stored `DimReduction$percent_var` as strings
such as `"85.43 %"`. `PlotPCA()` and `PlotScree()` parsed them back with
`gsub()`. `extractEigen()` did too, until the `Contribution` fix removed
that code. The stored values were also rounded to 2 decimals.

**Fix.**

- `R/PCA_Functions.R`: `percent_var` is now a named numeric vector,
  `100 * sdev^2 / sum(sdev^2)`. It is not rounded and sums to 100. The
  startup message prints it as `PC1: 85.43%`.
- `R/utils.R`: new internal helper `.formatPercent()` formats the numbers
  for labels and messages.
- `R/PlotPCA.R`: the axis labels use `.formatPercent()` and no longer use
  `gsub()`.
- `R/PlotScree.R`: plots the numeric values directly and no longer uses
  `gsub()`.
- `R/classes.R`: the validity check rejects a non-numeric `percent_var`.
  The error message tells you to re-run `RunPCA()`.
- Vignette: describes `percent_var` as numeric, and shows
  `round(head(pca$percent_var), 2)`.
- Tests (`test-RunPCA.R`, `test-extractEigen.R`): check that the type is
  numeric, that the names are correct, that the values sum to 100 and match
  `prcomp()`, that the validity check rejects strings, and that the
  `PlotPCA()` labels and `PlotScree()` data come from the numeric values.
  The `extractEigen()` fixture now uses numbers.

**User impact.** Code that parses `percent_var` with
`gsub()` must use the numbers directly. Objects saved with an older version
fail `validObject()` until you re-run `RunPCA()`.

### Consistent PCA terms: scores and loadings

**Problem.** `DimReduction$Loadings` was renamed to `Scores` in 2.0.1, but
some text still used the old name:

- The validity check in `R/classes.R` checked for `Scores` but its error
  message said `Loadings`.
- The `PlotEigenHeatmap()` error asked for `'Loadings'`.
- The `@return` docs of `DimReduction()` and `generatePCs()` listed
  `Loadings`.
- `ExtractPCA()` called the sample scores "loadings" in its title,
  description, details and return docs.
- The vignette called `$Eigenvectors` "gene contributions". That clashes
  with the new `Contribution` column.

**Fix.** The package now uses one set of terms everywhere:

- **Scores** are per-sample coordinates: `DimReduction$Scores`
  (`prcomp()$x`) and the output of `ExtractPCA()`.
- **Loadings** (rotation) are per-gene weights: `DimReduction$Eigenvectors`
  (`prcomp()$rotation`).
- **Contribution** is the per-gene percent share of the squared loadings
  (from `extractEigen()`).

The element names are unchanged. Only messages and docs were changed. The
`generatePCs()`, `DimReduction()` and `geneda` class docs now describe each
element. The `ExtractPCA()` Rd title is now "Extract PCA Scores and
Metadata".

### `extractEigen()`: `PctVar` renamed to `Contribution` — **BREAKING CHANGE**

> **Breaking change.** The `PctVar` column of `extractEigen()` and
> `PlotEigenHeatmap()$topGenes` is now `Contribution`. Code that reads
> `$PctVar` will break. The values are unchanged.

**Problem.** `extractEigen()` returned a column called `PctVar`. The docs,
the vignette and the `PlotEigenHeatmap()` row annotation described it as
"percent variation per gene". That description was wrong.

- The code computed `(v^2 * pcVar) / sum(v^2 * pcVar) * 100`.
- `pcVar` (the PC's percent variance) is one number, so it cancels out of
  the ratio.
- Eigenvectors have unit length (`sum(v^2) == 1`).
- So the value was `100 * v^2`: each gene's share of the PC's squared
  loadings (its *contribution* to the PC). It is not the variance the gene
  explains.

**Fix.**

- `R/PCA_Functions.R`: the column is now `Contribution`, computed as
  `100 * v^2 / sum(v^2)`. The numbers are the same as before. The unused
  `percent_var` lookup was removed. A new `@details` section gives the
  formula and says that the value is not the variance the gene explains.
- `R/PlotEigenHeatmap.R`: `topGenes` and the heatmap row annotation now use
  `Contribution`. The description and comments were updated.
- `man/extractEigen.Rd`, `man/PlotEigenHeatmap.Rd`: regenerated with
  roxygen2.
- `vignettes/Snail1_Vignette.Rmd` (and the pre-rendered `.md`): the text
  now defines the contribution correctly. The printed column header is now
  `Contribution`.
- `tests/testthat/test-extractEigen.R`: tests now check for the
  `Contribution` column, check it against `100 * v^2 / sum(v^2)`, check that
  `PctVar` is gone, and check that the result does not change when
  `percent_var` changes.

**User impact.** Code that reads `$PctVar` from
`extractEigen()` or `PlotEigenHeatmap()$topGenes` must use `$Contribution`.

## Bug fixes and deprecations

### `PlotDistances()`: plotted order did not match `order`; `reorder = FALSE` had no effect

**Problem.**

- The sample order was computed with `hclust(sampleDists)`. Then
  `pheatmap(..., cluster_rows = TRUE, cluster_cols = TRUE)` clustered again,
  this time on the Euclidean distances between *rows of the distance
  matrix*, not on the sample distances.
- So the dendrogram heights were wrong. On Snail1KO the old tree had
  heights of 33–139, but the real sample distances give 23–81. The plotted
  order could also differ from the returned `order`. For example, with 8
  noisy samples (seed 3) the returned order was `S7 S8 S3 S1 S2 S5 S4 S6`,
  but the plot showed `S8 S3 S1 S2 S5 S4 S7 S6`.
- `reorder = FALSE` still gave a clustered plot.
- The docs said it drew "a ggplot2 heatmap" and that "no annotation strips
  are drawn". Neither was true.

**Fix (`R/PlotDistances.R`).** It is still a pheatmap.

- `cluster_rows` and `cluster_cols` are set to `reorder`.
- `clustering_distance_rows = sampleDists` and
  `clustering_distance_cols = sampleDists`, with
  `clustering_method = "complete"`, so the dendrograms use the sample
  distances.
- pheatmap gets the distance matrix in the original sample order, because it
  applies the tree order itself. The returned `dist_matrix` is in the plotted
  order.
- The docs now say it is a `pheatmap` heatmap, that `meta_cols` are drawn as
  column annotation strips coloured by `palettes`, how the clustering is
  done, and what each returned element is.
- Tests (`test-PlotDistances.R`): check that the plotted row and column
  order equals `order` and `hclust(dist(...), "complete")`, that the
  dendrogram heights equal the sample-distance heights, that `method =
  "manhattan"` is followed, that `reorder = FALSE` gives no trees and the
  original order, that annotation strips and their legend are drawn, and the
  seed-3 regression case.

### `PlotOrdCorr()`: p-values adjusted for multiple testing — **BEHAVIOUR CHANGE**

> **Behaviour change.** As in `PlotEigenCorr()`, the stars now come from
> adjusted p-values (BH by default). Raw p-values are still returned in
> `pval_matrix`, but they are never used for the plot.

**Problem.** `PlotOrdCorr()` had the same problem as `PlotEigenCorr()`: the
p-values for NMDS axis x covariate tests were not corrected for multiple
testing.

**Fix (`R/PlotOrdCorr.R`, `R/utils.R`).**

- `PlotOrdCorr()` has the same `p_adjust_method` argument (default
  `"BH"`) and the same return elements (`padj_matrix`, `pval_matrix`,
  `stars`, `p_adjust_method`). It also has the `PAdj` plot column and the
  caption.
- The adjustment and star code is now in shared internal helpers in
  `R/utils.R`, `.adjustAndStar()` and `.starsCaption()`. `PlotEigenCorr()`
  and `PlotOrdCorr()` both use them, so the two functions cannot drift
  apart.
- The `@returns` docs list every element, including
  `variance_explained`. The vignette note on `PlotOrdCorr()` mentions the
  correction.
- New file `tests/testthat/test-PlotOrdCorr.R`, the first tests for this
  function. The tests skip when vegan is not installed. They check that
  `padj_matrix` equals `p.adjust(pval_matrix)`, that the stars and plot
  data follow padj, the caption, and that `"none"` and unknown methods give
  an error. In the test data, `noise2` vs NMDS1 has raw p = 0.033 (one
  star before this fix) and BH padj = 0.099 (no star now).


### `PlotEigenCorr()`: p-values adjusted for multiple testing — **BEHAVIOUR CHANGE**

> **Behaviour change.** The stars on the heatmap now come from adjusted
> p-values (Benjamini-Hochberg by default), so some cells that had stars
> before may lose them. `pval_matrix` still holds the raw p-values, but the
> plot always uses the adjusted values. There is no way to plot raw
> p-values.

**Problem.** Each PC x covariate cell is a separate test, but the p-values
and stars were not corrected for multiple testing. With `num_pcs = 10` and
several covariates there are dozens of tests, and the vignette data has only
n = 6 samples.

**Fix (`R/PlotEigenCorr.R`).**

- New argument `p_adjust_method` (default `"BH"`). It can be any
  `stats::p.adjust.methods` entry except `"none"`: both p-values and padj
  are always returned, and the plot always uses padj. The internal helper
  `.checkAdjustMethod()` (in `R/utils.R`) does the validation for both
  heatmap functions.
- `stats::p.adjust()` is applied once across all tested cells, after
  untestable covariates are dropped.
- New return element `padj_matrix`, used for the stars and the plot. The
  plot data has a new `PAdj` column, and a caption names the correction
  method. `pval_matrix` (raw) and the new `p_adjust_method` element are
  returned for reference.
- The docs and vignette explain the correction and warn that with n = 6
  the results are exploratory.
- Tests (`test-PlotEigenCorr.R`): check that `padj_matrix` equals
  `p.adjust(pval_matrix)`, that the stars and plot data follow padj (with a
  case that has raw p < 0.05), that the caption names the method, and that
  `"none"` and unknown methods give an error.
- On the vignette data, Genotype vs PC1 is still `***` (raw p = 1.4e-06,
  BH padj = 7.2e-06).

### `SetDEGs()`: "genes as rownames" check did nothing

**Problem.** The check was `is.null(rownames(deg_table))`, but a data.frame
always has row names. A table with automatic row names (`1`, `2`, ...),
such as one with gene IDs in a column, was accepted without an error.

**Fix (`R/DEG_Functions.R`).** The check now uses
`.row_names_info(deg_table, type = 1L) < 0L`, which is true only for
automatic row names. The error message tells the user to set
`rownames(deg_table)`. The `@param` docs now state the requirement. A test
in `test-DEGFunctions.R` checks that a table with gene IDs in a column is
rejected, and that the same table with gene row names is accepted.

### Gene variances use `MatrixGenerics::rowVars()`

**Problem.** Gene variances were computed with `apply(x, 1, var)`, which is
slow and is not the Bioconductor idiom.

**Fix.** `FindVariableFeatures()` (`R/PCA_Functions.R`), `PlotHVGVariance()`
(`R/PlotHVGVariance.R`), the `generatePCs()` example and `test-RunPCA.R`
now use `MatrixGenerics::rowVars()`. `MatrixGenerics` was added to
`Imports`. The results are the same, and the names are kept.

### `PlotMA()`: duplicate size scale, `padj` check, `NA` handling

**Problem.**

- The plot added `scale_size()` twice. ggplot2 kept the second one and
  printed "Scale for size is already present". The second scale's title,
  "-log10 P-value", was also wrong: the size shows the *adjusted* p-value.
- `padj` was used for colour, size and labels, but the required-column
  check only covered `log2FoldChange` and `baseMean`.
- `NA` values in `padj` (from DESeq2 independent filtering or outlier
  detection) gave `NA` sizes. ggplot2 then dropped those rows with a
  "Removed rows" warning. The label subsets also used `padj` with no `NA`
  handling.

**Fix (`R/PlotMA.R`).**

- There is now one `scale_size()`, titled "-log10 adj. P-value", with
  range 0.5 to 12 (the range that was in effect before).
- `padj` is now a required column. The error lists all missing columns.
- Genes with `NA` `padj` are not significant: they are black, drawn at the
  smallest size (`-log10(padj)` is set to 0) and never labelled. One
  significance flag is computed once and used for colour and both label
  layers.
- The legend now says `Adj. P <= alpha & |LFC| >= l2fc`, which matches the
  inclusive thresholds.
- Tests (`test-PlotMA.R`): check that a missing `padj` gives an error, that
  `NA` `padj` is black, size 0 and unlabelled, that the plot builds with no
  warnings, and that there is exactly one size scale.

### `PlotPCA()`: deprecated `aes_string()` replaced

**Problem.** `aes_string()` has been deprecated since ggplot2 3.0.0. It was
used for the main mapping and the `shape_by` mapping, and every call gave
a deprecation warning.

**Fix (`R/PlotPCA.R`).** Both mappings now use
`aes(x = .data[[x]], y = .data[[y]], color = .data[[color_by]])` and
`aes(shape = .data[[shape_by]])`. A stray empty argument in `labs()` was
also removed. A new test in `test-RunPCA.R` checks that the plot gives no
warnings, that the correct PCs are mapped, and that `shape_by` works.

### `GenSave()`: device left open on error

**Problem.** `GenSave()` opened the graphics device before it checked the
heatmap object. On an invalid object it stopped with an error and the
device stayed open, so later plots went to the wrong device. The device
was also not closed if drawing failed. It used
`"gtable" %in% class(x)`.

**Fix (`R/GenSave.R`).**

- All input is checked before a device is opened: the gtable (with
  `inherits(x, "gtable")`) and the file extension.
- After the device opens, `on.exit(grDevices::dev.off(dev), add = TRUE)`
  closes it on every exit path. `dev` is the device that `GenSave()`
  opened, so the user's other devices are not touched.
- The five-way `if` chain was replaced by a lookup table of devices.
- A plain `pheatmap` object is now accepted as well as an RGenEDA result
  list. The docs already said that it was.
- Tests (`test-GenSave.R`): check that `dev.list()` does not change after
  invalid input, an unsupported extension, or a drawing error, and that a
  plain `pheatmap` object works.

## Packaging and style fixes

### Scope narrowed to bulk RNA-seq; outdated vignette text and typos

**Problem.** The vignette overview said the package works "across multiple
omics data types", but every function and the example are bulk RNA-seq /
DESeq2. The `PlotOrdCorr()` note and docs used microbiome wording
("microbiome data", "beta diversity"). The vignette also had outdated text
(`plotHVGVariance()`) and several typos.

**Fix.**

- Claim narrowed to bulk RNA-seq in the vignette overview (`.Rmd` and the
  pre-rendered `.md`), the `DESCRIPTION` `Description:` field, `?RGenEDA`,
  the `geneda` class docs and the README tagline. The overview now says the
  package is built around the DESeq2 workflow (`DESeqDataSet`, `vst()` /
  `rlog()` output and `results()` tables).
- `PlotOrdCorr()`: the vignette note and the function description now
  describe NMDS on a sample dissimilarity matrix (Bray-Curtis by default)
  and no longer mention microbiome data.
- Outdated text: `plotHVGVariance()` changed to `PlotHVGVariance()`. The
  sentence was also fixed: "can be plotted with". The "`$percent_var` …
  up to PC5" text was already fixed with the numeric `percent_var` change.
- Typos: "distrubutions" changed to "distributions", "Eucliden" to
  "Euclidean" (headings, table of contents, example file names and the
  README change-log), "has be pre-wrangled" to "have been pre-wrangled",
  "overriden" to "overridden", and "was not ran" to "was not run".
- Not changed: the package `Title:` is still "Genomic Exploratory Data
  Analysis".

**Check.** The vignette rebuilds in `R CMD check`: 0 errors, 0 warnings,
and the same 2 environment and size notes.

### SummarizedExperiment interoperability: two-way conversion and direct input

**Problem.** `geneda` holds the same data as a `SummarizedExperiment`
(assays, colData, rowData, metadata), but `as.geneda()` only converted one
way, from an SE to a geneda. Users could not get RGenEDA results back into a
`SummarizedExperiment` for other Bioconductor tools. Every function also
needed a `geneda` object.

**Fix.** The class is kept, with three additions:

1. **`geneda` → `SummarizedExperiment` coercion** (`R/se-interop.R`,
   `?"geneda-SummarizedExperiment"`). `as(obj, "SummarizedExperiment")`
   uses this layout:
   - `assays`: `counts` (if its rows match) and `normalized`.
   - `colData`: the sample metadata.
   - `rowData`: `isHVG` and `hvgRank`.
   - `metadata(se)$RGenEDA`: `HVGs`, `DimReduction`, `DEGs`, `counts` (only
     when its rows differ from `normalized`) and `version`.

   `as.geneda()` now reads this layout back. It uses an existing
   `"normalized"` assay, restores HVGs, PCA and DEG tables, and keeps colData
   names unchanged (`optional = TRUE`). So
   `as.geneda(as(obj, "SummarizedExperiment"))` is `identical()` to `obj`.
   Stored results that no longer match the SE (after subsetting by features
   or samples) are dropped with a warning. Extra matrix attributes, such as
   the `betaPriorVar` and `intercept` attributes from `DESeq2::rlog()`, are
   kept.
2. **Main functions accept a `SummarizedExperiment` directly.** Every
   analysis and plot function goes through the internal `.toGeneda()`, which
   gives a clear error for other classes. The functions that add results
   (`FindVariableFeatures()`, `RunPCA()`, `SetDEGs()`, `FilterDEGs()`) use
   `.seRoundTrip()`. They return the *caller's own object*, of the same class
   (e.g. a `DESeqDataSet` with its design and assays), with the results
   added in the layout above. A `DESeqTransform` does not get an extra
   assay, because its own values are the normalized values. The
   `@param object` and `@return` docs say this.
3. **Vignette.** A new section, "Why a separate class?", maps each slot to
   its SE location. It explains the reasons for the class: one normalized
   matrix shared by every function, validated links between the parts, DEG
   tables of any size, and dense matrices for bulk-scale PCA. It shows
   `RunPCA()` run directly on the rlog `DESeqTransform`. A new section,
   "Moving results back to a SummarizedExperiment", shows the coercion and
   the round trip, which gives `TRUE` on the vignette object.

Other changes: `S4Vectors` is added to `Imports` (`metadata()`,
`DataFrame()`), and `NAMESPACE` has `exportMethods(coerce)` and
`importFrom(methods, coerce)`. The `as.geneda()` docs say that only the
counts and normalized assays are made dense. New file
`tests/testthat/test-se-interop.R` checks the SE layout, the lossless round
trip (also with counts on different rows and with extra matrix attributes),
the drop-with-warning after subsetting, that the modifying functions return
the input SE with the user's assays and colData unchanged, that SE and
geneda inputs give the same outputs for the read-only functions, the error
for other classes, and that DESeqDataSet and DESeqTransform inputs keep
their class.

**Check.** 273 tests pass. `R CMD check`: 0 errors, 0 warnings, 2 notes.
The notes are the environment's clock and an installed size of 5.2 MB,
which is the built vignette HTML (`doc` is 3.3 MB). The source tarball is
3.3 MB, under the Bioconductor limit of 5 MB. `BiocCheck`: 0 errors,
0 warnings, and the same 6 notes as before.

### Documentation clean-up, package man page and `inst/CITATION`

**Problem.**

- Accessor titles had developer notes in them, for example "Access counts
  matrix (avoid name clash with DESeq2::counts)".
- `generatePCs()` pointed to `RGenEDA::plotVariance`, which does not exist.
- `FindVariableFeatures()` pointed to `normalized(object)`, which does not
  exist.
- There was no proper `?RGenEDA` page, and there was no `inst/CITATION`.
  The README cites a protocols.io DOI.

**Fix.**

- `R/accessors.R`: new titles. `getCounts()` is "Access the raw counts
  matrix", `getNormalized()` is "Access the normalized expression matrix",
  `getMetadata()` is "Access the sample metadata", `HVGs()` is "Access the
  highly variable gene IDs", `DimReduction()` is "Access the PCA results",
  and `DEGs()` is "Access a stored DEG table".
- `R/PCA_Functions.R`: the `generatePCs()` `vars` argument now says it takes
  a named numeric vector of variances (e.g. from
  `MatrixGenerics::rowVars()`) or a character vector of gene names (e.g.
  `HVGs(object)`). `FindVariableFeatures()` now refers to
  `getNormalized(object)`.
- `R/RGenEDA-package.R`: `?RGenEDA` (also `?"RGenEDA-package"`) is now a
  real page. It is no longer `@keywords internal`. It gives an overview,
  the typical workflow with links to each function, the accessors, the
  example data, the vignette and how to cite the package.
- New `inst/CITATION` with two entries:
  - the protocols.io record (Martinez MJ, 2025, "RGenEDA (R package) v1",
    doi:10.17504/protocols.io.bp2l6z6rdgqe/v1). The title, author and date
    come from the DOI's Crossref metadata.
  - the package manual, with the version taken from `DESCRIPTION`.

**Check.** `R CMD check`: 0 errors, 0 warnings, 1 note (the clock note
from this environment). `BiocCheck(new-package = TRUE)`: 0 errors,
0 warnings, 6 notes. None of the notes come from this fix: a suggested
`GenomicVariation` biocView, an ORCID iD, a local
`vignettes/Snail1_Vignette.html` intermediate file, 13 functions over 50
lines, one vignette line over 80 characters (a URL), and the Bioc-Devel
mailing-list check.

### `Snail1KO` dataset: build script and documentation

**Problem.** `data/Snail1KO.rda` had no script showing how it was built from
E-MTAB-5244. `?Snail1KO` did not give the dimensions, organism or cell line,
gene ID type, or licence of the source data.

**Fix.**

- New script `inst/scripts/make-Snail1KO.R`. It rebuilds the object from
  the Expression Atlas downloads, and the result is `identical()` to the
  shipped object (checked 2026-11-22). It records:
  - Source: the Expression Atlas raw-count file for E-MTAB-5244, with its
    URL, download date, size and MD5. The script warns if the MD5 changes.
  - Processing by Atlas: iRAP 1.0.3, TopHat2 2.1.1 against GRCh38
    (Ensembl release 107), and HTSeq 0.11.0 counts.
  - Gene IDs: Ensembl gene IDs relabelled to gene symbols. For each symbol
    only the first gene in Ensembl ID order is kept. This removes 1,194
    duplicated symbols and 21,302 genes with no symbol. Result: 39,364 genes
    x 6 samples, integer counts, no other filtering.
  - Metadata: built from the Atlas experiment-design file. Genotype is
    recoded to `WT` / `Snai1_KO`.
  - Licence: Expression Atlas data are CC BY 4.0.
- `R/SnailKO.R` / `man/Snail1KO.Rd`: now give the dimensions (39,364 genes
  x 6 samples), the organism and cell line (human HS578T), the gene ID type
  (gene symbols from Ensembl release 107 IDs), how the data were processed,
  the source URL, a "License" section (CC BY 4.0), and how to find the
  script with `system.file()`.

**Known quirk (not changed).** Keeping the first row per symbol also keeps
one gene that has no symbol (ENSG00000083622, zero counts in all samples).
Its row name is the empty string `""`. The script reproduces this so that
it matches the shipped data.

### Imports: `importFrom()` only, `grDevices` added, `scales` and `vegan` removed

**Problem.**

- `grDevices` was used (`importFrom(grDevices, colorRampPalette)`,
  `grDevices::png()` and others) but was not in `Imports`.
- `import(scales)` was never used.
- dplyr, ggplot2, ggrepel, grid, patchwork, RColorBrewer, scales, tidyr
  and vegan were imported whole with `import()`.
- `vegan` was only called as `vegan::` in `PlotOrdCorr()`, but it was a
  hard dependency.

**Fix.**

- `DESCRIPTION`: `grDevices` added to `Imports`. `scales` removed. `vegan`
  moved to `Suggests`.
- All `@import` tags were removed. The new file `R/RGenEDA-package.R` has
  `@importFrom` tags for only the functions that the code calls without
  `pkg::` (found with `codetools::findGlobals()`):
  - dplyr: `left_join`, `mutate`, `n`
  - ggplot2: the geoms, scales, themes and helpers in use, plus `.data`
  - ggrepel: `geom_label_repel`
  - grid: `unit`
  - patchwork: `plot_annotation`, `plot_layout`
  - tidyr: `pivot_longer`

  RColorBrewer and grid are always called as `pkg::`, so they stay in
  `Imports` with no `importFrom()` (apart from `grid::unit`).
- `PlotOrdCorr()` checks `requireNamespace("vegan")` and gives a clear
  error if vegan is missing. Its example now uses
  `@examplesIf requireNamespace("vegan", quietly = TRUE)`.
- `R/zzz.R`: added `count` and `negLog10Padj` to `globalVariables()`.
  `count` was hidden before because `import(dplyr)` made `dplyr::count`
  visible.

**User impact.** `PlotOrdCorr()` now needs `vegan` to be installed
separately.

**Check.** `R CMD check`: 0 errors, 0 warnings, 1 note. The note is
"unable to verify current time", which comes from the environment.

## Packaging and style fixes (0.99.0)

These changes are also in `NEWS.md`.

- Version changed to 0.99.0 for Bioconductor submission (continues from
  2.1.0).
- Added `as.geneda()` to build a `geneda` object from a
  `SummarizedExperiment`, `DESeqDataSet` or `DESeqTransform`. It carries over
  counts, normalized values, `colData` and optional DEG tables.
- Added `biocViews`, `URL` and `BugReports` fields to `DESCRIPTION`.
- Removed `LazyData`. Example data now loads with `data(Snail1KO)`.
- Added runnable examples for all exported functions, using the bundled
  `Snail1KO` data.
- Documented return values for `FindHVDEGs()` and the `show()` method.
- Replaced `sapply()` with `vapply()`, and `1:n` with `seq_len()`.
- The vignette uses accessors instead of `@` slot access, and reports
  `sessionInfo()`.
- Restyled code to Bioconductor style (4-space indents, 80-character
  lines).
- Added PCA tests (`test-RunPCA.R`, `test-extractEigen.R`) and `as.geneda()`
  tests.
