#' RGenEDA: Genomic Exploratory Data Analysis
#'
#' RGenEDA collects the plots and summaries that are repeated in most bulk
#' RNA-seq exploratory data analyses (EDA), built around the DESeq2 workflow,
#' into one S4 object, `geneda`, and a set of functions that work on it.
#'
#' @section Typical workflow:
#' \enumerate{
#'   \item Build a `geneda` object with [GenEDA()] (from matrices) or
#'     [as.geneda()] (from a `SummarizedExperiment`, `DESeqDataSet` or
#'     `DESeqTransform`).
#'   \item Check the data with [PlotCountDist()], [PlotDistances()] and
#'     [PlotHVGVariance()].
#'   \item Select highly variable genes with [FindVariableFeatures()] and run
#'     PCA with [RunPCA()]. Explore it with [PlotScree()], [PlotPCA()],
#'     [extractEigen()], [PlotEigenHeatmap()] and [PlotEigenCorr()]; use
#'     [PlotOrdCorr()] for NMDS ordination instead of PCA.
#'   \item Store differential expression results with [SetDEGs()], filter them
#'     with [FilterDEGs()] and [FindHVDEGs()], and plot them with [PlotMA()],
#'     [PlotVolcano()], [PlotPValHist()] or [DEDashboard()].
#'   \item Save heatmaps with [GenSave()].
#' }
#'
#' Data are read with the accessors [getCounts()], [getNormalized()],
#' [getMetadata()], [HVGs()], [DimReduction()] and [DEGs()]. The example
#' dataset is [Snail1KO]. See `vignette("Snail1_Vignette", package =
#' "RGenEDA")` for a full walk-through.
#'
#' @section Citation:
#' Run `citation("RGenEDA")` to see how to cite the package.
#'
#' @seealso
#' \url{https://github.com/mikemartinez99/RGenEDA},
#' \url{https://mikemartinez99.github.io/RGenEDA/}
"_PACKAGE"

## Package-wide imports. Functions from other packages that are not listed
## here are called with `pkg::fun()`.
#' @importFrom dplyr left_join mutate n
#' @importFrom ggplot2 .data aes coord_cartesian element_blank element_text
#'   facet_wrap geom_boxplot geom_col geom_hline geom_line geom_point
#'   geom_text geom_tile geom_violin geom_vline ggplot ggtitle guide_legend
#'   guides labs scale_color_manual scale_fill_gradientn scale_fill_manual
#'   scale_size scale_x_continuous theme theme_bw theme_classic theme_minimal
#'   xlab xlim ylab ylim
#' @importFrom ggrepel geom_label_repel
#' @importFrom grid unit
#' @importFrom patchwork plot_annotation plot_layout
#' @importFrom tidyr pivot_longer
NULL
