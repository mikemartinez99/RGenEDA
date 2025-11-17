#' DEDashboard
#'
#' @description
#' Generate a multi-panel differential expression dashboard for a single assay in a `geneda` object.
#' The dashboard includes:
#'   1. P-value histogram
#'   2. MA plot
#'   3. Volcano plot
#'
#' @param object A `geneda` object containing DEGs in the DEGs slot.
#' @param assay Character. Name of the assay to visualize.
#' @param alpha Numeric. Significance threshold for highlighting p-values and DEGs (default 0.05).
#' @param l2fc Numeric. Log2 fold-change threshold for MA and volcano plots.
#' @param num Character. Name of the numerator group for volcano plot labeling.
#' @param den Character. Name of the denominator group for volcano plot labeling.
#'
#' @return A `patchwork` object combining plots.
#' @import ggplot2
#' @import patchwork
#' @export
DEDashboard <- function(object, assay, alpha = 0.05, l2fc = 1, num = NULL, den = NULL) {
  stopifnot(methods::is(object, "geneda"))

  if (!assay %in% names(object@DEGs)) {
    stop(paste("Assay", assay, "was not found in DEGs slot!"))
  }

  df <- DEGs(object, assay)
  if (nrow(df) == 0) stop("No differential expression results found in object@DEGs$DEG")

  requiredCols <- c("log2FoldChange", "baseMean", "pvalue", "padj")
  if (!all(requiredCols %in% colnames(df))) {
    stop("Columns missing! Required: log2FoldChange, baseMean, pvalue, padj")
  }

  # Generate plots
  phist <- PlotPValHist(object, assay, bins = 40, alpha)
  ma <- PlotMA(object, assay, alpha, l2fc)
  volc <- PlotVolcano(object, assay, alpha, l2fc, den, num)
  volc <- volc + guides(size = "none")

  # Assemble dashboard
  dashboard <- (phist) / (ma | volc) +
    plot_layout(heights = c(1, 2)) +
    plot_annotation(
      title = paste("Differential Expression Dashboard:", assay),
      theme = theme(plot.title = element_text(size = 20, face = "bold"))
    )

  return(dashboard)
}
