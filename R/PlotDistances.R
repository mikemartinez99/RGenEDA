#' Plot a Euclidean Distance Heatmap
#'
#' @description Compute a sample-by-sample distance matrix from a `geneda`
#' object and plot a ggplot2 heatmap. Optionally subset metadata columns and
#' reorder rows/columns by hierarchical clustering.
#'
#' @import ggplot2
#' @import scales
#' @import RColorBrewer
#' @import grid
#'
#' @param object A `geneda` object
#' @param method Distance method for stats::dist (default "euclidean")
#' @param reorder Logical; if TRUE, reorder rows/cols by hclust on the distance
#'   matrix. Default TRUE
#' @param meta_cols Optional character vector of metadata columns to display in
#'   the title/subtitle for quick reference (no annotation strips are drawn).
#' @param palettes Optional named list of color vectors that can be used
#'   downstream (returned) if the caller wants to annotate elsewhere.
#' @param return One of c("object", "plot"). If "plot", draws the heatmap
#'   immediately via grid and still returns the result list. Default "object".
#'
#' @returns A list with elements: dist_matrix, order (character vector of sample
#'   names), heatmap (pheatmap object), palettes (as passed through)
#' @export
#'
#' @examples
#' \donttest{
#' mock_norm <- matrix(rnorm(12, mean = 0, sd = 2), nrow = 4, ncol = 3)
#' colnames(mock_norm) <- paste0("Sample", 1:3)
#' rownames(mock_norm) <- paste0("Gene", 1:4)
#' mock_meta <- data.frame(condition = c("A","B","A"), row.names = colnames(mock_norm))
#' obj <- GenEDA(normalized = mock_norm, metadata = mock_meta)
#' colorList <- list(condition = c("A" = "red", "B" = "blue"))
#' PlotDistances(
#'     obj,
#'     meta_cols = c("condition"),
#'     palettes = colorList,
#'     return = "plot")
#'     }
PlotDistances <- function(object,
                          method = "euclidean",
                          reorder = TRUE,
                          meta_cols = NULL,
                          palettes = NULL,
                          return = c("object", "plot")) {
  stopifnot(methods::is(object, "geneda"))
  return <- match.arg(return)

  MAT <- object@normalized
  META <- object@metadata
  stopifnot(!is.null(colnames(MAT)), !is.null(rownames(META)))
  if (!identical(rownames(META), colnames(MAT))) {
    stop("Rownames of metadata must match colnames of normalized matrix.")
  }

  # Optional metadata subset for subtitle context
  if (!is.null(meta_cols)) {
    missing_cols <- setdiff(meta_cols, colnames(META))
    if (length(missing_cols) > 0L) {
      stop(paste0("The following metadata columns were not found: ", paste(missing_cols, collapse = ", ")))
    }
    META <- META[, meta_cols, drop = FALSE]
  }

  # Compute distances across samples (transpose to samples x features)
  sampleDists <- stats::dist(t(MAT), method = method)
  distMat <- as.matrix(sampleDists)
  diag(distMat) <- 0

  # Reorder by hierarchical clustering if requested
  ord <- colnames(MAT)
  if (isTRUE(reorder)) {
    hc <- stats::hclust(sampleDists, method = "complete")
    ord <- hc$labels[hc$order]
    distMat <- distMat[ord, ord, drop = FALSE]
  }

  # Optional annotations for columns using meta_cols and palettes
  annotation_col <- NULL
  annotation_colors <- NULL
  if (!is.null(meta_cols) && length(meta_cols) > 0L) {
    annotation_col <- META[, meta_cols, drop = FALSE]
    if (!is.null(palettes)) {
      # palettes should be a named list mapping each metadata column to a named vector of colors
      annotation_colors <- lapply(meta_cols, function(f) {
        pal <- palettes[[f]]
        if (is.null(pal)) return(NULL)
        # Ensure levels are factors to map colors correctly
        if (!is.factor(annotation_col[[f]])) annotation_col[[f]] <- as.factor(annotation_col[[f]])
        pal[levels(annotation_col[[f]])]
      })
      names(annotation_colors) <- meta_cols
    }
  }

  # Choose a reversed blue palette so darker = more similar (lower distance)
  heatmap_colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(255)

  # Build pheatmap (shows dendrograms by default)
  ph <- pheatmap::pheatmap(
    mat = distMat,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    color = heatmap_colors,
    legend = TRUE,
    silent = TRUE
  )

  if (identical(return, "plot")) {
    grid::grid.newpage(); grid::grid.draw(ph$gtable)
  }

  return(list(
    dist_matrix = distMat,
    order = ord,
    heatmap = ph,
    palettes = palettes
  ))
}

