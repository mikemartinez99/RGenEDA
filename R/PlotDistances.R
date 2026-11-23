#' Plot a Euclidean Distance Heatmap
#'
#' @description Compute a sample-by-sample distance matrix from a `geneda`
#' object and plot it as a `pheatmap` heatmap. Optionally add metadata columns
#' as column annotation strips, and order rows and columns by hierarchical
#' clustering of the sample distances.
#'
#' @details When `reorder = TRUE`, samples are clustered once with
#' `hclust(dist(t(normalized), method), method = "complete")`. The heatmap uses
#' these same sample distances for its dendrograms
#' (`clustering_distance_rows` / `clustering_distance_cols`), so the plotted
#' order is the same as the returned `order`. When `reorder = FALSE`, no
#' clustering is done, no dendrograms are drawn, and samples are shown in the
#' column order of the normalized matrix.
#'
#' @param object A `geneda` object (or a `SummarizedExperiment`; see
#'   [geneda-SummarizedExperiment]).
#' @param method Distance method for stats::dist (default "euclidean")
#' @param reorder Logical; if `TRUE` (default), order rows and columns by
#'   hierarchical clustering of the sample distances and draw dendrograms. If
#'   `FALSE`, keep the sample order of the normalized matrix.
#' @param meta_cols Optional character vector of metadata columns to draw as
#'   column annotation strips above the heatmap.
#' @param palettes Optional named list of color vectors for the annotation
#'   strips, one per `meta_cols` entry, each named by level (e.g.
#'   `list(Genotype = c(WT = "grey40", KO = "firebrick3"))`). Columns without
#'   a palette use pheatmap's default colors. Also returned unchanged.
#' @param return One of c("object", "plot"). If "plot", draws the heatmap
#'   immediately via grid and still returns the result list. Default "object".
#'
#' @returns A list with elements:
#'   \describe{
#'     \item{`dist_matrix`}{sample x sample distance matrix, in the plotted
#'       order.}
#'     \item{`order`}{character vector of sample names in the plotted
#'       order.}
#'     \item{`heatmap`}{the `pheatmap` object.}
#'     \item{`palettes`}{`palettes`, as passed in.}
#'   }
#' @export
#'
#' @examples
#' data(Snail1KO)
#' counts <- as.matrix(Snail1KO$rawCounts)
#' counts <- counts[rowSums(counts >= 10) >= 3, ]
#' meta <- Snail1KO$metadata
#' rownames(meta) <- meta$Sample
#' norm <- log2(t(t(counts) / colSums(counts)) * 1e6 + 1)
#' obj <- GenEDA(normalized = norm, metadata = meta, counts = counts)
#' colorList <- list(Genotype = c(WT = "grey40", Snai1_KO = "firebrick3"))
#' dists <- PlotDistances(obj,
#'     meta_cols = "Genotype", palettes = colorList,
#'     return = "plot"
#' )
PlotDistances <- function(
    object,
    method = "euclidean",
    reorder = TRUE,
    meta_cols = NULL,
    palettes = NULL,
    return = c("object", "plot")) {
    object <- .toGeneda(object)
    return <- match.arg(return)

    MAT <- object@normalized
    META <- object@metadata
    stopifnot(!is.null(colnames(MAT)), !is.null(rownames(META)))
    if (!identical(rownames(META), colnames(MAT))) {
        cli::cli_abort(
            "Rownames of metadata must match colnames of normalized matrix."
        )
    }

    # Optional metadata subset for subtitle context
    if (!is.null(meta_cols)) {
        missing_cols <- setdiff(meta_cols, colnames(META))
        if (length(missing_cols) > 0L) {
            missing_str <- paste(missing_cols, collapse = ", ")
            cli::cli_abort(
                "The following metadata columns were not found: {missing_str}"
            )
        }
        META <- META[, meta_cols, drop = FALSE]
    }

    # Compute distances across samples (transpose to samples x features)
    sampleDists <- stats::dist(t(MAT), method = method)
    distMat <- as.matrix(sampleDists)
    diag(distMat) <- 0

    # Cluster the samples once on their distances (if requested). The heatmap
    # below uses the same distances and linkage, so its order matches `ord`.
    reorder <- isTRUE(reorder)
    ord <- colnames(MAT)
    if (reorder) {
        hc <- stats::hclust(sampleDists, method = "complete")
        ord <- hc$labels[hc$order]
    }

    # Optional annotations for columns using meta_cols and palettes
    annotation_col <- NULL
    annotation_colors <- NULL
    if (!is.null(meta_cols) && length(meta_cols) > 0L) {
        annotation_col <- META[, meta_cols, drop = FALSE]
        if (!is.null(palettes)) {
            # palettes should be a named list mapping each metadata column to a
            # named vector of colors
            annotation_colors <- lapply(meta_cols, function(f) {
                pal <- palettes[[f]]
                if (is.null(pal)) {
                    return(NULL)
                }
                # Ensure levels are factors to map colors correctly
                if (!is.factor(annotation_col[[f]])) {
                    annotation_col[[f]] <- as.factor(annotation_col[[f]])
                }
                pal[levels(annotation_col[[f]])]
            })
            names(annotation_colors) <- meta_cols
        }
    }

    # Choose a reversed blue palette so darker = more similar (lower distance)
    heatmap_colors <- colorRampPalette(
        rev(RColorBrewer::brewer.pal(9, "Blues"))
    )(255)

    # Build pheatmap. `distMat` is passed in the original sample order,
    # because pheatmap applies the clustering order itself; it clusters on
    # the sample distances, not on Euclidean distances between rows of
    # `distMat`.
    ph <- pheatmap::pheatmap(
        mat = distMat,
        annotation_col = annotation_col,
        annotation_colors = annotation_colors,
        cluster_rows = reorder,
        cluster_cols = reorder,
        clustering_distance_rows = sampleDists,
        clustering_distance_cols = sampleDists,
        clustering_method = "complete",
        color = heatmap_colors,
        legend = TRUE,
        silent = TRUE
    )

    if (identical(return, "plot")) {
        grid::grid.newpage()
        grid::grid.draw(ph$gtable)
    }

    return(list(
        dist_matrix = distMat[ord, ord, drop = FALSE],
        order = ord,
        heatmap = ph,
        palettes = palettes
    ))
}
