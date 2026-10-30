#' PlotEigenHeatmap
#'
#' @description
#' Generates a pheatmap of the top genes contributing to a specific principal
#' component from a `geneda` object. Genes are selected by absolute eigenvector
#' loading for the chosen PC, and expression values are Z-score scaled per gene.
#' Optionally annotate columns with metadata variables and add row annotations
#' showing gene-level percent variance.
#'
#' @param object A `geneda` object containing normalized expression data and
#'   PCA in `DimReduction`.
#' @param pc Character. Principal component to visualize (e.g., "PC1").
#' @param n Integer. Number of genes to select by absolute loading (default:
#'   25).
#' @param annotate_by Character vector of metadata column names for column
#'   annotations (optional).
#' @param annotate_colors Named list of color vectors for metadata columns.
#'   Names should match `annotate_by`.
#' @param return One of c("object","plot"). If "plot", draws the heatmap via
#'   grid.
#'
#' @return A list with `topGenes` (data.frame with EigenVector and PctVar),
#'   `expression` (scaled matrix), and `heatmap` (pheatmap object)
#'
#' @importFrom pheatmap pheatmap
#' @importFrom grDevices colorRampPalette
#' @import grid
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
#' obj <- FindVariableFeatures(obj, nfeatures = 1000)
#' obj <- RunPCA(obj)
#' res <- PlotEigenHeatmap(obj,
#'     pc = "PC1", n = 25,
#'     annotate_by = "Genotype",
#'     annotate_colors = list(Genotype = c(
#'         WT = "grey40", Snai1_KO =
#'             "firebrick3"
#'     )),
#'     return = "plot"
#' )
PlotEigenHeatmap <- function(
    object,
    pc = "PC1",
    n = 25,
    annotate_by = NULL,
    annotate_colors = NULL,
    return = c("object", "plot")) {
    stopifnot(methods::is(object, "geneda"))
    return <- match.arg(return)

    if (
        !(
            length(object@DimReduction) > 0L && all(
                c("Eigenvectors", "Scores") %in% names(object@DimReduction)
            )
        )
    ) {
        cli::cli_abort(
            paste0(
                "DimReduction is missing required entries 'Eigenvectors' ",
                "and 'Loadings'. Run RunPCA()."
            )
        )
    }

    # Extract gene-level percent variance using extractEigen
    pc_df <- extractEigen(object, component = pc)
    pc_df$AbsLoading <- abs(pc_df$EigenVector)
    pc_df <- pc_df[order(pc_df$AbsLoading, decreasing = TRUE), ]
    pc_top <- pc_df[
        seq_len(min(n, nrow(pc_df))), c("Gene", "EigenVector", "PctVar")
    ]
    rownames(pc_top) <- pc_top$Gene

    # Expression submatrix and Z-score scaling by gene
    mat <- object@normalized
    common_genes <- intersect(rownames(mat), pc_top$Gene)
    if (length(common_genes) == 0L) {
        cli::cli_abort(
            paste0(
                "No overlap between eigenvector gene IDs and expression ",
                "matrix rownames. Check rownames(normalized) vs ",
                "rownames(DimReduction$Eigenvectors)."
            )
        )
    }
    top_ex <- mat[common_genes, , drop = FALSE]

    # Column annotations from metadata
    annotation_col <- NULL
    if (!is.null(annotate_by) && length(annotate_by) > 0L) {
        meta <- object@metadata
        missing_cols <- setdiff(annotate_by, colnames(meta))
        if (length(missing_cols) > 0L) {
            missing_str <- paste(missing_cols, collapse = ", ")
            cli::cli_abort("Missing metadata columns: {missing_str}")
        }
        annotation_col <- meta[, annotate_by, drop = FALSE]
        # Ensure factors for color mapping
        for (cn in colnames(annotation_col)) {
            if (!is.factor(annotation_col[[cn]])) {
                annotation_col[[cn]] <- as.factor(annotation_col[[cn]])
            }
        }
    }

    # Row annotations: Percent Variance
    annotation_row <- data.frame(
        PctVar = pc_top[rownames(top_ex), "PctVar", drop = TRUE]
    )
    rownames(annotation_row) <- rownames(top_ex)

    # Colors
    heatmap_colors <- colorRampPalette(
        c("#313695", "#74ADD1", "#FFFFFF", "#F46D43", "#A50026")
    )(255)

    # Build annotation_colors list for column annotations (if provided)
    annotation_colors_list <- annotate_colors
    # Note: PctVar is numeric/continuous, pheatmap will use default gradient for
    # continuous annotations

    ph <- pheatmap::pheatmap(
        mat = as.matrix(top_ex),
        scale = "row",
        annotation_col = annotation_col,
        annotation_row = annotation_row,
        annotation_colors = annotation_colors_list,
        color = heatmap_colors,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        legend = TRUE,
        silent = TRUE
    )

    if (identical(return, "plot")) {
        grid::grid.newpage()
        grid::grid.draw(ph$gtable)
    }

    return(list(
        topGenes = pc_top,
        expression = top_ex,
        heatmap = ph
    ))
}
