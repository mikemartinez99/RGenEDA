#' Plot Eigen-Correlation With Metadata Heatmap
#'
#' @description Calculate Eigen-correlations from a `geneda` object and plot a
#' publication-quality heatmap using ggplot2 with numeric labels and
#' significance stars. For continuous variables, Pearson correlation is used.
#'
#' Every PC x covariate cell is a separate test, so p-values are adjusted for
#' multiple testing with [stats::p.adjust()] across all tested cells (method
#' `p_adjust_method`, default Benjamini-Hochberg). The significance stars on
#' the heatmap use the adjusted p-values: padj < 0.001, < 0.01, < 0.05 = three,
#' two, one stars, respectively. With few samples (e.g. n = 6) these tests
#' have little power, so treat the results as exploratory.
#'
#' Two-level categorical variables are coded by factor level order, so the sign
#' of the correlation follows the level order. Categorical variables with more
#' than two levels are summarized by the multiple correlation R (the square root
#' of R^2 from a one-way linear model, score ~ variable). R is unsigned (0 to 1)
#' and its p-value comes from the model F-test, so results do not depend on how
#' levels are labelled. Variables that cannot be tested (a single level,
#' constant
#' values, or no residual degrees of freedom, e.g. a sample ID column) are
#' dropped.
#'
#' Requires PCA scores in `DimReduction(object)`; does not recompute PCA.
#'
#' @importFrom stats prcomp cor.test na.omit var as.formula
#'
#' @param object A `geneda` object (or a `SummarizedExperiment`; see
#'   [geneda-SummarizedExperiment]) containing `normalized` and `metadata`, and
#'   optionally `DimReduction` scores.
#' @param num_pcs Number of principal components to correlate.
#' @param meta_cols Optional character vector of metadata column names to
#'   include.
#'   Defaults to all metadata columns.
#' @param p_adjust_method Multiple-testing correction passed to
#'   [stats::p.adjust()]: one of `"BH"` (default), `"BY"`, `"fdr"`,
#'   `"holm"`, `"hochberg"`, `"hommel"` or `"bonferroni"`. `"none"` is not
#'   allowed: raw p-values are always returned in `pval_matrix`, and the plot
#'   always uses adjusted p-values.
#'
#' @returns A list with elements:
#'   \describe{
#'     \item{`cor_matrix`}{correlation (or multiple correlation R) per
#'       covariate x PC.}
#'     \item{`padj_matrix`}{adjusted p-values (`p_adjust_method`). These are
#'       the values used for the stars and the plot.}
#'     \item{`pval_matrix`}{unadjusted p-values, for reference.}
#'     \item{`stars`}{significance stars from `padj_matrix`.}
#'     \item{`p_adjust_method`}{the correction method used.}
#'     \item{`plot`}{the heatmap (`ggplot`).}
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
#' obj <- FindVariableFeatures(obj, nfeatures = 1000)
#' obj <- RunPCA(obj)
#' ec <- PlotEigenCorr(obj, num_pcs = 5, meta_cols = "Genotype")
#' ec$plot
PlotEigenCorr <- function(
    object,
    num_pcs = 10,
    meta_cols = NULL,
    p_adjust_method = "BH") {
    object <- .toGeneda(object)
    p_adjust_method <- .checkAdjustMethod(p_adjust_method)

    MAT <- object@normalized
    META <- object@metadata

    #----- Copy metadata (categorical columns are handled in .assocTest)
    META_num <- META

    #----- Optionally subset metadata columns
    if (!is.null(meta_cols)) {
        missing_cols <- setdiff(meta_cols, colnames(META_num))
        if (length(missing_cols) > 0L) {
            missing_str <- paste(missing_cols, collapse = ", ")
            cli::cli_abort(
                "The following metadata columns were not found: {missing_str}"
            )
        }
        META_num <- META_num[, meta_cols, drop = FALSE]
    }

    #----- Ensure META rownames match those of MAT
    if (!all(rownames(META_num) == colnames(MAT))) {
        cli::cli_abort(
            "Rownames of metadata must match the colnames of normalized matrix."
        )
    }

    #----- Obtain PCs from stored DimReduction (required)
    if (
        !(
            length(object@DimReduction) > 0L && "Scores" %in% names(
                object@DimReduction
            )
        )
    ) {
        cli::cli_abort(
            paste0(
                "DimReduction slot is empty. Please run RunPCA() before ",
                "calling PlotEigenCorr()."
            )
        )
    }
    pcs_mat <- as.data.frame(object@DimReduction[["Scores"]])
    # Align PCs to metadata order if necessary
    if (!identical(rownames(pcs_mat), rownames(META_num))) {
        pcs_mat <- pcs_mat[rownames(META_num), , drop = FALSE]
    }
    pcs <- pcs_mat[, seq_len(min(ncol(pcs_mat), num_pcs)), drop = FALSE]

    #----- Create empty matrices to store correlations and p-values
    cor_matrix <- matrix(NA, ncol = ncol(pcs), nrow = ncol(META_num))
    pval_matrix <- matrix(NA, ncol = ncol(pcs), nrow = ncol(META_num))
    rownames(cor_matrix) <- rownames(pval_matrix) <- colnames(META_num)
    colnames(cor_matrix) <- colnames(pval_matrix) <- colnames(pcs)

    #----- Loop through META variables and calculate correlation for each PC
    for (meta_var in colnames(META_num)) {
        for (pc in colnames(pcs)) {
            test <- .assocTest(pcs[[pc]], META_num[[meta_var]])
            cor_matrix[meta_var, pc] <- test[["estimate"]]
            pval_matrix[meta_var, pc] <- test[["p.value"]]
        }
    }

    #----- Remove variables that could not be tested (drop the same rows from
    #----- both)
    keep <- stats::complete.cases(cor_matrix) & stats::complete.cases(
        pval_matrix
    )
    cor_matrix <- cor_matrix[keep, , drop = FALSE]
    pval_matrix <- pval_matrix[keep, , drop = FALSE]

    #----- Adjust p-values across all tested PC x covariate cells, and
    #----- create significance stars from the adjusted values
    adj <- .adjustAndStar(pval_matrix, p_adjust_method)
    padj_matrix <- adj$padj
    stars <- adj$stars

    #----- Prepare data for ggplot (long format)
    df <- expand.grid(
        Meta = rownames(cor_matrix),
        PC = colnames(cor_matrix),
        stringsAsFactors = FALSE
    )
    # Element-wise matrix lookup using paired row/col indices (avoids outer
    # product)
    row_idx <- match(df$Meta, rownames(cor_matrix))
    col_idx <- match(df$PC, colnames(cor_matrix))
    df$Correlation <- cor_matrix[cbind(row_idx, col_idx)]
    df$PValue <- pval_matrix[cbind(row_idx, col_idx)]
    df$PAdj <- padj_matrix[cbind(row_idx, col_idx)]
    df$Stars <- stars[cbind(row_idx, col_idx)]

    #----- Set colors similar to RdBu reversed
    heatmap_colors <- colorRampPalette(
        rev(RColorBrewer::brewer.pal(9, "RdBu"))
    )(255)

    #----- Build heatmap with labels and stars
    p <- ggplot2::ggplot(df, aes(x = PC, y = Meta, fill = Correlation)) +
        geom_tile(color = NA) +
        scale_fill_gradientn(
            colors = heatmap_colors,
            limits = c(-1, 1),
            name = "Correlation"
        ) +
        geom_text(aes(label = sprintf("%.2f", Correlation)), size = 4) +
        geom_text(aes(label = Stars), nudge_y = -0.30, size = 5) +
        theme_minimal(base_size = 16) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
            axis.text.y = element_text(face = "bold"),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
        ) +
        labs(caption = .starsCaption(p_adjust_method))

    #----- Return eigencorrelation data and plot
    return(list(
        cor_matrix = cor_matrix,
        padj_matrix = padj_matrix,
        pval_matrix = pval_matrix,
        stars = stars,
        p_adjust_method = p_adjust_method,
        plot = p
    ))
}
