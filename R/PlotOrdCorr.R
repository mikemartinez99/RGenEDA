#' Plot Ordination Correlations With Metadata Heatmap
#'
#' @description Calculate correlations between NMDS ordination axes (from a
#' sample dissimilarity matrix, Bray-Curtis by default) and metadata, and plot a
#' ggplot2 heatmap with numeric labels and significance stars.
#'
#' For continuous variables, Pearson correlation is used.
#'
#' Every NMDS axis x covariate cell is a separate test, so p-values are
#' adjusted for multiple testing with [stats::p.adjust()] across all tested
#' cells (method `p_adjust_method`, default Benjamini-Hochberg). The
#' significance stars on the heatmap use the adjusted p-values: padj < 0.001,
#' < 0.01, < 0.05 = three, two, one stars, respectively. With few samples
#' these tests have little power, so treat the results as exploratory.
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
#' Requires the `vegan` package (in Suggests).
#'
#' Default num_mds is 10, but function will internally set max number of NMDS
#' to be n-1 where n is your number of samples.
#'
#' @importFrom stats as.dist cor dist
#'
#' @param object A `geneda` object (or a `SummarizedExperiment`; see
#'   [geneda-SummarizedExperiment]) containing `normalized` and `metadata`.
#' @param num_mds Number of NMDS axes to correlate. Default = 10
#' @param meta_cols Optional character vector of metadata column names to
#'   include.
#'   Defaults to all metadata columns.
#' @param distance Distance metric for `vegan::vegdist()` (default "bray").
#' @param p_adjust_method Multiple-testing correction passed to
#'   [stats::p.adjust()]: one of `"BH"` (default), `"BY"`, `"fdr"`,
#'   `"holm"`, `"hochberg"`, `"hommel"` or `"bonferroni"`. `"none"` is not
#'   allowed: raw p-values are always returned in `pval_matrix`, and the plot
#'   always uses adjusted p-values.
#'
#' @returns A list with elements:
#'   \describe{
#'     \item{`cor_matrix`}{correlation (or multiple correlation R) per
#'       covariate x NMDS axis.}
#'     \item{`padj_matrix`}{adjusted p-values (`p_adjust_method`). These are
#'       the values used for the stars and the plot.}
#'     \item{`pval_matrix`}{unadjusted p-values, for reference.}
#'     \item{`stars`}{significance stars from `padj_matrix`.}
#'     \item{`p_adjust_method`}{the correction method used.}
#'     \item{`plot`}{the heatmap (`ggplot`).}
#'     \item{`variance_explained`}{squared correlations between ordination
#'       distances and observed dissimilarities, per NMDS axis.}
#'   }
#' @export
#'
#' @examplesIf requireNamespace("vegan", quietly = TRUE)
#' data(Snail1KO)
#' counts <- as.matrix(Snail1KO$rawCounts)
#' counts <- counts[rowSums(counts >= 10) >= 3, ]
#' meta <- Snail1KO$metadata
#' rownames(meta) <- meta$Sample
#' norm <- log2(t(t(counts) / colSums(counts)) * 1e6 + 1)
#' obj <- GenEDA(normalized = norm, metadata = meta, counts = counts)
#' set.seed(42)
#' oc <- PlotOrdCorr(obj, num_mds = 3, meta_cols = "Genotype")
#' oc$plot
PlotOrdCorr <- function(
    object,
    num_mds = 10,
    meta_cols = NULL,
    distance = "bray",
    p_adjust_method = "BH") {
    object <- .toGeneda(object)
    p_adjust_method <- .checkAdjustMethod(p_adjust_method)

    if (!requireNamespace("vegan", quietly = TRUE)) {
        cli::cli_abort(
            paste0(
                "PlotOrdCorr() needs the 'vegan' package. Install it with ",
                "install.packages(\"vegan\")."
            )
        )
    }

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

    #----- Compute distance matrix across samples (rows must be samples)
    distmat <- vegan::vegdist(t(MAT), method = distance)

    #----- Adjust k if it's larger than n - 1
    max_k <- nrow(as.matrix(distmat)) - 1
    if (num_mds > max_k) {
        cli::cli_inform(
            paste0(
                "Reducing num_mds from {num_mds} to {max_k} (max allowed ",
                "for {max_k + 1} samples)."
            )
        )
        num_mds <- max_k
    }

    #----- NMDS ordination
    dte <- vegan::metaMDS(
        distmat,
        distance = distance,
        k = num_mds,
        trymax = 100,
        trace = FALSE
    )
    nmds_scores <- as.data.frame(vegan::scores(dte))

    #----- Approximate variance explained by each NMDS axis
    # Compute Euclidean distances between points in NMDS space for each k
    coords <- as.matrix(nmds_scores)
    original_dist <- as.dist(distmat)
    var_explained <- numeric(ncol(coords))

    for (i in seq_len(ncol(coords))) {
        reduced_dist <- dist(coords[, seq_len(i), drop = FALSE])
        r2 <- cor(original_dist, reduced_dist)^2
        var_explained[i] <- ifelse(
            i == 1,
            r2,
            r2 - sum(var_explained[seq_len(i - 1)])
        )
    }
    names(var_explained) <- colnames(nmds_scores)

    # Align to metadata order if necessary
    if (!identical(rownames(nmds_scores), rownames(META_num))) {
        nmds_scores <- nmds_scores[rownames(META_num), , drop = FALSE]
    }

    #----- Create empty matrices to store correlations and p-values
    cor_matrix <- matrix(NA, ncol = ncol(nmds_scores), nrow = ncol(META_num))
    pval_matrix <- matrix(NA, ncol = ncol(nmds_scores), nrow = ncol(META_num))
    rownames(cor_matrix) <- rownames(pval_matrix) <- colnames(META_num)
    colnames(cor_matrix) <- colnames(pval_matrix) <- colnames(nmds_scores)

    #----- Loop through META variables and calculate correlation for each NMDS
    #----- axis
    for (meta_var in colnames(META_num)) {
        for (axis in colnames(nmds_scores)) {
            test <- .assocTest(nmds_scores[[axis]], META_num[[meta_var]])
            cor_matrix[meta_var, axis] <- test[["estimate"]]
            pval_matrix[meta_var, axis] <- test[["p.value"]]
        }
    }

    #----- Remove variables that could not be tested (drop the same rows from
    #----- both)
    keep <- stats::complete.cases(cor_matrix) & stats::complete.cases(
        pval_matrix
    )
    cor_matrix <- cor_matrix[keep, , drop = FALSE]
    pval_matrix <- pval_matrix[keep, , drop = FALSE]

    #----- Adjust p-values across all tested axis x covariate cells, and
    #----- create significance stars from the adjusted values
    adj <- .adjustAndStar(pval_matrix, p_adjust_method)
    padj_matrix <- adj$padj
    stars <- adj$stars

    #----- Prepare data for ggplot (long format)
    df <- expand.grid(
        Meta = rownames(cor_matrix),
        Axis = colnames(cor_matrix),
        stringsAsFactors = FALSE
    )
    row_idx <- match(df$Meta, rownames(cor_matrix))
    col_idx <- match(df$Axis, colnames(cor_matrix))
    df$Correlation <- cor_matrix[cbind(row_idx, col_idx)]
    df$PValue <- pval_matrix[cbind(row_idx, col_idx)]
    df$PAdj <- padj_matrix[cbind(row_idx, col_idx)]
    df$Stars <- stars[cbind(row_idx, col_idx)]

    #----- Set colors similar to RdBu reversed
    heatmap_colors <- colorRampPalette(
        rev(RColorBrewer::brewer.pal(9, "RdBu"))
    )(255)

    #----- Build heatmap with labels and stars
    p <- ggplot2::ggplot(df, aes(x = Axis, y = Meta, fill = Correlation)) +
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
            axis.title.y = element_blank()
        ) +
        labs(caption = .starsCaption(p_adjust_method))

    #----- Return ord-correlation data and plot
    return(list(
        cor_matrix = cor_matrix,
        padj_matrix = padj_matrix,
        pval_matrix = pval_matrix,
        stars = stars,
        p_adjust_method = p_adjust_method,
        plot = p,
        variance_explained = var_explained
    ))
}
