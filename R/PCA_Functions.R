#' Find and store highly variable genes (HVGs)
#'
#' Computes per-feature variance on `getNormalized(object)`, ranks genes by
#' decreasing variance, and stores the top `nfeatures` gene IDs into the `HVGs`
#' slot. Returns the updated object.
#'
#' @param object A `geneda` object (or a `SummarizedExperiment`; see
#'   [geneda-SummarizedExperiment]).
#' @param nfeatures Number of HVGs to store
#' @return Updated `geneda` object.
#'   If `object` is a `SummarizedExperiment`, the same object with the
#'   results added (see [geneda-SummarizedExperiment]).
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
#' length(HVGs(obj))
FindVariableFeatures <- function(object, nfeatures) {
    if (methods::is(object, "SummarizedExperiment")) {
        return(.seRoundTrip(object, FindVariableFeatures, nfeatures))
    }
    object <- .toGeneda(object)
    nfeatures <- as.integer(nfeatures)

    vars <- MatrixGenerics::rowVars(object@normalized, na.rm = TRUE)
    if (is.null(names(vars))) names(vars) <- rownames(object@normalized)
    vars <- sort(vars, decreasing = TRUE)
    nfeatures <- max(1L, min(nfeatures, length(vars)))
    object@HVGs <- names(vars)[seq_len(nfeatures)]
    validObject(object)
    object
}

#' generatePCs
#'
#' @description Generate principal component analysis data that can be used
#'   in downstream analyses.
#'
#' @param mat A data matrix where rows are features and columns are samples
#' @param vars Either a named numeric vector of gene variances (for example
#'   from `MatrixGenerics::rowVars(mat)`), from which the top `nFeatures`
#'   genes are used, or a character vector of gene names to use directly
#'   (for example `HVGs(object)`).
#' @param nFeatures Number of top features to generate principal components on.
#'
#' @returns A list with three elements:
#'   \describe{
#'     \item{`Scores`}{data.frame of sample scores (samples x PCs), i.e.
#'       `prcomp()$x`.}
#'     \item{`Eigenvectors`}{data.frame of gene loadings (genes x PCs), i.e.
#'       the rotation matrix `prcomp()$rotation`.}
#'     \item{`percent_var`}{named numeric vector of the percent of total
#'       variance explained by each PC (0-100; sums to 100).}
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
#' vars <- MatrixGenerics::rowVars(norm)
#' pcs <- generatePCs(norm, vars, nFeatures = 1000)
#' pcs$percent_var
generatePCs <- function(mat, vars, nFeatures) {
    var_features_n <- nFeatures
    select <- if (is.character(vars)) {
        vars
    } else {
        order(vars, decreasing = TRUE)[seq_len(var_features_n)]
    }
    vsd_sub <- mat[select, ]
    vsd_sub <- t(vsd_sub)
    pca <- prcomp(vsd_sub)

    nPC <- ncol(pca$x)

    # Percent of total variance explained by each PC (numeric, 0-100)
    percentVar <- pca$sdev^2 / sum(pca$sdev^2) * 100
    names(percentVar) <- paste0("PC", seq_len(nPC))

    cli::cli_inform("Percent variance explained:")
    shownVar <- percentVar[seq_len(min(10, length(percentVar)))]
    cli::cli_inform(
        paste0(names(shownVar), ": ", .formatPercent(shownVar))
    )

    pca_res <- list()
    pca_df <- as.data.frame(pca$x[, seq_len(nPC), drop = FALSE])
    pca_eigenvecs <- as.data.frame(pca$rotation[, seq_len(nPC), drop = FALSE])
    pca_res[["Scores"]] <- pca_df
    pca_res[["Eigenvectors"]] <- pca_eigenvecs
    pca_res[["percent_var"]] <- percentVar
    return(pca_res)
}


#' Run PCA and store in `DimReduction`
#'
#' Uses `generatePCs()` under the hood.
#' If `HVGs` are empty, selects HVGs via `FindVariableFeatures(object,
#' nfeatures)`
#' with default `nfeatures = 2000`. If `HVGs` are present, the PCA uses
#' `NFEATURES = length(HVGs(object))`. This aligns the PCA feature count with
#' the
#' HVG selection while allowing an explicit override of `nfeatures` when empty.
#'
#' @param object A `geneda` object (or a `SummarizedExperiment`; see
#'   [geneda-SummarizedExperiment]).
#' @param nfeatures Number of features to use when HVGs are empty. Default =
#'   2000
#' @return Updated `geneda` object with `DimReduction` filled.
#'   If `object` is a `SummarizedExperiment`, the same object with the
#'   results added (see [geneda-SummarizedExperiment]).
#' @importFrom methods validObject
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
#' obj <- RunPCA(obj, nfeatures = 1000)
#' DimReduction(obj)$percent_var
RunPCA <- function(object, nfeatures = 2000) {
    if (methods::is(object, "SummarizedExperiment")) {
        return(.seRoundTrip(object, RunPCA, nfeatures))
    }
    object <- .toGeneda(object)

    if (length(object@HVGs) == 0L) {
        cli::cli_inform(
            paste0(
                "HVG slot is empty. Running FindVariableFeatures with top ",
                "{nfeatures} genes"
            )
        )
        object <- FindVariableFeatures(object, nfeatures = nfeatures)
    }
    nFeatUse <- length(object@HVGs)
    cli::cli_inform("Calculating principal components from top {nFeatUse} HVGs")
    object@DimReduction <- generatePCs(object@normalized, object@HVGs, nFeatUse)
    validObject(object)
    object
}

#' Extract PCA Scores and Metadata
#'
#' @description
#' This function extracts the PCA sample scores (`DimReduction$Scores`) of a
#' `geneda` object and combines them with the associated metadata. It ensures
#' that the metadata has valid rownames and aligns the PCA scores accordingly.
#'
#' @param object A `geneda` object (or a `SummarizedExperiment`; see
#'   [geneda-SummarizedExperiment]) containing PCA results in the `DimReduction`
#'   slot and sample-level metadata in the `metadata` slot.
#'
#' @details
#' The function performs several checks:
#' - Ensures the input object is of class `geneda`.
#' - Verifies that the `DimReduction` slot contains PCA scores.
#' - Confirms that the metadata has valid rownames.
#' - Reorders the PCA scores to match the order of metadata rows.
#'
#' If metadata rownames are missing or invalid, the function throws an error.
#'
#' @return
#' A `data.frame` combining PCA scores and sample metadata, where rows
#' correspond to samples and columns include principal component scores and
#' metadata fields.
#'
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
#' pca_df <- ExtractPCA(obj)
#' head(pca_df[, c("PC1", "PC2", "Genotype")])
ExtractPCA <- function(object) {
    object <- .toGeneda(object)

    if (length(object@DimReduction) == 0L) {
        cli::cli_abort("DimReduction slot is empty. Please use RunPCA.")
    }
    pcaRes <- object@DimReduction[["Scores"]]
    meta <- object@metadata
    if (is.null(rownames(meta)) || any(rownames(meta) == "")) {
        cli::cli_abort(
            paste0(
                "Metadata does not contain valid rownames. Please ensure ",
                "metadata rows are named."
            )
        )
    }
    order <- rownames(meta)
    pcaRes <- pcaRes[order, ]
    pcaRes <- cbind(pcaRes, meta)
    return(pcaRes)
}

#' extractEigen
#'
#' @description Extract gene-level eigenvector loadings for a PC of interest
#'   and calculate each gene's contribution to that PC.
#'
#' @details The contribution of gene \eqn{i} to a PC with eigenvector
#'   \eqn{v} is \eqn{100 \cdot v_i^2 / \sum_j v_j^2}, i.e. the gene's share
#'   of the PC's squared loadings. Because eigenvectors have unit length this
#'   equals \eqn{100 \cdot v_i^2}, and contributions sum to 100 within each
#'   PC. This is not the percent of total variance explained by the gene;
#'   for the variance explained by the PC itself see the `percent_var` entry
#'   of the `DimReduction` slot.
#'
#' @param object A `geneda` object (or a `SummarizedExperiment`; see
#'   [geneda-SummarizedExperiment]) containing PCA results in the `DimReduction`
#'   slot
#' @param component Principal component to use (i.e., PC1)
#'
#' @returns A data.frame with columns `Gene`, `EigenVector` (the gene's
#'   loading on the PC) and `Contribution` (the gene's percent contribution
#'   to the PC; sums to 100)
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
#' pc1 <- extractEigen(obj, component = "PC1")
#' head(pc1[order(-abs(pc1$EigenVector)), ])
extractEigen <- function(object, component) {
    object <- .toGeneda(object)

    if (length(object@DimReduction) == 0L) {
        cli::cli_abort("DimReduction slot is empty. Please use RunPCA.")
    }
    eigenvecs <- object@DimReduction[["Eigenvectors"]][[component]]
    geneNames <- rownames(object@DimReduction$Eigenvectors)

    # Gene contribution: share of the PC's squared loadings (sums to 100)
    contribution <- eigenvecs^2 / sum(eigenvecs^2) * 100

    # Save as dataframe
    df <- data.frame(
        Gene = geneNames,
        EigenVector = eigenvecs,
        Contribution = contribution
    )

    return(df)
}
