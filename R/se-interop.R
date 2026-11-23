## Interoperability between `geneda` and `SummarizedExperiment`
##
## Layout of RGenEDA results inside a SummarizedExperiment:
##   assays           "counts" (raw counts, when the rows match) and
##                    "normalized" (the matrix used by RGenEDA)
##   colData          sample metadata
##   rowData          `isHVG` (logical) and `hvgRank` (integer, NA if not HVG)
##   metadata()$RGenEDA
##                    list(HVGs, DimReduction, DEGs, counts, version)
##                    `counts` is only used when the counts matrix has
##                    different rows from the normalized matrix.
## as.geneda() reads this layout back, so geneda -> SE -> geneda is lossless.

#' Convert a geneda object to a SummarizedExperiment
#'
#' @description
#' `as(object, "SummarizedExperiment")` stores everything in a `geneda`
#' object in a standard `SummarizedExperiment`, so the results can be used
#' with other Bioconductor tools. [as.geneda()] reads the same layout back,
#' so converting to a `SummarizedExperiment` and back gives the same
#' `geneda` object.
#'
#' @details
#' The `SummarizedExperiment` holds:
#' \describe{
#'   \item{`assays`}{`normalized` (the normalized matrix) and `counts` (the
#'     raw counts), if they have the same features as `normalized`.}
#'   \item{`colData`}{the sample metadata.}
#'   \item{`rowData`}{`isHVG` (logical: is the gene a highly variable gene)
#'     and `hvgRank` (rank by variance, `NA` for genes that are not HVGs).}
#'   \item{`metadata(se)$RGenEDA`}{a list with `HVGs`, `DimReduction` (PCA
#'     scores, loadings and percent variance), `DEGs` (the stored DEG
#'     tables), `counts` (only if the counts have different features from
#'     `normalized`) and `version` (the RGenEDA version).}
#' }
#'
#' All main RGenEDA functions also accept a `SummarizedExperiment` directly.
#' Functions that add results ([FindVariableFeatures()], [RunPCA()],
#' [SetDEGs()], [FilterDEGs()]) then return the input `SummarizedExperiment`
#' (including subclasses such as `DESeqDataSet`) with the results added in
#' this layout.
#'
#' @name geneda-SummarizedExperiment
#' @aliases coerce,geneda,SummarizedExperiment-method
#' @exportMethod coerce
#' @importFrom methods coerce
#'
#' @seealso [as.geneda()] for the reverse conversion.
#'
#' @return A `SummarizedExperiment`.
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
#'
#' se <- as(obj, "SummarizedExperiment")
#' se
#' head(SummarizedExperiment::rowData(se))
#' names(S4Vectors::metadata(se)$RGenEDA)
#'
#' # Round trip back to a geneda object
#' obj2 <- as.geneda(se)
#' identical(DimReduction(obj2), DimReduction(obj))
NULL

setAs("geneda", "SummarizedExperiment", function(from) {
    assays <- list(normalized = from@normalized)
    if (.countsMatchRows(from)) {
        assays <- c(list(counts = from@counts), assays)
    }
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = assays,
        colData = S4Vectors::DataFrame(from@metadata, check.names = FALSE)
    )
    .storeGeneda(se, from)
})

#' TRUE if the counts slot can be stored as an assay next to `normalized`
#' @keywords internal
#' @noRd
.countsMatchRows <- function(obj) {
    !is.null(obj@counts) &&
        identical(rownames(obj@counts), rownames(obj@normalized)) &&
        identical(colnames(obj@counts), colnames(obj@normalized))
}

#' Write the RGenEDA results of `obj` into the SummarizedExperiment `se`
#'
#' Adds (or replaces) the `normalized` assay, the `isHVG` / `hvgRank`
#' rowData columns and `metadata(se)$RGenEDA`. Everything else in `se`
#' (other assays, rowRanges, colData, the class) is kept.
#'
#' @keywords internal
#' @noRd
.storeGeneda <- function(se, obj) {
    norm <- obj@normalized
    # Subset only when needed: subsetting drops extra matrix attributes
    # (e.g. those DESeq2::rlog() adds), which would break the round trip
    if (!identical(rownames(norm), rownames(se)) ||
        !identical(colnames(norm), colnames(se))) {
        norm <- norm[rownames(se), colnames(se), drop = FALSE]
    }
    # A DESeqTransform already holds the normalized values as its first
    # (unnamed) assay, and as.geneda() reads them from there: leave it as is
    isTransform <- methods::is(se, "DESeqTransform") &&
        isTRUE(all.equal(
            unname(as.matrix(SummarizedExperiment::assay(se, 1L))),
            unname(norm)
        ))
    if (!isTransform) {
        SummarizedExperiment::assay(se, "normalized", withDimnames = FALSE) <-
            .stripDimnames(norm)
    }

    rank <- match(rownames(se), obj@HVGs)
    rd <- SummarizedExperiment::rowData(se)
    rd$isHVG <- !is.na(rank)
    rd$hvgRank <- rank
    SummarizedExperiment::rowData(se) <- rd

    md <- S4Vectors::metadata(se)
    md$RGenEDA <- list(
        HVGs = obj@HVGs,
        DimReduction = obj@DimReduction,
        DEGs = obj@DEGs,
        counts = if (.countsMatchRows(obj)) NULL else obj@counts,
        version = as.character(utils::packageVersion("RGenEDA"))
    )
    S4Vectors::metadata(se) <- md
    se
}

#' Restore RGenEDA results stored by .storeGeneda() into a geneda object
#'
#' Results that no longer fit `obj` (for example after the
#' SummarizedExperiment was subset) are dropped with a warning.
#'
#' @keywords internal
#' @noRd
.restoreGeneda <- function(obj, stored) {
    if (is.null(stored)) {
        return(obj)
    }
    if (is.null(obj@counts) && !is.null(stored$counts)) {
        obj@counts <- stored$counts
    }

    hvgs <- stored$HVGs
    dr <- stored$DimReduction
    if (length(hvgs) > 0L &&
        !all(hvgs %in% rownames(obj@normalized))) {
        cli::cli_warn(paste0(
            "Stored HVGs and PCA results do not match the features of the ",
            "SummarizedExperiment and were dropped. Re-run ",
            "FindVariableFeatures() and RunPCA()."
        ))
        hvgs <- character(0L)
        dr <- list()
    }
    if (length(dr) > 0L &&
        !identical(rownames(dr$Scores), colnames(obj@normalized))) {
        cli::cli_warn(paste0(
            "Stored PCA results do not match the samples of the ",
            "SummarizedExperiment and were dropped. Re-run RunPCA()."
        ))
        dr <- list()
    }
    obj@HVGs <- if (is.null(hvgs)) character(0L) else hvgs
    obj@DimReduction <- if (is.null(dr)) list() else dr

    for (nm in names(stored$DEGs)) {
        obj@DEGs[[nm]] <- stored$DEGs[[nm]]
    }
    methods::validObject(obj)
    obj
}

#' Get a geneda object from a geneda or SummarizedExperiment input
#'
#' A SummarizedExperiment is converted with as.geneda() defaults: the
#' "normalized" assay if present, otherwise normalized values computed from
#' the "counts" assay.
#'
#' @keywords internal
#' @noRd
.toGeneda <- function(object) {
    if (methods::is(object, "geneda")) {
        return(object)
    }
    if (methods::is(object, "SummarizedExperiment")) {
        return(as.geneda(object))
    }
    cli::cli_abort(paste0(
        "{.arg object} must be a {.cls geneda} object or a ",
        "{.cls SummarizedExperiment}, not {.cls {class(object)}}."
    ))
}

#' Run a geneda-modifying function on a SummarizedExperiment
#'
#' Converts `se` to geneda, calls `FUN`, and writes the results back into
#' `se`, so the caller gets its own SummarizedExperiment (and subclass) back.
#'
#' @keywords internal
#' @noRd
.seRoundTrip <- function(se, FUN, ...) {
    .storeGeneda(se, FUN(as.geneda(se), ...))
}

#' Remove dimnames but keep all other attributes of a matrix
#' @keywords internal
#' @noRd
.stripDimnames <- function(m) {
    dimnames(m) <- NULL
    m
}
