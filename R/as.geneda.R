#' Convert a SummarizedExperiment to a geneda object
#'
#' @description
#' Builds a `geneda` object from a `SummarizedExperiment`, including
#' `DESeqDataSet` objects from DESeq2 and `DESeqTransform` objects returned by
#' `DESeq2::vst()` or `DESeq2::rlog()`. Raw counts, normalized values, and
#' sample metadata (`colData`) are carried over, so the object is ready for
#' `FindVariableFeatures()`, `RunPCA()`, and the plotting functions.
#'
#' @param x A `SummarizedExperiment` (or subclass such as `DESeqDataSet` or
#'   `DESeqTransform`) with features as rows and samples as columns.
#' @param normalized Where the normalized values come from. One of:
#'   * `NULL` (default): computed from the counts assay using `transform`.
#'     If `x` is a `DESeqTransform`, its transformed values are used instead.
#'   * A single assay name or index in `x`.
#'   * A `SummarizedExperiment` with the same features and samples as `x`,
#'     for example the output of `DESeq2::vst(x)`. Its first assay is used.
#'   * A matrix with the same row and column names as `x`.
#' @param counts Name or index of the raw counts assay in `x`, or `NULL` to
#'   leave the `counts` slot empty. Default `"counts"`. If the assay is not
#'   present, the `counts` slot is left empty, unless counts are needed to
#'   compute `normalized`.
#' @param transform Normalization used when `normalized = NULL`: `"vst"`
#'   (variance stabilizing transformation from DESeq2, `blind = TRUE`) or
#'   `"logcpm"` (log2 counts per million with a pseudocount of 1). Default
#'   `"vst"`, which requires the DESeq2 package and integer counts.
#' @param DEGs Optional named list of DEG tables (`data.frame` or
#'   `DESeqResults`) to store with `SetDEGs()`, using the list names as assay
#'   names.
#' @param ... Not used.
#'
#' @details
#' Sample metadata are taken from `colData(x)` and converted to a
#' `data.frame`. Columns added by DESeq2, such as `sizeFactor`, are kept and
#' can be tested with `PlotEigenCorr()` like any other covariate.
#'
#' Assays stored as sparse or delayed matrices are converted to base
#' matrices with `as.matrix()`.
#'
#' @return A `geneda` object.
#'
#' @seealso [GenEDA()] to build a `geneda` object from matrices.
#'
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment assay assayNames colData
#' @export
#' @rdname as.geneda
#'
#' @examples
#' data(Snail1KO)
#' counts <- as.matrix(Snail1KO$rawCounts)
#' counts <- counts[rowSums(counts >= 10) >= 3, ]
#' meta <- Snail1KO$metadata
#' rownames(meta) <- meta$Sample
#'
#' se <- SummarizedExperiment::SummarizedExperiment(
#'     assays = list(counts = counts),
#'     colData = meta
#' )
#' obj <- as.geneda(se, transform = "logcpm")
#' obj
#'
#' # From a DESeqDataSet, with VST-normalized values and DEGs
#' if (requireNamespace("DESeq2", quietly = TRUE)) {
#'     meta$Genotype <- factor(meta$Genotype, levels = c("WT", "Snai1_KO"))
#'     dds <- DESeq2::DESeqDataSetFromMatrix(counts, meta, design = ~Genotype)
#'     dds <- DESeq2::DESeq(dds, quiet = TRUE)
#'     obj <- as.geneda(dds,
#'         normalized = DESeq2::vst(dds),
#'         DEGs = list(KO_vs_WT = DESeq2::results(dds))
#'     )
#'     obj
#' }
setGeneric("as.geneda", function(x, ...) standardGeneric("as.geneda"))

#' @rdname as.geneda
setMethod("as.geneda", "SummarizedExperiment", function(
    x,
    normalized = NULL,
    counts = "counts",
    transform = c("vst", "logcpm"),
    DEGs = NULL,
    ...) {
    transform <- match.arg(transform)

    if (is.null(colnames(x))) {
        cli::cli_abort("{.arg x} must have column names (sample IDs).")
    }
    if (is.null(rownames(x))) {
        cli::cli_abort("{.arg x} must have row names (feature IDs).")
    }

    countMat <- .seCounts(x, counts)
    normMat <- .resolveNormalized(x, normalized, countMat, transform)

    metadata <- as.data.frame(colData(x))
    rownames(metadata) <- colnames(x)

    GenEDA(
        normalized = normMat,
        metadata = metadata,
        counts = countMat,
        DEGs = DEGs
    )
})

#' Choose or compute the normalized matrix for as.geneda()
#' @keywords internal
#' @noRd
.resolveNormalized <- function(x, normalized, countMat, transform) {
    if (!is.null(normalized)) {
        return(.seNormalized(x, normalized))
    }
    # A DESeqTransform already holds normalized values in its first assay
    if (methods::is(x, "DESeqTransform")) {
        return(as.matrix(assay(x, 1L)))
    }
    if (is.null(countMat)) {
        cli::cli_abort(c(
            "No counts assay found to compute normalized values from.",
            "i" = paste0(
                "Set {.arg counts} to the raw counts assay, or pass ",
                "{.arg normalized}."
            )
        ))
    }
    .normalizeCounts(countMat, transform)
}

#' Store a named list of DEG tables with SetDEGs()
#' @keywords internal
#' @noRd
.addDEGList <- function(obj, DEGs) {
    if (is.null(DEGs)) {
        return(obj)
    }
    if (!is.list(DEGs) || is.data.frame(DEGs) || is.null(names(DEGs)) ||
        any(names(DEGs) == "")) {
        cli::cli_abort("{.arg DEGs} must be a named list of DEG tables.")
    }
    for (nm in names(DEGs)) {
        obj <- SetDEGs(obj, as.data.frame(DEGs[[nm]]), assay = nm)
    }
    obj
}

#' Pull the raw counts assay from a SummarizedExperiment
#' @keywords internal
#' @noRd
.seCounts <- function(x, counts) {
    if (is.null(counts)) {
        return(NULL)
    }
    available <- if (is.character(counts)) {
        counts %in% assayNames(x)
    } else {
        counts <= length(SummarizedExperiment::assays(x))
    }
    if (!available) {
        return(NULL)
    }
    as.matrix(assay(x, counts))
}

#' Resolve the `normalized` argument of as.geneda() to a matrix
#' @keywords internal
#' @noRd
.seNormalized <- function(x, normalized) {
    if (methods::is(normalized, "SummarizedExperiment")) {
        normMat <- as.matrix(assay(normalized, 1L))
    } else if (is.matrix(normalized) || methods::is(normalized, "Matrix")) {
        normMat <- as.matrix(normalized)
    } else if (length(normalized) == 1L &&
        (is.character(normalized) || is.numeric(normalized))) {
        if (is.character(normalized) && !normalized %in% assayNames(x)) {
            cli::cli_abort("Assay {.val {normalized}} not found in {.arg x}.")
        }
        normMat <- as.matrix(assay(x, normalized))
    } else {
        cli::cli_abort(paste0(
            "{.arg normalized} must be NULL, an assay name, a ",
            "SummarizedExperiment, or a matrix."
        ))
    }

    if (!identical(colnames(normMat), colnames(x))) {
        cli::cli_abort(
            "Column names of {.arg normalized} must match those of {.arg x}."
        )
    }
    if (!identical(rownames(normMat), rownames(x))) {
        cli::cli_abort(
            "Row names of {.arg normalized} must match those of {.arg x}."
        )
    }
    normMat
}

#' Normalize a raw counts matrix for as.geneda()
#' @keywords internal
#' @noRd
.normalizeCounts <- function(countMat, transform) {
    if (transform == "logcpm") {
        cli::cli_inform("Normalizing counts with log2(CPM + 1).")
        libSize <- colSums(countMat)
        return(log2(t(t(countMat) / libSize) * 1e6 + 1))
    }

    if (!requireNamespace("DESeq2", quietly = TRUE)) {
        cli::cli_abort(c(
            "{.code transform = \"vst\"} requires the {.pkg DESeq2} package.",
            "i" = "Install DESeq2, or use {.code transform = \"logcpm\"}."
        ))
    }
    if (any(countMat != round(countMat), na.rm = TRUE)) {
        cli::cli_abort(c(
            "{.code transform = \"vst\"} requires integer counts.",
            "i" = paste0(
                "Round the counts first, or use ",
                "{.code transform = \"logcpm\"}."
            )
        ))
    }
    cli::cli_inform("Normalizing counts with DESeq2 vst (blind = TRUE).")
    storage.mode(countMat) <- "integer"
    # vst() subsamples 1000 genes to fit the trend; use the full
    # transformation for smaller matrices
    if (nrow(countMat) < 1000L) {
        DESeq2::varianceStabilizingTransformation(countMat, blind = TRUE)
    } else {
        DESeq2::vst(countMat, blind = TRUE)
    }
}
