# optional: setClassUnion remains separate
setClassUnion("matrixOrNULL", c("matrix", "NULL"))

#' geneda S4 class
#'
#' @description An S4 container for exploratory genomic data analysis.
#' Stores optional raw counts, normalized data, sample metadata,
#' highly variable genes (HVGs), and PCA dimensionality reduction results.
#'
#' @slot counts Optional counts matrix (features x samples). Can be NULL.
#' @slot normalized Normalized expression matrix (features x samples).
#' @slot metadata Sample-level metadata `data.frame` (rows = samples).
#' @slot HVGs Character vector of selected highly variable gene IDs (row names).
#' @slot DimReduction List for PCA results with `Loadings`, `Eigenvectors`, `percent_var`.
#' @slot DEGs List container for differential expression results. Contains `DEG`
#'   (data.frame) for unfiltered results, and optionally named slots for filtered
#'   results (e.g., `DEGs$assay1`, `DEGs$assay2`).
#'
#' @exportClass geneda

setClass(
  "geneda",
  slots = list(
    counts = "matrixOrNULL",
    normalized = "matrix",
    metadata = "data.frame",
    HVGs = "character",
    DimReduction = "list",
    DEGs = "list"
  )
)


.valid_geneda <- function(object) {
  errors <- character()

  # Check dimensions and naming consistency
  if (is.null(colnames(object@normalized))) {
    errors <- c(errors, "'normalized' must have column names (sample IDs).")
  }
  if (is.null(rownames(object@normalized))) {
    errors <- c(errors, "'normalized' must have row names (feature IDs).")
  }
  if (nrow(object@metadata) == 0L) {
    errors <- c(errors, "'metadata' must have rows (samples).")
  }
  if (!is.null(colnames(object@normalized)) && !is.null(rownames(object@metadata))) {
    if (!identical(rownames(object@metadata), colnames(object@normalized))) {
      errors <- c(errors, "Row names of 'metadata' must match column names of 'normalized' (same order).")
    }
  }
  if (!is.null(object@counts)) {
    if (!is.matrix(object@counts)) {
      errors <- c(errors, "'counts' must be a matrix or NULL.")
    } else {
      if (!identical(colnames(object@counts), colnames(object@normalized))) {
        errors <- c(errors, "Column names of 'counts' must match 'normalized'.")
      }
      if (is.null(rownames(object@counts))) {
        errors <- c(errors, "'counts' must have row names (feature IDs).")
      }
    }
  }
  if (length(object@HVGs) > 0L) {
    missing_feats <- setdiff(object@HVGs, rownames(object@normalized))
    if (length(missing_feats) > 0L) {
      errors <- c(errors, "All 'HVGs' must be present in row names of 'normalized'.")
    }
  }
  if (length(object@DimReduction) > 0L) {
    needed_pcs <- c("Loadings", "Eigenvectors", "percent_var")
    if (!all(needed_pcs %in% names(object@DimReduction))) {
      errors <- c(errors, "'DimReduction' must contain names: Loadings, Eigenvectors, percent_var.")
    }
  }

  if (length(object@DEGs) > 0L) {
    if (!is.null(object@DEGs$DEG) && !is.data.frame(object@DEGs$DEG)) {
      errors <- c(errors, "'DEGs$DEG' must be a data.frame or NULL.")
    }
  }

  if (length(errors)) errors else TRUE
}

setValidity("geneda", .valid_geneda)

#' Construct a geneda object
#'
#' @param normalized Normalized expression matrix (features x samples).
#' @param metadata Sample metadata `data.frame` with row names matching `colnames(normalized)`.
#' @param counts Optional counts matrix (features x samples).
#' @param DEGs Optional DEG results from `DESeq2`
#' Note: This constructor does not compute HVGs or PCA. Use
#' `FindVariableFeatures()` and `RunPCA()` afterward.
#'
#' @return A `geneda` object.
#' @importFrom methods new
#' @export
GenEDA <- function(normalized,
                   metadata,
                   counts = NULL,
                   DEGs = NULL) {
  stopifnot(is.matrix(normalized))
  stopifnot(is.data.frame(metadata))
  if (!is.null(counts)) stopifnot(is.matrix(counts))

  # Align metadata rows to normalized sample order if identical set
  if (!is.null(rownames(metadata)) && !is.null(colnames(normalized))) {
    if (setequal(rownames(metadata), colnames(normalized)) &&
        !identical(rownames(metadata), colnames(normalized))) {
      metadata <- metadata[colnames(normalized), , drop = FALSE]
    }
  }

  new("geneda",
      counts = counts,
      normalized = normalized,
      metadata = metadata,
      HVGs = character(0L),
      DimReduction = list(),
      DEGs = list(DEG = NULL))
}
