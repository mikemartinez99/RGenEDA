#' geneda S4 class and helpers
#'
#' @description An S4 container for exploratory genomic data analysis.
#' Stores optional raw counts, normalized data, sample metadata, a place for
#' selected highly variable genes (HVGs), and a place for PCA results.
#'
#' @slot counts Optional counts matrix (features x samples). Can be NULL.
#' @slot normalized Normalized expression matrix (features x samples).
#' @slot metadata Sample-level metadata `data.frame` (rows = samples).
#' @slot HVGs Character vector of selected highly variable gene IDs (row names).
#' @slot DimReduction List for PCA results with `Loadings`, `Eigenvectors`, `percent_var`.
#'
#' @exportClass geneda
setClassUnion("matrixOrNULL", c("matrix", "NULL"))

setClass(
  "geneda",
  slots = list(
    counts = "matrixOrNULL",
    normalized = "matrix",
    metadata = "data.frame",
    HVGs = "character",
    DimReduction = "list"
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

  if (length(errors)) errors else TRUE
}

setValidity("geneda", .valid_geneda)

#' Construct a geneda object
#'
#' @param normalized Normalized expression matrix (features x samples).
#' @param metadata Sample metadata `data.frame` with row names matching `colnames(normalized)`.
#' @param counts Optional counts matrix (features x samples).
#' Note: This constructor does not compute HVGs or PCA. Use
#' `FindVariableFeatures()` and `RunPCA()` afterward.
#'
#' @return A `geneda` object.
#' @export
GenEDA <- function(normalized,
                   metadata,
                   counts = NULL) {
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
      DimReduction = list())
}

#' Access counts matrix
#' @param object A `geneda` object
#' @return Matrix or NULL
#' @export
counts <- function(object) {
  object@counts
}

#' Access normalized matrix
#' @param object A `geneda` object
#' @return Matrix
#' @export
normalized <- function(object) {
  object@normalized
}

#' Access metadata
#' @param object A `geneda` object
#' @return data.frame
#' @export
metadata <- function(object) {
  object@metadata
}

#' Access HVG IDs
#' @param object A `geneda` object
#' @return Character vector
#' @export
HVGs <- function(object) {
  object@HVGs
}

#' Access dimensional reduction results list
#' @param object A `geneda` object
#' @return List with `Loadings`, `Eigenvectors`, `percent_var`
#' @export
DimReduction <- function(object) {
  object@DimReduction
}

#' Find and store highly variable genes (HVGs)
#'
#' Computes per-feature variance on `normalized(object)`, ranks genes by
#' decreasing variance, and stores the top `nfeatures` gene IDs into the `HVGs`
#' slot. Returns the updated object.
#'
#' @param object A `geneda` object
#' @param nfeatures Number of HVGs to store
#' @return Updated `geneda` object
#' @export
FindVariableFeatures <- function(object, nfeatures) {
  stopifnot(methods::is(object, "geneda"))
  nfeatures <- as.integer(nfeatures)

  vars <- apply(object@normalized, 1L, stats::var, na.rm = TRUE)
  if (is.null(names(vars))) names(vars) <- rownames(object@normalized)
  vars <- sort(vars, decreasing = TRUE)
  nfeatures <- max(1L, min(nfeatures, length(vars)))
  object@HVGs <- names(vars)[seq_len(nfeatures)]
  validObject(object)
  object
}

#' Run PCA and store in `DimReduction`
#'
#' Uses `generatePCs()` under the hood. If `nfeatures` is NULL and `HVGs` are
#' present, uses the length of `HVGs`. Otherwise uses the provided `nfeatures`.
#'
#' @param object A `geneda` object
#' @param nfeatures Optional number of features to use for PCA
#' @return Updated `geneda` object with `DimReduction` filled
#' @export
RunPCA <- function(object, nfeatures = NULL) {
  stopifnot(methods::is(object, "geneda"))
  vars <- apply(object@normalized, 1L, stats::var, na.rm = TRUE)
  if (is.null(names(vars))) names(vars) <- rownames(object@normalized)

  if (is.null(nfeatures)) {
    if (length(object@HVGs) > 0L) {
      nfeatures <- length(object@HVGs)
    } else {
      stop("Provide 'nfeatures' or run FindVariableFeatures() first.")
    }
  }
  nfeatures <- as.integer(nfeatures)
  nfeatures <- max(1L, min(nfeatures, length(vars)))

  object@DimReduction <- generatePCs(MAT = object@normalized, VARS = vars, NFEATURES = nfeatures)
  validObject(object)
  object
}

#' Show method for geneda
#' @param object A `geneda` object
setMethod(
  "show",
  signature = "geneda",
  definition = function(object) {
    cat("geneda object\n")
    cat(sprintf("  features: %d\n", nrow(object@normalized)))
    cat(sprintf("  samples:  %d\n", ncol(object@normalized)))
    cat(sprintf("  HVGs: %d\n", length(object@HVGs)))
    cat(sprintf("  counts: %s\n", if (is.null(object@counts)) "NULL" else "present"))
    invisible(object)
  }
)


