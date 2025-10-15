#' geneda S4 class and helpers
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
#' @slot DEGs List container for differential expression results with two entries:
#'   `unfiltered` (data.frame) and `filtered` (data.frame).
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
    needed_degs <- c("unfiltered", "filtered")
    if (!all(needed_degs %in% names(object@DEGs))) {
      errors <- c(errors, "'DEGs' must contain names: unfiltered, filtered.")
    }
    if (!is.null(object@DEGs$unfiltered) && !is.data.frame(object@DEGs$unfiltered)) {
      errors <- c(errors, "'DEGs$unfiltered' must be a data.frame or NULL.")
    }
    if (!is.null(object@DEGs$filtered) && !is.data.frame(object@DEGs$filtered)) {
      errors <- c(errors, "'DEGs$filtered' must be a data.frame or NULL.")
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
      DimReduction = list(),
      DEGs = list(unfiltered = NULL, filtered = NULL))
}

#' Access counts matrix (avoid name clash with DESeq2::counts)
#' @param object A `geneda` object
#' @return Matrix or NULL
#' @export
getCounts <- function(object) {
  object@counts
}

#' Access normalized matrix (avoid name clash with generics)
#' @param object A `geneda` object
#' @return Matrix
#' @export
getNormalized <- function(object) {
  object@normalized
}

#' Access metadata (avoid name clash with S4 generics)
#' @param object A `geneda` object
#' @return data.frame
#' @export
getMetadata <- function(object) {
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

#' Access DEGs container
#' @param object A `geneda` object
#' @return List with `unfiltered` and `filtered` data.frames
#' @export
DEGs <- function(object) {
  object@DEGs
}

#' Set unfiltered DEGs on the object
#'
#' @param object A `geneda` object
#' @param deg_table A data.frame of DESeq2-like results containing at least
#'   columns `log2FoldChange` and `padj`
#' @return Updated `geneda` object with `DEGs$unfiltered` set and `DEGs$filtered` cleared
#' @export
SetDEGs <- function(object, deg_table) {
  stopifnot(methods::is(object, "geneda"))
  stopifnot(is.data.frame(deg_table))
  req_cols <- c("log2FoldChange", "padj")
  missing_cols <- setdiff(req_cols, colnames(deg_table))
  if (length(missing_cols) > 0L) {
    stop(paste0("DEG table must contain columns: ", paste(req_cols, collapse = ", ")))
  }
  object@DEGs$unfiltered <- deg_table
  object@DEGs$filtered <- NULL
  validObject(object)
  object
}

#' Filter DEGs by padj and absolute log2FoldChange
#'
#' @param object A `geneda` object with `DEGs$unfiltered` set
#' @param padj_thresh Adjusted p-value threshold (<=)
#' @param log2FC_thresh Absolute log2 fold change threshold (>=)
#' @return Updated `geneda` object with `DEGs$filtered` set
#' @export
FilterDEGs <- function(object, padj_thresh = 0.05, log2FC_thresh = 1.0) {
  stopifnot(methods::is(object, "geneda"))
  if (is.null(object@DEGs$unfiltered)) {
    stop("DEGs$unfiltered is NULL. Use SetDEGs(object, deg_table) first.")
  }
  df <- object@DEGs$unfiltered
  req_cols <- c("log2FoldChange", "padj")
  missing_cols <- setdiff(req_cols, colnames(df))
  if (length(missing_cols) > 0L) {
    stop(paste0("DEG table must contain columns: ", paste(req_cols, collapse = ", ")))
  }
  filt <- stats::complete.cases(df$log2FoldChange, df$padj) &
    abs(df$log2FoldChange) >= log2FC_thresh & df$padj <= padj_thresh
  object@DEGs$filtered <- df[filt, , drop = FALSE]
  validObject(object)
  object
}

#' Run PCA and store in `DimReduction`
#'
#' Uses `generatePCs()` under the hood.
#' If `HVGs` are empty, selects HVGs via `FindVariableFeatures(object, nfeatures)`
#' with default `nfeatures = 2000`. If `HVGs` are present, the PCA uses
#' `NFEATURES = length(HVGs(object))`. This aligns the PCA feature count with the
#' HVG selection while allowing an explicit override of `nfeatures` when empty.
#'
#' @param object A `geneda` object
#' @param nfeatures Number of features to use when HVGs are empty. Default = 2000
#' @return Updated `geneda` object with `DimReduction` filled
#' @export
RunPCA <- function(object, nfeatures = 2000) {
  stopifnot(methods::is(object, "geneda"))

  if (length(object@HVGs) == 0L){
    message("HVG slot is empty. Running FindVariableFeatures with top 2000 genes")
    object <- FindVariableFeatures(object, nfeatures = nfeatures)
  }
  nFeatUse <- length(object@HVGs)
  message(paste("Calculating principal components from top", nFeatUse, "HVGs"))
  object@DimReduction <- generatePCs(object@normalized, object@HVGs, nFeatUse)
  validObject(object)
  object
}

#' Extract PCA Loadings and Metadata
#'
#' @description
#' This function extracts PCA loadings stored in the `DimReduction` slot of a
#' `geneda` object and combines them with the associated metadata. It ensures
#' that the metadata has valid rownames and aligns the PCA loadings accordingly.
#'
#' @param object A `geneda` object containing PCA results in the `DimReduction` slot
#'   and sample-level metadata in the `metadata` slot.
#'
#' @details
#' The function performs several checks:
#' - Ensures the input object is of class `geneda`.
#' - Verifies that the `DimReduction` slot contains PCA loadings.
#' - Confirms that the metadata has valid rownames.
#' - Reorders the PCA loadings to match the order of metadata rows.
#'
#' If metadata rownames are missing or invalid, the function throws an error.
#'
#' @return
#' A `data.frame` combining PCA loadings and sample metadata, where rows correspond
#' to samples and columns include principal component loadings and metadata fields.
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' pca_results <- ExtractPCA(my_geneda_object)
#' head(pca_results)
#' }
#'
#' @seealso [RunPCA()], [extractLoadings()]
#' @export
ExtractPCA <- function(object) {
  stopifnot(methods::is(object, "geneda"))

  if (length(object@DimReduction) == 0L) {
    message("DimReduction slot is empty. Please use RunPCA.")
  }
  pcaRes <- object@DimReduction[["Loadings"]]
  meta <- object@metadata
  if (is.null(rownames(meta)) || any(rownames(meta) == "")) {
    stop("Metadata does not contain valid rownames. Please ensure metadata rows are named.")
  }
  order <- rownames(meta)
  pcaRes <- pcaRes[order,]
  pcaRes <- cbind(pcaRes, meta)
  return(pcaRes)
}



#' plotHVGVariance
#'
#' @description Plot variance for all genes (features) on a `geneda` object using ggplot2.
#'
#' @param object A `geneda` object created by `GenEDA()`.
#' @param transform One of c("log") to log-transform variance, or NULL for none. Default = NULL
#' @param dropTopN Number of top-most variable genes to drop prior to plotting. Default = 0
#'
#' @returns A ggplot2 object visualizing the variance of HVGs
#' @export
#'
#' @examples
#' # obj <- GenEDA(vsd, meta)
#' # p <- plotHVGVariance(obj, transform = "log")
#' # print(p)
plotHVGVariance <- function(object, transform = NULL, dropTopN = 0) {
  stopifnot(methods::is(object, "geneda"))

  mat <- getNormalized(object)
  vars <- apply(mat, 1L, var)
  if (is.null(names(vars))) names(vars) <- rownames(mat)
  vars <- sort(vars, decreasing = TRUE)

  if (!is.null(dropTopN) && dropTopN > 0) {
    vars <- vars[-seq_len(min(dropTopN, length(vars)))]
  }

  df <- data.frame(
    Index = seq_along(vars),
    Variance = as.numeric(vars)
  )

  if (!is.null(transform)) {
    transform <- match.arg(transform, c("log"))
  }
  if (!is.null(transform) && transform == "log") {
    df$Variance <- log(df$Variance)
  }

  ymax <- max(df$Variance, na.rm = TRUE)

  lines_df <- data.frame(
    x = c(2000, 1000, 500, 250),
    label = factor(c("2000", "1000", "500", "250"), levels = c("2000", "1000", "500", "250"))
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = Index, y = Variance)) +
    ggplot2::geom_point(alpha = 0.6, size = 1.5) +
    ggplot2::geom_vline(data = lines_df, ggplot2::aes(xintercept = x, color = label), linetype = "dashed", show.legend = FALSE) +
    # Dummy horizontal segments for legend keys only (placed off-plot)
    ggplot2::geom_segment(
      data = lines_df,
      ggplot2::aes(x = 0, xend = 1, y = -Inf, yend = -Inf, color = label),
      inherit.aes = FALSE,
      show.legend = TRUE,
      linewidth = 0.8,
      linetype = "solid"
    ) +
    ggplot2::scale_color_manual(
      name = "Num. HVGs",
      values = c(
        "250" = "#7570B3",
        "500" = "#009E73",
        "1000" = "#E69F00",
        "2000" = "#D55E00"
      )
    ) +
    ggplot2::labs(
      title = "Gene Expression Variance",
      x = "nGenes",
      y = ifelse(!is.null(transform) && transform == "log", "Log(Variance)", "Variance")
    ) +
    ggplot2::theme_classic(base_size = 16) +
    ggplot2::theme(
      axis.title = ggplot2::element_text(face = "bold"),
      legend.position = c(0.98, 0.98),
      legend.justification = c(1, 1),
      legend.direction = "horizontal",
      legend.title = ggplot2::element_text(face = "bold"),
      legend.box.background = ggplot2::element_rect(color = "black", fill = ggplot2::alpha("white", 0.95)),
      legend.background = ggplot2::element_rect(fill = ggplot2::alpha("white", 0))
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(linetype = "solid", linewidth = 1.0)))

  return(p)
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


