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

#' Access DEGs container
#' @param object A `geneda` object
#' @param assay The name of the DEGs slot to return
#' @return List with `DEG` (unfiltered data.frame) and optionally named filtered result sets
#' @export
DEGs <- function(object, assay) {
  object@DEGs[[assay]]
}

#' Set unfiltered DEGs on the object
#'
#' @param object A `geneda` object
#' @param deg_table A data.frame of DESeq2-like results containing at least
#'   columns `log2FoldChange` and `padj`
#' @return Updated `geneda` object with `DEGs$DEG` set
#' @export
SetDEGs <- function(object, deg_table) {
  stopifnot(methods::is(object, "geneda"))
  stopifnot(is.data.frame(deg_table))
  req_cols <- c("log2FoldChange", "padj")
  missing_cols <- setdiff(req_cols, colnames(deg_table))
  if (length(missing_cols) > 0L) {
    stop(paste0("DEG table must contain columns: ", paste(req_cols, collapse = ", ")))
  }
  if (is.null(rownames(deg_tale))) {
    stop("DEG table does not have genes as rownames!")
  }
  object@DEGs$DEG <- deg_table
  validObject(object)
  object
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


#' Filter DEGs by padj and absolute log2FoldChange
#'
#' Filters the unfiltered DEGs in `DEGs$DEG` and stores the filtered results
#' in a new named slot `DEGs[[assayName]]`. Multiple filtered result sets can
#' be stored with different assay names.
#'
#' @param object A `geneda` object with `DEGs$DEG` set
#' @param padj_thresh Adjusted p-value threshold (<=)
#' @param log2FC_thresh Absolute log2 fold change threshold (>=)
#' @param assayName Character name for the filtered result set (e.g., "padj05_lfc1")
#' @return Updated `geneda` object with filtered results stored in `DEGs[[assayName]]`
#' @importFrom methods validObject
#' @export
FilterDEGs <- function(object, padj_thresh = 0.05, log2FC_thresh = 1.0, assayName) {
  stopifnot(methods::is(object, "geneda"))
  if (is.null(object@DEGs$DEG)) {
    stop("DEGs$DEG is NULL. Use SetDEGs(object, deg_table) first.")
  }
  if (missing(assayName) || is.null(assayName) || !is.character(assayName) || length(assayName) != 1L) {
    stop("assayName must be a single character string.")
  }
  if (assayName == "DEG") {
    stop("assayName cannot be 'DEG' (reserved for unfiltered results).")
  }
  df <- object@DEGs$DEG
  req_cols <- c("log2FoldChange", "padj")
  missing_cols <- setdiff(req_cols, colnames(df))
  if (length(missing_cols) > 0L) {
    stop(paste0("DEG table must contain columns: ", paste(req_cols, collapse = ", ")))
  }
  filt <- stats::complete.cases(df$log2FoldChange, df$padj) &
    abs(df$log2FoldChange) >= log2FC_thresh & df$padj <= padj_thresh
  object@DEGs[[assayName]] <- df[filt, , drop = FALSE]
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
#' @importFrom methods validObject
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
