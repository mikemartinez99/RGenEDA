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

#' generatePCs
#'
#' @description Generate principal component analysis data that can be used in downstream analyses.
#'
#' @param mat A data matrix where rows are features and columns are samples
#' @param vars A vector of gene variances (can calculate using RGenEDA::plotVariance)
#' @param nFeatures Number of top features to generate principal components on.
#'
#' @returns A list consisting of 3 slots: Loadings, Eigenvectors, and percent_var
#' @export
generatePCs <- function(mat, vars, nFeatures) {
  var_features_n <- nFeatures
  select <- order(vars, decreasing = TRUE)[1:var_features_n]
  vsd_sub <- mat[select,]
  vsd_sub <- t(vsd_sub)
  pca <- prcomp(vsd_sub)

  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  percentVar <- percentVar[1:10]

  message("Percent variations:")
  percentVar <- paste(round(percentVar*100,2), "%", sep = " ")
  names(percentVar) <- c("PC1", "PC2", "PC3", "PC4", "PC5",
                         "PC6", "PC7", "PC8", "PC9", "PC10")
  print(percentVar)

  pca_res <- list()
  pca_df <- as.data.frame(pca$x)
  pca_eigenvecs <- as.data.frame(pca$rotation[,1:10])
  pca_res[["Loadings"]] <- pca_df
  pca_res[["Eigenvectors"]] <- pca_eigenvecs
  pca_res[["percent_var"]] <- percentVar
  return(pca_res)
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

#' extractEigen
#'
#' @description Extract gene-level Eigen-vectors for PCs of interest and calculate percent variance.
#'
#' @param object A `geneda` object containing PCA results in the `DimReduction` slot
#' @param component Principal component to use (i.e., PC1)
#'
#' @returns A data-frame containing gene loading information and percent variation per gene
#' @export
#'
#' @examples
#' \donttest{
#' mock_norm <- matrix(rnorm(12, mean = 0, sd = 2), nrow = 4, ncol = 3)
#' colnames(mock_norm) <- paste0("Sample", 1:3)
#' rownames(mock_norm) <- paste0("Gene", 1:4)
#' mock_meta <- data.frame(condition = c("A","B","A"), row.names = colnames(mock_norm))
#' obj <- GenEDA(normalized = mock_norm, metadata = mock_meta)
#' extractEigen(object = obj, component = "PC1")
#'}
extractEigen <- function(object, component) {
  stopifnot(methods::is(object, "geneda"))

  if (length(object@DimReduction) == 0L) {
    message("DimReduction slot is empty. Please use RunPCA.")
  }
  eigenvecs <- object@DimReduction[["Eigenvectors"]][[component]]
  percentVar <- object@DimReduction[["percent_var"]][[component]]
  geneNames <- rownames(object@DimReduction$Eigenvectors)

  # Get gene-level variance
  pcVar <- as.numeric(gsub("%", "", percentVar))/100
  genePercent <- (eigenvecs^2 * pcVar) / sum(eigenvecs^2 * pcVar) * 100

  # Save as dataframe
  df <- data.frame(
    Gene = geneNames,
    EigenVector = eigenvecs,
    PctVar = genePercent
  )

  return(df)

}
