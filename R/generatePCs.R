#' generatePCs
#'
#' @description Generate principal component analysis data that can be used in downstream analyses.
#'
#' @param MAT A data matrix where rows are features and columns are samples
#' @param VARS A vector of gene variances (can calculate using DACExploreR::plotVariance)
#' @param NFEATURES Number of top features to generate principal components on.
#'
#' @returns A list consisting of 3 slots: Loadings, Eigenvectors, and percent_var
#'
#' @examples # Calculate variance
#' @examples vars <- DACExplorer::plotVariance(matrix, output_folder)
#' @examples pcs <- plotVariance(matrix, vars, 2000)

generatePCs <- function(MAT, VARS, NFEATURES) {
  #----- Set a variable for the number of genes (features) to be used for PCA and clustering
  var_features_n <- NFEATURES
  #----- Order variance and select the rows (genes) with the most variance
  select <- order(VARS, decreasing = TRUE)[1:var_features_n]
  #----- Subset vsd values for genes by top variance ranks
  vsd_sub <- MAT[select,]
  #----- Transpose the matrix
  vsd_sub <- t(vsd_sub)
  #----- Run principal component analysis
  message(paste0("Running PCA on ", NFEATURES, " most variable features..."))
  pca <- prcomp(vsd_sub)
  #----- extract the variance explained by each PC
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  #----- subset for first 5 elements
  percentVar <- percentVar[1:5]
  #----- add names to the percentVar vector
  message("Percent variations:")
  percentVar <- paste(round(percentVar*100,2), "%", sep = " ")
  names(percentVar) <- c("PC1", "PC2", "PC3", "PC4", "PC5")
  print(percentVar)
  # ---- Construct a data frame with PC loadings and add sample labels
  pca_res <- list()
  pca_df <- as.data.frame(pca$x)
  pca_eigenvecs <- as.data.frame(pca$rotation[,1:2])
  pca_res[["Loadings"]] <- pca_df
  pca_res[["Eigenvectors"]] <- pca_eigenvecs
  pca_res[["percent_var"]] <- percentVar
  #----- Return PCA values
  return(pca_res)
}
