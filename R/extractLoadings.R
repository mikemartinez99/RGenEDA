#' extractEigen
#'
#' @description Extract gene-level Eigen-vectors for PCs of interest and calculate percent variance.
#'
#' @param object A `geneda` object containing PCA results in the `DimReduction` slot
#' @param component Principal component to use (i.e., PC1)
#'
#' @returns A data-frame containing gene loading information and percent variation per gene
#' @export
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
