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
