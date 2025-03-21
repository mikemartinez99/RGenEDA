#' plotVariance
#'
#' @description Analyze gene variance to determine thresholds for PCA. Plot gene by variance plot.
#'
#' @param MAT A data matrix where rows are features and columns are samples
#' @param OUTPUT A directory path to where the output should be saved. Path must end in a "/"
#' @param LOG Log transform variance. Default = FALSE
#' @param nFeaturesDrop Number of top most variable features to drop. Useful if top genes are extremely variable. Default = NULL
#'
#' @returns A vector of variances named by gene. A gene by variance plot in tiff format at the specified output location.
#'
#' @examples output_folder <- c("Path/to/outputs/")
#' @examples # Plot variance as it
#' @examples plotVariance(matrix, output_folder)
#' @examples # Plot logged transformed variance
#' @examples plotVariance(matrix, output_folder, LOG = TRUE)
#' @examples # Plot variance with top 3 features dropped
#' @examples plotVariance(matrix, output_folder, nFeaturesDrop = 3)
#' @examples # If running PCA, save outputs to variable to access the features
#' @examples vars <- plotVariance(matrix, output_folder)

plotVariance <- function(MAT, OUTPUT, LOG = FALSE, nFeaturesDrop = NULL) {
  #----- Get dimensions of matrix
  rows <- nrow(MAT)
  columns <- ncol(MAT)

  #----- Extract the vsd assay as a matrix
  message(paste0("Calculating variance on ", rows, " x ", columns, " matrix..."))
  expMat <- MAT

  #----- Calculate variance in a row-wise fashion and sort in decreasing order
  var <- apply(expMat, 1, var)
  var <- sort(var, decreasing = TRUE)

  #----- Drop top nFeaturesDrop most variable features if specified and valid
  if (!is.null(nFeaturesDrop) && nFeaturesDrop > 0) {
    message(paste0("Dropping top ", nFeaturesDrop, " most variable features..."))
    var <- var[-seq_len(min(nFeaturesDrop, length(var)))]
  }

  #----- Reset plotting window to 1 row x 1 column
  par(mfrow = c(1,1))

  #----- Initialize a plot file
  tiff(paste0(OUTPUT, "Variable_Features.tiff"), width = 8, height = 8, units = "in", res = 200)

  #----- Plot variance for genes across samples
  plot(
    var,
    las = 1,
    main = "Sample Gene Expression Variance",
    xlab = "nGenes",
    cex.lab = 1.4,
    cex.axis = 1.1,
    font.lab = 2,
    ylab = ifelse(LOG, "Log(Variance)", "Variance"),
    log = if (LOG) "y" else ""  # Log-transform y-axis if LOG = TRUE
  )

  #----- Add vertical lines at specific gene number indexes
  abline(v = 1000, col = "#E69F00", lty = 2)
  abline(v = 500, col = "#009E73", lty = 2)
  abline(v = 250, col = "#7570B3", lty = 2)

  #----- Add labels to the lines
  text(
    x = c(1000, 500, 250),
    y = max(var) * 0.9,  # Adjust Y position slightly below max variance
    labels = c("1000", "500", "250"),
    col = c("#E69F00", "#009E73", "#7570B3"),
    pos = 4,  # Right-aligned text
    cex = 0.8,
    srt = 45# Adjust font size
  )

  #----- Close the graphics device
  dev.off()

  #----- Return the variances
  return(var)
}
