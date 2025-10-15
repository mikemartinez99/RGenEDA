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
#' @export
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

  mat <- normalized(object)
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
