#' Plot PCA results from a GenEDA object
#'
#' @description
#' Visualize principal component analysis (PCA) results stored within a
#' `GenEDA` object. This function extracts PCA data via \code{\link{ExtractPCA}}
#' and provides a flexible ggplot2-based visualization interface.
#'
#' @param object A `GenEDA` object containing PCA results in the `DimReduction` slot.
#' @param x Numeric or character. The principal component to plot on the x-axis (e.g., 1 or "PC1").
#' @param y Numeric or character. The principal component to plot on the y-axis (e.g., 2 or "PC2").
#' @param color_by Character. Column name in the metadata used to color points.
#' @param colors Vector. Custom colors to use for plotting
#' @param split_by Character (optional). Column name in metadata used for faceting (creates separate panels).
#' @param shape_by Character (optional). Column name in metadata used to control point shape.
#' @param return_data Logical (default = FALSE), whether or not to return pca dataframe for more custom plotting.
#'
#' @return A `ggplot` object displaying the PCA scatter plot, or a list of pca_df and plot if `return_data = TRUE`
#' @examples
#' \dontrun{
#' p <- PlotPCA(obj, x = 1, y = 2, color_by = "condition",
#'      colors = c("untreated" = "red", "treated" = "blue"),
#'      split_by = "library")
#' p
#' }
#' @export
PlotPCA <- function(object,
                    x = 1,
                    y = 2,
                    color_by,
                    colors = NULL,
                    split_by = NULL,
                    shape_by = NULL,
                    return_data = FALSE) {

  # Check input type
  stopifnot(methods::is(object, "geneda"))

  # Extract PCA data
  pca_df <- ExtractPCA(object)

  # Ensure components exist
  if (is.numeric(x)) x <- paste0("PC", x)
  if (is.numeric(y)) y <- paste0("PC", y)

  if (!(x %in% colnames(pca_df))) stop(paste("Column", x, "not found in PCA data."))
  if (!(y %in% colnames(pca_df))) stop(paste("Column", y, "not found in PCA data."))

  # Retrieve variance explained
  percent_var <- object@DimReduction$percent_var
  xlab <- percent_var[[x]]
  xlab <- gsub(" ", "", xlab)
  ylab <- percent_var[[y]]
  ylab <- gsub(" ", "", ylab)

  # Check color_by
  if (!(color_by %in% colnames(pca_df))) stop(paste("Column", color_by, "not found in PCA data."))

  # Check custom colors
  if (!is.null(colors)) {
    if (!(color_by %in% colnames(pca_df))) stop(paste("Column", color_by, "not found in PCA data."))
    # Ensure all levels are mapped
    unique_vals <- unique(pca_df[[color_by]])
    if (!all(unique_vals %in% names(colors))) {
      stop(paste0(
        "Color vector does not include all unique values of '", color_by, "'. Missing: ",
        paste(setdiff(unique_vals, names(colors)), collapse = ", ")
      ))
    }
  }

  # Construct ggplot
  p <- ggplot(pca_df, aes_string(x = x, y = y, color = color_by)) +
    geom_point(size = 3, alpha = 0.8) +
    theme_classic(base_size = 16) +
    theme(
      axis.title = element_text(face = "bold"),
      legend.position = "right"
    ) +
    labs(
      x = paste0(x, ": ", xlab),
      y = paste0(y, ": ", ylab),,
      color = color_by
    )

  # Apply custom colors if provided
  if (!is.null(colors)) {
    p <- p + scale_color_manual(values = colors)
  }

  # Add shape mapping if specified
  if (!is.null(shape_by)) {
    if (!(shape_by %in% colnames(pca_df))) stop(paste("Column", shape_by, "not found in PCA data."))
    p <- p + aes_string(shape = shape_by)
  }

  # Add faceting if requested
  if (!is.null(split_by)) {
    if (!(split_by %in% colnames(pca_df))) stop(paste("Column", split_by, "not found in PCA data."))
    p <- p + facet_wrap(as.formula(paste("~", split_by))) +
      theme(strip.text = element_text(face = "bold"))
  }

  # Return either plot or both plot + data
  if (return_data) {
    return(list(plot = p, data = pca_df))
  } else {
    return(p)
  }
}
