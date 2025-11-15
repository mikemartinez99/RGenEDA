#' Save a pheatmap Objects generated from RGenEDA to File
#'
#' Saves a `pheatmap` derived `RGenEDA` heatmap object to a file in various formats
#' including PNG, JPEG, TIFF, BMP, PDF, and SVG.
#'
#' @param pheatmap_obj A `pheatmap` object returned by an RGenEDA or pheatmap function.
#' @param filename Character string specifying the path and filename where the heatmap should be saved.
#' @param width Numeric, the width of the output image. Default is 8.
#' @param height Numeric, the height of the output image. Default is 6.
#' @param units Character, units for width and height when saving raster images. Default is "in".
#' @param res Numeric, the resolution (in dpi) for raster images. Default is 300.
#' @param ... Additional arguments passed to the graphics device function.
#'
#' @return Invisibly returns the filename of the saved heatmap.
#'
#' @details
#' This function is similar in spirit to `ggsave()` but works specifically with `pheatmap` objects.
#'
#' @examples
#' \dontrun{
#' library(pheatmap)
#' mat <- matrix(rnorm(100), 10, 10)
#' hm <- pheatmap(mat)
#' GenSave(hm, "heatmap.png")
#' }
#'
#' @export
GenSave <- function(pheatmap_obj, filename, width = 8, height = 6, units = "in", res = 300, ...) {

  pheatmap_obj <- pheatmap_obj$heatmap

  # Determine file extension
  ext <- tools::file_ext(filename)
  ext <- tolower(ext)

  # Open graphics device based on extension
  if (ext %in% c("png", "jpeg", "jpg", "tiff", "bmp")) {
    grDevices::png(filename, width = width, height = height, units = units, res = res, ...)
  } else if (ext %in% c("pdf")) {
    grDevices::pdf(filename, width = width, height = height, ...)
  } else if (ext %in% c("svg")) {
    grDevices::svg(filename, width = width, height = height, ...)
  } else {
    stop("Unsupported file format: ", ext)
  }

  # Draw the heatmap onto the current device
  if ("gtable" %in% class(pheatmap_obj$gtable)) {
    grid::grid.newpage()
    grid::grid.draw(pheatmap_obj$gtable)
  } else {
    stop("pheatmap object does not contain a gtable. Did you use pheatmap(..., silent = TRUE)?")
  }

  # Close the device
  grDevices::dev.off()

  invisible(filename)
}
