#' Save a pheatmap Objects generated from RGenEDA to File
#'
#' Saves a `pheatmap` derived `RGenEDA` heatmap object to a file in various
#' formats
#' including PNG, JPEG, TIFF, BMP, PDF, and SVG.
#'
#' @param pheatmap_obj A `pheatmap` object returned by an RGenEDA or pheatmap
#'   function.
#' @param filename Character string specifying the path and filename where
#'   the heatmap should be saved.
#' @param width Numeric, the width of the output image. Default is 8.
#' @param height Numeric, the height of the output image. Default is 6.
#' @param units Character, units for width and height: one of "in", "cm",
#'   "mm", or "px".
#'   Default is "in". Values are converted to inches for PDF and SVG output.
#' @param res Numeric, the resolution (in dpi) for raster images. Default is
#'   300.
#' @param ... Additional arguments passed to the graphics device function.
#'
#' @return Invisibly returns the filename of the saved heatmap.
#'
#' @details
#' This function is similar in spirit to `ggsave()` but works specifically
#' with `pheatmap` objects.
#'
#' @export
#'
#' @examples
#' data(Snail1KO)
#' counts <- as.matrix(Snail1KO$rawCounts)
#' counts <- counts[rowSums(counts >= 10) >= 3, ]
#' meta <- Snail1KO$metadata
#' rownames(meta) <- meta$Sample
#' norm <- log2(t(t(counts) / colSums(counts)) * 1e6 + 1)
#' obj <- GenEDA(normalized = norm, metadata = meta, counts = counts)
#' dists <- PlotDistances(obj, meta_cols = "Genotype")
#' out <- tempfile(fileext = ".png")
#' GenSave(dists, out, width = 6, height = 5)
#' file.exists(out)
GenSave <- function(
    pheatmap_obj,
    filename,
    width = 8,
    height = 6,
    units = "in",
    res = 300,
    ...) {
    pheatmap_obj <- pheatmap_obj$heatmap

    # Determine file extension
    ext <- tools::file_ext(filename)
    ext <- tolower(ext)

    # Vector devices only take inches
    units <- match.arg(units, c("in", "cm", "mm", "px"))
    toInches <- switch(units,
        "in" = 1,
        "cm" = 1 / 2.54,
        "mm" = 1 / 25.4,
        "px" = 1 / res
    )

    # Open graphics device based on extension
    if (ext == "png") {
        grDevices::png(
            filename,
            width = width,
            height = height,
            units = units,
            res = res,
            ...
        )
    } else if (ext %in% c("jpeg", "jpg")) {
        grDevices::jpeg(
            filename,
            width = width,
            height = height,
            units = units,
            res = res,
            ...
        )
    } else if (ext %in% c("tiff", "tif")) {
        grDevices::tiff(
            filename,
            width = width,
            height = height,
            units = units,
            res = res,
            ...
        )
    } else if (ext == "bmp") {
        grDevices::bmp(
            filename,
            width = width,
            height = height,
            units = units,
            res = res,
            ...
        )
    } else if (ext %in% c("pdf")) {
        grDevices::pdf(
            filename,
            width = width * toInches,
            height = height * toInches,
            ...
        )
    } else if (ext %in% c("svg")) {
        grDevices::svg(
            filename,
            width = width * toInches,
            height = height * toInches,
            ...
        )
    } else {
        cli::cli_abort("Unsupported file format: {ext}")
    }

    # Draw the heatmap onto the current device
    if ("gtable" %in% class(pheatmap_obj$gtable)) {
        grid::grid.newpage()
        grid::grid.draw(pheatmap_obj$gtable)
    } else {
        cli::cli_abort(
            paste0(
                "pheatmap object does not contain a gtable. Did you use ",
                "pheatmap(..., silent = TRUE)?"
            )
        )
    }

    # Close the device
    grDevices::dev.off()

    invisible(filename)
}
