#' PlotCountDist
#'
#' @description
#' Generates violin plots with embedded boxplots showing the distribution
#' of normalized or VST-transformed RNA-Seq counts for each sample in a
#' `geneda` object. Optionally, plots can be faceted by a metadata variable.
#'
#' @param object A `geneda` object containing normalized expression values in `@normalized`
#'               and sample metadata in `@metadata`.
#' @param split_by Character. Optional column name from `object@metadata` used
#'                 for faceting (default: NULL, no faceting).
#'
#' @return A ggplot2 object displaying the distribution of counts per sample.
#'
#' @examples
#' \dontrun{
#' # Basic plot without faceting
#' PlotCountDist(obj)
#'
#' # Facet by a metadata variable "condition"
#' PlotCountDist(obj, split_by = "condition")
#' }
#'
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @export
#'
#' @examples
#' \donttest{
#' plotCountDist(object)
#' }
PlotCountDist <- function(object, split_by = NULL) {
  stopifnot(methods::is(object, "geneda"))

  # Convert normalized matrix to long format
  counts_long <- as.data.frame(object@normalized) |>
    mutate(Gene = rownames(object@normalized)) |>
    pivot_longer(
      cols = -Gene,
      names_to = "Sample",
      values_to = "Count"
    )

  # Join metadata
  meta <- as.data.frame(object@metadata)
  meta$Sample <- rownames(meta)  # ensure sample column matches
  counts_long <- counts_long |>
    left_join(meta, by = "Sample")

  # Base ggplot
  p <- ggplot(counts_long, aes(x = Sample, y = Count, fill = Sample)) +
    geom_violin(trim = FALSE, alpha = 0.5) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.5) +
    theme_minimal(base_size = 16) +
    labs(
      x = "",
      y = "Normalized Count",
      title = ""
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      legend.position = "none"
    )

  # Add faceting if split_by is provided
  if (!is.null(split_by) && split_by %in% colnames(meta)) {
    p <- p + facet_wrap(as.formula(paste0("~", split_by)), scales = "free_x")
  } else if (!is.null(split_by)) {
    warning(paste("Column", split_by, "not found in metadata. Ignoring faceting."))
  }

  return(p)
}
