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
