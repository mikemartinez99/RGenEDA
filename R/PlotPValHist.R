#' PlotPValueHistogram
#'
#' @description Generate a p-value histogram from a DEG table inside a `geneda` object.
#' This function uses raw p-values only and highlights p-values below a given
#' threshold (default 0.05).
#'
#' @param object A `geneda` object containing DEGs in the DEGs slot.
#' @param assay The name of the DEG table to use (must contain a "pvalue" column).
#' @param bins Number of histogram bins (default 40).
#' @param alpha Significance threshold (default 0.05). Bars representing
#'        p-values < alpha will be highlighted.
#' @param title Optional title for the plot.
#'
#' @return A ggplot2 object
#' @export
PlotPValueHistogram <- function(object, assay, bins = 40, alpha = 0.05, title = NULL) {
  stopifnot(methods::is(object, "geneda"))

  #----- Pull DEG table
  if (!assay %in% names(object@DEGs)) {
    stop(paste("Assay", assay, "was not found in DEGs slot!"))
  }
  df <- DEGs(object, assay)

  #----- Input validation
  if (!"pvalue" %in% colnames(df)) {
    stop("The DEG table must contain a 'pvalue' column.")
  }

  #----- Drop missing values
  df <- df[!is.na(df$pvalue), , drop = FALSE]

  #----- Annotate group for highlighting
  df$p_group <- ifelse(df$pvalue < alpha,
                       paste0("p < ", alpha),
                       paste0("p >= ", alpha))

  #----- Set GenEDA-style colors
  group1 <- paste0("p < ", alpha)
  group2 <- paste0("p >= ", alpha)
  groupColors <- c(
    "firebrick3",
    "darkgrey"
  )

  #----- Plot
  plt <- ggplot(df, aes(x = pvalue, fill = p_group)) +
    geom_histogram(bins = bins, color = "black", alpha = 0.85) +
    theme_classic(base_size = 16) +
    scale_fill_manual(values = groupColors) +
    xlab("Raw p-value") +
    ylab("Count") +
    theme(
      axis.text.x = element_text(colour="black", size = 14),
      axis.text.y = element_text(colour="black", size = 14),
      axis.title.x = element_text(colour="black", size = 16, face="bold"),
      axis.title.y = element_text(colour="black", size = 16, face="bold"),
      legend.title = element_blank(),
      legend.text = element_text(colour="black", size = 14),
      title = element_text(colour="black", size = 20)
    ) +
    scale_x_continuous(limits = c(0, 1))

  if (!is.null(title)) {
    plt <- plt + ggtitle(title)
  }

  return(plt)
}
