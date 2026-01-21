#' PlotScree
#'
#' @description Generate a scree plot
#'
#' @param object A `geneda` object containing PCA results in the DimReduction slot.
#'
#' @return A ggplot2 object
#' @export
PlotScree <- function(object) {
  stopifnot(methods::is(object, "geneda"))

  #----- Pull PCA results
  if (length(object@DimReduction) == 0) {
    stop("No dimensionality reduction results were found.")
  }

  #----- Grab the variance results
  pcVars <- data.frame(
    PC = factor(
      names(object@DimReduction$percent_var),
      levels = names(object@DimReduction$percent_var)
    ),
    Variance = as.numeric(
      gsub(" %", "", object@DimReduction$percent_var)
    )
  )

  #----- Plot
  scree <- ggplot(pcVars, aes(x = PC, y = Variance)) +
    geom_point(shape = 21, color = "black", fill = "steelblue", size = 3, stroke = 0.5) +
    geom_line(aes(group = 1), color = "steelblue", linewidth = 1) +
    theme_bw(base_size = 16) +
    theme(axis.title = element_text(face = "bold"),
          axis.text.x = element_text(face = "bold", angle = 45, hjust = 1),
          strip.text = element_text(face = "bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none") +
    labs(x = "",
         y = "Variance Explained (%)")

  return(scree)
}
