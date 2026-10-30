#' PlotScree
#'
#' @description Generate a scree plot
#'
#' @param object A `geneda` object containing PCA results in the DimReduction
#'   slot.
#'
#' @return A ggplot2 object
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
#' obj <- FindVariableFeatures(obj, nfeatures = 1000)
#' obj <- RunPCA(obj)
#' PlotScree(obj)
PlotScree <- function(object) {
    stopifnot(methods::is(object, "geneda"))

    #----- Pull PCA results
    if (length(object@DimReduction) == 0) {
        cli::cli_abort("No dimensionality reduction results were found.")
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
        geom_point(
            shape = 21,
            color = "black",
            fill = "steelblue",
            size = 3,
            stroke = 0.5
        ) +
        geom_line(aes(group = 1), color = "steelblue", linewidth = 1) +
        theme_bw(base_size = 16) +
        theme(
            axis.title = element_text(face = "bold"),
            axis.text.x = element_text(face = "bold", angle = 45, hjust = 1),
            strip.text = element_text(face = "bold"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = "none"
        ) +
        labs(
            x = "",
            y = "Variance Explained (%)"
        )

    return(scree)
}
