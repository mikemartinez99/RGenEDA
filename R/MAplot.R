#' MAplot
#'
#' @description Generate publicatrion quality MA plot for bulk RNA-Seq differential expression results
#'
#' @import ggplot2
#' @import ggblend
#' @import ggrepel
#'
#' @param results A dataframe of bulk RNA-Seq differential expression results from DESeq2
#' @param numerator Numerator level name
#' @param refLevel Reference level name
#' @param log2FC_thresh Absolute value of log2FC magnitude threshold
#' @param padj_thresh Adjusted p-value significance threshold
#' @param title Plot title
#' @param figDir Output directory ending in "/"
#' @param figName Output figure name with .png extension
#'
#' @returns An output figure and ggplot object


MAplot <- function(results,
                   numerator,
                   refLevel,
                   log2FC_thresh,
                   padj_thresh,
                   title,
                   figDir,
                   figName) {
  #----- Set colors
  colors <- c("firebrick",
              "steelblue",
              "#b0b0b0")

  #----- Set group thresholds
  denThresh <- log2FC_thresh * -1

  #----- Add group column to results
  results$Group <- ifelse(results$log2FoldChange < denThresh & results$padj <= padj_thresh, paste0("Enriched in ", refLevel),
                          ifelse(results$log2FoldChange > numThresh & results$padj <= padj_thresh, paste0("Enriched in ", num), "ns"))

  #----- Filter to just rows with top 15% base mean values (for labelling)
  threshold <- quantile(results$baseMean, 0.85, na.rm = TRUE)
  filtered <- subset(results, baseMean > threshold)

  #----- Get the top 5 downregulated and top 5 upregulated genes based on l2FC
  top_hitsB <- filtered[order(filtered$log2FoldChange, decreasing = FALSE), ][1:5, ]
  top_hitsA <- filtered[order(filtered$log2FoldChange, decreasing = TRUE), ][1:5, ]

  #----- Combine the hits
  top_hits <- rbind(top_hitsB, top_hitsA)
  top_hits$feature <- rownames(top_hits)

  #----- Trace hits
  trace_filter <- abs(results$log2FoldChange) > log2FC_thresh & results$padj <= padj_thresh

  #----- Plot
  MAFigure <- ggplot(results, aes(x = log(baseMean), y = log2FoldChange, fill = Group, size = -log10(padj))) +
    geom_point(data = results, aes(fill = Group),
               color = "white", fill = "grey85", shape = 21, alpha = 0.4, stroke = 0.3) +
    geom_point(data = results[trace_filter, ],
               aes(fill = Group), color = "black", shape = 21, alpha = 0.75, stroke = 0.8) +
    geom_text_repel(data = top_hits, aes(label = feature), size = 4, max.overlaps = 110) +
    annotate("text", x = -Inf, y = 1.2, label = paste0("Enriched in ", numerator), hjust = -0.1, vjust = 0, fontface = 4) +
    annotate("text", x = -Inf, y = -1.2, label = paste0("Enriched in ", refLevel), hjust = -0.1, vjust = 1, fontface = 4) +
    geom_hline(yintercept = -1, linetype = "dashed") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_color_hue(direction = 1) +
    labs(
      y = "Log2 Fold Change",
      x = "Log10 Base Mean",
      size = "-log10 padj",
      color = "",
      fill = ""
    ) +
    scale_size(range = c(1, 10)) +
    theme_classic(base_size = 16) +
    theme(axis.title = element_text(face = "bold")) +
    guides(fill = "none")

  ggsave(paste0(figDir, figName), MAFigure, width = 8, height = 8)
  message(paste0("Plotted ", paste0(figDir, figName)))
  return(MAFigure)
}
