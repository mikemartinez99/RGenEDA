#' PlotVolcano
#'
#' @description Plot volcano plot from a `geneda` object. Requires that you have ran
#' `SetDEGs` on your object following differential expression in `DESeq2`.
#'
#' @import ggplot2
#' @import ggrepel
#' @import RColorBrewer
#'
#' @param object A `geneda` object containing `DEGs` from `SetDEGs` method
#' @param alpha Threshold for adjusted p-values (padj column from `DESeq2`)
#' @param fc Absolute value log2Fold-change magnitude threshold (log2FoldChange column)
#' @param den Denominator (reference level for comparison, num vs. den)
#' @param num Numerator (numerator level for comparison, num vs. den)
#' @param title Optional character vector of what plot should be titled.
#'
#' @returns A `ggplot2` object
#'
#' @examples
#' \donttest{
#' mock_norm <- matrix(rnorm(5000 * 6, mean = 0, sd = 2), nrow = 5000, ncol = 6)
#' colnames(mock_norm) <- paste0("Sample", 1:6)
#' rownames(mock_norm) <- paste0("Gene", 1:5000)
#'
#' # Sample metadata
#' mock_meta <- data.frame(condition = c("A","B","A","B","A","B"),
#' row.names = colnames(mock_norm))
#'
#' # Construct GenEDA object
#' obj <- GenEDA(normalized = mock_norm, metadata = mock_meta)
#'
#' # make dummy DEG table
#' n <- 1000
#' genes <- paste0("Gene", 1:n)
#' log2FoldChange <- runif(n, min = -5, max = 5)
#' padj <- runif(n, min = 0, max = 1)
#' df <- data.frame(
#'  log2FoldChange = log2FoldChange,
#'  padj = padj,
#'  stringsAsFactors = FALSE)
#' rownames(df) <- genes
#' obj <- SetDEGs(obj, df)
#'
#' PlotVolcano(obj, 0.05, 1, "Denominator Group", "Numerator Group", "Test")
#' }
#' @export

PlotVolcano <- function(object, alpha, fc, den, num, title = NULL) {
  stopifnot(methods::is(object, "geneda"))

  df <- DEGs(object, "DEG")
  if (nrow(df) == 0) {
    stop("No differential expression results found in object@DEGs$DEG")
  }

  requiredCols <- c("log2FoldChange", "padj")
  if (!all(requiredCols %in% colnames(df))) {
    stop("log2FoldChange or padj or both are missing from DEG table!")
  }

  df$Group <- NA
  df$Group[df$padj > alpha & abs(df$log2FoldChange) < fc] <- "ns"
  df$Group[df$padj > alpha & abs(df$log2FoldChange) > fc] <- "ns"
  df$Group[df$padj < alpha & abs(df$log2FoldChange) < fc] <- paste("padj <", alpha)
  df$Group[df$padj < alpha & df$log2FoldChange > fc] <- paste("Upregulated in", num)
  df$Group[df$padj < alpha & df$log2FoldChange < -fc] <- paste("Upregulated in", den)

  df$Group <- factor(df$Group, levels = c(
    paste("Upregulated in", den),
    paste("Upregulated in", num),
    paste("padj <", alpha),
    "ns"
  ))

  groupLevels <- c(
    paste("Upregulated in", num),
    paste("Upregulated in", den),
    paste("padj <", alpha),
    "ns")
  colors <- c(
    "forestgreen",
    "firebrick",
    "steelblue3",
    "darkgrey")
  groupColors <- setNames(colors, groupLevels)

  df$gene <- rownames(df)

  cap_pos_xvalue <- max(abs(df$log2FoldChange)) + 0.4
  cap_neg_xvalue <- cap_pos_xvalue/-1

  resPlot <- ggplot(df, aes(x = log2FoldChange, y = -log10(padj), fill = Group)) +
    geom_point(alpha = 0.5, shape = 21)  +
    theme_classic(base_size = 16) +
    ylab("-log10(padj)") +
    xlab("Log2FoldChange") +
    geom_hline(yintercept = -log10(alpha), color = "black", linetype = "dashed", size = 0.4) +
    labs(fill = "Group") +
    scale_fill_manual(values = groupColors) +
    scale_size(name = "-log10 P-value", range = c(0.5, 12)) +
    xlim(cap_neg_xvalue, cap_pos_xvalue) +
    geom_label_repel(data = subset(df, log2FoldChange > fc & padj < alpha), aes(label = gene),
                     box.padding   = 0.35,
                     nudge_x = 0.05,
                     nudge_y = 0.04,
                     point.padding = 0.5,
                     label.size = 0.15,
                     max.overlaps = 20,
                     segment.size = 0.3, fill = "grey90",
                     segment.color = 'grey50', size = 3.5) +
    geom_label_repel(data = subset(df, log2FoldChange < -fc & padj < alpha), aes(label = gene),
                     box.padding   = 0.35,
                     nudge_x = 0.05,
                     nudge_y = -0.04,
                     point.padding = 0.5,
                     label.size = 0.15,
                     segment.size = 0.3, fill = "grey90",
                     max.overlaps = 20,
                     segment.color = 'grey50', size = 3.5) +
    geom_vline(xintercept = fc, colour = "black", linetype="dotted") +
    geom_vline(xintercept = -fc, colour = "black", linetype="dotted") +
    theme(legend.key.size = unit(1, "cm"),
          title = element_text(colour="black", size = 20),
          axis.text.x=element_text(colour="black", size = 14),
          axis.text.y=element_text(colour="black", size = 14),
          axis.title.x=element_text(colour="black", size = 16, face = "bold"),
          axis.title.y=element_text(colour="black", size = 16, face = "bold"),
          legend.text = element_text(colour="black", size = 14),
          legend.title = element_text(colour="black", size = 14)) +
    guides(size=guide_legend(override.aes = list(fill="black", alpha=1)))

  if (!is.null(title)) {
    resPlot + ggtitle(title)
  }

  return(resPlot)
}
