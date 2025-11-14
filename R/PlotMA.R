#' PlotMA
#'
#' @description Plot MA plot from a `geneda` object. Requires that you have ran
#' `SetDEGs` on your object following differential expression in `DESeq2`.
#'
#'
#' @import ggplot2
#' @import ggrepel
#' @import RColorBrewer
#'
#' @param object A `geneda` object containing `DEGs` from `SetDEGs` method
#' @param assay The DEG slot to use for visualization
#' @param alpha Threshold for adjusted p-values (padj column from `DESeq2`)
#' @param fc Absolute value log2Fold-change magnitude threshold (log2FoldChange column)
#' @param title Optional character vector of what plot should be titled.
#'
#' @returns A `ggplot2` object
#' @export
PlotMA <- function(object, assay, alpha, fc, title = NULL) {
  stopifnot(methods::is(object, "geneda"))

  if (!assay %in% names(object@DEGs)) {
    stop(paste("Assay", assay, "was not found in DEGs slot!"))
  } else {
    df <- DEGs(object, assay)
    if (nrow(df) == 0) {
      stop("No differential expression results found in object@DEGs$DEG")
    }
  }

  requiredCols <- c("log2FoldChange", "baseMean")
  if (!all(requiredCols %in% colnames(df))) {
    stop("log2FoldChange or baseMean or both are missing from DEG table!")
  }

  df$col <- NA
  df$col[df$padj < alpha & abs(df$log2FoldChange) > fc] <- "red"
  df$col[is.na(df$col)] <- "black"

  df$gene <- rownames(df)

  cap_pos_yvalue <- max(abs(df$log2FoldChange)) + 0.4
  cap_neg_yvalue <- cap_pos_yvalue/-1

  resPlot <- ggplot(df, aes(x = baseMean, y = log2FoldChange, fill=col, size = -log10(padj))) +
    geom_point(alpha = 0.5, shape = 21)  +
    theme_classic() +
    scale_x_continuous(trans='log2') +
    ylab("Log2 fold change") +
    xlab("Mean normalized counts (Log2, SCT)") +
    geom_hline(yintercept = 0, color = "blue", linetype = "dashed", size = 0.4) +
    labs(size = expression("-log"[10]*" adj. P-value)"), fill = "sig") +
    scale_size(range = c(1, 12)) +
    scale_fill_manual(name = " ",
                      values = c("black", "indianred"),
                      labels = c("Not. sig.", paste0("Adj. P<0.05 & |LFC| > ", fc))) +
    scale_size(name = "-log10 P-value", range = c(0.5, 12)) +
    ylim(cap_neg_yvalue, cap_pos_yvalue) +
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
    geom_hline(yintercept = fc, colour = "black", linetype="dotted") +
    geom_hline(yintercept = -fc, colour = "black", linetype="dotted") +
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


# ggplot2
# ggrepel
