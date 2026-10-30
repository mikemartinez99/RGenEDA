#' PlotVolcano
#'
#' @description Plot volcano plot from a `geneda` object. Requires that you
#'   have ran
#' `SetDEGs` on your object following differential expression in `DESeq2`.
#'
#' @import ggplot2
#' @import ggrepel
#' @import RColorBrewer
#' @importFrom stats setNames
#'
#' @param object A `geneda` object containing `DEGs` from `SetDEGs` method
#' @param assay The DEG slot to use for visualization
#' @param alpha Threshold for adjusted p-values (padj column from `DESeq2`)
#' @param l2fc Absolute value log2Fold-change magnitude threshold
#'   (log2FoldChange column)
#' @param den Denominator (reference level for comparison, num vs. den).
#'   If `NULL`, the legend uses "denominator".
#' @param num Numerator (numerator level for comparison, num vs. den).
#'   If `NULL`, the legend uses "numerator".
#' @param title Optional character vector of what plot should be titled.
#'
#' @details Thresholds are inclusive (`padj <= alpha`, `|log2FoldChange| >=
#'   l2fc`),
#' matching `FilterDEGs()`. Genes with a missing `log2FoldChange` or `padj`
#' are not plotted.
#'
#' @returns A `ggplot2` object
#'
#' @export
#'
#' @examplesIf requireNamespace("DESeq2", quietly = TRUE)
#' data(Snail1KO)
#' counts <- as.matrix(Snail1KO$rawCounts)
#' counts <- counts[rowSums(counts >= 10) >= 3, ]
#' meta <- Snail1KO$metadata
#' rownames(meta) <- meta$Sample
#' meta$Genotype <- factor(meta$Genotype, levels = c("WT", "Snai1_KO"))
#' norm <- log2(t(t(counts) / colSums(counts)) * 1e6 + 1)
#' obj <- GenEDA(normalized = norm, metadata = meta, counts = counts)
#'
#' # Differential expression with DESeq2 (Snai1_KO vs WT)
#' dds <- DESeq2::DESeqDataSetFromMatrix(counts, meta, design = ~Genotype)
#' dds <- DESeq2::DESeq(dds, quiet = TRUE)
#' res <- as.data.frame(DESeq2::results(dds))
#' obj <- SetDEGs(obj, res, assay = "KO_vs_WT")
#' PlotVolcano(obj,
#'     assay = "KO_vs_WT", alpha = 0.05, l2fc = 1,
#'     den = "WT", num = "Snai1_KO"
#' )
PlotVolcano <- function(
    object,
    assay,
    alpha,
    l2fc,
    den = NULL,
    num = NULL,
    title = NULL) {
    stopifnot(methods::is(object, "geneda"))

    if (!assay %in% names(object@DEGs)) {
        cli::cli_abort("Assay {assay} was not found in DEGs slot!")
    } else {
        df <- DEGs(object, assay)
        if (nrow(df) == 0) {
            cli::cli_abort(
                "No differential expression results found in Assay {assay}"
            )
        }
    }

    requiredCols <- c("log2FoldChange", "padj")
    if (!all(requiredCols %in% colnames(df))) {
        cli::cli_abort(
            "log2FoldChange or padj or both are missing from DEG table!"
        )
    }

    if (is.null(num)) num <- "numerator"
    if (is.null(den)) den <- "denominator"
    if (identical(num, den)) {
        cli::cli_abort("num and den must be different group names.")
    }

    # Only drop genes missing the columns used for plotting
    df <- df[stats::complete.cases(df$log2FoldChange, df$padj), , drop = FALSE]

    sig <- df$padj <= alpha
    df$Group <- "ns"
    df$Group[sig] <- paste("padj <", alpha)
    df$Group[sig & df$log2FoldChange >= l2fc] <- paste("Upregulated in", num)
    df$Group[sig & df$log2FoldChange <= -l2fc] <- paste("Upregulated in", den)

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
        "ns"
    )
    colors <- c(
        "forestgreen",
        "firebrick",
        "steelblue3",
        "darkgrey"
    )
    groupColors <- stats::setNames(colors, groupLevels)

    df$gene <- rownames(df)

    cap_pos_xvalue <- max(abs(df$log2FoldChange)) + 0.4
    cap_neg_xvalue <- cap_pos_xvalue / -1

    resPlot <- ggplot(
        df,
        aes(x = log2FoldChange, y = -log10(padj), fill = Group)
    ) +
        geom_point(alpha = 0.5, shape = 21) +
        theme_classic(base_size = 16) +
        ylab("-log10(padj)") +
        xlab("Log2FoldChange") +
        geom_hline(
            yintercept = -log10(alpha),
            color = "black",
            linetype = "dashed",
            linewidth = 0.4
        ) +
        labs(fill = "Group") +
        scale_fill_manual(values = groupColors) +
        scale_size(name = "-log10 P-value", range = c(0.5, 12)) +
        xlim(cap_neg_xvalue, cap_pos_xvalue) +
        geom_label_repel(
            data = subset(df, log2FoldChange >= l2fc & padj <= alpha), aes(
                label = gene
            ),
            box.padding = 0.35,
            nudge_x = 0.05,
            nudge_y = 0.04,
            point.padding = 0.5,
            label.size = 0.15,
            max.overlaps = 20,
            segment.size = 0.3, fill = "grey90",
            segment.color = "grey50", size = 3.5
        ) +
        geom_label_repel(
            data = subset(df, log2FoldChange <= -l2fc & padj <= alpha), aes(
                label = gene
            ),
            box.padding = 0.35,
            nudge_x = 0.05,
            nudge_y = -0.04,
            point.padding = 0.5,
            label.size = 0.15,
            segment.size = 0.3, fill = "grey90",
            max.overlaps = 20,
            segment.color = "grey50", size = 3.5
        ) +
        geom_vline(xintercept = l2fc, colour = "black", linetype = "dotted") +
        geom_vline(xintercept = -l2fc, colour = "black", linetype = "dotted") +
        theme(
            legend.key.size = unit(1, "cm"),
            title = element_text(colour = "black", size = 20),
            axis.text.x = element_text(colour = "black", size = 14),
            axis.text.y = element_text(colour = "black", size = 14),
            axis.title.x = element_text(
                colour = "black",
                size = 16,
                face = "bold"
            ),
            axis.title.y = element_text(
                colour = "black",
                size = 16,
                face = "bold"
            ),
            legend.text = element_text(colour = "black", size = 14),
            legend.title = element_text(colour = "black", size = 14)
        ) +
        guides(
            size = guide_legend(override.aes = list(fill = "black", alpha = 1))
        )

    if (!is.null(title)) {
        resPlot <- resPlot + ggtitle(title)
        return(resPlot)
    } else {
        return(resPlot)
    }
}
