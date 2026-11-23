library(testthat)
library(RGenEDA)

# Helper: geneda object with HVGs, PCA and one DEG table
create_full_geneda <- function() {
    set.seed(11)
    counts <- matrix(rpois(60 * 6, lambda = 50), nrow = 60)
    counts[1:10, 4:6] <- counts[1:10, 4:6] * 4L
    dimnames(counts) <- list(paste0("Gene", 1:60), paste0("S", 1:6))
    meta <- data.frame(
        condition = factor(rep(c("A", "B"), each = 3)),
        `batch id` = rep(c("x", "y"), 3),
        row.names = colnames(counts),
        check.names = FALSE
    )
    norm <- log2(counts + 1)
    obj <- GenEDA(normalized = norm, metadata = meta, counts = counts)
    obj <- suppressMessages(RunPCA(obj, nfeatures = 20))
    deg <- data.frame(
        baseMean = rowMeans(counts),
        log2FoldChange = c(rep(2, 10), rep(0, 50)),
        padj = c(rep(0.001, 10), rep(0.9, 50)),
        row.names = rownames(counts)
    )
    SetDEGs(obj, deg, "B_vs_A")
}

#----- geneda -> SummarizedExperiment
test_that("as(geneda, 'SummarizedExperiment') stores every slot", {
    obj <- create_full_geneda()
    se <- as(obj, "SummarizedExperiment")

    expect_s4_class(se, "SummarizedExperiment")
    expect_equal(SummarizedExperiment::assayNames(se),
        c("counts", "normalized"))
    expect_equal(SummarizedExperiment::assay(se, "normalized"),
        getNormalized(obj))
    expect_equal(SummarizedExperiment::assay(se, "counts"), getCounts(obj))
    expect_equal(colnames(SummarizedExperiment::colData(se)),
        colnames(getMetadata(obj)))

    rd <- SummarizedExperiment::rowData(se)
    expect_setequal(rownames(se)[rd$isHVG], HVGs(obj))
    expect_equal(rownames(se)[order(rd$hvgRank, na.last = NA)], HVGs(obj))

    stored <- S4Vectors::metadata(se)$RGenEDA
    expect_identical(stored$HVGs, HVGs(obj))
    expect_identical(stored$DimReduction, DimReduction(obj))
    expect_identical(stored$DEGs, obj@DEGs)
    expect_null(stored$counts)
})

#----- Round trip
test_that("geneda -> SummarizedExperiment -> geneda is lossless", {
    obj <- create_full_geneda()
    obj2 <- as.geneda(as(obj, "SummarizedExperiment"))
    expect_identical(obj2, obj)
})

test_that("counts with different features are kept in metadata", {
    obj <- create_full_geneda()
    obj@counts <- obj@counts[1:50, ]
    se <- as(obj, "SummarizedExperiment")

    expect_equal(SummarizedExperiment::assayNames(se), "normalized")
    expect_identical(S4Vectors::metadata(se)$RGenEDA$counts, obj@counts)
    expect_identical(as.geneda(se), obj)
})

test_that("as.geneda drops stored results that no longer match", {
    se <- as(create_full_geneda(), "SummarizedExperiment")
    expect_warning(sub <- as.geneda(se[1:5, ]), "do not match the features")
    expect_length(HVGs(sub), 0L)
    expect_length(DimReduction(sub), 0L)

    expect_warning(sub2 <- as.geneda(se[, 1:4]), "do not match the samples")
    expect_length(DimReduction(sub2), 0L)
})

#----- Main functions accept a SummarizedExperiment
create_se <- function() {
    obj <- create_full_geneda()
    SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = getCounts(obj),
            normalized = getNormalized(obj)),
        colData = getMetadata(obj)
    )
}

test_that("modifying functions return the input SummarizedExperiment", {
    se <- create_se()
    se2 <- FindVariableFeatures(se, nfeatures = 15)
    expect_s4_class(se2, "SummarizedExperiment")
    expect_equal(sum(SummarizedExperiment::rowData(se2)$isHVG), 15L)

    se3 <- suppressMessages(RunPCA(se2))
    pca <- S4Vectors::metadata(se3)$RGenEDA$DimReduction
    expect_equal(nrow(pca$Eigenvectors), 15L)

    deg <- DEGs(create_full_geneda(), "B_vs_A")
    se4 <- SetDEGs(se3, deg, "B_vs_A")
    se5 <- FilterDEGs(se4, "B_vs_A", alpha = 0.05, l2fc = 1,
        saveAssay = "B_vs_A_sig")
    expect_equal(names(S4Vectors::metadata(se5)$RGenEDA$DEGs),
        c("B_vs_A", "B_vs_A_sig"))
    expect_equal(nrow(S4Vectors::metadata(se5)$RGenEDA$DEGs$B_vs_A_sig), 10L)

    # The user's own assays and colData are kept
    expect_identical(SummarizedExperiment::assay(se5, "counts"),
        SummarizedExperiment::assay(se, "counts"))
    expect_identical(SummarizedExperiment::colData(se5),
        SummarizedExperiment::colData(se))
})

test_that("SummarizedExperiment and geneda inputs give the same results", {
    obj <- create_full_geneda()
    se <- as(obj, "SummarizedExperiment")

    expect_identical(ExtractPCA(se), ExtractPCA(obj))
    expect_identical(extractEigen(se, "PC1"), extractEigen(obj, "PC1"))
    expect_identical(
        PlotEigenCorr(se, num_pcs = 2)$padj_matrix,
        PlotEigenCorr(obj, num_pcs = 2)$padj_matrix
    )
    expect_identical(
        PlotPCA(se, color_by = "condition")$data,
        PlotPCA(obj, color_by = "condition")$data
    )
    expect_identical(
        FindHVDEGs(se, "B_vs_A", direction = "both"),
        FindHVDEGs(obj, "B_vs_A", direction = "both")
    )
    expect_identical(SummarizeDEGs(se), SummarizeDEGs(obj))
    expect_s3_class(PlotMA(se, "B_vs_A", 0.05, 1), "ggplot")
    expect_s3_class(PlotVolcano(se, "B_vs_A", 0.05, 1), "ggplot")
    expect_s3_class(PlotScree(se), "ggplot")
})

test_that("other input classes give a clear error", {
    expect_error(RunPCA(matrix(1:4, 2)), "must be a")
    expect_error(PlotScree(list()), "must be a")
})

#----- DESeq2 objects keep their class
test_that("DESeqDataSet input is returned as a DESeqDataSet", {
    skip_if_not_installed("DESeq2")
    obj <- create_full_geneda()
    dds <- DESeq2::DESeqDataSetFromMatrix(
        getCounts(obj), getMetadata(obj), design = ~condition
    )
    dds2 <- suppressMessages(RunPCA(dds, nfeatures = 20))
    expect_s4_class(dds2, "DESeqDataSet")
    expect_identical(DESeq2::design(dds2), DESeq2::design(dds))
    expect_true("normalized" %in% SummarizedExperiment::assayNames(dds2))

    # The second call reuses the stored normalized values and PCA
    expect_identical(ExtractPCA(dds2)[, "PC1"],
        S4Vectors::metadata(dds2)$RGenEDA$DimReduction$Scores$PC1)
})

test_that("DESeqTransform input is not given an extra assay", {
    skip_if_not_installed("DESeq2")
    obj <- create_full_geneda()
    dds <- DESeq2::DESeqDataSetFromMatrix(
        getCounts(obj), getMetadata(obj), design = ~condition
    )
    vsd <- DESeq2::varianceStabilizingTransformation(dds, blind = TRUE)
    vsd2 <- suppressMessages(RunPCA(vsd, nfeatures = 20))
    expect_s4_class(vsd2, "DESeqTransform")
    expect_identical(SummarizedExperiment::assayNames(vsd2),
        SummarizedExperiment::assayNames(vsd))
    expect_length(HVGs(as.geneda(vsd2)), 20L)
})

#----- Extra matrix attributes (e.g. from DESeq2::rlog()) survive the trip
test_that("round trip keeps extra attributes of the normalized matrix", {
    obj <- create_full_geneda()
    attr(obj@normalized, "intercept") <- seq_len(nrow(obj@normalized))
    obj2 <- as.geneda(as(obj, "SummarizedExperiment"))
    expect_identical(obj2, obj)
    expect_identical(attr(obj2@normalized, "intercept"),
        seq_len(nrow(obj@normalized)))
})
