library(testthat)
library(methods)

#-------------------------------
# Helper: Create a test SummarizedExperiment
#-------------------------------
create_test_se <- function(nGenes = 50, extraAssay = FALSE) {
    set.seed(1)
    cts <- matrix(rpois(nGenes * 4, lambda = 50), nrow = nGenes)
    rownames(cts) <- paste0("Gene", seq_len(nGenes))
    colnames(cts) <- paste0("S", 1:4)
    meta <- data.frame(
        condition = factor(c("A", "A", "B", "B")),
        batch = c("X", "Y", "X", "Y"),
        row.names = colnames(cts)
    )
    assays <- list(counts = cts)
    if (extraAssay) assays$logcounts <- log2(cts + 1)
    SummarizedExperiment::SummarizedExperiment(assays = assays, colData = meta)
}

#-------------------------------
# Tests for as.geneda()
#-------------------------------
test_that("as.geneda carries over counts and colData", {
    se <- create_test_se()
    obj <- suppressMessages(as.geneda(se, transform = "logcpm"))

    expect_s4_class(obj, "geneda")
    expect_equal(getCounts(obj), SummarizedExperiment::assay(se, "counts"))
    expect_equal(rownames(getMetadata(obj)), colnames(se))
    expect_equal(getMetadata(obj)$condition, se$condition)
    expect_equal(colnames(getNormalized(obj)), colnames(se))
})

test_that("as.geneda logcpm matches log2(CPM + 1)", {
    se <- create_test_se()
    cts <- SummarizedExperiment::assay(se, "counts")
    expected <- log2(t(t(cts) / colSums(cts)) * 1e6 + 1)

    expect_message(
        obj <- as.geneda(se, transform = "logcpm"),
        "log2\\(CPM \\+ 1\\)"
    )
    expect_equal(getNormalized(obj), expected)
})

test_that("as.geneda uses a named assay for normalized values", {
    se <- create_test_se(extraAssay = TRUE)
    obj <- as.geneda(se, normalized = "logcounts")

    expect_equal(
        getNormalized(obj),
        SummarizedExperiment::assay(se, "logcounts")
    )
})

test_that("as.geneda accepts a matrix for normalized values", {
    se <- create_test_se()
    norm <- log2(SummarizedExperiment::assay(se, "counts") + 1)
    obj <- as.geneda(se, normalized = norm)

    expect_equal(getNormalized(obj), norm)
})

test_that("as.geneda rejects normalized values that do not match x", {
    se <- create_test_se()
    norm <- log2(SummarizedExperiment::assay(se, "counts") + 1)

    expect_error(
        as.geneda(se, normalized = norm[, 4:1]),
        "Column names"
    )
    expect_error(
        as.geneda(se, normalized = norm[10:1, ]),
        "Row names"
    )
    expect_error(as.geneda(se, normalized = "missing"), "not found")
    expect_error(as.geneda(se, normalized = list()), "must be NULL")
})

test_that("as.geneda handles counts = NULL and a missing counts assay", {
    se <- create_test_se(extraAssay = TRUE)
    obj <- as.geneda(se, normalized = "logcounts", counts = NULL)
    expect_null(getCounts(obj))

    obj <- as.geneda(se, normalized = "logcounts", counts = "raw")
    expect_null(getCounts(obj))

    expect_error(
        as.geneda(se, counts = "raw", transform = "logcpm"),
        "No counts assay found"
    )
})

test_that("as.geneda stores a named list of DEG tables", {
    se <- create_test_se()
    deg <- data.frame(
        log2FoldChange = rnorm(50),
        padj = runif(50),
        row.names = rownames(se)
    )
    obj <- suppressMessages(as.geneda(se,
        transform = "logcpm",
        DEGs = list(A_vs_B = deg)
    ))

    expect_equal(DEGs(obj, "A_vs_B"), deg)

    expect_error(
        suppressMessages(as.geneda(se, transform = "logcpm", DEGs = deg)),
        "named list"
    )
    expect_error(
        suppressMessages(
            as.geneda(se, transform = "logcpm", DEGs = list(deg))
        ),
        "named list"
    )
})

test_that("as.geneda output works with RunPCA", {
    se <- create_test_se()
    obj <- suppressMessages(as.geneda(se, transform = "logcpm"))
    obj <- suppressMessages(RunPCA(obj, nfeatures = 20))

    expect_equal(rownames(DimReduction(obj)$Scores), colnames(se))
})

test_that("as.geneda vst requires integer counts", {
    skip_if_not_installed("DESeq2")
    se <- create_test_se()
    SummarizedExperiment::assay(se, "counts") <-
        SummarizedExperiment::assay(se, "counts") + 0.5

    expect_error(as.geneda(se), "integer counts")
})

test_that("as.geneda vst matches DESeq2 on a small matrix", {
    skip_if_not_installed("DESeq2")
    se <- create_test_se()
    cts <- SummarizedExperiment::assay(se, "counts")
    storage.mode(cts) <- "integer"
    expected <- suppressMessages(
        DESeq2::varianceStabilizingTransformation(cts, blind = TRUE)
    )

    msgs <- capture_messages(obj <- as.geneda(se))
    expect_match(msgs, "vst", all = FALSE)
    expect_equal(getNormalized(obj), expected)
})

test_that("as.geneda works with DESeqDataSet and DESeqTransform", {
    skip_if_not_installed("DESeq2")
    se <- create_test_se()
    dds <- suppressMessages(
        DESeq2::DESeqDataSet(se, design = ~condition)
    )
    vsd <- suppressMessages(
        DESeq2::varianceStabilizingTransformation(dds, blind = TRUE)
    )

    # DESeqDataSet with a DESeqTransform for normalized values
    obj <- as.geneda(dds, normalized = vsd)
    expect_equal(getNormalized(obj), SummarizedExperiment::assay(vsd))
    expect_equal(getCounts(obj), DESeq2::counts(dds))

    # A DESeqTransform on its own supplies normalized values only
    obj <- as.geneda(vsd)
    expect_equal(getNormalized(obj), SummarizedExperiment::assay(vsd))
    expect_null(getCounts(obj))
})
