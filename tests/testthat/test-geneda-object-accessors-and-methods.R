library(testthat)
library(ggplot2)
library(methods)

#-------------------------------
# Helper: Create a test geneda object
#-------------------------------
create_test_geneda <- function() {
    # 6 features x 4 samples for testing
    mat <- matrix(runif(24), nrow = 6)
    rownames(mat) <- paste0("Gene", 1:6)
    colnames(mat) <- paste0("S", 1:4)
    meta <- data.frame(
        condition = c("A", "A", "B", "B"),
        batch = c("X", "Y", "X", "Y"),
        row.names = colnames(mat),
        stringsAsFactors = FALSE
    )
    GenEDA(normalized = mat, metadata = meta)
}

#-------------------------------
# Tests for geneda constructor & accessors
#-------------------------------
test_that("GenEDA constructor creates valid geneda object", {
    obj <- create_test_geneda()
    expect_s4_class(obj, "geneda")
    expect_equal(getNormalized(obj), obj@normalized)
    expect_equal(getMetadata(obj), obj@metadata)
    expect_equal(getCounts(obj), NULL)
    expect_equal(HVGs(obj), character(0))
    expect_equal(DimReduction(obj), list())
})

test_that("GenEDA handles counts argument", {
    obj <- create_test_geneda()
    counts <- obj@normalized * 10
    obj2 <- GenEDA(
        normalized = obj@normalized,
        metadata = obj@metadata,
        counts = counts
    )
    expect_equal(getCounts(obj2), counts)
})
test_that("GenEDA stores a named list of DEG tables", {
    obj <- create_test_geneda()
    deg <- data.frame(
        log2FoldChange = c(1, -2, 0.5, 3, -1, 0),
        padj = c(0.01, 0.2, 0.04, 0.001, 0.5, 1),
        row.names = rownames(obj@normalized)
    )
    obj2 <- GenEDA(
        normalized = obj@normalized,
        metadata = obj@metadata,
        DEGs = list(A_vs_B = deg)
    )
    expect_equal(names(obj2@DEGs), "A_vs_B")
    expect_equal(DEGs(obj2, "A_vs_B"), deg)

    expect_error(
        GenEDA(obj@normalized, obj@metadata, DEGs = deg),
        "named list"
    )
    expect_error(
        GenEDA(obj@normalized, obj@metadata, DEGs = list(deg)),
        "named list"
    )
    expect_error(
        GenEDA(obj@normalized, obj@metadata, DEGs = list(x = deg["padj"])),
        "must contain columns"
    )
})

#-------------------------------
# Tests for HVG selection
#-------------------------------
test_that("FindVariableFeatures selects correct number of HVGs", {
    obj <- create_test_geneda()
    obj <- FindVariableFeatures(obj, nfeatures = 3)
    expect_equal(length(HVGs(obj)), 3)
    expect_true(all(HVGs(obj) %in% rownames(obj@normalized)))
})

test_that("FindVariableFeatures clamps nfeatures appropriately", {
    obj <- create_test_geneda()
    obj <- FindVariableFeatures(obj, nfeatures = 10) # more than genes
    expect_equal(length(HVGs(obj)), nrow(obj@normalized))
    obj <- FindVariableFeatures(obj, nfeatures = 0) # zero
    expect_equal(length(HVGs(obj)), 1)
})

library(testthat)
library(ggplot2)
library(methods)

#-------------------------------
# Helper: Create a large test geneda object
#-------------------------------
create_test_geneda_large <- function(n_genes = 5000, n_samples = 50) {
    mat <- matrix(rnorm(n_genes * n_samples, mean = 10, sd = 2),
        nrow = n_genes, ncol = n_samples
    )
    rownames(mat) <- paste0("Gene", seq_len(n_genes))
    colnames(mat) <- paste0("S", seq_len(n_samples))

    meta <- data.frame(
        condition = rep(c("A", "B"), length.out = n_samples),
        batch = rep(c("X", "Y"), length.out = n_samples),
        row.names = colnames(mat),
        stringsAsFactors = FALSE
    )

    GenEDA(normalized = mat, metadata = meta)
}

#-------------------------------
# PCA Tests
#-------------------------------
test_that("RunPCA fills DimReduction slot for large matrix", {
    obj <- create_test_geneda_large()
    obj <- RunPCA(obj, nfeatures = 2000) # realistic top HVGs
    dr <- DimReduction(obj)

    expect_true(all(c("Scores", "Eigenvectors", "percent_var") %in% names(dr)))
    expect_s3_class(dr$Scores, "data.frame")
    expect_s3_class(dr$Eigenvectors, "data.frame")
    expect_equal(nrow(dr$Scores), ncol(obj@normalized))
    expect_equal(nrow(dr$Eigenvectors), length(HVGs(obj)))
})

test_that("ExtractPCA returns combined data.frame for large matrix", {
    obj <- create_test_geneda_large()
    obj <- RunPCA(obj, nfeatures = 2000)
    df <- ExtractPCA(obj)

    expect_s3_class(df, "data.frame")
    expect_equal(nrow(df), ncol(obj@normalized))
    expect_true(all(colnames(obj@metadata) %in% colnames(df)))
})


test_that("PlotPCA returns ggplot object for large matrix", {
    obj <- create_test_geneda_large()
    obj <- RunPCA(obj, nfeatures = 2000)
    p <- PlotPCA(obj, x = 1, y = 2, color_by = "condition")
    expect_s3_class(p, "ggplot")
})

test_that("PlotPCA errors if color_by column missing in metadata", {
    obj <- create_test_geneda_large()
    obj <- RunPCA(obj, nfeatures = 2000)

    # Use a column that definitely doesn't exist
    expect_error(
        PlotPCA(obj, x = 1, y = 2, color_by = "nonexistent_column"),
        "Column nonexistent_column not found in PCA data"
    )
})

test_that("PlotPCA with custom colors works", {
    obj <- create_test_geneda_large()
    obj <- RunPCA(obj, nfeatures = 2000)
    cols <- c("A" = "red", "B" = "blue")
    p <- PlotPCA(obj, x = 1, y = 2, color_by = "condition", colors = cols)
    expect_s3_class(p, "ggplot")
})

#----- GenEDA() starts with an empty DEGs slot, so "DEG" is a valid assay name
test_that("SetDEGs accepts assay name 'DEG' on a freshly constructed object", {
    mat <- matrix(
        rnorm(20),
        nrow = 5,
        dimnames = list(paste0("Gene", 1:5), paste0("S", 1:4))
    )
    meta <- data.frame(
        condition = c("A", "A", "B", "B"),
        row.names = colnames(mat)
    )
    obj <- GenEDA(normalized = mat, metadata = meta)
    expect_length(obj@DEGs, 0)

    deg <- data.frame(
        log2FoldChange = rnorm(5),
        padj = runif(5),
        row.names = rownames(mat)
    )
    obj <- SetDEGs(obj, deg, "DEG")
    expect_equal(DEGs(obj, "DEG"), deg)
})

#----- show() reports DEG assays on one line
test_that("show() reports no DEGs when fresh and lists assays on one line", {
    mat <- matrix(
        rnorm(20),
        nrow = 5,
        dimnames = list(paste0("Gene", 1:5), paste0("S", 1:4))
    )
    meta <- data.frame(
        condition = c("A", "A", "B", "B"),
        row.names = colnames(mat)
    )
    obj <- GenEDA(normalized = mat, metadata = meta)

    out <- capture.output(show(obj))
    expect_true(any(grepl("^  DEGs: none$", out)))

    deg <- data.frame(
        log2FoldChange = rnorm(5),
        padj = runif(5),
        row.names = rownames(mat)
    )
    obj <- SetDEGs(SetDEGs(obj, deg, "a"), deg, "b")
    out <- capture.output(show(obj))
    expect_equal(sum(grepl("DEGs:", out)), 1L)
    expect_true(any(grepl("^  DEGs: a, b$", out)))
})
