library(testthat)

context("Testing SetDEGs, FilterDEGs, and FindHVDEGs")

# Helper function to create a dummy GenEDA object
create_test_geneda <- function() {
    obj <- methods::new("geneda")
    obj@normalized <- matrix(rnorm(1000), nrow = 100, ncol = 10)
    rownames(obj@normalized) <- paste0("Gene", 1:100)
    colnames(obj@normalized) <- paste0("Sample", 1:10)

    obj@metadata <- data.frame(Sample = colnames(obj@normalized))
    rownames(obj@metadata) <- colnames(obj@normalized) # <-- THIS LINE FIXES IT

    obj@DEGs <- list()
    obj@HVGs <- character(0)
    obj
}

# Create a dummy DESeq2-like DEG table
dummy_deg <- data.frame(
    log2FoldChange = rnorm(100),
    padj = runif(100, 0, 1)
)
rownames(dummy_deg) <- paste0("Gene", 1:100)

# ===== SetDEGs =====#

#----- SetDEGs sets DEG slot correctly
test_that("SetDEGs sets DEG slot correctly", {
    obj <- create_test_geneda()

    # Happy path
    obj2 <- SetDEGs(obj, dummy_deg, assay = "DEG")
    expect_true("DEG" %in% names(obj2@DEGs))
    expect_equal(obj2@DEGs$DEG, dummy_deg)
})

#----- SetDEGs errors if log2FoldChange column is missing
test_that("SetDEGs error if DEG table is missing log2FoldChange column", {
    obj <- create_test_geneda()
    deg_missing <- dummy_deg[, "log2FoldChange", drop = FALSE]
    expect_error(
        SetDEGs(obj, deg_missing, "DEG"),
        "DEG table must contain columns"
    )
})

#----- SetDEGs errors if padj column is missing
test_that("SetDEGs error if DEG table is missing padj column", {
    obj <- create_test_geneda()
    deg_missing <- dummy_deg[, "padj", drop = FALSE]
    expect_error(
        SetDEGs(obj, deg_missing, "DEG"),
        "DEG table must contain columns"
    )
})

#----- SetDEGs errors if assay already exists
test_that("SetDEGs errors if assay already exists", {
    obj <- create_test_geneda()
    obj@DEGs <- list(DEG = dummy_deg)
    expect_error(SetDEGs(obj, dummy_deg, "DEG"), "Assay DEG already exists!")
})


# ===== FilterDEGs =====#

#----- Test that FilterDEGs filters DEGs correctly
test_that("FilterDEGs filters DEGs correctly", {
    obj <- create_test_geneda()
    obj <- SetDEGs(obj, dummy_deg, "DEG")

    # Happy path
    filtered_obj <- FilterDEGs(
        obj,
        assay = "DEG",
        alpha = 0.5,
        l2fc = 0.5,
        saveAssay = "filtered"
    )
    expect_true("filtered" %in% names(filtered_obj@DEGs))
    expect_true(all(abs(filtered_obj@DEGs$filtered$log2FoldChange) >= 0.5))
    expect_true(all(filtered_obj@DEGs$filtered$padj <= 0.5))
})

#----- Test that FilterDEGs errors if assay does not exist
test_that("FilterDEGs errors if assay does not exist", {
    obj <- create_test_geneda()
    obj <- SetDEGs(obj, dummy_deg, "DEG")
    expect_error(
        FilterDEGs(obj, "degs", saveAssay = "filtered"),
        "Assay degs is NULL. Use SetDEGs\\(\\)\\ first."
    )
})

#----- Test that FilterDEGs errors if saveAssay argument is not provided
test_that("FilterDEGs errors if saveAssay argument is not provided", {
    obj <- create_test_geneda()
    obj <- SetDEGs(obj, dummy_deg, "DEG")
    expect_error(
        FilterDEGs(obj, "DEG", 0.05, 0.5),
        "saveAssay must be a single character string."
    )
})

# ===== FindHVDEGs =====#

#----- Test that FindHVDEGs returns correct overlaps
test_that("FindHVDEGs returns correct overlaps", {
    obj <- create_test_geneda()
    obj <- SetDEGs(obj, dummy_deg, "DEG")
    obj@HVGs <- paste0("Gene", 1:5)

    # Error if DEG slot missing
    expect_error(
        FindHVDEGs(obj, "NONEXISTENT"),
        paste0(
            "Assay NONEXISTENT is NULL. Use SetDEGs\\(\\)\\ or ",
            "FilterDEGs\\(\\)\\ first!"
        )
    )

    # Error if HVGs empty
    obj2 <- create_test_geneda()
    obj2 <- SetDEGs(obj2, dummy_deg, "DEG")
    expect_error(FindHVDEGs(obj2, "DEG"), "HVG Slot is empty")

    # Positive direction
    pos_genes <- dummy_deg[dummy_deg$log2FoldChange > 0, ]
    expected_pos <- intersect(rownames(pos_genes), obj@HVGs)
    res_pos <- FindHVDEGs(obj, "DEG", "positive")
    expect_equal(res_pos, expected_pos)

    # Negative direction
    neg_genes <- dummy_deg[dummy_deg$log2FoldChange < 0, ]
    expected_neg <- intersect(rownames(neg_genes), obj@HVGs)
    res_neg <- FindHVDEGs(obj, "DEG", "negative")
    expect_equal(res_neg, expected_neg)

    # Both
    res_both <- FindHVDEGs(obj, "DEG", "both")
    expect_true(all(names(res_both) %in% c("positive", "negative", "both")))
    expect_equal(res_both$positive, expected_pos)
    expect_equal(res_both$negative, expected_neg)
    expect_equal(res_both$both, c(expected_pos, expected_neg))
})

#----- FindHVDEGs direction argument validation works
test_that("FindHVDEGs direction argument validation works", {
    obj <- create_test_geneda()
    obj <- SetDEGs(obj, dummy_deg, "DEG")
    obj@HVGs <- paste0("Gene", 1:5)

    expect_error(
        FindHVDEGs(obj, "DEG", "invalid"),
        "'arg' should be one of .positive., .negative., .both."
    )
})


# ===== SummarizeDEGs =====#

#----- SummarizeDEGs ignores NA values in columns it does not use
test_that("SummarizeDEGs counts genes with NA in unrelated columns", {
    obj <- create_test_geneda()
    deg <- data.frame(
        log2FoldChange = c(3, 3, 3, -3, 0.1),
        padj = c(1e-4, 1e-4, 1e-4, 1e-4, 0.9),
        symbol = c("a", NA, NA, "d", "e"),
        row.names = paste0("Gene", 1:5)
    )
    obj <- SetDEGs(obj, deg, "DEG")
    res <- SummarizeDEGs(obj, assays = "DEG")

    expect_equal(res[1, "DEG"], 3)
    expect_equal(res[2, "DEG"], 1)
})

#----- SummarizeDEGs drops genes missing log2FoldChange or padj
test_that("SummarizeDEGs skips genes with NA log2FoldChange or padj", {
    obj <- create_test_geneda()
    deg <- data.frame(
        log2FoldChange = c(3, NA, 3),
        padj = c(1e-4, 1e-4, NA),
        row.names = paste0("Gene", 1:3)
    )
    obj <- SetDEGs(obj, deg, "DEG")
    expect_equal(SummarizeDEGs(obj, assays = "DEG")[1, "DEG"], 1)
})

#----- SummarizeDEGs and FilterDEGs agree at the thresholds
test_that("SummarizeDEGs thresholds match FilterDEGs at boundaries", {
    obj <- create_test_geneda()
    deg <- data.frame(
        log2FoldChange = c(1, -1, 2, 0.5, 1.5),
        padj = c(0.05, 0.01, 0.05, 0.01, 0.2),
        row.names = paste0("Gene", 1:5)
    )
    obj <- SetDEGs(obj, deg, "DEG")
    obj <- FilterDEGs(obj, "DEG", alpha = 0.05, l2fc = 1, saveAssay = "filt")
    res <- SummarizeDEGs(obj, alpha = 0.05, lfc1 = 1, lfc2 = 2, assays = "DEG")

    expect_equal(res[1, "DEG"] + res[2, "DEG"], nrow(DEGs(obj, "filt")))
    expect_equal(res[3, "DEG"], 1)
})

#----- SummarizeDEGs warns about unknown assays
test_that("SummarizeDEGs warns and skips assays that do not exist", {
    obj <- create_test_geneda()
    obj <- SetDEGs(obj, dummy_deg, "DEG")
    expect_warning(
        SummarizeDEGs(obj, assays = c("DEG", "nope")),
        "Assay nope not found in DEGs slot; skipping."
    )
})
