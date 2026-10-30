library(testthat)
library(RGenEDA)

# Helper: geneda object with a full DESeq2-like DEG table
create_test_geneda <- function() {
    set.seed(1)
    mat <- matrix(
        rnorm(60),
        nrow = 10,
        dimnames = list(paste0("Gene", 1:10), paste0("S", 1:6))
    )
    meta <- data.frame(
        condition = rep(c("A", "B"), each = 3),
        row.names = colnames(mat)
    )
    obj <- GenEDA(normalized = mat, metadata = meta)
    deg <- data.frame(
        baseMean = runif(10, 10, 1000),
        log2FoldChange = c(3, -3, rnorm(8, 0, 0.2)),
        pvalue = c(1e-6, 1e-6, runif(8)),
        padj = c(1e-4, 1e-4, runif(8, 0.5, 1)),
        row.names = rownames(mat)
    )
    SetDEGs(obj, deg, "DEG")
}

#----- DEDashboard works with default num/den
test_that("DEDashboard works when num and den are left as NULL", {
    db <- DEDashboard(create_test_geneda(), "DEG")
    expect_s3_class(db, "patchwork")
})

#----- DEDashboard passes group names through to the volcano plot
test_that("DEDashboard passes num and den to the volcano plot", {
    db <- DEDashboard(create_test_geneda(), "DEG", num = "KO", den = "WT")
    expect_s3_class(db, "patchwork")
})

#----- DEDashboard errors when required columns are missing
test_that("DEDashboard errors when required columns are missing", {
    obj <- create_test_geneda()
    obj@DEGs$DEG$baseMean <- NULL
    expect_error(
        DEDashboard(obj, "DEG"),
        "Columns missing! Required: log2FoldChange, baseMean, pvalue, padj"
    )
})
