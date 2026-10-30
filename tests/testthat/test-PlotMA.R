library(testthat)
library(RGenEDA)

# Helper: geneda object with a DESeq2-like DEG table
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
        log2FoldChange = c(3, -3, 1, rnorm(7, 0, 0.2)),
        padj = c(1e-4, 1e-4, 0.01, runif(7, 0.5, 1)),
        row.names = rownames(mat)
    )
    SetDEGs(obj, deg, "DEG")
}

#----- PlotMA applies the title
test_that("PlotMA applies the title argument", {
    p <- PlotMA(create_test_geneda(), "DEG", 0.05, 1, title = "My MA Plot")
    expect_equal(p$labels$title, "My MA Plot")
})

#----- PlotMA has no title by default
test_that("PlotMA has no title when title is NULL", {
    p <- PlotMA(create_test_geneda(), "DEG", 0.05, 1)
    expect_null(p$labels$title)
})

#----- PlotMA legend reflects alpha and uses inclusive thresholds
test_that("PlotMA legend uses alpha and highlights genes on the thresholds", {
    p <- PlotMA(create_test_geneda(), "DEG", 0.01, 1)
    fillScale <- p$scales$get_scales("fill")
    expect_true(any(grepl("Adj. P<0.01", fillScale$labels)))

    # Gene3 has padj == alpha and |LFC| == l2fc
    expect_equal(p$data["Gene3", "col"], "red")
})

#----- PlotMA errors for missing assay
test_that("PlotMA errors when assay is missing", {
    expect_error(
        PlotMA(create_test_geneda(), "nope", 0.05, 1),
        "Assay nope was not found in DEGs slot!"
    )
})
