library(testthat)
library(RGenEDA)

# Helper: low-variance genes come first in the matrix and high-variance genes
# last, so selecting the first n rows (instead of the HVGs) would be caught
create_test_geneda <- function() {
    set.seed(42)
    nLow <- 50
    nHigh <- 10
    low <- matrix(rnorm(nLow * 6, sd = 0.01), nrow = nLow)
    high <- matrix(rnorm(nHigh * 6, sd = 5), nrow = nHigh)
    mat <- rbind(low, high)
    rownames(mat) <- c(
        paste0("Low", seq_len(nLow)),
        paste0("High", seq_len(nHigh))
    )
    colnames(mat) <- paste0("Sample", 1:6)

    meta <- data.frame(condition = rep(c("A", "B"), each = 3))
    rownames(meta) <- colnames(mat)

    GenEDA(normalized = mat, metadata = meta)
}

#----- RunPCA uses the genes stored in the HVGs slot
test_that("RunPCA computes PCA on the HVGs, not the first rows of the matrix", {
    obj <- create_test_geneda()
    obj <- FindVariableFeatures(obj, nfeatures = 10)
    expect_true(all(grepl("^High", HVGs(obj))))

    obj <- suppressMessages(RunPCA(obj))
    pcaGenes <- rownames(DimReduction(obj)$Eigenvectors)
    expect_identical(pcaGenes, HVGs(obj))
})

#----- RunPCA scores match prcomp on the HVG submatrix
test_that("RunPCA scores match prcomp on the HVG submatrix", {
    obj <- create_test_geneda()
    obj <- FindVariableFeatures(obj, nfeatures = 10)
    obj <- suppressMessages(RunPCA(obj))

    ref <- prcomp(t(getNormalized(obj)[HVGs(obj), ]))
    expect_equal(as.matrix(DimReduction(obj)$Scores), ref$x, ignore_attr = TRUE)
})

#----- RunPCA respects a user-chosen HVG set
test_that("RunPCA uses a manually set HVG vector", {
    obj <- create_test_geneda()
    obj@HVGs <- c("Low5", "High3", "Low40", "High9")

    obj <- suppressMessages(RunPCA(obj))
    expect_setequal(rownames(DimReduction(obj)$Eigenvectors), obj@HVGs)
})

#----- generatePCs still selects the top features from numeric variances
test_that("generatePCs selects top features when given numeric variances", {
    obj <- create_test_geneda()
    mat <- getNormalized(obj)
    vars <- MatrixGenerics::rowVars(mat)

    res <- suppressMessages(capture.output(pcs <- generatePCs(mat, vars, 10)))
    expect_setequal(rownames(pcs$Eigenvectors), paste0("High", 1:10))
})

#----- percent_var is stored as numbers, not strings
test_that("RunPCA stores percent_var as a named numeric percent", {
    obj <- create_test_geneda()
    obj <- suppressMessages(RunPCA(obj, nfeatures = 20))
    pv <- DimReduction(obj)$percent_var
    pca <- stats::prcomp(t(getNormalized(obj)[HVGs(obj), ]))

    expect_type(pv, "double")
    expect_equal(names(pv), paste0("PC", seq_along(pv)))
    expect_equal(sum(pv), 100)
    expect_equal(unname(pv), pca$sdev^2 / sum(pca$sdev^2) * 100)
})

test_that("validity rejects character percent_var", {
    obj <- create_test_geneda()
    obj <- suppressMessages(RunPCA(obj, nfeatures = 20))
    obj@DimReduction$percent_var <- c(PC1 = "50 %")

    expect_error(validObject(obj), "percent_var' must be numeric")
})

test_that("PlotPCA and PlotScree use numeric percent_var", {
    obj <- create_test_geneda()
    obj <- suppressMessages(RunPCA(obj, nfeatures = 20))
    pv <- DimReduction(obj)$percent_var

    p <- PlotPCA(obj, x = "PC1", y = "PC2", color_by = "condition")
    expect_equal(p$labels$x, sprintf("PC1: %.2f%%", pv[["PC1"]]))

    s <- PlotScree(obj)
    expect_equal(s$data$Variance, unname(pv))
})

#----- PlotPCA uses tidy evaluation, not aes_string()
test_that("PlotPCA maps columns without deprecation warnings", {
    obj <- create_test_geneda()
    obj <- suppressMessages(RunPCA(obj, nfeatures = 20))
    meta <- getMetadata(obj)
    meta$batch <- rep(c("b1", "b2"), 3)
    obj@metadata <- meta

    expect_no_warning(
        p <- PlotPCA(obj, x = 1, y = 3, color_by = "condition",
            shape_by = "batch")
    )
    built <- ggplot2::ggplot_build(p)
    expect_equal(built$data[[1]]$x, DimReduction(obj)$Scores$PC1,
        ignore_attr = TRUE)
    expect_equal(built$data[[1]]$y, DimReduction(obj)$Scores$PC3,
        ignore_attr = TRUE)
    expect_equal(length(unique(built$data[[1]]$shape)), 2L)
})
