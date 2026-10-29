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
  rownames(mat) <- c(paste0("Low", seq_len(nLow)), paste0("High", seq_len(nHigh)))
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
  vars <- apply(mat, 1L, stats::var)

  res <- suppressMessages(capture.output(pcs <- generatePCs(mat, vars, 10)))
  expect_setequal(rownames(pcs$Eigenvectors), paste0("High", 1:10))
})
