library(testthat)
library(RGenEDA)

# Helper to create geneda object with DimReduction including Eigenvectors
create_test_geneda <- function() {
  # normalized matrix (genes x samples)
  mat <- matrix(1:12, nrow = 3)
  colnames(mat) <- c("S1", "S2", "S3", "S4")
  rownames(mat) <- c("Gene1", "Gene2", "Gene3")

  # metadata
  meta <- data.frame(condition = c("A", "A", "B", "B"))
  rownames(meta) <- colnames(mat)

  # create geneda object
  obj <- GenEDA(normalized = mat, metadata = meta)

  # Add DimReduction slot
  eigenvectors <- data.frame(
    PC1 = c(0.5, -0.5, 0.7),
    PC2 = c(-0.3, 0.6, 0.2)
  )
  rownames(eigenvectors) <- rownames(mat)
  percent_var <- c(PC1 = "50%", PC2 = "30%")

  obj@DimReduction <- list(
    Eigenvectors = eigenvectors,
    percent_var = percent_var
  )

  return(obj)
}

test_that("extractEigen returns a data.frame", {
  obj <- create_test_geneda()

  res <- extractEigen(obj, "PC1")

  expect_s3_class(res, "data.frame")
  expect_true(all(c("Gene", "EigenVector", "PctVar") %in% colnames(res)))
})

test_that("extractEigen computes PctVar correctly", {
  obj <- create_test_geneda()

  res <- extractEigen(obj, "PC1")

  expect_equal(sum(res$PctVar), 100)
})

test_that("extractEigen preserves gene names", {
  obj <- create_test_geneda()

  res <- extractEigen(obj, "PC1")

  expect_equal(res$Gene, rownames(obj@DimReduction$Eigenvectors))
})

test_that("extractEigen returns correct EigenVector values", {
  obj <- create_test_geneda()

  res <- extractEigen(obj, "PC2")

  expect_equal(res$EigenVector, obj@DimReduction$Eigenvectors$PC2)
})

test_that("extractEigen warns when DimReduction slot is empty", {
  obj <- create_test_geneda()
  obj@DimReduction <- list()

  expect_message(res <- extractEigen(obj, "PC1"),
                 "DimReduction slot is empty. Please use RunPCA")
})
