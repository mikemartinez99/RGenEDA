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
  obj2 <- GenEDA(normalized = obj@normalized, metadata = obj@metadata, counts = counts)
  expect_equal(getCounts(obj2), counts)
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
  obj <- FindVariableFeatures(obj, nfeatures = 10)  # more than genes
  expect_equal(length(HVGs(obj)), nrow(obj@normalized))
  obj <- FindVariableFeatures(obj, nfeatures = 0)  # zero
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
                nrow = n_genes, ncol = n_samples)
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

  expect_true(all(c("Loadings", "Eigenvectors", "percent_var") %in% names(dr)))
  expect_s3_class(dr$Loadings, "data.frame")
  expect_s3_class(dr$Eigenvectors, "data.frame")
  expect_equal(nrow(dr$Loadings), ncol(obj@normalized))
  expect_equal(nrow(dr$Eigenvectors), length(HVGs(obj)))
  expect_equal(length(dr$percent_var), min(5, length(HVGs(obj))))
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
  expect_error(PlotPCA(obj, x = 1, y = 2, color_by = "nonexistent_column"),
               "Column nonexistent_column not found in PCA data")
})

test_that("PlotPCA with custom colors works", {
  obj <- create_test_geneda_large()
  obj <- RunPCA(obj, nfeatures = 2000)
  cols <- c("A" = "red", "B" = "blue")
  p <- PlotPCA(obj, x = 1, y = 2, color_by = "condition", colors = cols)
  expect_s3_class(p, "ggplot")
})
