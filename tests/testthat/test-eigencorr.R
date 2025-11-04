library(testthat)
library(RGenEDA)
library(ggplot2)

# Helper function to create a minimal geneda object with DimReduction
create_test_geneda <- function() {
  # normalized matrix (genes x samples)
  mat <- matrix(1:12, nrow = 3)
  colnames(mat) <- c("S1", "S2", "S3", "S4")
  rownames(mat) <- c("Gene1", "Gene2", "Gene3")

  # metadata (samples x metadata)
  meta <- data.frame(condition = c("A", "A", "B", "B"),
                     batch = c("X", "X", "Y", "Y"))
  rownames(meta) <- colnames(mat)

  # create geneda object
  obj <- GenEDA(normalized = mat, metadata = meta)

  # add DimReduction slot after creation
  loadings <- data.frame(PC1 = 1:4, PC2 = c(2,1,4,3))
  rownames(loadings) <- colnames(mat)
  obj@DimReduction <- list(Loadings = loadings)

  return(obj)
}

test_that("eigencorr returns correct structure", {
  obj <- create_test_geneda()

  res <- eigencorr(obj, num_pcs = 2)

  expect_type(res, "list")
  expect_true(all(c("cor_matrix", "pval_matrix", "stars", "plot") %in% names(res)))
  expect_equal(dim(res$cor_matrix), c(ncol(obj@metadata), 2))
  expect_equal(dim(res$pval_matrix), c(ncol(obj@metadata), 2))
  expect_s3_class(res$plot, "gg")
})

test_that("eigencorr subsets metadata columns", {
  obj <- create_test_geneda()

  res <- eigencorr(obj, num_pcs = 2, meta_cols = "batch")

  expect_equal(rownames(res$cor_matrix), "batch")
  expect_equal(rownames(res$pval_matrix), "batch")
})

test_that("eigencorr errors when metadata column missing", {
  obj <- create_test_geneda()

  expect_error(
    eigencorr(obj, meta_cols = "nonexistent"),
    "The following metadata columns were not found: nonexistent"
  )
})

test_that("eigencorr errors when DimReduction slot missing", {
  obj <- create_test_geneda()
  obj@DimReduction <- list()  # empty

  expect_error(
    eigencorr(obj),
    "DimReduction slot is empty. Please run RunPCA\\(\\) before calling eigencorr\\(\\)\\."
  )
})

test_that("eigencorr computes stars correctly", {
  obj <- create_test_geneda()

  res <- eigencorr(obj, num_pcs = 2)

  expect_equal(dim(res$stars), dim(res$cor_matrix))
  expect_true(all(res$stars %in% c("", "*", "**", "***")))
})

test_that("eigencorr plot aesthetics", {
  obj <- create_test_geneda()

  res <- eigencorr(obj, num_pcs = 2)

  p <- res$plot
  expect_s3_class(p, "gg")

  # check that x and y aesthetics are PC and Meta
  aes_mapping <- ggplot2::ggplot_build(p)$plot$mapping
  expect_true("x" %in% names(aes_mapping))
  expect_true("y" %in% names(aes_mapping))
})
