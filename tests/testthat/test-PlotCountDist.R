library(testthat)
library(RGenEDA)
library(ggplot2)

test_that("PlotCountDist returns a ggplot object", {
  # Minimal geneda object
  mat <- matrix(1:12, nrow = 3)
  colnames(mat) <- c("S1", "S2", "S3", "S4")
  meta <- data.frame(condition = c("A", "A", "B", "B"))
  rownames(meta) <- colnames(mat)
  rownames(mat) <- c("Gene1", "Gene2", "Gene3")

  obj <- GenEDA(normalized = mat, metadata = meta)

  p <- PlotCountDist(obj)

  expect_s3_class(p, "ggplot")
})

test_that("PlotCountDist handles faceting", {
  mat <- matrix(1:12, nrow = 3)
  colnames(mat) <- c("S1", "S2", "S3", "S4")
  meta <- data.frame(condition = c("A", "A", "B", "B"))
  rownames(meta) <- colnames(mat)
  rownames(mat) <- c("Gene1", "Gene2", "Gene3")


  obj <- GenEDA(normalized = mat, metadata = meta)

  p <- PlotCountDist(obj, split_by = "condition")

  expect_s3_class(p, "ggplot")

  # Check that facetting is applied
  expect_true("FacetWrap" %in% class(p$facet))
})

test_that("PlotCountDist warns for invalid split_by column", {
  mat <- matrix(1:12, nrow = 3)
  colnames(mat) <- c("S1", "S2", "S3", "S4")
  meta <- data.frame(condition = c("A", "A", "B", "B"))
  rownames(meta) <- colnames(mat)
  rownames(mat) <- c("Gene1", "Gene2", "Gene3")


  obj <- GenEDA(normalized = mat, metadata = meta)

  expect_warning(PlotCountDist(obj, split_by = "nonexistent"),
                 "Column nonexistent not found in metadata")
})
