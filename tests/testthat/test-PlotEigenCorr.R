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
  scores <- data.frame(PC1 = 1:4, PC2 = c(2,1,4,3))
  rownames(scores) <- colnames(mat)
  obj@DimReduction <- list(Scores = scores)

  return(obj)
}

test_that("PlotEigenCorr returns correct structure", {
  obj <- create_test_geneda()

  res <- PlotEigenCorr(obj, num_pcs = 2)

  expect_type(res, "list")
  expect_true(all(c("cor_matrix", "pval_matrix", "stars", "plot") %in% names(res)))
  expect_equal(dim(res$cor_matrix), c(ncol(obj@metadata), 2))
  expect_equal(dim(res$pval_matrix), c(ncol(obj@metadata), 2))
  expect_s3_class(res$plot, "gg")
})

test_that("PlotEigenCorr subsets metadata columns", {
  obj <- create_test_geneda()

  res <- PlotEigenCorr(obj, num_pcs = 2, meta_cols = "batch")

  expect_equal(rownames(res$cor_matrix), "batch")
  expect_equal(rownames(res$pval_matrix), "batch")
})

test_that("PlotEigenCorr errors when metadata column missing", {
  obj <- create_test_geneda()

  expect_error(
    PlotEigenCorr(obj, meta_cols = "nonexistent"),
    "The following metadata columns were not found: nonexistent"
  )
})

test_that("PlotEigenCorr errors when DimReduction slot missing", {
  obj <- create_test_geneda()
  obj@DimReduction <- list()  # empty

  expect_error(
    PlotEigenCorr(obj),
    "DimReduction slot is empty\\. Please run RunPCA\\(\\) before calling PlotEigenCorr\\(\\)\\."
  )
})

test_that("PlotEigenCorr computes stars correctly", {
  obj <- create_test_geneda()

  res <- PlotEigenCorr(obj, num_pcs = 2)

  expect_equal(dim(res$stars), dim(res$cor_matrix))
  expect_true(all(res$stars %in% c("", "*", "**", "***")))
})

test_that("PlotEigenCorr plot aesthetics", {
  obj <- create_test_geneda()

  res <- PlotEigenCorr(obj, num_pcs = 2)

  p <- res$plot
  expect_s3_class(p, "gg")

  # check that x and y aesthetics are PC and Meta
  aes_mapping <- ggplot2::ggplot_build(p)$plot$mapping
  expect_true("x" %in% names(aes_mapping))
  expect_true("y" %in% names(aes_mapping))
})

# Helper with a 3-level categorical variable for association tests
create_multilevel_geneda <- function(batch) {
  mat <- matrix(rnorm(6 * 5), nrow = 5)
  colnames(mat) <- paste0("S", 1:6)
  rownames(mat) <- paste0("Gene", 1:5)
  meta <- data.frame(batch = batch, sample_id = colnames(mat), row.names = colnames(mat))
  obj <- GenEDA(normalized = mat, metadata = meta)
  scores <- data.frame(PC1 = c(1, 1.2, 5, 5.3, 9, 9.1), PC2 = c(3, -1, 2, 0, -2, 1))
  rownames(scores) <- colnames(mat)
  obj@DimReduction <- list(Scores = scores)
  obj
}

#----- Results for multi-level factors do not depend on level labels
test_that("PlotEigenCorr is invariant to relabelling a multi-level factor", {
  res1 <- PlotEigenCorr(create_multilevel_geneda(c("a", "a", "b", "b", "c", "c")), num_pcs = 2)
  res2 <- PlotEigenCorr(create_multilevel_geneda(c("z", "z", "b", "b", "m", "m")), num_pcs = 2)

  expect_equal(res1$cor_matrix["batch", ], res2$cor_matrix["batch", ])
  expect_equal(res1$pval_matrix["batch", ], res2$pval_matrix["batch", ])
})

#----- Multi-level factors use the unsigned multiple correlation from lm()
test_that("PlotEigenCorr uses sqrt(R^2) and the F-test for multi-level factors", {
  obj <- create_multilevel_geneda(c("a", "a", "b", "b", "c", "c"))
  res <- PlotEigenCorr(obj, num_pcs = 2)

  fit <- summary(lm(obj@DimReduction$Scores$PC1 ~ factor(obj@metadata$batch)))
  fstat <- fit$fstatistic
  expect_equal(res$cor_matrix["batch", "PC1"], sqrt(fit$r.squared))
  expect_equal(res$pval_matrix["batch", "PC1"],
               unname(pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE)))
  expect_true(all(res$cor_matrix["batch", ] >= 0))
})

#----- Untestable variables are dropped from both matrices together
test_that("PlotEigenCorr drops untestable variables from cor and pval matrices", {
  obj <- create_multilevel_geneda(c("a", "a", "b", "b", "c", "c"))
  res <- PlotEigenCorr(obj, num_pcs = 2)

  # sample_id has one level per sample, so there are no residual degrees of freedom
  expect_false("sample_id" %in% rownames(res$cor_matrix))
  expect_identical(rownames(res$cor_matrix), rownames(res$pval_matrix))
  expect_identical(rownames(res$stars), rownames(res$cor_matrix))
})

#----- Two-level factors still use Pearson correlation
test_that("PlotEigenCorr keeps Pearson correlation for two-level factors", {
  obj <- create_test_geneda()
  res <- PlotEigenCorr(obj, num_pcs = 1, meta_cols = "condition")
  expected <- cor(obj@DimReduction$Scores$PC1, as.numeric(factor(obj@metadata$condition)))
  expect_equal(res$cor_matrix["condition", "PC1"], expected)
})
