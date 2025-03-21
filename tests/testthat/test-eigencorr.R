test_that("eigencorr produces correct output structure", {
  # Sample input data
  MAT <- matrix(rnorm(100), nrow = 10)
  rownames(MAT) <- paste0("Gene", 1:10)
  colnames(MAT) <- paste0("Sample", 1:10)

  # Output
  output_path <- "/users/mike/Desktop/"

  META <- data.frame(Group = rep(c("A", "B"), each = 5))
  rownames(META) <- colnames(MAT)

  OUTPUT <- tempdir()  # Temporary output directory

  # Run function
  result <- eigencorr(MAT, META, NUM_PCS = 5, output_path)

  # Check if result is a list
  expect_type(result, "list")

  # Check if required elements exist
  expect_true(all(c("cor_matrix", "pval_matrix", "stars") %in% names(result)))

  # Check matrix dimensions
  expect_equal(dim(result$cor_matrix), c(ncol(META), 5))
  expect_equal(dim(result$pval_matrix), c(ncol(META), 5))
})
