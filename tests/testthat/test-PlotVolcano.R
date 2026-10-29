library(testthat)
library(ggplot2)

context("PlotVolcano tests")

# Helper function to create a dummy GenEDA object
create_test_geneda <- function() {
  obj <- methods::new("geneda")
  obj@normalized <- matrix(rnorm(1000), nrow = 100, ncol = 10)
  rownames(obj@normalized) <- paste0("Gene", 1:100)
  colnames(obj@normalized) <- paste0("Sample", 1:10)

  obj@metadata <- data.frame(Sample = colnames(obj@normalized))
  rownames(obj@metadata) <- colnames(obj@normalized)  # <-- THIS LINE FIXES IT

  obj@DEGs <- list()
  obj@HVGs <- character(0)
  obj
}

# Create a dummy DESeq2-like DEG table
dummy_deg <- data.frame(
  log2FoldChange = rnorm(100),
  padj = runif(100, 0, 1)
)
rownames(dummy_deg) <- paste0("Gene", 1:100)

#----- PlotVolcano errors if object is not GenEDA
test_that("PlotVolcano errors if object is not geneda", {
  dummy_obj <- list()
  expect_error(
    PlotVolcano(dummy_obj, "DEG", 0.05, 1, "Den", "Num"),
    "methods::is\\(\\object, \\\"geneda\\\"\\) is not TRUE"
  )
})

#----- PlotVolcano errors is assay is missing
test_that("PlotVolcano errors if assay missing", {
  obj <- create_test_geneda()
  expect_error(
    PlotVolcano(obj, "NONEXISTENT", 0.05, 1, "Den", "Num"),
    "Assay NONEXISTENT was not found in DEGs slot!"
  )
})

#----- PlotVolcano errors if DEG table is empty
test_that("PlotVolcano errors if DEG table empty", {
  obj <- create_test_geneda()
  obj <- SetDEGs(obj, dummy_deg, "DEG")
  obj@DEGs$DEG <- obj@DEGs$DEG[0, ]  # empty
  expect_error(
    PlotVolcano(obj, "DEG", 0.05, 1, "Den", "Num"),
    "No differential expression results found in Assay DEG"
  )
})

#----- PlotVolcano returns a ggplot object
test_that("PlotVolcano returns a ggplot object", {
  obj <- create_test_geneda()
  obj <- SetDEGs(obj, dummy_deg, "DEG")
  p <- PlotVolcano(obj, "DEG", 0.05, 1, "Den", "Num")
  expect_s3_class(p, "ggplot")
})

#----- PlotVolcano correctly assigns Groups
test_that("PlotVolcano correctly assigns Groups", {
  obj <- create_test_geneda()
  df <- data.frame(
    log2FoldChange = c(2, -0.5, -2, 0.2, 3),
    padj = c(0.01, 0.2, 0.03, 0.5, 0.001),
    stringsAsFactors = FALSE
  )
  rownames(df) <- paste0("Gene", 1:5)
  obj <- SetDEGs(obj, df, "DEG")

  p <- PlotVolcano(obj, "DEG", alpha = 0.05, l2fc = 1, "Den", "Num")

  groups <- p$data$Group
  expect_true(all(groups %in% c("Upregulated in Den",
                                "Upregulated in Num",
                                paste("padj <", 0.05),
                                "ns")))
})

#----- PlotVolcano respects title argument
test_that("PlotVolcano respects title argument", {
  obj <- create_test_geneda()
  obj <- SetDEGs(obj, dummy_deg, "DEG")
  p <- PlotVolcano(obj, "DEG", 0.05, 1, "Den", "Num", title = "My Volcano Plot")
  # Check that ggtitle has been added
  expect_true("My Volcano Plot" %in% p$labels$title)
})

#----- PlotVolcano works without num/den
test_that("PlotVolcano uses default group labels when num and den are NULL", {
  obj <- create_test_geneda()
  obj <- SetDEGs(obj, dummy_deg, "DEG")
  p <- PlotVolcano(obj, "DEG", 0.05, 1)
  expect_s3_class(p, "ggplot")
  expect_true(all(c("Upregulated in numerator", "Upregulated in denominator") %in% levels(p$data$Group)))
})

#----- PlotVolcano errors when num and den are the same
test_that("PlotVolcano errors when num and den are identical", {
  obj <- create_test_geneda()
  obj <- SetDEGs(obj, dummy_deg, "DEG")
  expect_error(PlotVolcano(obj, "DEG", 0.05, 1, "A", "A"),
               "num and den must be different group names.")
})

#----- Genes exactly on a threshold are kept and grouped like FilterDEGs
test_that("PlotVolcano keeps genes that sit exactly on the thresholds", {
  obj <- create_test_geneda()
  deg <- data.frame(
    log2FoldChange = c(1, -1, 0.5, 2),
    padj = c(0.01, 0.05, 0.05, 0.5),
    row.names = paste0("Gene", 1:4)
  )
  obj <- SetDEGs(obj, deg, "DEG")
  p <- PlotVolcano(obj, "DEG", 0.05, 1, "Den", "Num")
  groups <- setNames(as.character(p$data$Group), p$data$gene)

  expect_equal(nrow(p$data), 4)
  expect_equal(groups[["Gene1"]], "Upregulated in Num")
  expect_equal(groups[["Gene2"]], "Upregulated in Den")
  expect_equal(groups[["Gene3"]], "padj < 0.05")
  expect_equal(groups[["Gene4"]], "ns")
})

#----- NA values in unrelated columns do not drop genes
test_that("PlotVolcano keeps genes with NA in columns it does not use", {
  obj <- create_test_geneda()
  deg <- dummy_deg
  deg$symbol <- NA
  obj <- SetDEGs(obj, deg, "DEG")
  p <- PlotVolcano(obj, "DEG", 0.05, 1, "Den", "Num")
  expect_equal(nrow(p$data), nrow(deg))
})
