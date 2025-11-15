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

  p <- PlotVolcano(obj, "DEG", alpha = 0.05, fc = 1, "Den", "Num")

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
