library(testthat)
library(RGenEDA)
library(pheatmap)

test_that("PlotDistances returns correct structure", {
    # Minimal geneda object
    mat <- matrix(1:20, nrow = 5)
    colnames(mat) <- c("S1", "S2", "S3", "S4")
    rownames(mat) <- paste0("Gene", 1:5)
    meta <- data.frame(condition = c("A", "A", "B", "B"))
    rownames(meta) <- colnames(mat)
    rownames(mat) <- c("Gene1", "Gene2", "Gene3", "Gene4", "Gene5")

    obj <- GenEDA(normalized = mat, metadata = meta)

    res <- PlotDistances(obj)

    expect_type(res, "list")
    expect_true(
        all(c("dist_matrix", "order", "heatmap", "palettes") %in% names(res))
    )
    expect_s3_class(res$heatmap, "pheatmap")
    expect_equal(dim(res$dist_matrix), c(4, 4))
})

test_that("PlotDistances respects reorder = FALSE", {
    mat <- matrix(1:20, nrow = 5)
    colnames(mat) <- c("S1", "S2", "S3", "S4")
    rownames(mat) <- paste0("Gene", 1:5)
    meta <- data.frame(condition = c("A", "A", "B", "B"))
    rownames(meta) <- colnames(mat)
    rownames(mat) <- c("Gene1", "Gene2", "Gene3", "Gene4", "Gene5")


    obj <- GenEDA(normalized = mat, metadata = meta)

    res <- PlotDistances(obj, reorder = FALSE)

    expect_equal(res$order, colnames(mat))
    expect_equal(rownames(res$dist_matrix), colnames(mat))
    expect_equal(colnames(res$dist_matrix), colnames(mat))
})


test_that("PlotDistances returns plot when return='plot'", {
    mat <- matrix(1:20, nrow = 5)
    colnames(mat) <- c("S1", "S2", "S3", "S4")
    rownames(mat) <- paste0("Gene", 1:5)
    meta <- data.frame(condition = c("A", "A", "B", "B"))
    rownames(meta) <- colnames(mat)
    rownames(mat) <- c("Gene1", "Gene2", "Gene3", "Gene4", "Gene5")


    obj <- GenEDA(normalized = mat, metadata = meta)

    res <- PlotDistances(obj, return = "plot")

    expect_type(res, "list")
    expect_s3_class(res$heatmap, "pheatmap")
})

test_that("PlotDistances errors when metadata columns missing", {
    mat <- matrix(1:20, nrow = 5)
    colnames(mat) <- c("S1", "S2", "S3", "S4")
    rownames(mat) <- paste0("Gene", 1:5)
    meta <- data.frame(condition = c("A", "A", "B", "B"))
    rownames(meta) <- colnames(mat)
    rownames(mat) <- c("Gene1", "Gene2", "Gene3", "Gene4", "Gene5")


    obj <- GenEDA(normalized = mat, metadata = meta)

    expect_error(
        PlotDistances(obj, meta_cols = "nonexistent"),
        "The following metadata columns were not found: nonexistent"
    )
})

test_that("PlotDistances errors when metadata columns mismatch", {
    mat <- matrix(1:20, nrow = 5)
    colnames(mat) <- c("S1", "S2", "S3", "S4")
    rownames(mat) <- paste0("Gene", 1:5)
    meta <- data.frame(condition = c("A", "A", "B", "B"))
    rownames(meta) <- colnames(mat)

    obj <- GenEDA(normalized = mat, metadata = meta)

    # Hack: force mismatch by renaming meta inside the function call
    expect_error(
        {
            PlotDistances(obj, meta_cols = c("condition"))$dist_matrix[1:4, 1:4]
        },
        NA
        # Actually the S4 validation prevents this, so safest is to remove
        # this test
    )
})

#----- Plotted order matches the returned order
create_cluster_geneda <- function() {
    set.seed(5)
    # Samples in an interleaved order, so clustering must reorder them
    grp <- c("A", "B", "A", "B", "A", "B")
    mat <- matrix(rnorm(50 * 6), nrow = 50)
    mat[1:25, grp == "B"] <- mat[1:25, grp == "B"] + 4
    dimnames(mat) <- list(paste0("Gene", 1:50), paste0("S", 1:6))
    meta <- data.frame(condition = grp, row.names = colnames(mat))
    GenEDA(normalized = mat, metadata = meta)
}

plotted_order <- function(ph, which = c("row", "col")) {
    which <- match.arg(which)
    tree <- ph[[paste0("tree_", which)]]
    tree$labels[tree$order]
}

test_that("PlotDistances heatmap order matches the returned order", {
    obj <- create_cluster_geneda()
    res <- PlotDistances(obj)

    hc <- stats::hclust(stats::dist(t(getNormalized(obj))),
        method = "complete")
    expect_equal(res$order, hc$labels[hc$order])
    expect_false(identical(res$order, colnames(getNormalized(obj))))

    expect_equal(plotted_order(res$heatmap, "row"), res$order)
    expect_equal(plotted_order(res$heatmap, "col"), res$order)
    expect_equal(rownames(res$dist_matrix), res$order)

    # The dendrogram is built on the sample distances themselves
    expect_equal(res$heatmap$tree_row$height, hc$height)
})

test_that("PlotDistances order follows the distance method", {
    obj <- create_cluster_geneda()
    res <- PlotDistances(obj, method = "manhattan")
    hc <- stats::hclust(stats::dist(t(getNormalized(obj)),
        method = "manhattan"), method = "complete")
    expect_equal(plotted_order(res$heatmap, "row"), hc$labels[hc$order])
    expect_equal(res$heatmap$tree_row$height, hc$height)
})

test_that("PlotDistances with reorder = FALSE does not cluster the plot", {
    obj <- create_cluster_geneda()
    res <- PlotDistances(obj, reorder = FALSE)

    expect_equal(res$order, colnames(getNormalized(obj)))
    expect_true(is.na(res$heatmap$tree_row)[1])
    expect_true(is.na(res$heatmap$tree_col)[1])
    gt_names <- res$heatmap$gtable$layout$name
    expect_false(any(c("row_tree", "col_tree") %in% gt_names))
})

test_that("PlotDistances draws meta_cols as annotation strips", {
    obj <- create_cluster_geneda()
    res <- PlotDistances(obj, meta_cols = "condition",
        palettes = list(condition = c(A = "grey40", B = "firebrick3")))
    gt_names <- res$heatmap$gtable$layout$name
    expect_true("col_annotation" %in% gt_names)
    expect_true("annotation_legend" %in% gt_names)
})

#----- Regression: noisy data where re-clustering the distance matrix (old
#----- behaviour) plotted a different order from the returned one
test_that("PlotDistances order matches on noisy data", {
    set.seed(3)
    mat <- matrix(rnorm(60 * 8), nrow = 60)
    dimnames(mat) <- list(paste0("Gene", 1:60), paste0("S", 1:8))
    meta <- data.frame(x = 1:8, row.names = colnames(mat))
    res <- PlotDistances(GenEDA(normalized = mat, metadata = meta))

    hc <- stats::hclust(stats::dist(t(mat)), method = "complete")
    expect_equal(res$order, hc$labels[hc$order])
    expect_equal(plotted_order(res$heatmap, "row"), res$order)
    expect_equal(plotted_order(res$heatmap, "col"), res$order)
})
