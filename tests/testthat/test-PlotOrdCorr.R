library(testthat)
library(RGenEDA)

# Helper: positive "expression" values (Bray-Curtis needs non-negative data)
# with one strong group effect, so NMDS axis 1 tracks the condition
create_ord_geneda <- function() {
    set.seed(3)
    n <- 12
    grp <- rep(c("A", "B"), each = n / 2)
    mat <- matrix(rpois(200 * n, lambda = 20), nrow = 200)
    mat[1:20, grp == "B"] <- mat[1:20, grp == "B"] + 15
    colnames(mat) <- paste0("S", seq_len(n))
    rownames(mat) <- paste0("Gene", 1:200)
    meta <- data.frame(
        condition = grp,
        noise1 = rnorm(n),
        noise2 = rnorm(n),
        row.names = colnames(mat)
    )
    GenEDA(normalized = mat, metadata = meta)
}

test_that("PlotOrdCorr adjusts p-values across all axis x covariate cells", {
    skip_if_not_installed("vegan")
    set.seed(1)
    res <- PlotOrdCorr(create_ord_geneda(), num_mds = 2)

    expect_true(all(c("cor_matrix", "padj_matrix", "pval_matrix", "stars",
        "p_adjust_method", "plot", "variance_explained") %in% names(res)))
    expect_equal(res$p_adjust_method, "BH")
    expect_equal(dimnames(res$padj_matrix), dimnames(res$pval_matrix))
    expect_equal(
        as.vector(res$padj_matrix),
        stats::p.adjust(as.vector(res$pval_matrix), method = "BH")
    )
})

test_that("PlotOrdCorr stars and plot use adjusted p-values", {
    skip_if_not_installed("vegan")
    set.seed(1)
    res <- PlotOrdCorr(create_ord_geneda(), num_mds = 2,
        p_adjust_method = "bonferroni")
    expected <- ifelse(res$padj_matrix < 0.001, "***",
        ifelse(res$padj_matrix < 0.01, "**",
            ifelse(res$padj_matrix < 0.05, "*", "")))
    expect_equal(res$stars, expected)

    pd <- res$plot$data
    expect_equal(
        pd$PAdj,
        res$padj_matrix[cbind(
            match(pd$Meta, rownames(res$padj_matrix)),
            match(pd$Axis, colnames(res$padj_matrix))
        )]
    )
    expect_match(res$plot$labels$caption, "bonferroni")
})

test_that("PlotOrdCorr does not allow p_adjust_method = 'none'", {
    skip_if_not_installed("vegan")
    expect_error(
        PlotOrdCorr(create_ord_geneda(), p_adjust_method = "none"),
        "must be one of"
    )
})

test_that("PlotOrdCorr rejects unknown p_adjust_method", {
    skip_if_not_installed("vegan")
    expect_error(
        PlotOrdCorr(create_ord_geneda(), p_adjust_method = "magic")
    )
})
