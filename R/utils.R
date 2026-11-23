#' Test association between a score vector and one metadata variable
#'
#' Numeric variables and two-level categorical variables use Pearson correlation
#' (two-level variables are coded by factor level order). Categorical variables
#' with more than two levels use the multiple correlation R (square root of R^2
#' from `lm(score ~ variable)`) with the overall F-test p-value, so the result
#' does not depend on how levels are labelled or ordered.
#'
#' @param score Numeric vector of PC or ordination scores
#' @param var A metadata column aligned to `score`
#' @return Named numeric vector with `estimate` and `p.value` (NA when
#'   untestable)
#' @keywords internal
#' @noRd
.assocTest <- function(score, var) {
    untestable <- c(estimate = NA_real_, p.value = NA_real_)

    if (is.numeric(var)) {
        x <- var
    } else {
        f <- droplevels(as.factor(var))
        if (nlevels(f) < 2L) {
            return(untestable)
        }

        if (nlevels(f) > 2L) {
            ok <- stats::complete.cases(score, f)
            fit <- summary(stats::lm(score[ok] ~ f[ok]))
            fstat <- fit$fstatistic
            if (is.null(fstat) || !is.finite(fstat[["value"]]) ||
                fstat[["dendf"]] < 1) {
                return(untestable)
            }
            pval <- stats::pf(
                fstat[["value"]],
                fstat[["numdf"]],
                fstat[["dendf"]],
                lower.tail = FALSE
            )
            return(c(estimate = sqrt(fit$r.squared), p.value = unname(pval)))
        }
        x <- as.numeric(f)
    }

    ok <- stats::complete.cases(score, x)
    if (sum(ok) < 3L || stats::sd(x[ok]) == 0 || stats::sd(score[ok]) == 0) {
        return(untestable)
    }
    test <- stats::cor.test(score[ok], x[ok], method = "pearson")
    c(estimate = unname(test$estimate), p.value = test$p.value)
}

#' Format a numeric percent for labels and messages
#'
#' @param x Numeric percent (0-100), e.g. from `DimReduction$percent_var`
#' @param digits Number of decimal places
#' @return Character vector such as `"77.87%"`
#' @keywords internal
#' @noRd
.formatPercent <- function(x, digits = 2L) {
    paste0(formatC(x, format = "f", digits = digits), "%")
}

#' Validate a p-value adjustment method
#'
#' Allows every [stats::p.adjust.methods] entry except `"none"`: heatmaps must
#' always be based on adjusted p-values.
#'
#' @param method Character scalar
#' @return `method`, if valid
#' @keywords internal
#' @noRd
.checkAdjustMethod <- function(method) {
    allowed <- setdiff(stats::p.adjust.methods, "none")
    if (!is.character(method) || length(method) != 1L ||
        !method %in% allowed) {
        allowed_str <- paste(allowed, collapse = ", ")
        cli::cli_abort(
            "{.arg p_adjust_method} must be one of: {allowed_str}."
        )
    }
    method
}

#' Adjust a matrix of p-values and convert them to significance stars
#'
#' Applies [stats::p.adjust()] once across all cells of `pval_matrix`, so the
#' correction covers every axis x covariate test in the heatmap.
#'
#' @param pval_matrix Numeric matrix of raw p-values (no `NA`)
#' @param method A method accepted by `.checkAdjustMethod()`
#' @return List with `padj` (matrix, same dimnames as `pval_matrix`) and
#'   `stars` (character matrix: `"***"` < 0.001, `"**"` < 0.01, `"*"` < 0.05)
#' @keywords internal
#' @noRd
.adjustAndStar <- function(pval_matrix, method) {
    padj <- pval_matrix
    padj[] <- stats::p.adjust(as.vector(pval_matrix), method = method)
    stars <- ifelse(padj < 0.001, "***",
        ifelse(padj < 0.01, "**",
            ifelse(padj < 0.05, "*", "")
        )
    )
    list(padj = padj, stars = stars)
}

#' Caption describing the significance stars on a correlation heatmap
#' @param method The p-value adjustment method used
#' @return Character scalar
#' @keywords internal
#' @noRd
.starsCaption <- function(method) {
    paste0(
        "* padj < 0.05, ** < 0.01, *** < 0.001 (", method, " adjusted)"
    )
}
