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
#' @return Named numeric vector with `estimate` and `p.value` (NA when untestable)
#' @keywords internal
#' @noRd
.assocTest <- function(score, var) {
  untestable <- c(estimate = NA_real_, p.value = NA_real_)

  if (is.numeric(var)) {
    x <- var
  } else {
    f <- droplevels(as.factor(var))
    if (nlevels(f) < 2L) return(untestable)

    if (nlevels(f) > 2L) {
      ok <- stats::complete.cases(score, f)
      fit <- summary(stats::lm(score[ok] ~ f[ok]))
      fstat <- fit$fstatistic
      if (is.null(fstat) || !is.finite(fstat[["value"]]) || fstat[["dendf"]] < 1) return(untestable)
      pval <- stats::pf(fstat[["value"]], fstat[["numdf"]], fstat[["dendf"]], lower.tail = FALSE)
      return(c(estimate = sqrt(fit$r.squared), p.value = unname(pval)))
    }
    x <- as.numeric(f)
  }

  ok <- stats::complete.cases(score, x)
  if (sum(ok) < 3L || stats::sd(x[ok]) == 0 || stats::sd(score[ok]) == 0) return(untestable)
  test <- stats::cor.test(score[ok], x[ok], method = "pearson")
  c(estimate = unname(test$estimate), p.value = test$p.value)
}
