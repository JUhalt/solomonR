#' Hedges g for a two-group contrast with CI
#'
#' Returns Hedges' bias-corrected standardized mean difference and a
#' confidence interval for the population standardized mean difference from
#' the noncentral t distribution of the pooled-variance t statistic
#' (Cumming & Finch, 2001; Kelley, 2007).
#' @references
#' Cumming, G., & Finch, S. (2001). A primer on the understanding, use, and
#' calculation of confidence intervals that are based on central and
#' noncentral distributions. *Educational and Psychological Measurement,
#' 61*(4), 532-574.
#'
#' Hedges, L. V. (1981). Distribution theory for Glass's estimator of effect
#' size and related estimators. *Journal of Educational Statistics, 6*(2),
#' 107-128.
#'
#' Kelley, K. (2007). Confidence intervals for standardized effect sizes:
#' Theory, application, and implementation. *Journal of Statistical
#' Software, 20*(8), 1-24.
#' @keywords internal
hedges_g_ci <- function(m1, m0, s1, s0, n1, n0, conf = 0.95) {
  sp2 <- ((n1-1)*s1^2 + (n0-1)*s0^2) / (n1 + n0 - 2)
  sp  <- sqrt(sp2)
  d   <- (m1 - m0) / sp
  J   <- 1 - 3/(4*(n1+n0)-9)
  g   <- J * d
  t   <- d / sqrt(1/n1 + 1/n0)
  ci  <- .smd_ci(t, n1, n0, conf)
  c(lower = unname(ci[["lower"]]), g = g, upper = unname(ci[["upper"]]))
}
