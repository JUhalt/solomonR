# Internal helper for a Wald-based partial R-squared
#
# For an ordinary Gaussian linear model with conventional covariance,
# t^2 / (t^2 + df_residual) is the usual partial R-squared associated
# with a one-df contrast, and its confidence interval follows from the
# noncentral F distribution (Steiger, 2004).
#
# With a robust covariance matrix (for example HC3), the same
# transformation is retained as a descriptive Wald-based approximation.
# It should not be interpreted as an exact model-decomposition R-squared,
# and no confidence interval is reported because the noncentral F pivot
# does not apply.
#
# For non-Gaussian GLMs, this quantity is not reported.
contrast_r2_ci <- function(
    model,
    L,
    vcov,
    conf = 0.95,
    conventional = FALSE
) {

  unavailable <- function(type) {
    list(
      r2 = NA_real_,
      r2_lo = NA_real_,
      r2_hi = NA_real_,
      type = type
    )
  }

  fam <- stats::family(model)$family

  if (!identical(fam, "gaussian")) {
    return(unavailable("not_available"))
  }

  b <- stats::coef(model)

  L_full <- stats::setNames(
    numeric(length(b)),
    names(b)
  )

  common <- intersect(
    names(L),
    names(b)
  )

  L_full[common] <- L[common]

  V <- as.matrix(vcov)

  V <- V[
    names(b),
    names(b),
    drop = FALSE
  ]

  estimate <- sum(
    L_full * b
  )

  variance <- as.numeric(
    t(L_full) %*% V %*% L_full
  )

  if (!is.finite(variance) || variance <= 0) {
    return(unavailable("wald_partial"))
  }

  statistic <- estimate / sqrt(variance)

  df_residual <- stats::df.residual(model)

  if (!is.finite(df_residual) || df_residual <= 0) {
    return(unavailable("wald_partial"))
  }

  r2 <- statistic^2 /
    (statistic^2 + df_residual)

  if (isTRUE(conventional)) {
    ci <- .partial_r2_ci(statistic^2, 1, df_residual, conf)
    return(
      list(
        r2 = unname(r2),
        r2_lo = unname(ci[["lower"]]),
        r2_hi = unname(ci[["upper"]]),
        type = "partial"
      )
    )
  }

  list(
    r2 = unname(r2),
    r2_lo = NA_real_,
    r2_hi = NA_real_,
    type = "wald_partial"
  )
}
