# Internal helper for a Wald-based partial R-squared
#
# For an ordinary Gaussian linear model with conventional covariance,
# t^2 / (t^2 + df_residual) is the usual partial R-squared associated
# with a one-df contrast.
#
# With a robust covariance matrix (for example HC3), the same
# transformation is retained as a descriptive Wald-based approximation.
# It should not be interpreted as an exact model-decomposition R-squared.
#
# Confidence intervals are intentionally not produced here. Earlier
# versions transformed the endpoints of a confidence interval for the
# signed Wald statistic; that procedure can yield invalid R-squared
# intervals when the interval crosses zero.
#
# For non-Gaussian GLMs, this quantity is not reported.
contrast_r2_ci <- function(
    model,
    L,
    vcov,
    conf = 0.95
) {

  # Retain conf in the signature for backward compatibility.
  invisible(conf)

  fam <- stats::family(model)$family

  if (!identical(fam, "gaussian")) {
    return(
      list(
        r2 = NA_real_,
        r2_lo = NA_real_,
        r2_hi = NA_real_,
        type = "not_available"
      )
    )
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
    return(
      list(
        r2 = NA_real_,
        r2_lo = NA_real_,
        r2_hi = NA_real_,
        type = "wald_partial"
      )
    )
  }

  statistic <- estimate / sqrt(variance)

  df_residual <- stats::df.residual(model)

  if (!is.finite(df_residual) || df_residual <= 0) {
    return(
      list(
        r2 = NA_real_,
        r2_lo = NA_real_,
        r2_hi = NA_real_,
        type = "wald_partial"
      )
    )
  }

  r2 <- statistic^2 /
    (statistic^2 + df_residual)

  list(
    r2 = unname(r2),

    # Intentionally unavailable until a defensible CI procedure
    # is implemented.
    r2_lo = NA_real_,
    r2_hi = NA_real_,

    type = "wald_partial"
  )
}
