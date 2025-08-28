# Wald-test-based partial/semipartial R^2 for a 1-df linear contrast L %*% beta = 0
# Returns lower, r2, upper (approx via z CI transformed to R^2)
contrast_r2_ci <- function(model, L, vcov, conf = 0.95) {
  b <- stats::coef(model)
  est <- as.numeric(drop(L %*% b))
  se2 <- as.numeric(drop(L %*% vcov %*% t(L)))
  if (!is.finite(se2) || se2 <= 0) return(c(lower = NA, r2 = NA, upper = NA))
  z <- est / sqrt(se2)
  df <- model$df.residual
  r2 <- z^2 / (z^2 + df)
  zcrit <- stats::qnorm(1 - (1-conf)/2)
  zlo <- z - zcrit; zhi <- z + zcrit
  r2lo <- zlo^2/(zlo^2 + df); r2hi <- zhi^2/(zhi^2 + df)
  c(lower = max(0, r2lo), r2 = max(0, r2), upper = min(1, r2hi))
}
