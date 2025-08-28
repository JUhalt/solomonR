# Wald-test-based semi-partial R^2 for a 1-df linear contrast L %*% beta = 0
# L can be a named numeric vector (preferred) or a full-length numeric vector.
contrast_r2_ci <- function(model, L, vcov, conf = 0.95) {
  b <- stats::coef(model)
  p <- length(b)

  # Align L to coefficient vector by NAME; pad with 0 for absent terms
  if (is.null(names(L))) {
    if (length(L) != p) stop("Contrast vector L has no names and wrong length.")
    Lvec <- as.numeric(L)
    names(Lvec) <- names(b)
  } else {
    Lvec <- rep(0, p)
    names(Lvec) <- names(b)
    common <- intersect(names(L), names(b))
    Lvec[common] <- L[common]
  }

  # Ensure VCOV matches coef order
  V <- as.matrix(vcov)[names(b), names(b), drop = FALSE]

  # Variance of the contrast
  Lm <- matrix(Lvec, nrow = 1)
  est <- as.numeric(Lm %*% b)
  se2 <- as.numeric(Lm %*% V %*% t(Lm))
  if (!is.finite(se2) || se2 <= 0) return(c(lower = NA, r2 = NA, upper = NA))

  z  <- est / sqrt(se2)
  df <- model$df.residual
  r2 <- z^2 / (z^2 + df)

  zcrit <- stats::qnorm(1 - (1 - conf) / 2)
  zlo   <- z - zcrit
  zhi   <- z + zcrit
  r2lo  <- zlo^2 / (zlo^2 + df)
  r2hi  <- zhi^2 / (zhi^2 + df)

  c(lower = max(0, r2lo), r2 = max(0, r2), upper = min(1, r2hi))
}
