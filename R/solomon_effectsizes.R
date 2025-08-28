#' Hedges g for a two-group contrast with CI
#' @keywords internal
hedges_g_ci <- function(m1, m0, s1, s0, n1, n0, conf = 0.95) {
  sp2 <- ((n1-1)*s1^2 + (n0-1)*s0^2) / (n1 + n0 - 2)
  sp  <- sqrt(sp2)
  d   <- (m1 - m0) / sp
  J   <- 1 - 3/(4*(n1+n0)-9)
  g   <- J * d
  se  <- sqrt((n1+n0)/(n1*n0) + g^2/(2*(n1+n0-2)))
  z   <- stats::qnorm(1 - (1-conf)/2)
  c(lower = g - z*se, g = g, upper = g + z*se)
}

#' Semi-partial R^2 and CI for a named coefficient using robust vcov
#' @keywords internal
spr2_ci <- function(model, term, vcov, conf = 0.95) {
  # Based on Wald z and relation to partial correlation
  coefs <- stats::coef(model)
  if (!term %in% names(coefs)) return(c(lower = NA, r2 = NA, upper = NA))
  b  <- unname(coefs[term]); se <- sqrt(vcov[term, term])
  z  <- b/se
  r2 <- z^2/(z^2 + model$df.residual)  # approx semi-partial R^2
  # CI via z CI transformed to r2 (delta approx)
  zcrit <- stats::qnorm(1 - (1-conf)/2)
  zlo <- z - zcrit; zhi <- z + zcrit
  r2lo <- zlo^2/(zlo^2 + model$df.residual)
  r2hi <- zhi^2/(zhi^2 + model$df.residual)
  c(lower = max(0, r2lo), r2 = max(0, r2), upper = min(1, r2hi))
}
