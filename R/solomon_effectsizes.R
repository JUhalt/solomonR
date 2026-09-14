#' Hedges g for a two-group contrast with CI
#' @references
#' Hedges, L. V. (1981). Distribution theory for Glass's estimator of effect
#' size and related estimators. *Journal of Educational Statistics, 6*(2),
#' 107-128.
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
