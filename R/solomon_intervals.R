# Internal confidence-interval helpers.
#
# Each interval follows the inferential framework of the method that
# produced the effect summary:
# - Wald intervals use the same reference distribution as the reported test
#   (t with the reported degrees of freedom, or the normal when df = Inf).
# - Standardized mean differences use the noncentral t method (Cumming &
#   Finch, 2001; Kelley, 2007).
# - Partial R-squared for a one-df effect in a conventional Gaussian linear
#   model uses the noncentral F method (Steiger, 2004).

.check_conf_level <- function(conf_level) {
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      is.na(conf_level) || conf_level <= 0 || conf_level >= 1) {
    stop("`conf_level` must be a single number between 0 and 1.", call. = FALSE)
  }
  invisible(conf_level)
}


# Wald interval with a t (df finite) or normal (df = Inf) critical value.
.wald_ci <- function(estimate, std.error, df = Inf, conf_level = 0.95) {
  crit <- stats::qt(1 - (1 - conf_level) / 2, df = df)
  cbind(
    conf.low = estimate - crit * std.error,
    conf.high = estimate + crit * std.error
  )
}


# Noncentrality parameter at which cdf(ncp) equals `prob`. The cumulative
# probability of an observed t or F statistic decreases as ncp increases.
.ncp_solve <- function(cdf, prob, lower, upper) {
  f <- function(ncp) cdf(ncp) - prob
  while (f(lower) < 0) lower <- lower - max(1, abs(lower))
  while (f(upper) > 0) upper <- upper + max(1, abs(upper))
  stats::uniroot(f, c(lower, upper), tol = 1e-10)$root
}


# Confidence limits for the noncentrality parameter of a t statistic.
.ncp_t_limits <- function(t, df, conf_level = 0.95) {
  alpha <- 1 - conf_level
  cdf <- function(ncp) suppressWarnings(stats::pt(t, df = df, ncp = ncp))
  width <- 10 + abs(t)
  c(
    lower = .ncp_solve(cdf, 1 - alpha / 2, t - width, t + width),
    upper = .ncp_solve(cdf, alpha / 2, t - width, t + width)
  )
}


# Confidence interval for the population standardized mean difference of
# two independent groups, from the noncentral t distribution.
.smd_ci <- function(t, n1, n0, conf_level = 0.95) {
  .ncp_t_limits(t, df = n1 + n0 - 2, conf_level = conf_level) *
    sqrt(1 / n1 + 1 / n0)
}


# Confidence interval for partial R-squared of an effect with `df1`
# numerator degrees of freedom in a conventional Gaussian linear model:
# limits for the noncentrality parameter lambda are converted with
# lambda / (lambda + df1 + df2 + 1).
.partial_r2_ci <- function(F, df1, df2, conf_level = 0.95) {
  alpha <- 1 - conf_level
  cdf <- function(ncp) suppressWarnings(stats::pf(F, df1, df2, ncp = ncp))

  limit <- function(prob) {
    if (cdf(0) <= prob) {
      return(0)
    }
    .ncp_solve(cdf, prob, 0, max(10, 2 * F * df1))
  }

  lambda <- c(lower = limit(1 - alpha / 2), upper = limit(alpha / 2))
  lambda / (lambda + df1 + df2 + 1)
}
