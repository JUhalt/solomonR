# Power and design planning for Solomon four-group designs.
#
# The data-generating mechanism follows the protocol posted on issue #18:
# every participant has a latent baseline X ~ N(0, 1) that only pretested
# participants observe, so structural pretest absence is part of the design.

.solomon_power_cells <- function(n) {
  if (is.list(n)) n <- unlist(n)
  if (length(n) == 1L) n <- rep(n, 4L)
  if (length(n) != 4L) {
    stop("n must be a single cell size or four cell sizes (n1, n2, n3, n4).")
  }
  if (is.null(names(n)) || !all(c("n1", "n2", "n3", "n4") %in% names(n))) {
    names(n) <- c("n1", "n2", "n3", "n4")
  }
  n <- n[c("n1", "n2", "n3", "n4")]
  if (any(!is.finite(n)) || any(n != round(n)) || any(n < 2)) {
    stop("Each Solomon cell needs at least 2 participants, given as whole numbers.")
  }
  storage.mode(n) <- "integer"
  n
}

.check_power_inputs <- function(delta, rho, sens, sigma, alpha, sims) {
  if (!is.finite(delta) || !is.finite(sens)) {
    stop("delta and sens must be finite numbers.")
  }
  if (!is.finite(rho) || abs(rho) > 1) {
    stop("rho must be between -1 and 1.")
  }
  if (!is.finite(sigma) || sigma <= 0) {
    stop("sigma must be positive.")
  }
  if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be between 0 and 1.")
  }
  if (!is.finite(sims) || sims < 1) {
    stop("sims must be a positive number of replications.")
  }
  invisible(TRUE)
}

.simulate_solomon_power <- function(cells, delta, rho, sens, sigma) {
  treat <- rep(c(1L, 0L, 1L, 0L), times = cells)
  pretested <- rep(c(1L, 1L, 0L, 0L), times = cells)
  n_total <- sum(cells)
  x <- stats::rnorm(n_total)
  e <- stats::rnorm(n_total)
  y <- delta * treat + sens * treat * pretested +
    sigma * (rho * x + sqrt(1 - rho^2) * e)
  data.frame(
    y_post = y,
    treat = treat,
    pretested = pretested,
    y_pre = ifelse(pretested == 1L, x, NA_real_)
  )
}

# Two-sided power of a t test with the given effect, standard error, and df.
.power_noncentral_t <- function(effect, se, df, alpha) {
  if (!is.finite(se) || se <= 0 || !is.finite(df) || df < 1) return(NA_real_)
  ncp <- effect / se
  crit <- stats::qt(1 - alpha / 2, df)
  stats::pt(-crit, df, ncp) + stats::pt(crit, df, ncp, lower.tail = FALSE)
}

# Normal-theory power for the four Solomon contrasts. The pretested contrast is
# an ANCOVA comparison, so its residual variance is reduced by the squared
# pretest-posttest correlation; the contrasts that combine pretest conditions
# use Welch-Satterthwaite degrees of freedom.
.solomon_power_analytic <- function(n = 50, delta = 0.3, rho = 0.5, sens = 0,
                                    sigma = 1, alpha = 0.05) {
  cells <- .solomon_power_cells(n)
  .check_power_inputs(delta, rho, sens, sigma, alpha, 1)

  var_pre <- sigma^2 * (1 - rho^2) * (1 / cells[["n1"]] + 1 / cells[["n2"]])
  var_un <- sigma^2 * (1 / cells[["n3"]] + 1 / cells[["n4"]])
  df_pre <- cells[["n1"]] + cells[["n2"]] - 3
  df_un <- cells[["n3"]] + cells[["n4"]] - 2
  df_combined <- (var_pre + var_un)^2 /
    (var_pre^2 / max(df_pre, 1) + var_un^2 / max(df_un, 1))

  effects <- c(
    "ATE (avg over pretest)" = delta + sens / 2,
    "Pretest x Treatment" = sens,
    "Treatment | pretested" = delta + sens,
    "Treatment | unpretested" = delta
  )
  ses <- c(
    sqrt(var_pre + var_un) / 2,
    sqrt(var_pre + var_un),
    sqrt(var_pre),
    sqrt(var_un)
  )
  dfs <- c(df_combined, df_combined, df_pre, df_un)

  data.frame(
    estimand = names(effects),
    true_effect = unname(effects),
    std.error = ses,
    df = dfs,
    power = mapply(.power_noncentral_t, unname(effects), ses, dfs,
                   MoreArgs = list(alpha = alpha)),
    row.names = NULL
  )
}

#' Power simulation for Solomon designs
#'
#' Simulates normally distributed Solomon four-group data and reports the
#' rejection rate of each Solomon test at the chosen `alpha`, with its Monte
#' Carlo standard error.
#'
#' Every participant has a latent baseline drawn from a standard normal
#' distribution, which only pretested participants observe, so structural
#' pretest absence is part of the design rather than missing data. The
#' posttest residual standard deviation is `sigma` in every cell, and the
#' pretest-posttest correlation among pretested participants is `rho`. The
#' treatment effect is `delta` among unpretested participants and
#' `delta + sens` among pretested participants, so the equal-weighted average
#' treatment effect is `delta + sens / 2`.
#'
#' The reported `power` is a rejection rate: when `true_effect` is zero it
#' estimates the Type I error rather than power.
#'
#' Tests are computed from the same functions users would call, so reported
#' power reflects the package's default inference: `fit_solomon_glm()` with
#' HC3 standard errors and t reference distributions, the 2x2 ANOVA
#' interaction, and the historical one-tailed Test I (Walton Braver & Braver, 1988).
#'
#' @section Validation:
#' A pre-specified simulation study of 126 scenarios and 315,000 replications
#' checked this function against normal-theory benchmarks. The protocol and
#' its amendment were posted to issue #18 before any results were examined,
#' and the article "Validating power_solomon()" on the package website reports
#' the study in full.
#'
#' - The 2x2 ANOVA interaction agreed with its analytic benchmark in every
#'   scenario, with mean differences of 0.002 or less.
#' - Test I held its nominal size under the complete null, rejecting between
#'   4.4% and 5.6% of the time at alpha = .05.
#' - Rejection rates were invariant to the residual scale, and no fit failed.
#' - Rejection rates from the unified GLM are conservative in small samples,
#'   because the HC3 standard errors the package uses by default are
#'   conservative there. Type I error averaged 0.041 with 10 participants per
#'   cell and 0.049 with 100, and simulated power fell below the normal-theory
#'   benchmark by 0.030 on average with 10 per cell (up to 0.083 in one
#'   scenario), 0.017 with 20, and 0.003 with 100.
#'
#' With 20 or fewer participants per cell, treat GLM-based power as a
#' conservative figure rather than an exact one.
#'
#' Only continuous outcomes are simulated. Binary and count outcomes are not
#' supported; for binary outcomes see [marginal_solomon()] and its simulation
#' validation on issue #43.
#'
#' @param n Cell sizes: a single number used for all four cells, or four
#'   sizes given as a list or vector with elements `n1` (pretested treatment),
#'   `n2` (pretested control), `n3` (unpretested treatment), and `n4`
#'   (unpretested control).
#' @param delta Treatment effect among unpretested participants, on the
#'   posttest scale.
#' @param rho Pretest-posttest correlation among pretested participants.
#' @param sens Sensitization: the additional treatment effect among pretested
#'   participants (0 = none).
#' @param sigma Posttest residual standard deviation in all cells.
#' @param sims Number of Monte Carlo replications.
#' @param stouffer Logical; if `TRUE`, also report the historical Test I
#'   rejection rate, evaluated one-tailed (treatment > control) as in
#'   [fit_solomon_classic()].
#' @param alpha Significance level. Default is 0.05.
#' @param seed Optional integer seed. The global random number state is
#'   restored afterwards.
#'
#' @return A data frame with one row per test, giving the estimand, the test
#'   used, the true effect, the rejection rate and its Monte Carlo standard
#'   error, the replications that produced a usable fit, the number of
#'   failures, and `alpha`. Test I rows report `NA` for `true_effect`, because
#'   the combination targets a directional hypothesis rather than a single
#'   contrast, and all Test I values are `NA` when `stouffer = FALSE`.
#'
#' @references
#' Morris, T. P., White, I. R., & Crowther, M. J. (2019). Using simulation
#' studies to evaluate statistical methods. *Statistics in Medicine, 38*(11),
#' 2074–2102. https://doi.org/10.1002/sim.8086
#'
#' Walton Braver, M. C., & Braver, S. L. (1988). Statistical treatment of the
#' Solomon four-group design: A meta-analytic approach. *Psychological Bulletin,
#' 104*(1), 150–154. https://doi.org/10.1037/0033-2909.104.1.150
#'
#' @examples
#' power_solomon(n = 30, delta = 0.5, sens = 0.2, sims = 50, seed = 1)
#'
#' @export
power_solomon <- function(n = 50,
                          delta = 0.3,
                          rho = 0.5,
                          sens = 0,
                          sigma = 1,
                          sims = 2000,
                          stouffer = TRUE,
                          alpha = 0.05,
                          seed = NULL) {

  cells <- .solomon_power_cells(n)
  sims <- as.integer(sims)
  .check_power_inputs(delta, rho, sens, sigma, alpha, sims)

  contrasts <- c(
    "ATE (avg over pretest)",
    "Pretest x Treatment",
    "Treatment | pretested",
    "Treatment | unpretested"
  )

  one_replicate <- function() {
    d <- .simulate_solomon_power(cells, delta, rho, sens, sigma)

    glm_fit <- fit_solomon_glm(d$y_post, d$treat, d$pretested, d$y_pre,
                               robust = "HC3")
    glm_p <- glm_fit$effects$p.value[match(contrasts, glm_fit$effects$contrast)]

    aov_tab <- stats::anova(
      stats::lm(y_post ~ factor(treat) * factor(pretested), data = d)
    )
    aov_p <- aov_tab[["Pr(>F)"]][match("factor(treat):factor(pretested)",
                                       rownames(aov_tab))]

    stouffer_p <- NA_real_
    if (isTRUE(stouffer)) {
      pre <- d[d$pretested == 1L, ]
      un <- d[d$pretested == 0L, ]
      ancova <- stats::lm(y_post ~ treat + y_pre, data = pre)
      posttest <- stats::lm(y_post ~ treat, data = un)
      t_pre <- summary(ancova)$coefficients["treat", "t value"]
      t_un <- summary(posttest)$coefficients["treat", "t value"]
      p_one <- c(
        stats::pt(t_pre, df = stats::df.residual(ancova), lower.tail = FALSE),
        stats::pt(t_un, df = stats::df.residual(posttest), lower.tail = FALSE)
      )
      eps <- .Machine$double.eps
      p_one <- pmin(pmax(p_one, eps), 1 - eps)
      stouffer_p <- stouffer_solomon(p_one)$p_meta_one_tailed
    }

    c(glm_p, aov_p, stouffer_p)
  }

  run <- function() {
    rejections <- matrix(NA_real_, nrow = sims, ncol = 6L)
    failures <- 0L
    for (s in seq_len(sims)) {
      p <- tryCatch(one_replicate(), error = function(e) NULL)
      if (is.null(p)) {
        failures <- failures + 1L
        next
      }
      rejections[s, ] <- as.numeric(p < alpha)
    }
    list(rejections = rejections, failures = failures)
  }

  res <- if (is.null(seed)) run() else withr::with_seed(seed, run())

  successful <- colSums(!is.na(res$rejections))
  power <- colMeans(res$rejections, na.rm = TRUE)
  power[successful == 0L] <- NA_real_
  mcse <- sqrt(power * (1 - power) / successful)

  data.frame(
    estimand = c(contrasts, "Pretest x Treatment", "Treatment (one-sided)"),
    test = c(rep("GLM (HC3, t)", 4L), "2x2 ANOVA interaction",
             "Test I (Walton Braver & Braver, 1988)"),
    true_effect = c(delta + sens / 2, sens, delta + sens, delta, sens, NA_real_),
    power = unname(power),
    mcse = unname(mcse),
    sims = unname(successful),
    failures = res$failures,
    alpha = alpha,
    row.names = NULL
  )
}
