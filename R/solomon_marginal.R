# Marginal Solomon contrasts for binary outcomes (issue #43).

.marginal_scales <- c(
  difference = "Risk difference",
  ratio = "Risk ratio",
  odds_ratio = "Odds ratio"
)

# Marginal (standardized) risks for the four Solomon cells from logistic
# coefficients `b`: each participant's risk is predicted under treatment and
# under control, and averaged within their pretest condition (Daniel et al.,
# 2021; Localio et al., 2007). Returns r(t, p) as a named vector.
.marginal_risks <- function(b, X, pretested) {
  b[is.na(b)] <- 0
  b_t <- b[["treat"]]
  b_tp <- if ("treat:pretested" %in% names(b)) b[["treat:pretested"]] else 0
  # Linear predictor with the treatment terms removed; setting treatment to t
  # then adds t * b_t, plus t * b_tp in the pretested condition.
  base <- drop(X %*% b) - X[, "treat"] * b_t -
    if (b_tp != 0) X[, "treat:pretested"] * b_tp else 0
  pre <- pretested == 1L
  c(
    t1p1 = mean(stats::plogis(base[pre] + b_t + b_tp)),
    t0p1 = mean(stats::plogis(base[pre])),
    t1p0 = mean(stats::plogis(base[!pre] + b_t)),
    t0p0 = mean(stats::plogis(base[!pre]))
  )
}

# The four Solomon contrasts on one scale from the four marginal risks. The
# average treatment effect compares risks averaged over equal numbers of
# pretested and unpretested participants. Ratio scales are undefined when a
# risk is 0 or 1 (within 1e-8, as under separation).
.marginal_contrasts <- function(r, scale) {
  g <- switch(scale, difference = identity, ratio = log, odds_ratio = stats::qlogis)
  if (scale != "difference" && .degenerate_risk(r)) {
    return(stats::setNames(rep(NA_real_, 4), .solomon_contrast_order))
  }
  pre <- g(r[["t1p1"]]) - g(r[["t0p1"]])
  un <- g(r[["t1p0"]]) - g(r[["t0p0"]])
  ate <- g((r[["t1p1"]] + r[["t1p0"]]) / 2) - g((r[["t0p1"]] + r[["t0p0"]]) / 2)
  stats::setNames(c(ate, pre - un, pre, un), .solomon_contrast_order)
}

.marginal_all <- function(b, X, pretested, scales) {
  r <- .marginal_risks(b, X, pretested)
  unlist(lapply(scales, function(s) .marginal_contrasts(r, s)), use.names = FALSE)
}

.degenerate_risk <- function(r) any(r < 1e-8 | r > 1 - 1e-8)

# Logistic regression by iteratively reweighted least squares, with the
# convergence rule of stats::glm.fit() (relative change in deviance below
# `epsilon`). A lean version for bootstrap refits, warm-started at `start`.
# Fitted risks that reach 0 or 1 in floating point (near separation) make the
# working response or the deviance non-finite; the fit is then reported as
# not converged, which the bootstrap counts as a failed resample.
.logistic_irls <- function(X, y, start, epsilon = 1e-8, maxit = 25L) {
  failed <- list(converged = FALSE)
  events <- y == 1
  deviance <- function(mu) -2 * (sum(log(mu[events])) + sum(log1p(-mu[!events])))
  b <- start
  eta <- drop(X %*% b)
  mu <- stats::plogis(eta)
  dev <- deviance(mu)
  if (!is.finite(dev)) return(failed)
  for (iter in seq_len(maxit)) {
    w <- mu * (1 - mu)
    z <- eta + (y - mu) / w
    if (!all(is.finite(z))) return(failed)
    b <- tryCatch(solve(crossprod(X, w * X), crossprod(X, w * z))[, 1],
                  error = function(e) NULL)
    if (is.null(b) || !all(is.finite(b))) return(failed)
    eta <- drop(X %*% b)
    mu <- stats::plogis(eta)
    dev_new <- deviance(mu)
    if (!is.finite(dev_new)) return(failed)
    if (abs(dev_new - dev) / (abs(dev_new) + 0.1) < epsilon) {
      return(list(coefficients = b, fitted.values = mu, converged = TRUE))
    }
    dev <- dev_new
  }
  list(coefficients = b, fitted.values = mu, converged = FALSE)
}

# A Solomon cell whose outcomes are all 0 or all 1 has no maximum-likelihood
# estimate on the logit scale (the fit stops at a risk near, not at, 0 or 1).
.cell_separated <- function(y, treat, pretested) {
  any(tapply(y, list(treat, pretested), function(v) all(v == 0) || all(v == 1)))
}

# A failed logistic fit: no convergence, a separated Solomon cell, or a fitted
# probability within 1e-8 of 0 or 1.
.logistic_failed <- function(fit, y, treat, pretested) {
  !isTRUE(fit$converged) || .cell_separated(y, treat, pretested) ||
    any(fit$fitted.values < 1e-8 | fit$fitted.values > 1 - 1e-8)
}


#' Marginal Solomon contrasts for binary outcomes
#'
#' Estimates the Solomon contrasts for a binary outcome as risk differences,
#' risk ratios, or odds ratios, comparing marginal risks in every cell.
#'
#' A logistic model that adjusts for the pretest estimates, among pretested
#' participants, a treatment effect conditional on the pretest, but among
#' unpretested participants a marginal effect, because they have no pretest.
#' Odds ratios are noncollapsible: when the pretest predicts the outcome, a
#' conditional and a marginal odds ratio differ even without confounding
#' (Daniel et al., 2021). The logit-scale Pretest x Treatment contrast of
#' [fit_solomon_glm()] can therefore be nonzero when there is no
#' sensitization. This function compares like with like.
#'
#' Marginal risks are estimated by standardization (Daniel et al., 2021;
#' Localio et al., 2007). In the pretested cells, each pretested participant's
#' risk is predicted under treatment and under control from the fitted
#' logistic model and averaged over all pretested participants. In the
#' unpretested cells the same is done over the unpretested participants;
#' without covariates, these are the observed cell proportions. The contrasts
#' are then:
#'
#' - `Treatment | pretested` and `Treatment | unpretested`: the treatment
#'   effect within each pretest condition;
#' - `Pretest x Treatment`: their difference on the chosen scale (a difference
#'   of risk differences, or a ratio of risk ratios or of odds ratios);
#' - `ATE (avg over pretest)`: the effect in a population with equal numbers
#'   of pretested and unpretested participants, computed from the averaged
#'   risks.
#'
#' Sensitization depends on the scale: an effect can be modified on the
#' risk-difference scale and not on the ratio scale, or the reverse. Report
#' the scale with every result.
#'
#' Intervals use a nonparametric bootstrap by default, resampling participants
#' within each of the four Solomon cells and reporting percentile intervals;
#' Daniel et al. (2021) and Localio et al. (2007) found the bootstrap performed
#' better than the delta method, which is available as `method = "delta"`
#' (using the fit's covariance matrix, with covariates held fixed). p-values
#' use the bootstrap or delta-method standard error on the analysis scale:
#' risk differences, log risk ratios, or log odds ratios. Bootstrap resamples
#' in which the logistic fit fails (no convergence, a Solomon cell with only
#' events or only non-events, or a fitted probability within 1e-8 of 0 or 1)
#' are excluded and counted; when more than 10% fail, intervals are not
#' reported.
#'
#' The package's simulation validation of this function is described on
#' issue #43.
#'
#' @param fit A fit from [fit_solomon_glm()] with `family = binomial()` (logit
#'   link). Cluster-robust (CR2) fits are not yet supported.
#' @param scale One or more of `"difference"`, `"ratio"`, and `"odds_ratio"`.
#' @param method `"bootstrap"` (default) or `"delta"`.
#' @param R Number of bootstrap resamples (at least 99). Default 999.
#' @param seed Optional integer seed for the bootstrap. The global random number
#'   state is restored afterwards.
#' @param conf_level Confidence level. Defaults to the fit's.
#'
#' @return An object of class `solomon_marginal` with `effects` (one row per
#'   scale and contrast: estimate and interval on the reporting scale, standard
#'   error on the analysis scale, p-value), `risks` (the four marginal risks),
#'   and the settings, including the number of failed bootstrap resamples.
#'
#' @references
#' Daniel, R., Zhang, J., & Farewell, D. (2021). Making apples from oranges:
#' Comparing noncollapsible effect estimators and their standard errors after
#' adjustment for different covariate sets. *Biometrical Journal, 63*(3),
#' 528–557. https://doi.org/10.1002/bimj.201900297
#'
#' Localio, A. R., Margolis, D. J., & Berlin, J. A. (2007). Relative risks and
#' confidence intervals were easily computed indirectly from multivariable
#' logistic regression. *Journal of Clinical Epidemiology, 60*(9), 874–882.
#' https://doi.org/10.1016/j.jclinepi.2006.12.001
#'
#' @seealso [fit_solomon_glm()], [fisher_solomon()]
#'
#' @examples
#' d <- solomon_example
#' d$passed <- as.integer(d$y_post > 55)
#' fit <- with(d, fit_solomon_glm(passed, treat, pretested, y_pre,
#'                                family = binomial()))
#' marginal_solomon(fit, method = "delta")
#' \donttest{
#' marginal_solomon(fit, R = 499, seed = 1)
#' }
#'
#' @export
marginal_solomon <- function(fit, scale = c("difference", "ratio", "odds_ratio"),
                             method = c("bootstrap", "delta"), R = 999, seed = NULL,
                             conf_level = fit$conf_level) {

  if (!inherits(fit, "solomon_glm")) {
    stop("`fit` must come from fit_solomon_glm().", call. = FALSE)
  }
  family <- stats::family(fit$model)
  if (!identical(family$family, "binomial") || !identical(family$link, "logit")) {
    stop("`fit` must use family = binomial() with the logit link.", call. = FALSE)
  }
  if (identical(fit$robust, "CR2")) {
    stop("Cluster-robust (CR2) fits are not yet supported; cluster-level inference ",
         "for binary outcomes is planned with issues #19 and #46.", call. = FALSE)
  }
  scale <- match.arg(scale, several.ok = TRUE)
  method <- match.arg(method)
  .check_conf_level(conf_level)
  if (method == "bootstrap" && (!is.numeric(R) || length(R) != 1L || R < 99)) {
    stop("`R` must be a single number of at least 99.", call. = FALSE)
  }

  model <- fit$model
  X <- stats::model.matrix(model)
  y <- model$y
  pretested <- as.integer(X[, "pretested"])
  treat <- as.integer(X[, "treat"])
  b <- stats::coef(model)

  estimate <- .marginal_all(b, X, pretested, scale)
  risks <- .marginal_risks(b, X, pretested)
  separated <- .cell_separated(y, treat, pretested) || .degenerate_risk(risks)
  if (separated) {
    estimate[rep(scale, each = 4) != "difference"] <- NA_real_
  }
  if (any(scale != "difference") && separated) {
    warning(structure(
      class = c("solomonR_sparse_cell_warning", "warning", "condition"),
      list(message = paste0(
        "A Solomon cell has only events or only non-events, so ratio-scale ",
        "contrasts are undefined and reported as NA."
      ), call = NULL)
    ))
  }

  alpha <- 1 - conf_level
  failures <- NA_integer_

  if (method == "bootstrap") {
    cells <- split(seq_along(y), interaction(treat, pretested, drop = TRUE))
    sizes <- lengths(cells)
    ends <- cumsum(sizes)
    run <- function() {
      t(vapply(seq_len(R), function(i) {
        idx <- unlist(lapply(cells, function(rows) rows[sample.int(length(rows), replace = TRUE)]),
                      use.names = FALSE)
        # Resampling is within cells, so the cells occupy consecutive blocks.
        events <- diff(c(0, cumsum(y[idx])[ends]))
        if (any(events == 0 | events == sizes)) return(rep(NA_real_, length(estimate)))
        Xi <- X[idx, , drop = FALSE]
        refit <- .logistic_irls(Xi, y[idx], start = b)
        if (!isTRUE(refit$converged) ||
            any(refit$fitted.values < 1e-8 | refit$fitted.values > 1 - 1e-8)) {
          return(rep(NA_real_, length(estimate)))
        }
        .marginal_all(refit$coefficients, Xi, pretested[idx], scale)
      }, numeric(length(estimate))))
    }
    boot <- if (is.null(seed)) run() else withr::with_seed(seed, run())
    failed <- apply(boot, 1, function(row) all(is.na(row)))
    failures <- sum(failed)
    boot <- boot[!failed, , drop = FALSE]
    std.error <- apply(boot, 2, stats::sd, na.rm = TRUE)
    ci <- t(apply(boot, 2, stats::quantile, probs = c(alpha / 2, 1 - alpha / 2),
                  na.rm = TRUE, names = FALSE))
    if (failures > 0.1 * R) {
      ci[] <- NA_real_
      warning(structure(
        class = c("solomonR_bootstrap_warning", "warning", "condition"),
        list(message = sprintf(
          paste0("The logistic fit failed in %d of %d bootstrap resamples (more ",
                 "than 10%%), usually because of sparse cells; intervals are not ",
                 "reported."), failures, R), call = NULL)
      ))
    }
  } else {
    V <- as.matrix(fit$vcov)
    grad <- vapply(seq_along(b), function(j) {
      h <- 1e-6 * max(1, abs(b[[j]]))
      up <- b; up[j] <- up[j] + h
      down <- b; down[j] <- down[j] - h
      (.marginal_all(up, X, pretested, scale) - .marginal_all(down, X, pretested, scale)) / (2 * h)
    }, numeric(length(estimate)))
    grad <- matrix(grad, nrow = length(estimate))
    std.error <- sqrt(rowSums((grad %*% V) * grad))
    z <- stats::qnorm(1 - alpha / 2)
    ci <- cbind(estimate - z * std.error, estimate + z * std.error)
  }

  p.value <- 2 * stats::pnorm(-abs(estimate / std.error))
  scale_col <- rep(scale, each = 4)
  ratio <- scale_col != "difference"
  report <- function(x) ifelse(ratio, exp(x), x)

  effects <- data.frame(
    scale = unname(.marginal_scales[scale_col]),
    contrast = rep(.solomon_contrast_order, length(scale)),
    estimate = report(estimate),
    conf.low = report(ci[, 1]),
    conf.high = report(ci[, 2]),
    std.error = std.error,
    p.value = p.value,
    stringsAsFactors = FALSE
  )

  structure(
    list(
      effects = effects,
      risks = data.frame(
        cell = c("Pretested, treatment", "Pretested, control",
                 "Unpretested, treatment", "Unpretested, control"),
        treat = c(1L, 0L, 1L, 0L),
        pretested = c(1L, 1L, 0L, 0L),
        risk = unname(risks),
        stringsAsFactors = FALSE
      ),
      method = method,
      R = if (method == "bootstrap") R else NA_integer_,
      failures = failures,
      conf_level = conf_level,
      vcov = if (method == "delta") .solomon_vcov_label(fit) else NA_character_,
      pretest_adjusted = "pre_obs" %in% colnames(X)
    ),
    class = "solomon_marginal"
  )
}


#' @export
print.solomon_marginal <- function(x, digits = 3, ...) {
  cat("Marginal Solomon contrasts for a binary outcome\n")
  cat(if (x$pretest_adjusted) {
    "Risks standardized over the pretest among pretested participants.\n"
  } else {
    "Risks from the fitted cell probabilities (no pretest adjustment).\n"
  })
  cat(if (x$method == "bootstrap") {
    sprintf("%s%% percentile intervals from %d cell-stratified bootstrap resamples (%d failed).\n",
            format(100 * x$conf_level), x$R, x$failures)
  } else {
    sprintf("%s%% delta-method intervals (%s covariance).\n", format(100 * x$conf_level), x$vcov)
  })
  cat("\nMarginal risks:\n")
  risks <- x$risks[, c("cell", "risk")]
  risks$risk <- round(risks$risk, digits)
  print(risks, row.names = FALSE)
  for (s in unique(x$effects$scale)) {
    cat("\n", s, ":\n", sep = "")
    tab <- x$effects[x$effects$scale == s, c("contrast", "estimate", "conf.low", "conf.high", "p.value")]
    tab[c("estimate", "conf.low", "conf.high")] <- lapply(tab[c("estimate", "conf.low", "conf.high")],
                                                          round, digits)
    tab$p.value <- format.pval(tab$p.value, digits = 2)
    print(tab, row.names = FALSE)
  }
  cat("\nSensitization depends on the scale; report the scale with every result.\n")
  invisible(x)
}


#' Historical categorical analysis of a binary Solomon outcome
#'
#' Reproduces the categorical path described by El Karkri et al. (2025b) for
#' qualitative outcomes: the treatment comparison is tested separately among
#' pretested and unpretested participants with Fisher's exact test (Pearson's
#' chi-square test is also reported), and pretest sensitization is "judged to
#' exist if a significant effect is observed for pre-tested groups but not for
#' non-pre-tested groups" (p. 7).
#'
#' This is provided for teaching and replication. The rule compares
#' significance, not effects: when the treatment effect is the same in both
#' pretest conditions, it still declares sensitization whenever the pretested
#' comparison reaches significance and the unpretested one does not. As Gelman
#' and Stern (2006) put it, "even large changes in significance levels can
#' correspond to small, nonsignificant changes in the underlying quantities"
#' (p. 328). A test of sensitization compares the effects themselves; see
#' [marginal_solomon()].
#'
#' @param y Binary outcome coded 0/1 (or logical).
#' @param treat Treatment indicator coded 0/1 (or logical).
#' @param pretested Pretest indicator coded 0/1 (or logical).
#' @param alpha Significance level for the historical rule. Default 0.05.
#'
#' @return An object of class `solomon_fisher` with `tests` (one row per pretest
#'   condition and for both combined: counts, proportions, the uncorrected
#'   Pearson chi-square, and Fisher's exact p-value) and `sensitization` (the
#'   historical rule's verdict).
#'
#' @references
#' El Karkri, M., Quesada, A., & Romero-Ariza, M. (2025b). Methodological
#' aspects of the Solomon four-group design: Detecting pre-test sensitisation
#' and analysing qualitative and quantitative variables in education research.
#' *Review of Education, 13*(1), Article e70050.
#' https://doi.org/10.1002/rev3.70050
#'
#' Gelman, A., & Stern, H. (2006). The difference between "significant" and
#' "not significant" is not itself statistically significant. *The American
#' Statistician, 60*(4), 328–331. https://doi.org/10.1198/000313006X152649
#'
#' @seealso [marginal_solomon()]
#'
#' @examples
#' # Kvalem et al. (1996): condom use at most recent intercourse, 6 months.
#' kvalem <- data.frame(
#'   treat = c(1, 1, 0, 0), pretested = c(1, 0, 1, 0),
#'   events = c(51, 21, 76, 69), n = c(73, 49, 148, 133)
#' )
#' d <- kvalem[rep(1:4, kvalem$n), c("treat", "pretested")]
#' d$y <- unlist(lapply(1:4, function(i) {
#'   rep(c(1, 0), c(kvalem$events[i], kvalem$n[i] - kvalem$events[i]))
#' }))
#' with(d, fisher_solomon(y, treat, pretested))
#'
#' @export
fisher_solomon <- function(y, treat, pretested, alpha = 0.05) {
  y <- .solomon_indicator(y, "y")
  treat <- .solomon_indicator(treat, "treat")
  pretested <- .solomon_indicator(pretested, "pretested")
  .solomon_check_lengths(y = y, treat = treat, pretested = pretested)

  ok <- !is.na(y) & !is.na(treat) & !is.na(pretested)
  test_rows <- function(rows, label) {
    tab <- table(factor(treat[rows], levels = c(1, 0)), factor(y[rows], levels = c(1, 0)))
    chi <- suppressWarnings(stats::chisq.test(tab, correct = FALSE))
    data.frame(
      condition = label,
      events_treatment = tab[1, 1], n_treatment = sum(tab[1, ]),
      events_control = tab[2, 1], n_control = sum(tab[2, ]),
      risk_treatment = tab[1, 1] / sum(tab[1, ]),
      risk_control = tab[2, 1] / sum(tab[2, ]),
      chisq = unname(chi$statistic),
      chisq_p = chi$p.value,
      fisher_p = stats::fisher.test(tab)$p.value,
      stringsAsFactors = FALSE
    )
  }
  tests <- rbind(
    test_rows(ok & pretested == 1L, "Pretested"),
    test_rows(ok & pretested == 0L, "Unpretested"),
    test_rows(ok, "Combined")
  )
  rownames(tests) <- NULL

  sig_pre <- tests$fisher_p[1] < alpha
  sig_un <- tests$fisher_p[2] < alpha
  structure(
    list(
      tests = tests,
      sensitization = sig_pre && !sig_un,
      alpha = alpha
    ),
    class = "solomon_fisher"
  )
}


#' @export
print.solomon_fisher <- function(x, digits = 3, ...) {
  cat("Historical categorical Solomon analysis (El Karkri et al., 2025b)\n\n")
  tab <- x$tests[, c("condition", "events_treatment", "n_treatment", "events_control",
                     "n_control", "chisq", "fisher_p")]
  tab$chisq <- round(tab$chisq, 2)
  tab$fisher_p <- format.pval(tab$fisher_p, digits = 2)
  names(tab) <- c("Condition", "Events (T)", "n (T)", "Events (C)", "n (C)",
                  "Chi-square", "Fisher p")
  print(tab, row.names = FALSE)
  cat(sprintf("\nHistorical rule (significant among pretested, not among unpretested, alpha = %s): %s\n",
              format(x$alpha), if (x$sensitization) "sensitization declared" else "not declared"))
  cat(.wrap_lines(paste(
    "Caution: this rule compares significance, not effects, and is not a test of",
    "the Pretest x Treatment interaction. See marginal_solomon()."
  )), "\n", sep = "")
  invisible(x)
}
