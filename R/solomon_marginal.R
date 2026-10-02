# Marginal Solomon contrasts for binary (issue #43) and count (issue #44)
# outcomes.

.marginal_scales <- c(
  difference = "Risk difference",
  ratio = "Risk ratio",
  odds_ratio = "Odds ratio"
)
.marginal_rate_scales <- c(
  difference = "Rate difference",
  ratio = "Rate ratio"
)

# Marginal (standardized) risks for the four Solomon cells from logistic
# coefficients `b`: each participant's risk is predicted under treatment and
# under control, and averaged within their pretest condition (Daniel et al.,
# 2021; Localio et al., 2007). Returns r(t, p) as a named vector.
.marginal_risks <- function(b, X, pretested, linkinv = stats::plogis) {
  b[is.na(b)] <- 0
  b_t <- b[["treat"]]
  b_tp <- if ("treat:pretested" %in% names(b)) b[["treat:pretested"]] else 0
  # Linear predictor with the treatment terms removed; setting treatment to t
  # then adds t * b_t, plus t * b_tp in the pretested condition.
  base <- drop(X %*% b) - X[, "treat"] * b_t -
    if (b_tp != 0) X[, "treat:pretested"] * b_tp else 0
  pre <- pretested == 1L
  c(
    t1p1 = mean(linkinv(base[pre] + b_t + b_tp)),
    t0p1 = mean(linkinv(base[pre])),
    t1p0 = mean(linkinv(base[!pre] + b_t)),
    t0p0 = mean(linkinv(base[!pre]))
  )
}

# The four Solomon contrasts on one scale from the four marginal risks. The
# average treatment effect compares risks averaged over equal numbers of
# pretested and unpretested participants. Ratio scales are undefined when a
# risk is 0 or 1 (within 1e-8, as under separation), or when a count
# outcome's rate is 0.
.marginal_contrasts <- function(r, scale, count = FALSE) {
  g <- switch(scale, difference = identity, ratio = log, odds_ratio = stats::qlogis)
  degenerate <- if (count) any(r < 1e-8) else .degenerate_risk(r)
  if (scale != "difference" && degenerate) {
    return(stats::setNames(rep(NA_real_, 4), .solomon_contrast_order))
  }
  pre <- g(r[["t1p1"]]) - g(r[["t0p1"]])
  un <- g(r[["t1p0"]]) - g(r[["t0p0"]])
  ate <- g((r[["t1p1"]] + r[["t1p0"]]) / 2) - g((r[["t0p1"]] + r[["t0p0"]]) / 2)
  stats::setNames(c(ate, pre - un, pre, un), .solomon_contrast_order)
}

.marginal_all <- function(b, X, pretested, scales, count = FALSE) {
  r <- .marginal_risks(b, X, pretested, linkinv = if (count) exp else stats::plogis)
  unlist(lapply(scales, function(s) .marginal_contrasts(r, s, count)), use.names = FALSE)
}

# A Solomon cell with no counts at all has no rate ratio (the log-link
# estimate diverges).
.cell_empty <- function(y, treat, pretested) {
  any(tapply(y, list(treat, pretested), function(v) all(v == 0)))
}

# Classed warning for link-scale Solomon contrasts from a noncollapsible link
# with the pretest as a covariate (#43).
.warn_noncollapsible <- function(link) {
  warning(structure(
    class = c("solomonR_noncollapsible_warning", "warning", "condition"),
    list(
      message = paste0(
        "With the ", link, " link and a pretest covariate, the Pretest x Treatment ",
        "contrast compares a treatment effect conditional on the pretest ",
        "(pretested participants) with a marginal one (unpretested participants). ",
        "These differ whenever the pretest predicts the outcome, even without ",
        "sensitization (Daniel et al., 2021). Use marginal_solomon() to compare ",
        "the effects on a common scale."
      ),
      call = NULL
    )
  ))
}

# Stop when a function that takes the fit of a four-group design receives
# the fit of a design with several treatments (issue #45). Used by
# marginal_solomon() and perm_solomon().
.stop_ngroup_fit <- function(fit, fun) {
  if (inherits(fit, "solomon_ngroup")) {
    stop(structure(
      class = c("solomonR_ngroup_unsupported", "error", "condition"),
      list(
        message = paste0(
          "`", fun, "()` takes the fit of a Solomon four-group design: one ",
          "treatment and a control, each with and without a pretest. `fit` ",
          "comes from a design with several treatments. For this analysis, ",
          "subset the data to one treatment and the control, and fit the ",
          "subset with `fit_solomon_glm()`."
        ),
        call = NULL
      )
    ))
  }
  invisible(NULL)
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


#' Marginal Solomon contrasts for binary and count outcomes
#'
#' `r lifecycle::badge("stable")`
#' Estimates the Solomon contrasts for a binary outcome as risk differences,
#' risk ratios, or odds ratios, and for a count outcome as rate differences or
#' rate ratios, comparing marginal risks or rates in every cell.
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
#' Count outcomes: for a fit with `family = poisson()` or
#' `family = "negative_binomial"`, rates per unit of
#' exposure are standardized in the same way and compared as rate differences
#' or rate ratios. Log-link rate ratios are collapsible (Daniel et al., 2021),
#' so they agree with the fitted model's coefficients when there are no other
#' covariates; rate differences depend on the covariate distribution. For
#' counts, intervals use the delta method with the fit's (by default robust
#' HC3) covariance, which Cameron and Trivedi (2013) recommend under
#' overdispersion; a bootstrap for counts has not been evaluated and is not
#' offered. The fit's Pearson dispersion statistic, and for negative-binomial
#' fits the estimated theta, are printed for description. No published
#' Solomon study with a count outcome has been identified, so this use of
#' count-data methods is a solomonR extension.
#'
#' Clustered fits: for binary outcomes fitted with `robust = "CR2"`, the
#' delta-method standard errors use the CR2 cluster-robust covariance (Bell &
#' McCaffrey, 2002), and intervals
#' and tests use the t distribution with Satterthwaite degrees of freedom for
#' the delta method's linear approximation (Pustejovsky & Tipton, 2018).
#' Applying those degrees of freedom to the linearized contrast is a solomonR
#' extension. Only the delta method is offered, because the bootstrap
#' resamples participants rather than clusters, and a classed warning flags
#' degrees of freedom below 4 (Tipton, 2015).
#'
#' In the package's simulation of cluster-randomized Solomon designs (issue
#' #64; 4 to 47 clusters per cell or arm, intracluster correlations of 0.02
#' and 0.10), these intervals met the coverage and Type I error tolerances
#' for risk differences in 140 of 144 scenario-contrasts and for odds ratios
#' in 137 of 144. The four risk-difference shortfalls, with coverage from
#' 0.935 to 0.961, all had an intracluster correlation of 0.10: 4 or 10
#' clusters per arm with pretesting assigned within clusters, 15 and 47
#' clusters per arm with pretesting within clusters, and 15, 15, 47, and 47
#' clusters per cell. Risk-ratio intervals met the tolerances in 111 of 144
#' and were conservative with 4 clusters per cell or arm (coverage up to
#' 0.973). A normal reference distribution undercovered (coverage down to
#' 0.880), so it is not used for clustered fits.
#'
#' In the three shortfalls with pretesting assigned within clusters, the
#' unweighted comparison of cluster-level summaries (Hayes & Moulton, 2017,
#' pp. 211–215) met the tolerances, so a cluster-level analysis is the
#' recommended alternative for risk differences in such designs when
#' clustering is strong. `method = "cluster_summary"` computes it, as the
#' study did: the unweighted mean of the cluster proportions in each arm,
#' compared with a t interval that uses separate variances and Satterthwaite
#' degrees of freedom. With pretesting assigned within clusters, each cluster
#' contributes its pretested and unpretested proportions, and treated and
#' control clusters are compared. It warns when an arm has fewer than four
#' clusters, the minimum Hayes and Moulton (2017, p. 128) recommend. The differences were small, and the cluster-level comparison met
#' the tolerances less often across all scenarios (119 of 144), mostly by
#' covering more than 96% of the time with 4 clusters per cell or arm.
#' [perm_solomon()] gives a cluster-level randomization test. Clustered count
#' fits have not been validated and are refused; for them, the log rate-ratio
#' contrasts of [fit_solomon_glm()] with CR2 covariance and [perm_solomon()]
#' account for clustering.
#'
#' The package's simulation validation of this function is described on
#' issue #43 (binary outcomes), issue #44 (count outcomes), issue #62
#' (negative-binomial fits), and issue #64 (clustered fits).
#'
#' @param fit A fit from [fit_solomon_glm()] with `family = binomial()` (logit
#'   link), `family = poisson()` (log link), or
#'   `family = "negative_binomial"`, with HC3 or model-based covariance, or,
#'   for binary outcomes, CR2 covariance. Designs with several treatments are
#'   not supported; see [fit_solomon_glm()].
#' @param scale One or more of `"difference"`, `"ratio"`, and, for binary
#'   outcomes, `"odds_ratio"`. For count outcomes the default is
#'   `c("difference", "ratio")`.
#' @param method `"bootstrap"` (default for binary outcomes), `"delta"`
#'   (the only method for count outcomes, and the default for clustered
#'   fits), or `"cluster_summary"` (clustered binary fits, risk differences
#'   only; see Clustered fits).
#' @param R Number of bootstrap resamples (at least 99). Default 999.
#' @param seed Optional integer seed for the bootstrap. The global random number
#'   state is restored afterwards.
#' @param conf_level Confidence level. Defaults to the fit's.
#'
#' @return An object of class `solomon_marginal` with `effects` (one row per
#'   scale and contrast: estimate and interval on the reporting scale, standard
#'   error on the analysis scale, p-value), `risks` (binary) or `rates`
#'   (counts, per unit of exposure) for the four cells, and the settings,
#'   including the number of failed bootstrap resamples.
#'
#' @references
#' Bell, R. M., & McCaffrey, D. F. (2002). Bias reduction in standard errors
#' for linear regression with multi-stage samples. *Survey Methodology, 28*(2),
#' 169–181.
#'
#' Cameron, A. C., & Trivedi, P. K. (2013). *Regression analysis of count data*
#' (2nd ed.). Cambridge University Press.
#' https://doi.org/10.1017/CBO9781139013567
#'
#' Daniel, R., Zhang, J., & Farewell, D. (2021). Making apples from oranges:
#' Comparing noncollapsible effect estimators and their standard errors after
#' adjustment for different covariate sets. *Biometrical Journal, 63*(3),
#' 528–557. https://doi.org/10.1002/bimj.201900297
#'
#' Hayes, R. J., & Moulton, L. H. (2017). *Cluster randomised trials* (2nd
#' ed.). Chapman and Hall/CRC. https://doi.org/10.4324/9781315370286
#'
#' Localio, A. R., Margolis, D. J., & Berlin, J. A. (2007). Relative risks and
#' confidence intervals were easily computed indirectly from multivariable
#' logistic regression. *Journal of Clinical Epidemiology, 60*(9), 874–882.
#' https://doi.org/10.1016/j.jclinepi.2006.12.001
#'
#' Pustejovsky, J. E., & Tipton, E. (2018). Small-sample methods for
#' cluster-robust variance estimation and hypothesis testing in fixed effects
#' models. *Journal of Business & Economic Statistics, 36*(4), 672–683.
#' https://doi.org/10.1080/07350015.2016.1247004
#'
#' Tipton, E. (2015). Small sample adjustments for robust variance estimation
#' with meta-regression. *Psychological Methods, 20*(3), 375–393.
#' https://doi.org/10.1037/met0000011
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
#' # A count outcome observed over different exposure times.
#' set.seed(2)
#' d$days <- runif(nrow(d), 5, 15)
#' d$visits <- rpois(nrow(d), d$days * exp(-1.5 + 0.3 * d$treat))
#' counts <- with(d, fit_solomon_glm(visits, treat, pretested, y_pre,
#'                                   family = poisson(), exposure = days))
#' marginal_solomon(counts)
#'
#' @export
marginal_solomon <- function(fit, scale = c("difference", "ratio", "odds_ratio"),
                             method = c("bootstrap", "delta", "cluster_summary"), R = 999,
                             seed = NULL, conf_level = fit$conf_level) {

  .stop_ngroup_fit(fit, "marginal_solomon")
  if (!inherits(fit, "solomon_glm")) {
    stop("`fit` must come from fit_solomon_glm().", call. = FALSE)
  }
  family <- stats::family(fit$model)
  binary <- identical(family$family, "binomial") && identical(family$link, "logit")
  count <- (identical(family$family, "poisson") || .is_negbin(family)) &&
    identical(family$link, "log")
  if (!binary && !count) {
    stop("`fit` must use family = binomial() (logit link), poisson() (log link), ",
         "or \"negative_binomial\".", call. = FALSE)
  }
  cr2 <- identical(fit$robust, "CR2")
  method_missing <- missing(method)
  if (cr2 && count) {
    stop("Clustered (CR2) fits are supported for binary outcomes, the case ",
         "validated in the package's simulation (issue #64). For counts, the ",
         "log rate-ratio contrasts of fit_solomon_glm(robust = \"CR2\") and ",
         "perm_solomon() account for clustering.", call. = FALSE)
  }
  if (!method_missing && identical(match.arg(method), "cluster_summary")) {
    if (!cr2 || !binary) {
      stop("method = \"cluster_summary\" needs a binary fit with robust = \"CR2\" ",
           "and its clusters.", call. = FALSE)
    }
    if (!missing(scale) && !identical(unique(scale), "difference")) {
      stop("Cluster-level summaries give risk differences only; use scale = \"difference\".",
           call. = FALSE)
    }
    .check_conf_level(conf_level)
    return(.marginal_cluster_summary(fit, conf_level))
  }
  if (cr2) {
    if (!method_missing && identical(match.arg(method), "bootstrap")) {
      stop("For clustered (CR2) fits the bootstrap is not offered; it resamples ",
           "participants, not clusters. Use method = \"delta\" or \"cluster_summary\".",
           call. = FALSE)
    }
    method <- "delta"
  }
  if (count && missing(scale)) scale <- c("difference", "ratio")
  scale <- match.arg(scale, several.ok = TRUE)
  method <- match.arg(method)
  if (count) {
    if ("odds_ratio" %in% scale) {
      stop("The odds-ratio scale applies only to binary outcomes.", call. = FALSE)
    }
    if (!method_missing && method == "bootstrap") {
      stop("For count outcomes only method = \"delta\" is offered; a bootstrap ",
           "for counts has not been evaluated.", call. = FALSE)
    }
    method <- "delta"
  }
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

  estimate <- .marginal_all(b, X, pretested, scale, count)
  risks <- .marginal_risks(b, X, pretested, linkinv = if (count) exp else stats::plogis)
  separated <- if (count) {
    .cell_empty(y, treat, pretested) || any(risks < 1e-8)
  } else {
    .cell_separated(y, treat, pretested) || .degenerate_risk(risks)
  }
  if (separated) {
    estimate[rep(scale, each = 4) != "difference"] <- NA_real_
  }
  if (any(scale != "difference") && separated) {
    warning(structure(
      class = c("solomonR_sparse_cell_warning", "warning", "condition"),
      list(message = if (count) {
        "A Solomon cell has no counts, so rate ratios are undefined and reported as NA."
      } else {
        paste0("A Solomon cell has only events or only non-events, so ratio-scale ",
               "contrasts are undefined and reported as NA.")
      }, call = NULL)
    ))
  }

  alpha <- 1 - conf_level
  failures <- NA_integer_
  df <- rep(Inf, length(estimate))

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
                 "reported. method = \"delta\" does not resample; in the package's ",
                 "simulation study (issue #43), its risk-difference intervals met ",
                 "the registered tolerances in most designs in which the bootstrap ",
                 "failed this often."), failures, R), call = NULL)
      ))
    }
  } else {
    V <- as.matrix(fit$vcov)
    grad <- vapply(seq_along(b), function(j) {
      h <- 1e-6 * max(1, abs(b[[j]]))
      up <- b; up[j] <- up[j] + h
      down <- b; down[j] <- down[j] - h
      (.marginal_all(up, X, pretested, scale, count) -
         .marginal_all(down, X, pretested, scale, count)) / (2 * h)
    }, numeric(length(estimate)))
    grad <- matrix(grad, nrow = length(estimate))
    std.error <- sqrt(rowSums((grad %*% V) * grad))
    if (cr2) {
      # Satterthwaite degrees of freedom for the delta method's linear
      # approximation (Pustejovsky & Tipton, 2018), from the complete-case
      # refit that the CR2 covariance was computed on.
      colnames(grad) <- names(b)
      fit_cr <- stats::glm(stats::formula(model), data = stats::model.frame(model),
                           family = fit$family)
      ok <- is.finite(estimate) & is.finite(std.error) & std.error > 0
      df[ok] <- vapply(which(ok), function(i) {
        as.data.frame(clubSandwich::linear_contrast(
          fit_cr, vcov = fit$vcov, contrasts = grad[i, , drop = FALSE], test = "Satterthwaite"
        ))$df
      }, numeric(1))
    }
    q <- ifelse(is.finite(df), stats::qt(1 - alpha / 2, pmax(df, 1e-8)), stats::qnorm(1 - alpha / 2))
    ci <- cbind(estimate - q * std.error, estimate + q * std.error)
  }

  p.value <- 2 * stats::pt(-abs(estimate / std.error), df)
  scale_col <- rep(scale, each = 4)
  ratio <- scale_col != "difference"
  report <- function(x) ifelse(ratio, exp(x), x)

  effects <- data.frame(
    scale = unname((if (count) .marginal_rate_scales else .marginal_scales)[scale_col]),
    contrast = rep(.solomon_contrast_order, length(scale)),
    estimate = report(estimate),
    conf.low = report(ci[, 1]),
    conf.high = report(ci[, 2]),
    std.error = std.error,
    df = df,
    p.value = p.value,
    stringsAsFactors = FALSE
  )

  if (cr2) {
    small_df <- is.finite(df) & df < 4
    if (any(small_df)) {
      .warn_cr2_small_df(paste(effects$scale[small_df], effects$contrast[small_df], sep = ": "),
                         df[small_df])
    }
  }

  cells <- data.frame(
    cell = c("Pretested, treatment", "Pretested, control",
             "Unpretested, treatment", "Unpretested, control"),
    treat = c(1L, 0L, 1L, 0L),
    pretested = c(1L, 1L, 0L, 0L),
    value = unname(risks),
    stringsAsFactors = FALSE
  )
  names(cells)[4] <- if (count) "rate" else "risk"

  structure(
    list(
      effects = effects,
      risks = if (!count) cells,
      rates = if (count) cells,
      outcome = if (count) "count" else "binary",
      dispersion = fit$dispersion,
      theta = fit$theta,
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
  count <- identical(x$outcome, "count")
  what <- if (count) "Rates" else "Risks"
  cat("Marginal Solomon contrasts for a", if (count) "count" else "binary", "outcome\n")
  cat(if (x$pretest_adjusted) {
    paste(what, "standardized over the pretest among pretested participants.\n")
  } else {
    paste(what, "from the fitted cells (no pretest adjustment).\n")
  })
  if (count && !is.null(x$theta)) {
    cat(sprintf("Negative-binomial (NB2) fit: theta = %.3g (alpha = 1/theta = %.3g).\n",
                x$theta[["theta"]], x$theta[["alpha"]]))
  } else if (count && is.finite(x$dispersion)) {
    cat(sprintf("Pearson dispersion: %.2f (values well above 1 indicate overdispersion).\n",
                x$dispersion))
  }
  cat(if (x$method == "bootstrap") {
    sprintf("%s%% percentile intervals from %d cell-stratified bootstrap resamples (%d failed).\n",
            format(100 * x$conf_level), x$R, x$failures)
  } else if (x$method == "cluster_summary") {
    sprintf(paste0("%s%% t intervals from unweighted cluster-level summaries, separate variances ",
                   "(Hayes & Moulton, 2017); design: %s; clusters: %s.\n"),
            format(100 * x$conf_level), x$design,
            paste(names(x$clusters), x$clusters, sep = " = ", collapse = ", "))
  } else {
    sprintf("%s%% delta-method intervals (%s covariance).\n", format(100 * x$conf_level), x$vcov)
  })
  if (count) {
    cat("\nMarginal rates per unit of exposure:\n")
    cells <- x$rates[, c("cell", "rate")]
    cells$rate <- round(cells$rate, digits)
  } else {
    cat("\nMarginal risks:\n")
    cells <- x$risks[, c("cell", "risk")]
    cells$risk <- round(cells$risk, digits)
  }
  print(cells, row.names = FALSE)
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
#' `r lifecycle::badge("stable")`
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
#' @param y_post Binary posttest outcome coded 0/1 (or logical).
#' @param treat Treatment indicator coded 0/1 (or logical). Designs with
#'   several treatments are not supported; see [fit_solomon_glm()].
#' @param pretested Pretest indicator coded 0/1 (or logical).
#' @param alpha Significance level for the historical rule. Default 0.05.
#' @param data Optional data frame. When supplied, the other data arguments
#'   are looked up in it first, as bare column names (`y_post = post`) or as
#'   strings (`y_post = "post"`).
#' @param y `r lifecycle::badge("deprecated")` Use `y_post`.
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
#' fisher_solomon(y, treat, pretested, data = d)
#'
#' @export
fisher_solomon <- function(y_post, treat, pretested, alpha = 0.05,
                           data = NULL,
                           y = deprecated()) {
  .solomon_data_args(
    data, c("y_post", "treat", "pretested", "y"),
    environment(), parent.frame()
  )
  .stop_ngroup_unsupported(treat, "fisher_solomon")
  if (lifecycle::is_present(y)) {
    .renamed_arg(!missing(y_post), "y", "y_post", "fisher_solomon")
    y_post <- y
  }
  y_post <- .solomon_indicator(y_post, "y_post")
  treat <- .solomon_indicator(treat, "treat")
  pretested <- .solomon_indicator(pretested, "pretested")
  .solomon_check_lengths(y_post = y_post, treat = treat, pretested = pretested)

  ok <- !is.na(y_post) & !is.na(treat) & !is.na(pretested)
  test_rows <- function(rows, label) {
    tab <- table(factor(treat[rows], levels = c(1, 0)), factor(y_post[rows], levels = c(1, 0)))
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
