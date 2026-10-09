# Longitudinal Solomon designs: a mixed model for repeated measures
# (issue #57).

.mmrm_fallbacks <- c(toeph = "heterogeneous Toeplitz", ar1h = "heterogeneous AR(1)",
                     csh = "heterogeneous compound symmetry")

# Fit the MMRM to the observed rows and estimate the Solomon contrasts at
# each occasion. `grouped` estimates the covariance separately for pretested
# and unpretested participants. Used by fit_solomon_mmrm() and by the
# validation study.
#
# `center`, when given, holds the variance of the estimated mean pretest at
# which the pretest is centered and its degrees of freedom (issue #104). The
# pretest effects at an occasion move by that occasion's pretest slope for
# each unit the mean moves, so they gain slope^2 times that variance, with
# Welch-Satterthwaite degrees of freedom; under the model the mean is
# independent of the regression estimates.
.mmrm_fit_contrasts <- function(observed, occasions, pre, grouped, df, conf_level,
                                center = NULL) {
  fixed <- if (pre) "y ~ occ * treat * pretested + occ:pre_obs" else "y ~ occ * treat * pretested"
  method <- if (df == "kenward-roger") "Kenward-Roger" else "Satterthwaite"
  fit_one <- function(type) {
    f <- stats::as.formula(sprintf("%s + %s(occ | %sid)", fixed, type,
                                   if (grouped) "pgrp / " else ""))
    # Kenward-Roger with the covariance parameterized by its elements (the
    # "linear" computation), which reproduces SAS PROC MIXED; see
    # ?fit_solomon_mmrm.
    fit <- tryCatch(
      if (df == "kenward-roger") {
        mmrm::mmrm(f, data = observed, reml = TRUE, method = method,
                   vcov = "Kenward-Roger-Linear")
      } else {
        mmrm::mmrm(f, data = observed, reml = TRUE, method = method)
      },
      error = function(e) NULL
    )
    if (!is.null(fit) && !identical(as.integer(mmrm::component(fit, "convergence")), 0L)) fit <- NULL
    fit
  }

  structure_used <- "unstructured"
  fit <- fit_one("us")
  if (is.null(fit)) {
    # Mallinckrodt et al. (2008, p. 312): if the unstructured model does not
    # converge, fit the listed alternatives and keep the converged one with
    # the smallest AIC.
    fits <- lapply(names(.mmrm_fallbacks), fit_one)
    ok <- !vapply(fits, is.null, logical(1))
    if (!any(ok)) {
      stop("The model did not converge with an unstructured covariance or any of the ",
           "fallback structures (", paste(.mmrm_fallbacks, collapse = ", "), ").", call. = FALSE)
    }
    aic <- vapply(fits[ok], stats::AIC, numeric(1))
    best <- which(ok)[which.min(aic)]
    fit <- fits[[best]]
    structure_used <- unname(.mmrm_fallbacks[best])
  }

  beta <- mmrm::component(fit, "beta_est")
  tt <- stats::delete.response(stats::terms(stats::as.formula(fixed)))
  point <- function(t, tr, pr, pre_obs = 0) {
    nd <- data.frame(occ = factor(occasions[t], levels = occasions), treat = tr, pretested = pr,
                     pre_obs = pre_obs)
    stats::model.matrix(tt, nd)[1, ]
  }
  weights <- function(v) {
    missing_terms <- setdiff(names(v)[abs(v) > 0], names(beta))
    if (length(missing_terms)) {
      stop("A Solomon contrast cannot be estimated because these terms are not estimable: ",
           paste(missing_terms, collapse = ", "), ". Check that every group is observed at ",
           "every occasion.", call. = FALSE)
    }
    out <- stats::setNames(numeric(length(beta)), names(beta))
    out[intersect(names(v), names(beta))] <- v[intersect(names(v), names(beta))]
    out
  }
  # `slope`: the pretest slope at the occasion of a pretest effect, which
  # adds the variance of the estimated center.
  one <- function(L, occasion, contrast, slope = 0) {
    r <- mmrm::df_1d(fit, weights(L))
    est <- r$est
    se <- r$se
    dfc <- r$df
    stat <- r$t_stat
    p <- r$p_val
    if (!is.null(center) && slope != 0) {
      v_center <- slope^2 * center$variance
      dfc <- .welch_df(c(se^2, v_center), c(dfc, center$df))
      se <- sqrt(se^2 + v_center)
      stat <- est / se
      p <- 2 * stats::pt(-abs(stat), dfc)
    }
    crit <- stats::qt(1 - (1 - conf_level) / 2, dfc)
    data.frame(occasion = occasion, contrast = contrast, estimate = est, std.error = se,
               statistic = stat, df = dfc, p.value = p,
               conf.low = est - crit * se, conf.high = est + crit * se,
               stringsAsFactors = FALSE)
  }
  sens <- list()
  rows <- list()
  for (t in seq_along(occasions)) {
    pre_t <- point(t, 1, 1) - point(t, 0, 1)
    un_t <- point(t, 1, 0) - point(t, 0, 0)
    sens[[t]] <- pre_t - un_t
    # The pretest effects (issue #104): pretested minus unpretested
    # participants, with the pretest at pre_obs = 0, its center; `slope` is
    # the pretest slope at this occasion.
    control_t <- point(t, 0, 1) - point(t, 0, 0)
    treated_t <- point(t, 1, 1) - point(t, 1, 0)
    slope <- if (pre) sum(weights(point(t, 0, 1, pre_obs = 1) - point(t, 0, 1)) * beta) else 0
    rows <- c(rows, list(
      one((pre_t + un_t) / 2, occasions[t], "ATE (avg over pretest)"),
      one(sens[[t]], occasions[t], "Pretest x Treatment"),
      one(pre_t, occasions[t], "Treatment | pretested"),
      one(un_t, occasions[t], "Treatment | unpretested"),
      one(control_t, occasions[t], "Pretest effect | control", slope),
      one(treated_t, occasions[t], "Pretest effect | treated", slope),
      one((control_t + treated_t) / 2, occasions[t], "Pretest main effect", slope)
    ))
  }
  last <- length(occasions)
  rows <- c(rows, list(one(sens[[last]] - sens[[1]],
                           sprintf("%s vs %s", occasions[last], occasions[1]),
                           "Change in Pretest x Treatment")))
  list(effects = .effects_table(do.call(rbind, rows), keys = "occasion"), model = fit,
       covariance = structure_used)
}

#' Mixed model for repeated measures in a longitudinal Solomon design
#'
#' `r lifecycle::badge("experimental")`
#' Analyzes a Solomon four-group design with several posttest occasions by
#' a mixed model for repeated measures (MMRM), and estimates the four
#' Solomon contrasts and the pretest effects at each occasion, and the change
#' in pretest sensitization from the first occasion to the last. The model is
#' likelihood-based, so it is valid when posttests are missing at random,
#' for example when participants drop out depending on their earlier
#' scores (Fitzmaurice et al., 2011, pp. 497, 505). Per-occasion analyses of
#' the participants still observed, and unweighted generalized estimating
#' equations, require the stronger assumption that the posttests are missing
#' completely at random (Liang & Zeger, 1986; Fitzmaurice et al., 2011, pp.
#' 358, 498).
#'
#' @section Model:
#' The specification follows the example of Mallinckrodt et al. (2008, p.
#' 312), adapted to the Solomon design:
#' - **Fixed effects.** Occasion x Treatment x Pretested, all categorical.
#' - **Pretest adjustment.** The pretest enters as in [fit_solomon_glm()]:
#'   in the pretested groups, the score minus the mean pretest of the
#'   pretested participants in the model (each counted once; returned as
#'   `pretest_mean`), and 0 in the unpretested groups, with a separate slope
#'   at each occasion ("a full interaction of the covariate with time", p.
#'   312). Pretests absent by design are never imputed.
#' - **Pretest effects.** At each occasion, the pretest effects among
#'   controls, among treated participants, and their average compare
#'   pretested and unpretested participants at that mean pretest, as
#'   described in the section "The pretest effect" of [fit_solomon_glm()].
#'   Pretest effects may fade with time (Entwisle, 1961, p. 610), and these
#'   estimates show whether they do. Their standard errors include the
#'   sampling variance of the mean pretest, s^2 / n for n pretested
#'   participants, times the square of that occasion's pretest slope, and
#'   their degrees of freedom combine those of the contrast with n - 1 by the
#'   Welch-Satterthwaite formula (Satterthwaite, 1946; Welch, 1947). The
#'   pretest effects were not part of the validation study below.
#' - **Covariance.** Unstructured within participant, estimated by
#'   restricted maximum likelihood (Laird & Ware, 1982). With a pretest, it
#'   is estimated separately for pretested and unpretested participants: the
#'   pretest adjustment makes the residual covariance conditional on the
#'   pretest in the pretested groups only, so a single shared matrix would
#'   be misspecified.
#' - **Fallback.** If the unstructured model does not converge, the
#'   heterogeneous Toeplitz, heterogeneous AR(1), and heterogeneous compound
#'   symmetry structures are fitted, and the converged one with the smallest
#'   AIC is used (Mallinckrodt et al., 2008, p. 312). The output names the
#'   structure used.
#' - **Inference.** t tests and intervals with Kenward-Roger degrees of
#'   freedom and adjusted covariance (Kenward & Roger, 1997, as cited in
#'   Fitzmaurice et al., 2011, p. 101), as Mallinckrodt et al. (2008, p.
#'   312) specify, or with Satterthwaite (1946) degrees of freedom. The
#'   Kenward-Roger adjustment depends on how the covariance is
#'   parameterized; it is computed with the matrix's own elements as the
#'   parameters, the form that reproduces SAS PROC MIXED (mmrm's
#'   "Kenward-Roger-Linear"; Sabanes Bove et al., 2026).
#'   Fitzmaurice et al. (2011, p. 101) note that the small-sample properties
#'   of both approximations in longitudinal models "have not been
#'   extensively studied"; the simulation study posted on issue #57 compares
#'   them in Solomon designs.
#'
#' The model is fitted with the mmrm package (Sabanes Bove et al., 2026),
#' which must be installed.
#'
#' @section Validation:
#' A simulation study under a protocol posted on issue #57 before any run
#' (12 scenarios, 2,000 replications each; see the article "Longitudinal
#' Designs: Validating the Repeated-Measures Analysis") generated monotone
#' dropout that depended on the last posttest. With 30 to 120 participants
#' per group:
#' - **Coverage and Type I error.** With Kenward-Roger degrees of freedom,
#'   95% intervals covered from 0.937 to 0.962 of the time, and Type I
#'   error where the true value was zero ranged from 0.041 to 0.0615.
#'   Satterthwaite degrees of freedom performed alike.
#' - **Bias.** The contrasts were unbiased; the largest bias was 0.029 SD.
#' - **Normal reference distribution.** It gave coverage as low as 0.935
#'   at 30 per group.
#' - **Shared covariance.** A covariance shared by all four groups misstated
#'   the standard errors by up to 12% for the unpretested contrasts and 17%
#'   for the pretested ones.
#' - **Complete-case analyses.** Per-occasion analyses of the participants
#'   still observed were biased by up to 0.09 SD at the last occasion.
#'
#' @section Lifecycle:
#' Experimental. The study's pre-specified rule for choosing the default
#' degrees of freedom, every tolerance met in 90% of cells and in every cell
#' with 60 or more per group, was not met by either approximation (145 of
#' 156 cells each; 97 and 96 of 104), although Monte Carlo error alone makes
#' the second part unlikely to be met even by a correctly calibrated method.
#' Kenward-Roger remains the default, as Mallinckrodt et al. (2008)
#' specify.
#'
#' @param y_post Numeric posttest scores, one per participant and occasion
#'   (long format), with `NA` for missing posttests.
#' @param treat Treatment indicator coded 0/1 (or logical), the same in every
#'   row of a participant. Designs with several treatments are not supported;
#'   see [fit_solomon_glm()].
#' @param pretested Pretest indicator coded 0/1 (or logical), the same in
#'   every row of a participant.
#' @param id Participant identifier.
#' @param occasion Posttest occasion. A factor keeps its level order;
#'   otherwise the sorted unique values are used.
#' @param y_pre Optional pretest score, the same in every row of a
#'   participant and missing by design for unpretested participants.
#' @param df Degrees of freedom for tests and intervals: `"kenward-roger"`
#'   (default) or `"satterthwaite"`.
#' @param conf_level Confidence level. Default 0.95.
#' @param data Optional data frame. When supplied, the other data arguments
#'   are looked up in it first, as bare column names (`y_post = post`) or as
#'   strings (`y_post = "post"`).
#'
#' @return An object of class `solomon_mmrm` with `effects` (the four
#'   Solomon contrasts and the three pretest effects at each occasion, then
#'   the change in sensitization, in the columns `occasion`, `contrast`,
#'   `estimate`, `std.error`, `statistic`, `df`, `p.value`, `conf.low`, and
#'   `conf.high`), `model` (the mmrm fit), `covariance` (the
#'   structure used), `observed` (observed posttests by group and occasion),
#'   `pretest_mean` (the center of the pretest; `NA` without `y_pre`), and the
#'   settings used.
#'
#' @references
#' Entwisle, D. R. (1961). Interactive effects of pretesting. *Educational and
#' Psychological Measurement, 21*(3), 607–620.
#' https://doi.org/10.1177/001316446102100307
#'
#' Fitzmaurice, G. M., Laird, N. M., & Ware, J. H. (2011). *Applied
#' longitudinal analysis* (2nd ed.). Wiley. https://doi.org/10.1002/9781119513469
#'
#' Laird, N. M., & Ware, J. H. (1982). Random-effects models for longitudinal
#' data. *Biometrics, 38*(4), 963–974. https://doi.org/10.2307/2529876
#'
#' Liang, K.-Y., & Zeger, S. L. (1986). Longitudinal data analysis using
#' generalized linear models. *Biometrika, 73*(1), 13–22.
#' https://doi.org/10.1093/biomet/73.1.13
#'
#' Mallinckrodt, C. H., Lane, P. W., Schnell, D., Peng, Y., & Mancuso, J. P.
#' (2008). Recommendations for the primary analysis of continuous endpoints
#' in longitudinal clinical trials. *Drug Information Journal, 42*(4),
#' 303–319. https://doi.org/10.1177/009286150804200402
#'
#' Sabanes Bove, D., Li, L., Dedic, J., Kelkhoff, D., Kunzmann, K., Lang, B.
#' M., Stock, C., Wang, Y., James, D., Sidi, J., Leibovitz, D., Sjoberg, D.
#' D., Krieger, N. I., Panagos, A., & Jones, J. (2026). *mmrm: Mixed models
#' for repeated measures* (Version 0.3.18) \[R package\].
#' https://doi.org/10.32614/CRAN.package.mmrm
#'
#' Satterthwaite, F. E. (1946). An approximate distribution of estimates of
#' variance components. *Biometrics Bulletin, 2*(6), 110–114.
#' https://doi.org/10.2307/3002019
#'
#' Welch, B. L. (1947). The generalization of "Student's" problem when several
#' different population variances are involved. *Biometrika, 34*(1–2), 28–35.
#' https://doi.org/10.1093/biomet/34.1-2.28
#'
#' @seealso [fit_solomon_glm()] for a single posttest occasion; [jordaan2014]
#'   and the article "Worked Example: Repeated Posttests" for a published
#'   study with three posttest occasions.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("mmrm", quietly = TRUE)) {
#'   # Three posttest occasions for the example participants, with dropout.
#'   set.seed(57)
#'   d <- solomon_example
#'   d$id <- seq_len(nrow(d))
#'   long <- do.call(rbind, lapply(1:3, function(t) {
#'     x <- d
#'     x$occasion <- t
#'     x$y_post <- d$y_post - 2 * (t - 1) + rnorm(nrow(d), 0, 4)
#'     x
#'   }))
#'   long$y_post[long$occasion == 3 & runif(nrow(long)) < 0.3] <- NA
#'   fit_solomon_mmrm(y_post, treat, pretested, id, occasion, y_pre, data = long)
#' }
#' }
#'
#' @export
fit_solomon_mmrm <- function(y_post, treat, pretested, id, occasion, y_pre = NULL,
                             df = c("kenward-roger", "satterthwaite"), conf_level = 0.95,
                             data = NULL) {
  .solomon_data_args(
    data, c("y_post", "treat", "pretested", "id", "occasion", "y_pre"),
    environment(), parent.frame()
  )
  .stop_ngroup_unsupported(treat, "fit_solomon_mmrm")
  if (!requireNamespace("mmrm", quietly = TRUE)) {
    stop("fit_solomon_mmrm() needs the mmrm package: install.packages(\"mmrm\").", call. = FALSE)
  }
  df <- match.arg(df)
  .check_conf_level(conf_level)
  treat <- .solomon_indicator(treat, "treat")
  pretested <- .solomon_indicator(pretested, "pretested")
  .solomon_check_lengths(y_post = y_post, treat = treat, pretested = pretested, id = id,
                         occasion = occasion, y_pre = y_pre)
  if (!is.numeric(y_post)) stop("`y_post` must be numeric.", call. = FALSE)
  if (anyNA(id) || anyNA(occasion)) stop("`id` and `occasion` must not be missing.", call. = FALSE)
  occ <- if (is.factor(occasion)) droplevels(occasion) else
    factor(occasion, levels = sort(unique(occasion)))
  occasions <- levels(occ)
  if (length(occasions) < 2L) {
    stop("`occasion` must have at least two posttest occasions; for one, use fit_solomon_glm().",
         call. = FALSE)
  }
  if (anyDuplicated(data.frame(id, occ))) {
    stop("Each participant can have only one row per occasion.", call. = FALSE)
  }
  constant <- function(x, name) {
    varies <- tapply(x, id, function(v) length(unique(v[!is.na(v)])) > 1L)
    if (any(varies)) {
      stop("`", name, "` must be the same in every row of a participant.", call. = FALSE)
    }
  }
  constant(treat, "treat")
  constant(pretested, "pretested")

  keep <- !is.na(treat) & !is.na(pretested)
  pre <- !is.null(y_pre)
  excluded <- 0L
  if (pre) {
    if (!is.numeric(y_pre)) stop("`y_pre` must be numeric.", call. = FALSE)
    constant(y_pre, "y_pre")
    if (any(pretested == 0L & !is.na(y_pre), na.rm = TRUE)) {
      warning("Observed `y_pre` values were supplied for unpretested participants; ",
              "these values are ignored.", call. = FALSE)
    }
    y_pre[pretested %in% 0L] <- NA
    no_pre <- keep & pretested == 1L & is.na(y_pre)
    if (any(no_pre)) {
      excluded <- length(unique(id[no_pre]))
      warning(excluded, " pretested participant(s) have missing pretest scores and are ",
              "excluded. Missing pretests are not imputed; see check_solomon_missing().",
              call. = FALSE)
      keep <- keep & !(id %in% id[no_pre])
    }
  }
  d <- data.frame(y = y_post, treat = treat, pretested = pretested, occ = occ,
                  id = factor(id), stringsAsFactors = FALSE)
  if (pre) d$pre_obs <- ifelse(pretested == 1L, y_pre, 0)
  d <- d[keep, , drop = FALSE]
  d$pgrp <- factor(ifelse(d$pretested == 1L, "pretested", "unpretested"))
  observed <- d[!is.na(d$y), , drop = FALSE]

  # Center the pretest at its mean among the pretested participants in the
  # model, each counted once, as in fit_solomon_glm() (issue #104).
  pretest_mean <- NA_real_
  center <- NULL
  if (pre) {
    first <- !duplicated(observed$id) & observed$pretested == 1L
    if (any(first)) pretest_mean <- mean(observed$pre_obs[first])
    if (is.finite(pretest_mean)) {
      observed$pre_obs <- ifelse(observed$pretested == 1L, observed$pre_obs - pretest_mean, 0)
    }
    # The variance of that mean, s^2 / n, which the pretest effects carry.
    n_pre <- sum(first)
    if (n_pre > 1L) {
      center <- list(variance = stats::var(observed$pre_obs[first]) / n_pre, df = n_pre - 1)
    }
  }

  res <- .mmrm_fit_contrasts(observed, occasions, pre, grouped = pre, df = df,
                             conf_level = conf_level, center = center)
  if (res$covariance != "unstructured") {
    warning(structure(
      class = c("solomonR_mmrm_fallback_warning", "warning", "condition"),
      list(message = sprintf(paste0(
        "The unstructured covariance did not converge; the %s structure, the converged ",
        "fallback with the smallest AIC, was used (Mallinckrodt et al., 2008, p. 312)."
      ), res$covariance), call = NULL)
    ))
  }

  groups <- factor(.solomon_group_labels[.solomon_group(observed$treat, observed$pretested)],
                   levels = .solomon_group_labels)
  counts <- as.data.frame.matrix(table(groups, observed$occ))
  participants <- unique(d[, c("id", "treat", "pretested")])
  structure(
    list(
      effects = res$effects,
      model = res$model,
      covariance = res$covariance,
      grouped = pre,
      pretest = pre,
      pretest_mean = pretest_mean,
      df_method = df,
      conf_level = conf_level,
      occasions = occasions,
      observed = counts,
      excluded = excluded,
      data = data.frame(treat = participants$treat, pretested = participants$pretested)
    ),
    class = "solomon_mmrm"
  )
}

#' @export
print.solomon_mmrm <- function(x, digits = 3, ...) {
  cat("Solomon MMRM: mixed model for repeated measures (Mallinckrodt et al., 2008)\n")
  cat("Occasions: ", paste(x$occasions, collapse = ", "), "\n", sep = "")
  cat("Covariance: ", x$covariance,
      if (x$grouped) ", separately for pretested and unpretested participants" else "",
      " (REML)\n", sep = "")
  cat("Degrees of freedom: ",
      if (x$df_method == "kenward-roger") "Kenward-Roger" else "Satterthwaite", "\n", sep = "")
  if (is.finite(.null_na(x$pretest_mean))) {
    cat(sprintf("Pretest centered at the pretested participants' mean: %.*f\n",
                digits, x$pretest_mean))
    cat("Pretest effects: standard errors include the sampling variance of that mean.\n")
  }
  if (x$excluded > 0L) cat("Excluded participants: ", x$excluded, "\n", sep = "")
  cat("\nObserved posttests by group and occasion:\n")
  print(x$observed)
  cat("\n")
  eff <- x$effects
  tab <- data.frame(
    Occasion = eff$occasion,
    Contrast = eff$contrast,
    Estimate = round(eff$estimate, digits),
    SE = round(eff$std.error, digits),
    df = round(eff$df, 1),
    t = round(eff$statistic, 2),
    p = ifelse(eff$p.value < .001, "<.001", sprintf("%.3f", eff$p.value)),
    CI = sprintf("[%s, %s]", format(round(eff$conf.low, digits), nsmall = digits),
                 format(round(eff$conf.high, digits), nsmall = digits)),
    check.names = FALSE
  )
  names(tab)[8] <- sprintf("%s%% CI", format(100 * x$conf_level))
  print(tab, row.names = FALSE, right = FALSE)
  invisible(x)
}
