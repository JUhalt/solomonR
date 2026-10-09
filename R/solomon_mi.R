# Multiple imputation of missing posttests, with delta-adjusted sensitivity
# analysis and a tipping-point summary (issue #82).

.solomon_group_labels <- c("pretested treatment", "pretested control",
                           "unpretested treatment", "unpretested control")

# Solomon group of each participant: 1 pretested treatment, 2 pretested
# control, 3 unpretested treatment, 4 unpretested control.
.solomon_group <- function(treat, pretested) {
  1L + (1L - treat) + 2L * (1L - pretested)
}

.mi_delta <- function(delta) {
  if (!is.numeric(delta) || !length(delta) %in% c(1L, 4L) || any(!is.finite(delta))) {
    stop("`delta` must be one number, or four numbers for Groups 1-4 (pretested ",
         "treatment, pretested control, unpretested treatment, unpretested control).",
         call. = FALSE)
  }
  stats::setNames(rep_len(as.numeric(delta), 4L), .solomon_group_labels)
}

# Participants analyzed, their groups, and those excluded.
.mi_prepare <- function(y_post, treat, pretested, y_pre) {
  treat <- .solomon_indicator(treat, "treat")
  pretested <- .solomon_indicator(pretested, "pretested")
  .solomon_check_lengths(y_post = y_post, treat = treat, pretested = pretested, y_pre = y_pre)
  if (!is.numeric(y_post)) {
    stop("`y_post` must be numeric: fit_solomon_mi() imputes continuous posttests.",
         call. = FALSE)
  }
  keep <- !is.na(treat) & !is.na(pretested)
  if (!is.null(y_pre)) {
    if (!is.numeric(y_pre)) stop("`y_pre` must be numeric.", call. = FALSE)
    if (any(pretested == 0L & !is.na(y_pre), na.rm = TRUE)) {
      warning("Observed `y_pre` values were supplied for unpretested participants; ",
              "these values are ignored.", call. = FALSE)
    }
    y_pre[pretested %in% 0L] <- NA
    no_pre <- keep & pretested == 1L & is.na(y_pre)
    if (any(no_pre)) {
      warning(sum(no_pre), " pretested participant(s) have missing pretest scores and ",
              "are excluded. Missing pretests are not imputed; see check_solomon_missing().",
              call. = FALSE)
      keep <- keep & !no_pre
    }
  }
  list(
    y_post = y_post[keep], treat = treat[keep], pretested = pretested[keep],
    y_pre = if (!is.null(y_pre)) y_pre[keep],
    group = .solomon_group(treat[keep], pretested[keep]),
    excluded = sum(!keep)
  )
}

# The imputation model of each group, fitted to its observed posttests: the
# posttest on the pretest in the pretested groups, and on a constant in the
# unpretested groups (Carpenter et al., 2023, section 3.1).
.mi_models <- function(d) {
  lapply(1:4, function(j) {
    rows <- which(d$group == j)
    obs <- !is.na(d$y_post[rows])
    if (all(obs)) return(NULL)
    use_pre <- !is.null(d$y_pre) && j <= 2L
    X <- if (use_pre) cbind(1, d$y_pre[rows]) else matrix(1, length(rows), 1L)
    p <- ncol(X)
    if (sum(obs) < p + 1L) {
      stop(sprintf(
        "Group %d (%s) has %d observed posttest(s); at least %d are needed to impute its missing posttests.",
        j, .solomon_group_labels[j], sum(obs), p + 1L
      ), call. = FALSE)
    }
    q <- qr(X[obs, , drop = FALSE])
    if (q$rank < p) {
      stop(sprintf("The pretest is constant among the observed posttests of Group %d (%s).",
                   j, .solomon_group_labels[j]), call. = FALSE)
    }
    yo <- d$y_post[rows][obs]
    beta <- qr.coef(q, yo)
    df <- sum(obs) - p
    list(
      rows_mis = rows[!obs],
      X_mis = X[!obs, , drop = FALSE],
      beta = beta,
      s2 = sum(qr.resid(q, yo)^2) / df,
      df = df,
      R = chol(chol2inv(qr.R(q)))
    )
  })
}

# One proper draw of the missing posttests of a group: the variance and the
# coefficients from their posterior, then the posttests (Carpenter et al.,
# 2023, pp. 81-83).
.mi_draw <- function(model) {
  s2 <- model$s2 * model$df / stats::rchisq(1L, model$df)
  beta <- model$beta + sqrt(s2) * drop(crossprod(model$R, stats::rnorm(length(model$beta))))
  drop(model$X_mis %*% beta) + stats::rnorm(nrow(model$X_mis), 0, sqrt(s2))
}

# The completed-data analysis of fit_solomon_glm() (linear model with the
# pretest adjustment, HC3 or model-based covariance), computed directly.
# Every completed data set has the same design matrix, so the contrast
# weights A = L (X'X)^-1 X' and the leverages are computed once: the
# estimates are A y and the HC3 variances are (A * A) w, with
# w = (e / (1 - h))^2 (MacKinnon & White, 1985). test-mi.R checks that the
# results equal fit_solomon_glm()'s.
.mi_design <- function(d, robust) {
  pre <- !is.null(d$y_pre)
  X <- cbind(1, d$treat, d$pretested,
             if (pre) ifelse(d$pretested == 1L, d$y_pre, 0),
             d$treat * d$pretested)
  L <- rbind(
    c(0, 1, 0, if (pre) 0, 0.5),  # ATE (avg over pretest)
    c(0, 0, 0, if (pre) 0, 1),    # Pretest x Treatment
    c(0, 1, 0, if (pre) 0, 1),    # Treatment | pretested
    c(0, 1, 0, if (pre) 0, 0)     # Treatment | unpretested
  )
  xtx_inv <- solve(crossprod(X))
  B <- xtx_inv %*% t(X)
  list(
    X = X, B = B, A = L %*% B,
    h = rowSums((X %*% xtx_inv) * X),
    lxl = rowSums((L %*% xtx_inv) * L),
    df = nrow(X) - ncol(X),
    robust = robust
  )
}

.mi_analyze <- function(des, y) {
  e <- y - drop(des$X %*% (des$B %*% y))
  variance <- if (des$robust == "HC3") {
    drop((des$A^2) %*% (e / (1 - des$h))^2)
  } else {
    des$lxl * sum(e^2) / des$df
  }
  list(estimate = drop(des$A %*% y), std.error = sqrt(variance))
}

# Rubin's rules (Carpenter et al., 2023, Eqs. 2.16 and 2.26) with the
# small-sample degrees of freedom of Barnard and Rubin (1999), as given by
# van Buuren (2018, Eqs. 2.30-2.32).
.rubin_pool <- function(est, se, df_com, conf_level, any_missing = TRUE) {
  m <- nrow(est)
  qbar <- colMeans(est)
  ubar <- colMeans(se^2)
  b <- if (any_missing) apply(est, 2, stats::var) else rep(0, ncol(est))
  total <- ubar + (1 + 1 / m) * b
  lambda <- (1 + 1 / m) * b / total
  riv <- (1 + 1 / m) * b / ubar
  df <- if (!any_missing) {
    df_com
  } else {
    df_old <- (m - 1) / pmax(lambda, 1e-10)^2
    df_obs <- (df_com + 1) / (df_com + 3) * df_com * (1 - lambda)
    ifelse(is.finite(df_com), df_old * df_obs / (df_old + df_obs), df_old)
  }
  se_total <- sqrt(total)
  crit <- stats::qt(1 - (1 - conf_level) / 2, df)
  statistic <- qbar / se_total
  data.frame(
    contrast = colnames(est),
    estimate = unname(qbar),
    std.error = unname(se_total),
    statistic = unname(statistic),
    df = unname(df),
    p.value = unname(2 * stats::pt(-abs(statistic), df)),
    conf.low = unname(qbar - crit * se_total),
    conf.high = unname(qbar + crit * se_total),
    fmi = if (any_missing) unname((riv + 2 / (df + 3)) / (riv + 1)) else 0,
    mc_se = unname(sqrt(b / m)),
    stringsAsFactors = FALSE
  )
}

#' Solomon analysis with multiply imputed posttests
#'
#' `r lifecycle::badge("experimental")`
#' Multiply imputes missing posttests, analyzes each completed data set with
#' [fit_solomon_glm()], and combines the four Solomon contrasts with Rubin's
#' rules. With `delta = 0` the imputations assume the posttests are missing
#' at random (MAR). A nonzero `delta` shifts the imputed posttests of each
#' group by a fixed amount, the delta-adjusted pattern-mixture sensitivity
#' analysis of Carpenter et al. (2023, section 10.3). [tipping_point_solomon()]
#' repeats it over a range of offsets.
#'
#' @section Method:
#' 1. **Imputation model.** A normal linear regression is fitted separately
#'    in each Solomon group to the participants with an observed posttest:
#'    on the pretest in the pretested groups, and on a constant in the
#'    unpretested groups, which have no pretest by design. The analysis
#'    contains the Pretest x Treatment interaction of two fully observed
#'    indicators, and imputing separately in the groups they define is the
#'    simplest approach to such interactions (Carpenter et al., 2023,
#'    section 6.3.5, p. 149).
#' 2. **Proper imputation.** For each imputation the residual variance and
#'    the coefficients are drawn from their posterior distribution, and the
#'    missing posttests are then drawn from the model (Carpenter et al.,
#'    2023, pp. 81-83).
#' 3. **Offsets.** Each imputed posttest in group j is shifted by
#'    `delta[j]`: those with a missing posttest are assumed to differ from
#'    those observed by an offset in each group, as in the pattern-mixture
#'    analysis of Little et al. (2012, p. 1358), and "a clinically plausible
#'    amount" is added to the imputed outcomes (White et al., 2011, "Perform
#'    Sensitivity Analyses" section, para. 1).
#' 4. **Analysis and pooling.** Each completed data set is analyzed with the
#'    model of [fit_solomon_glm()] (computed directly, since the design is
#'    the same in every completed data set), and the contrasts are combined
#'    with Rubin's rules:
#'    the mean estimate, with variance W + (1 + 1/m)B (Carpenter et al.,
#'    2023, Eq. 2.16). Tests and intervals use t with the small-sample degrees
#'    of freedom of Barnard and Rubin (1999, as cited in van Buuren, 2018,
#'    Eqs. 2.30-2.32), taking the complete-data residual degrees of freedom
#'    as their starting point.
#'
#' With fixed offsets, multiple imputation with Rubin's variance is
#' information-anchored: the sensitivity analysis neither adds nor removes
#' information relative to the MAR analysis (Cro et al., 2019; Carpenter et
#' al., 2023, p. 278).
#'
#' **Structural and incidental missingness.** Pretests absent by design are
#' never imputed. Pretested participants with a missing pretest are excluded
#' with a warning, as in [fit_solomon_glm()]; [check_solomon_missing()]
#' describes the options for them.
#'
#' **Number of imputations.** The default `m = 100` follows Carpenter et al.
#' (2023, p. 56), who note that p-values accurate to about .005 need at
#' least 100 imputations and advise erring "towards too many imputations
#' rather than too few". With much missing information, more are needed: in
#' the worked example on `mai2020`, where the fraction of missing information
#' is about 0.4, p-values varied by about .01 from one seed to another with
#' 100 imputations. The `mc_se` column gives the Monte Carlo standard error
#' of each estimate due to the finite number of imputations.
#'
#' @section Validation:
#' A simulation study under a protocol posted on issue #82 before any run
#' (24 scenarios, 2,000 replications each; see the article "Missing
#' Posttests: Validating the Sensitivity Analysis") found:
#' - **With 60 or more participants per group,** coverage of 95% intervals
#'   from 0.9415 to 0.9595 and Type I error within 0.040 to 0.060, under
#'   missing at random and under the pattern-mixture departures studied
#'   when the offsets were right. Bias exceeded 2 Monte Carlo standard
#'   errors in 2 of 64 contrasts, by at most 0.011 SD.
#' - **With 30 per group,** conservative intervals: coverage up to 0.967
#'   and model standard errors 3.9% above the empirical ones on average.
#' - **Missing at random:** agreement with the complete-case analysis, as
#'   theory predicts (Carpenter et al., 2023, p. 256).
#' - **Departures that differ between the pretested groups:** the
#'   analyses that assume missing at random biased the sensitization
#'   contrast by 0.10 to 0.16 SD.
#'
#' @section Lifecycle:
#' Experimental. The study's pre-specified rule for validation, every
#' tolerance met in 90% of cells and in every cell with 60 or more per group,
#' was not met: 84 of 96 cells and 62 of 64.
#'
#' @param y_post Numeric posttest scores, with `NA` for missing posttests.
#' @param treat Treatment indicator coded 0/1 (or logical). Designs with
#'   several treatments are not supported; see [fit_solomon_glm()].
#' @param pretested Pretest indicator coded 0/1 (or logical).
#' @param y_pre Optional numeric pretest scores, missing by design for
#'   unpretested participants.
#' @param delta Offsets added to the imputed posttests, on the posttest
#'   scale: one number for all four groups, or four numbers for Groups 1-4
#'   (pretested treatment, pretested control, unpretested treatment,
#'   unpretested control). The default, 0, is the MAR analysis.
#' @param m Number of imputations. Default 100.
#' @param robust Covariance for each completed-data analysis: `"HC3"`
#'   (default) or `"none"`; see [fit_solomon_glm()].
#' @param conf_level Confidence level. Default 0.95.
#' @param seed Optional random-number seed. The global random number state is
#'   restored afterwards.
#' @param data Optional data frame. When supplied, the other data arguments
#'   are looked up in it first, as bare column names (`y_post = post`) or as
#'   strings (`y_post = "post"`).
#'
#' @return An object of class `solomon_mi`, a list with:
#'   - `effects`: the pooled treatment contrasts, `ATE (avg over pretest)`,
#'     `Pretest x Treatment`, `Treatment | pretested`, and
#'     `Treatment | unpretested`, in the columns `contrast`, `estimate`,
#'     `std.error`, `statistic` (t), `df` (Barnard and Rubin's), `p.value`,
#'     `conf.low`, and `conf.high`, followed by the fraction of missing
#'     information `fmi` and the Monte Carlo standard error `mc_se`.
#'   - `conf_level`: the confidence level of the intervals.
#'   - `delta`: the offset added in each of the four groups.
#'   - `m`: the number of imputations.
#'   - `missing`: the posttests missing in each group (`group`, `n`,
#'     `missing`, and `proportion`).
#'   - `excluded`: the number of pretested participants left out because
#'     their pretest was missing.
#'   - `estimates` and `std_errors`: the estimates and standard errors of
#'     each imputation, one row for each imputation and one column for each
#'     contrast.
#'   - `df_com`: the complete-data degrees of freedom of each contrast.
#'   - `robust`: the covariance of each completed-data analysis.
#'   - `pretest` (whether pretest scores were supplied) and `data` (the
#'     design indicators), which solomonR's own functions use.
#'
#'   The `effects` table, `conf_level`, and [`tidy()`][solomon_output],
#'   which returns the table, are the stable interface of the result; see
#'   [solomon_output].
#'
#' @references
#' Carpenter, J. R., Bartlett, J. W., Morris, T. P., Wood, A. M., Quartagno,
#' M., & Kenward, M. G. (2023). *Multiple imputation and its application*
#' (2nd ed.). Wiley. https://doi.org/10.1002/9781119756118
#'
#' Cro, S., Carpenter, J. R., & Kenward, M. G. (2019). Information-anchored
#' sensitivity analysis: Theory and application. *Journal of the Royal
#' Statistical Society Series A: Statistics in Society, 182*(2), 623–645.
#' https://doi.org/10.1111/rssa.12423
#'
#' Little, R. J., D'Agostino, R., Cohen, M. L., Dickersin, K., Emerson, S.
#' S., Farrar, J. T., Frangakis, C., Hogan, J. W., Molenberghs, G., Murphy,
#' S. A., Neaton, J. D., Rotnitzky, A., Scharfstein, D., Shih, W. J.,
#' Siegel, J. P., & Stern, H. (2012). The prevention and treatment of
#' missing data in clinical trials. *The New England Journal of Medicine,
#' 367*(14), 1355–1360. https://doi.org/10.1056/NEJMsr1203730
#'
#' van Buuren, S. (2018). *Flexible imputation of missing data* (2nd ed.).
#' CRC Press. https://doi.org/10.1201/9780429492259
#'
#' White, I. R., Horton, N. J., Carpenter, J., & Pocock, S. J. (2011).
#' Strategy for intention to treat analysis in randomised trials with missing
#' outcome data. *BMJ, 342*, Article d40. https://doi.org/10.1136/bmj.d40
#'
#' @seealso [tipping_point_solomon()], [check_solomon_missing()]
#'
#' @examples
#' d <- solomon_example
#' set.seed(82)
#' d$y_post[sample(nrow(d), 20)] <- NA
#'
#' # Missing at random.
#' fit_solomon_mi(y_post, treat, pretested, y_pre, m = 20, seed = 1, data = d)
#'
#' # Missing posttests in the treatment groups 3 points lower than MAR predicts.
#' fit_solomon_mi(y_post, treat, pretested, y_pre, delta = c(-3, 0, -3, 0),
#'                m = 20, seed = 1, data = d)
#'
#' @export
fit_solomon_mi <- function(y_post, treat, pretested, y_pre = NULL, delta = 0, m = 100,
                           robust = c("HC3", "none"), conf_level = 0.95, seed = NULL,
                           data = NULL) {
  .solomon_data_args(
    data, c("y_post", "treat", "pretested", "y_pre"),
    environment(), parent.frame()
  )
  .stop_ngroup_unsupported(treat, "fit_solomon_mi")
  robust <- match.arg(robust)
  .check_conf_level(conf_level)
  delta <- .mi_delta(delta)
  m <- as.integer(m)
  if (length(m) != 1L || is.na(m) || m < 2L) {
    stop("`m` must be a whole number of at least 2.", call. = FALSE)
  }
  d <- .mi_prepare(y_post, treat, pretested, y_pre)
  if (!is.null(seed)) withr::local_seed(seed)

  models <- .mi_models(d)
  any_missing <- any(is.na(d$y_post))
  contrasts <- c("ATE (avg over pretest)", "Pretest x Treatment",
                 "Treatment | pretested", "Treatment | unpretested")
  est <- se <- matrix(NA_real_, m, 4L, dimnames = list(NULL, contrasts))
  des <- .mi_design(d, robust)
  df_com <- rep(des$df, 4L)
  for (k in seq_len(m)) {
    y_k <- d$y_post
    for (j in 1:4) {
      if (!is.null(models[[j]])) {
        y_k[models[[j]]$rows_mis] <- .mi_draw(models[[j]]) + delta[[j]]
      }
    }
    fit_k <- .mi_analyze(des, y_k)
    est[k, ] <- fit_k$estimate
    se[k, ] <- fit_k$std.error
  }

  structure(
    list(
      effects = .effects_table(.rubin_pool(est, se, df_com, conf_level, any_missing)),
      delta = delta,
      m = m,
      missing = .mi_missing_table(d),
      excluded = d$excluded,
      pretest = !is.null(d$y_pre),
      estimates = est,
      std_errors = se,
      df_com = stats::setNames(df_com, contrasts),
      robust = robust,
      conf_level = conf_level,
      data = data.frame(treat = d$treat, pretested = d$pretested)
    ),
    class = "solomon_mi"
  )
}

#' @export
print.solomon_mi <- function(x, digits = 3, ...) {
  cat("Solomon analysis with multiply imputed posttests (m = ", x$m, ")\n", sep = "")
  cat("Imputation: normal linear model in each group",
      if (x$pretest) ", on the pretest in the pretested groups" else "", "\n", sep = "")
  if (all(x$delta == 0)) {
    cat("Assumption: missing at random (delta = 0)\n")
  } else {
    cat("Offsets added to imputed posttests:\n")
    cat(sprintf("  %-22s %s\n", paste0(names(x$delta), ":"), format(x$delta, digits = digits)),
        sep = "")
  }
  miss <- x$missing
  cat("Missing posttests: ",
      paste(sprintf("%d of %d (%s)", miss$missing, miss$n, miss$group), collapse = "; "),
      "\n", sep = "")
  if (x$excluded > 0L) cat("Excluded participants: ", x$excluded, "\n", sep = "")
  cat("Pooling: Rubin's rules; Barnard-Rubin degrees of freedom; ",
      if (x$robust == "HC3") "HC3 standard errors" else "model-based standard errors",
      "\n\n", sep = "")
  eff <- x$effects
  tab <- data.frame(
    Contrast = eff$contrast,
    Estimate = round(eff$estimate, digits),
    SE = round(eff$std.error, digits),
    df = round(eff$df, 1),
    t = round(eff$statistic, 2),
    p = ifelse(eff$p.value < .001, "<.001", sprintf("%.3f", eff$p.value)),
    CI = sprintf("[%s, %s]", format(round(eff$conf.low, digits), nsmall = digits),
                 format(round(eff$conf.high, digits), nsmall = digits)),
    FMI = round(eff$fmi, 2),
    check.names = FALSE
  )
  names(tab)[7] <- sprintf("%s%% CI", format(100 * x$conf_level))
  print(tab, row.names = FALSE, right = FALSE)
  cat(sprintf("\nLargest Monte Carlo SE from the finite m: %s\n",
              format(max(eff$mc_se), digits = 2)))
  invisible(x)
}


# ---- Tipping point ----------------------------------------------------------------

.tipping_groups <- list(
  all = 1:4, treatment = c(1L, 3L), control = c(2L, 4L),
  pretested = 1:2, unpretested = 3:4
)

#' Tipping-point analysis for missing posttests
#'
#' `r lifecycle::badge("experimental")`
#' Repeats [fit_solomon_mi()] over a range of offsets added to the imputed
#' posttests of chosen Solomon groups, and reports the smallest offset in
#' each direction at which the conclusion about a contrast changes. White et
#' al. (2011, "Perform Sensitivity Analyses" section, para. 1) suggest
#' reporting "how large an amount should be added to or
#' subtracted from imputed outcomes" without changing the interpretation, and
#' Little et al. (2012, p. 1358) call a finding robust if it holds over the
#' plausible offsets.
#'
#' **Which groups.** Offsets confined to some groups test different
#' departures from missing at random. Offsets in the treatment groups
#' (`groups = "treatment"`) bear on the treatment effect. Offsets that differ
#' between the pretested and unpretested groups, for example in the
#' pretested treatment group alone (`groups = 1`), bear on the sensitization
#' contrast.
#'
#' **Common random numbers.** Every offset uses the same imputation draws,
#' so the estimates change smoothly with the offset. The location of the
#' tipping point still carries Monte Carlo error from the imputations: in
#' the worked example on `mai2020`, the tipping point for sensitization
#' ranged from 0.1 to 0.5 standard deviations across seeds with 100
#' imputations and was 0.3 with 2,000. Before reporting a tipping point,
#' increase `m` or compare a few seeds.
#'
#' @section Lifecycle:
#' Experimental, with [fit_solomon_mi()], whose validation study it shares.
#'
#' @inheritParams fit_solomon_mi
#' @param contrast The Solomon contrast to follow. Default is the average
#'   treatment effect.
#' @param groups The groups whose imputed posttests are shifted: `"all"`
#'   (default), `"treatment"`, `"control"`, `"pretested"`, `"unpretested"`,
#'   or group numbers (1 pretested treatment, 2 pretested control, 3
#'   unpretested treatment, 4 unpretested control).
#' @param deltas Offsets to try, on the posttest scale. The default is 21
#'   values from -1 to 1 pooled within-group standard deviations of the
#'   observed posttests. Zero is always included.
#' @param alpha Significance level that defines the conclusion. Default 0.05.
#' @param m Number of imputations for each offset. Default 100.
#' @param seed Optional random-number seed. When `NULL`, one is drawn so that
#'   every offset uses the same imputations.
#'
#' @return An object of class `solomon_tipping`, a list with:
#'   - `results`: one row for each offset, in the columns `delta` (the
#'     offset in posttest units), `delta_sd` (the offset in standard
#'     deviations), and then the pooled result in the columns of
#'     `fit_solomon_mi()$effects`: `contrast`, `estimate`, `std.error`,
#'     `statistic`, `df`, `p.value`, `conf.low`, and `conf.high`, followed
#'     by `significant` (whether `p.value` is below `alpha`).
#'   - `tipping`: the smallest negative and positive offsets at which the
#'     conclusion differs from the one under MAR, `NA` when it does not
#'     change within the range; and `tipping_sd`, the same in standard
#'     deviations.
#'   - `conf_level`: the confidence level of the intervals, 1 - `alpha`.
#'   - `contrast`, `groups` (the groups shifted), `alpha`, `m`, `seed`, and
#'     `robust`: the settings used.
#'   - `sd`: the pooled within-group standard deviation of the observed
#'     posttests.
#'   - `missing`: the posttests missing in each group, as in
#'     [fit_solomon_mi()].
#'   - `pretest` (whether pretest scores were supplied) and `data` (the
#'     design indicators), which solomonR's own functions use.
#'
#'   `results` has the columns and the labels of an effects table, and
#'   [`tidy()`][solomon_output] returns it. See [solomon_output] for the
#'   columns, the labels, and the parts of a result that are stable.
#'
#' @references
#' Little, R. J., D'Agostino, R., Cohen, M. L., Dickersin, K., Emerson, S.
#' S., Farrar, J. T., Frangakis, C., Hogan, J. W., Molenberghs, G., Murphy,
#' S. A., Neaton, J. D., Rotnitzky, A., Scharfstein, D., Shih, W. J.,
#' Siegel, J. P., & Stern, H. (2012). The prevention and treatment of
#' missing data in clinical trials. *The New England Journal of Medicine,
#' 367*(14), 1355–1360. https://doi.org/10.1056/NEJMsr1203730
#'
#' White, I. R., Horton, N. J., Carpenter, J., & Pocock, S. J. (2011).
#' Strategy for intention to treat analysis in randomised trials with missing
#' outcome data. *BMJ, 342*, Article d40. https://doi.org/10.1136/bmj.d40
#'
#' @seealso [fit_solomon_mi()], [plot_tipping_point()]
#'
#' @examples
#' d <- solomon_example
#' set.seed(82)
#' d$y_post[sample(nrow(d), 20)] <- NA
#' tipping_point_solomon(y_post, treat, pretested, y_pre, groups = "treatment",
#'                       deltas = seq(-10, 10, by = 2.5), m = 10, seed = 1, data = d)
#'
#' @export
tipping_point_solomon <- function(y_post, treat, pretested, y_pre = NULL,
                                  contrast = "ATE (avg over pretest)", groups = "all",
                                  deltas = NULL, alpha = 0.05, m = 100,
                                  robust = c("HC3", "none"), seed = NULL, data = NULL) {
  .solomon_data_args(
    data, c("y_post", "treat", "pretested", "y_pre"),
    environment(), parent.frame()
  )
  .stop_ngroup_unsupported(treat, "tipping_point_solomon")
  robust <- match.arg(robust)
  valid <- c("ATE (avg over pretest)", "Pretest x Treatment",
             "Treatment | pretested", "Treatment | unpretested")
  if (length(contrast) != 1L || !contrast %in% valid) {
    stop("`contrast` must be one of: ", paste(sprintf("\"%s\"", valid), collapse = ", "), ".",
         call. = FALSE)
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 0.5) {
    stop("`alpha` must be a single number between 0 and 0.5.", call. = FALSE)
  }
  shifted <- if (is.character(groups)) {
    if (length(groups) != 1L || !groups %in% names(.tipping_groups)) {
      stop("`groups` must be one of ", paste(sprintf("\"%s\"", names(.tipping_groups)), collapse = ", "),
           ", or group numbers from 1 to 4.", call. = FALSE)
    }
    .tipping_groups[[groups]]
  } else {
    if (!is.numeric(groups) || !length(groups) || any(!groups %in% 1:4) || anyDuplicated(groups)) {
      stop("Group numbers in `groups` must be distinct values from 1 to 4.", call. = FALSE)
    }
    sort(as.integer(groups))
  }
  pattern <- as.numeric(1:4 %in% shifted)

  # The scale for reporting offsets: the pooled within-group SD of the
  # observed posttests.
  d <- .mi_prepare(y_post, treat, pretested, y_pre)
  obs <- !is.na(d$y_post)
  resid <- d$y_post[obs] - stats::ave(d$y_post[obs], d$group[obs])
  sd_pooled <- sqrt(sum(resid^2) / (sum(obs) - length(unique(d$group[obs]))))

  if (is.null(deltas)) deltas <- seq(-1, 1, by = 0.1) * sd_pooled
  if (!is.numeric(deltas) || !length(deltas) || any(!is.finite(deltas))) {
    stop("`deltas` must be finite numbers.", call. = FALSE)
  }
  deltas <- sort(unique(c(0, deltas)))
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1L)

  rows <- lapply(deltas, function(delta) {
    fit <- fit_solomon_mi(d$y_post, d$treat, d$pretested, d$y_pre, delta = delta * pattern,
                          m = m, robust = robust, conf_level = 1 - alpha, seed = seed)
    fit$effects[fit$effects$contrast == contrast, ]
  })
  res <- do.call(rbind, rows)
  # The columns of an effects table (issue #110), indexed by the offset.
  res <- data.frame(delta = deltas, delta_sd = deltas / sd_pooled,
                    res[, .solomon_effect_columns],
                    significant = res$p.value < alpha, row.names = NULL,
                    stringsAsFactors = FALSE)
  res <- .effects_table(res, keys = c("delta", "delta_sd"))

  base <- res$significant[res$delta == 0]
  first_change <- function(side) {
    r <- res[side, , drop = FALSE]
    r <- r[order(abs(r$delta)), , drop = FALSE]
    hit <- which(r$significant != base)
    if (length(hit)) r$delta[hit[1L]] else NA_real_
  }
  tipping <- c(negative = first_change(res$delta < 0), positive = first_change(res$delta > 0))

  structure(
    list(
      results = res,
      tipping = tipping,
      tipping_sd = tipping / sd_pooled,
      contrast = contrast,
      groups = .solomon_group_labels[shifted],
      alpha = alpha,
      # The level of the intervals in `results`, under the name every
      # result with a table of contrasts uses (issue #110).
      conf_level = 1 - alpha,
      m = m,
      seed = seed,
      sd = sd_pooled,
      robust = robust,
      pretest = !is.null(d$y_pre),
      data = data.frame(treat = d$treat, pretested = d$pretested),
      missing = .mi_missing_table(d)
    ),
    class = "solomon_tipping"
  )
}

.mi_missing_table <- function(d) {
  n <- tabulate(d$group, 4L)
  n_missing <- tabulate(d$group[is.na(d$y_post)], 4L)
  data.frame(group = .solomon_group_labels, n = n, missing = n_missing,
             proportion = n_missing / n, stringsAsFactors = FALSE)
}

#' @export
print.solomon_tipping <- function(x, digits = 3, ...) {
  cat("Tipping-point analysis for missing posttests (m = ", x$m, " per offset)\n", sep = "")
  cat("Contrast: ", x$contrast, "\n", sep = "")
  cat("Offsets added to imputed posttests in: ", paste(x$groups, collapse = ", "), "\n", sep = "")
  cat(sprintf("Offset scale: pooled within-group SD of observed posttests = %s\n\n",
              format(x$sd, digits = digits)))
  r <- x$results
  tab <- data.frame(
    Offset = round(r$delta, digits),
    `Offset (SD)` = round(r$delta_sd, 2),
    Estimate = round(r$estimate, digits),
    CI = sprintf("[%s, %s]", format(round(r$conf.low, digits), nsmall = digits),
                 format(round(r$conf.high, digits), nsmall = digits)),
    p = ifelse(r$p.value < .001, "<.001", sprintf("%.3f", r$p.value)),
    check.names = FALSE
  )
  names(tab)[4] <- sprintf("%s%% CI", format(100 * (1 - x$alpha)))
  print(tab, row.names = FALSE, right = FALSE)
  base <- r[r$delta == 0, ]
  cat(sprintf("\nUnder MAR (offset 0): p %s, %s at alpha = %s.\n",
              if (base$p.value < .001) "< .001" else sprintf("= %.3f", base$p.value),
              if (base$significant) "significant" else "not significant", format(x$alpha)))
  for (side in c("negative", "positive")) {
    tp <- x$tipping[[side]]
    cat(if (is.na(tp)) {
      sprintf("No change in the conclusion over the %s offsets tried.\n", side)
    } else {
      sprintf("The conclusion changes at an offset of %s (%s SD).\n",
              format(round(tp, digits)), format(round(x$tipping_sd[[side]], 2)))
    })
  }
  invisible(x)
}

#' Plot a tipping-point analysis
#'
#' `r lifecycle::badge("experimental")`
#' Plots the pooled estimate of a Solomon contrast, with its confidence
#' interval, against the offset added to the imputed posttests, from
#' [tipping_point_solomon()]. Dashed lines mark the offsets at which the
#' conclusion changes.
#'
#' @param tipping A result from [tipping_point_solomon()].
#' @return A ggplot object.
#' @seealso [tipping_point_solomon()]
#' @examples
#' d <- solomon_example
#' set.seed(82)
#' d$y_post[sample(nrow(d), 20)] <- NA
#' tp <- tipping_point_solomon(y_post, treat, pretested, y_pre, groups = "treatment",
#'                             deltas = seq(-10, 10, by = 2.5), m = 10, seed = 1, data = d)
#' plot_tipping_point(tp)
#' @export
plot_tipping_point <- function(tipping) {
  if (!inherits(tipping, "solomon_tipping")) {
    stop("`tipping` must come from tipping_point_solomon().", call. = FALSE)
  }
  r <- tipping$results
  marks <- tipping$tipping[!is.na(tipping$tipping)]
  p <- ggplot2::ggplot(r, ggplot2::aes(x = delta, y = estimate)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
    ggplot2::geom_line() +
    ggplot2::geom_point(ggplot2::aes(shape = significant)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey40") +
    ggplot2::scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 19),
                                labels = c(`FALSE` = "No", `TRUE` = "Yes"),
                                name = sprintf("p < %s", format(tipping$alpha))) +
    ggplot2::labs(
      x = sprintf("Offset added to imputed posttests (%s)", paste(tipping$groups, collapse = ", ")),
      y = tipping$contrast, title = "Tipping-point analysis",
      subtitle = sprintf("Pooled estimates with %s%% intervals; m = %d per offset",
                         format(100 * (1 - tipping$alpha)), tipping$m)
    ) +
    ggplot2::theme_minimal(base_size = 12)
  if (length(marks)) {
    p <- p + ggplot2::geom_vline(xintercept = marks, linetype = "dashed")
  }
  p
}
