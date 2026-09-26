#' solomonR: Analyze Solomon Four-Group Designs
#'
#' Provides tools for analyzing, teaching, and studying the Solomon
#' four-group experimental design. The package includes the historical
#' Tests A-I workflow, unified generalized linear models with robust
#' covariance estimation, stratified randomization inference,
#' Solomon-specific full-information maximum likelihood, and observed-
#' and latent-variable structural equation models.
#'
#' The package distinguishes historically important procedures from
#' contemporary recommendations and emphasizes explicit Solomon-specific
#' estimands, structural pretest missingness, reproducible diagnostics,
#' and modern inference.
#'
#' Key historical references include Huck & Sandler (1973),
#' Walton Braver & Braver (1988), the Sawilowsky and Markman methodological
#' exchanges, and van Engelenburg (1999).
#'
#' @references
#' Campbell, D. T., & Stanley, J. C. (1963). *Experimental and
#' quasi-experimental designs for research*. Rand McNally.
#'
#' Huck, S. W., & Sandler, H. M. (1973). A note on the Solomon 4-group design:
#' Appropriate statistical analyses. *The Journal of Experimental Education,
#' 42*(2), 54–55. https://doi.org/10.1080/00220973.1973.11011460
#'
#' Sawilowsky, S. S., Kelley, D. L., Blair, R. C., & Markman, B. S. (1994).
#' Meta-analysis and the Solomon four-group design. *The Journal of Experimental
#' Education, 62*(4), 361–376. https://doi.org/10.1080/00220973.1994.9944140
#'
#' Solomon, R. L. (1949). An extension of control group design. *Psychological
#' Bulletin, 46*(2), 137–150. https://doi.org/10.1037/h0062958
#'
#' van Engelenburg, G. (1999). *Statistical analysis for the Solomon four-group
#' design* (Research Report 99-06). University of Twente.
#'
#' Walton Braver, M. C., & Braver, S. L. (1988). Statistical treatment of the
#' Solomon four-group design: A meta-analytic approach. *Psychological Bulletin,
#' 104*(1), 150–154. https://doi.org/10.1037/0033-2909.104.1.150
#'
#' @keywords Solomon four-group ANCOVA permutation maximum-likelihood SEM
#' @name solomonR
NULL

# ---- helpers ----

#' Convert p-value to Z (one-tailed) for Stouffer's method
#' @param p numeric vector of p-values assumed one-tailed and aligned in the same direction
#' @return numeric Z-scores
#' @references
#' Stouffer, S. A., Suchman, E. A., DeVinney, L. C., Star, S. A., & Williams, R.
#' M., Jr. (1949). *The American soldier: Adjustment during army life* (Vol. 1).
#' Princeton University Press.
#' @export
p_to_z <- function(p) {
  stats::qnorm(1 - p)
}

#' Stouffer's Z combiner (Walton Braver & Braver, 1988, Test I)
#'
#' Combine one-tailed p-values that test the *same directional* hypothesis into
#' a single Z. This is provided to reproduce the Walton Braver & Braver (1988)
#' meta-analytic option (Test I) for the Solomon four-group design. Use cautiously and document assumptions about homogeneity; see
#' the 1988–1990 exchanges for caveats.
#' @param p numeric vector of one-tailed p-values (same direction)
#' @return list with z_meta and p_meta (one-tailed)
#' @references
#' Sawilowsky, S. S., Kelley, D. L., Blair, R. C., & Markman, B. S. (1994).
#' Meta-analysis and the Solomon four-group design. *The Journal of Experimental
#' Education, 62*(4), 361–376. https://doi.org/10.1080/00220973.1994.9944140
#'
#' Stouffer, S. A., Suchman, E. A., DeVinney, L. C., Star, S. A., & Williams, R.
#' M., Jr. (1949). *The American soldier: Adjustment during army life* (Vol. 1).
#' Princeton University Press.
#'
#' Walton Braver, M. C., & Braver, S. L. (1988). Statistical treatment of the
#' Solomon four-group design: A meta-analytic approach. *Psychological Bulletin,
#' 104*(1), 150–154. https://doi.org/10.1037/0033-2909.104.1.150
#' @examples
#' stouffer_solomon(c(0.10, 0.11))
#' @export
stouffer_solomon <- function(p) {
  stopifnot(all(is.finite(p)), all(p > 0 & p < 1))
  z <- p_to_z(p)
  z_meta <- sum(z) / sqrt(length(z))
  list(z_meta = z_meta, p_meta_one_tailed = 1 - stats::pnorm(z_meta))
}

#' Fit the unified GLM for a Solomon Four-Group design
#'
#' Fits one generalized linear model to all four Solomon groups and reports
#' four Solomon contrasts: the equal-weighted average treatment effect, the
#' Pretest x Treatment (sensitization) contrast, and the treatment effect
#' within each pretesting condition.
#'
#' When `pretest_score` is supplied, it enters the model as `pre_obs`, equal
#' to the pretest score in the pretested groups and 0 in the unpretested
#' groups, so the structurally absent pretests do not remove Groups 3 and 4.
#' Regression adjustment for baseline covariates in randomized experiments,
#' and the case for pairing it with heteroskedasticity-robust standard errors,
#' is discussed by Lin (2013). Pretested participants with a missing pretest
#' score are excluded with a warning; incidental missingness is never imputed.
#'
#' @section Inference:
#' - `robust = "HC3"` (default): heteroskedasticity-consistent covariance
#'   (MacKinnon & White, 1985), recommended for routine use and particularly
#'   below about 250 observations (Long & Ervin, 2000; Hayes & Cai, 2007).
#'   Heteroskedasticity is expected in the unified Solomon model: when the
#'   pretest predicts the posttest, adjusting for it reduces residual variance
#'   only in the pretested groups.
#' - `robust = "none"`: conventional model-based covariance.
#' - `robust = "CR2"`: bias-reduced cluster-robust covariance (Bell &
#'   McCaffrey, 2002) with Satterthwaite degrees of freedom (Pustejovsky &
#'   Tipton, 2018), computed with the clubSandwich package. Use this when
#'   participants are nested in clusters such as classrooms or sites.
#'
#' Tests and confidence intervals use the same reference distribution. For
#' `"HC3"` and `"none"` it is the t distribution with residual degrees of
#' freedom when the dispersion is estimated, as in Gaussian models, which is
#' the conventional choice when t approximations are used with robust
#' standard errors (Imbens & Kolesár, 2016; Rajh-Weber et al., 2025). Families
#' with a fixed dispersion (binomial, Poisson) use the normal distribution.
#' With a noncollapsible link (such as the logit) and `pretest_score`, the
#' link-scale contrasts compare a treatment effect conditional on the pretest
#' among pretested participants with a marginal effect among unpretested
#' participants, who have no pretest. When the pretest predicts the outcome,
#' the Pretest x Treatment contrast is then nonzero even without
#' sensitization (Daniel et al., 2021): in the package's simulation study
#' (issue #43) it averaged 0.08 to 0.21 on the log-odds scale with no
#' sensitization present. Such fits give the classed warning
#' `solomonR_noncollapsible_warning`; for binary outcomes, estimate the
#' Solomon contrasts on a common scale with [marginal_solomon()]. Identity and
#' log links are collapsible and are not affected.
#'
#' The degrees of freedom are returned in the `df` columns (`Inf` for normal
#' reference distributions). Imbens and Kolesár (2016) further recommend
#' Bell-McCaffrey degrees of freedom for heteroskedasticity-robust intervals;
#' that refinement is under evaluation. For non-identity links, the contrasts
#' are on the link scale.
#'
#' In the package's simulation validation (issues #10 and #22), HC3 intervals
#' for the Solomon contrasts were conservative with 10 or fewer participants
#' per cell (mean coverage of nominal 95% intervals was 0.961 with 6 per cell
#' and 0.958 with 10, and the Pretest x Treatment test had a Type I error of
#' 0.034 with 6 per cell) and close to nominal with 20 or more (mean coverage
#' 0.951 to 0.954).
#'
#' Confidence intervals for the Wald partial R-squared use the noncentral F
#' method (Steiger, 2004) and are reported only for conventional Gaussian
#' fits; no corresponding interval is available with robust covariance.
#'
#' @param y numeric posttest vector
#' @param treat 0/1 (or logical) treatment indicator (1 = treatment)
#' @param pretested 0/1 (or logical) pretest indicator (1 = group received pretest)
#' @param pretest_score numeric vector for those pretested; NA for others
#' @param covariates optional data.frame of additional covariates
#' @param robust character: "HC3" (default), "none", or "CR2" (cluster-robust; requires `cluster`)
#' @param cluster optional clustering id (e.g., class/site), one value per
#'   participant. CR2 fits refuse designs in which a Solomon cell contains a
#'   single cluster, because cluster and condition are then confounded; see
#'   [validate_solomon()]. A classed warning (`solomonR_small_df_warning`)
#'   flags Solomon contrasts whose Satterthwaite degrees of freedom are below
#'   4, where Tipton (2015) advises that p-values not be trusted.
#' @param family model family (default gaussian())
#' @param conf_level confidence level for intervals (default 0.95)
#' @return An object of class `solomon_glm`: a list with the fitted model,
#'   coefficient and contrast tables (including degrees of freedom and
#'   confidence limits `conf.low` and `conf.high`), the covariance matrix,
#'   and the settings used.
#' @references
#' Bell, R. M., & McCaffrey, D. F. (2002). Bias reduction in standard errors for
#' linear regression with multi-stage samples. *Survey Methodology, 28*(2),
#' 169–181.
#'
#' Daniel, R., Zhang, J., & Farewell, D. (2021). Making apples from oranges:
#' Comparing noncollapsible effect estimators and their standard errors after
#' adjustment for different covariate sets. *Biometrical Journal, 63*(3),
#' 528–557. https://doi.org/10.1002/bimj.201900297
#'
#' Hayes, A. F., & Cai, L. (2007). Using heteroskedasticity-consistent standard
#' error estimators in OLS regression: An introduction and software
#' implementation. *Behavior Research Methods, 39*(4), 709–722.
#' https://doi.org/10.3758/BF03192961
#'
#' Imbens, G. W., & Kolesár, M. (2016). Robust standard errors in small samples:
#' Some practical advice. *The Review of Economics and Statistics, 98*(4),
#' 701–712. https://doi.org/10.1162/REST_a_00552
#'
#' Lin, W. (2013). Agnostic notes on regression adjustments to experimental
#' data: Reexamining Freedman's critique. *The Annals of Applied Statistics,
#' 7*(1), 295–318. https://doi.org/10.1214/12-AOAS583
#'
#' Long, J. S., & Ervin, L. H. (2000). Using heteroscedasticity consistent
#' standard errors in the linear regression model. *The American Statistician,
#' 54*(3), 217–224. https://doi.org/10.1080/00031305.2000.10474549
#'
#' MacKinnon, J. G., & White, H. (1985). Some heteroskedasticity-consistent
#' covariance matrix estimators with improved finite sample properties. *Journal
#' of Econometrics, 29*(3), 305–325.
#' https://doi.org/10.1016/0304-4076(85)90158-7
#'
#' Pustejovsky, J. E., & Tipton, E. (2018). Small-sample methods for
#' cluster-robust variance estimation and hypothesis testing in fixed effects
#' models. *Journal of Business & Economic Statistics, 36*(4), 672–683.
#' https://doi.org/10.1080/07350015.2016.1247004
#'
#' Rajh-Weber, H., Huber, S. E., & Arendasy, M. (2025). A practice-oriented
#' guide to statistical inference in linear modeling for non-normal or
#' heteroskedastic error distributions. *Behavior Research Methods, 57*(12),
#' Article 338. https://doi.org/10.3758/s13428-025-02801-4
#'
#' Solomon, R. L. (1949). An extension of control group design. *Psychological
#' Bulletin, 46*(2), 137–150. https://doi.org/10.1037/h0062958
#'
#' Steiger, J. H. (2004). Beyond the F test: Effect size confidence intervals
#' and tests of close fit in the analysis of variance and contrast analysis.
#' *Psychological Methods, 9*(2), 164–182.
#' https://doi.org/10.1037/1082-989X.9.2.164
#'
#' Tipton, E. (2015). Small sample adjustments for robust variance estimation
#' with meta-regression. *Psychological Methods, 20*(3), 375–393.
#' https://doi.org/10.1037/met0000011
#' @export
fit_solomon_glm <- function(y, treat, pretested, pretest_score = NULL,
                            covariates = NULL, robust = c("HC3", "none", "CR2"),
                            cluster = NULL, family = stats::gaussian(),
                            conf_level = 0.95) {
  robust <- match.arg(robust)
  .check_conf_level(conf_level)

  treat <- .solomon_indicator(treat, "treat")
  pretested <- .solomon_indicator(pretested, "pretested")

  .solomon_check_lengths(
    y = y,
    treat = treat,
    pretested = pretested,
    pretest_score = pretest_score,
    covariates = covariates,
    cluster = cluster
  )

  df <- data.frame(
    y = y,
    treat = treat,
    pretested = pretested
  )

  if (!is.null(pretest_score)) {

    n_incidental <- sum(pretested == 1L & is.na(pretest_score), na.rm = TRUE)

    if (n_incidental > 0L) {
      warning(
        n_incidental, " pretested participant(s) have missing pretest scores ",
        "and are excluded from the model. Incidental missingness is not imputed.",
        call. = FALSE
      )
    }

    if (any(pretested == 0L & !is.na(pretest_score), na.rm = TRUE)) {
      warning(
        "Observed pretest_score values were supplied for unpretested ",
        "participants; these values are ignored.",
        call. = FALSE
      )
    }

    # Safe pretest covariate: equals pretest_score in pretested rows, 0 otherwise
    df$pre_obs <- ifelse(df$pretested == 1, pretest_score, 0)
  }

  if (!is.null(covariates)) df <- cbind(df, covariates)

  # One model: treatment + pretest indicator + their interaction + (optional) pre_obs
  rhs <- c("treat*pretested",
           if (!is.null(pretest_score)) "pre_obs" else NULL,
           if (!is.null(covariates)) names(covariates) else NULL)
  fml <- stats::as.formula(paste("y ~", paste(rhs, collapse = " + ")))

  fit <- stats::glm(fml, data = df, family = family, na.action = stats::na.exclude)

  if (!is.null(pretest_score) && !stats::family(fit)$link %in% c("identity", "log")) {
    .warn_noncollapsible(stats::family(fit)$link)
  }

  # Robust VCOV
  fit_cr <- fit
  if (robust == "HC3") {
    vcovM <- sandwich::vcovHC(fit, type = "HC3")
  } else if (robust == "CR2") {
    if (is.null(cluster)) stop("CR2 requested but 'cluster' is NULL.", call. = FALSE)

    # Align the clustering variable with the rows retained by the model.
    used <- rep(TRUE, nrow(df))
    if (!is.null(fit$na.action)) {
      used[fit$na.action] <- FALSE
    }
    cluster_fit <- cluster[used]
    if (anyNA(cluster_fit)) {
      stop("`cluster` is missing for participants included in the model.", call. = FALSE)
    }
    single <- .cluster_structure(df$treat[used], df$pretested[used], cluster_fit)$single_cluster
    if (length(single)) {
      .stop_confounded_clusters(single)
    }

    # clubSandwich cannot use the NA-padded residuals of an na.exclude fit,
    # so CR2 quantities come from a complete-case refit of the same rows
    # (identical coefficients).
    fit_cr <- stats::glm(fml, data = df[used, , drop = FALSE], family = family)
    vcovM <- clubSandwich::vcovCR(fit_cr, cluster = cluster_fit, type = "CR2")
  } else {
    vcovM <- stats::vcov(fit)
  }

  # --- reference distribution ---
  # CR2 tests use Satterthwaite degrees of freedom. Otherwise, as in
  # summary.glm(), tests use t with residual df when the dispersion is
  # estimated (e.g., Gaussian models) and the normal distribution when it is
  # fixed (binomial, Poisson).
  dispersion_fixed <- stats::family(fit)$family %in% c("binomial", "poisson")
  df_model <- if (dispersion_fixed) Inf else stats::df.residual(fit)

  # --- tidy coefficients with chosen vcov ---
  tidy <- broom::tidy(fit)
  se_vec <- sqrt(diag(vcovM))
  tidy$std.error <- unname(se_vec[match(tidy$term, names(se_vec))])
  tidy$statistic <- tidy$estimate / tidy$std.error
  if (robust == "CR2") {
    ct <- clubSandwich::coef_test(fit_cr, vcov = vcovM, test = "Satterthwaite")
    tidy$df <- ct$df_Satt[match(tidy$term, ct$Coef)]
  } else {
    tidy$df <- df_model
  }
  tidy$p.value <- 2 * stats::pt(-abs(tidy$statistic), df = tidy$df)
  coef_ci <- .wald_ci(tidy$estimate, tidy$std.error, tidy$df, conf_level)
  tidy$conf.low <- unname(coef_ci[, "conf.low"])
  tidy$conf.high <- unname(coef_ci[, "conf.high"])

  # --- linear contrasts (one row each) ---
  cf <- stats::coef(fit)
  cn <- names(cf)
  Z <- function() { v <- numeric(length(cn)); names(v) <- cn; v }

  lin_contrast <- function(L) {
    L <- matrix(L, nrow = 1)
    est <- as.numeric(L %*% cf)
    se  <- sqrt(as.numeric(L %*% vcovM %*% t(L)))
    z   <- est / se
    dfc <- if (robust == "CR2") {
      as.data.frame(
        clubSandwich::linear_contrast(
          fit_cr,
          vcov = vcovM,
          contrasts = L,
          test = "Satterthwaite"
        )
      )$df
    } else {
      df_model
    }
    p   <- 2 * stats::pt(-abs(z), df = dfc)
    ci  <- .wald_ci(est, se, dfc, conf_level)
    c(estimate = est, std.error = se, statistic = z, p.value = p, df = dfc,
      conf.low = unname(ci[, "conf.low"]), conf.high = unname(ci[, "conf.high"]))
  }

  # Equal-weighted average treatment effect across pretest conditions:
  # beta_treat + 0.5 * beta_treat:pretested
  L_ate <- Z()
  L_ate["treat"] <- 1
  if ("treat:pretested" %in% cn) {
    L_ate["treat:pretested"] <- 0.5
  }

  # Pretest x Treatment interaction:
  # difference between the treatment effect in pretested versus
  # unpretested participants
  L_int <- Z()
  if ("treat:pretested" %in% cn) {
    L_int["treat:pretested"] <- 1
  }

  # Treatment effect among pretested participants:
  # beta_treat + beta_treat:pretested
  L_pre <- Z()
  L_pre["treat"] <- 1
  if ("treat:pretested" %in% cn) {
    L_pre["treat:pretested"] <- 1
  }

  # Treatment effect among unpretested participants:
  # beta_treat
  L_un <- Z()
  L_un["treat"] <- 1

  c_ate <- lin_contrast(L_ate)
  c_int <- lin_contrast(L_int)
  c_pre <- lin_contrast(L_pre)
  c_un  <- lin_contrast(L_un)

  # Wald-based partial R^2 for each contrast (intervals only for
  # conventional covariance)
  conventional <- robust == "none"
  r2_ate <- contrast_r2_ci(fit, L_ate, vcovM, conf_level, conventional)
  r2_int <- contrast_r2_ci(fit, L_int, vcovM, conf_level, conventional)
  r2_pre <- contrast_r2_ci(fit, L_pre, vcovM, conf_level, conventional)
  r2_un  <- contrast_r2_ci(fit, L_un,  vcovM, conf_level, conventional)

  effects <- data.frame(
    contrast = c(
      "ATE (avg over pretest)",
      "Pretest x Treatment",
      "Treatment | pretested",
      "Treatment | unpretested"
    ),

    estimate = c(
      unname(c_ate["estimate"]),
      unname(c_int["estimate"]),
      unname(c_pre["estimate"]),
      unname(c_un["estimate"])
    ),

    std.error = c(
      unname(c_ate["std.error"]),
      unname(c_int["std.error"]),
      unname(c_pre["std.error"]),
      unname(c_un["std.error"])
    ),

    statistic = c(
      unname(c_ate["statistic"]),
      unname(c_int["statistic"]),
      unname(c_pre["statistic"]),
      unname(c_un["statistic"])
    ),

    p.value = c(
      unname(c_ate["p.value"]),
      unname(c_int["p.value"]),
      unname(c_pre["p.value"]),
      unname(c_un["p.value"])
    ),

    df = c(
      unname(c_ate["df"]),
      unname(c_int["df"]),
      unname(c_pre["df"]),
      unname(c_un["df"])
    ),

    conf.low = c(
      unname(c_ate["conf.low"]),
      unname(c_int["conf.low"]),
      unname(c_pre["conf.low"]),
      unname(c_un["conf.low"])
    ),

    conf.high = c(
      unname(c_ate["conf.high"]),
      unname(c_int["conf.high"]),
      unname(c_pre["conf.high"]),
      unname(c_un["conf.high"])
    ),

    r2 = c(
      r2_ate$r2,
      r2_int$r2,
      r2_pre$r2,
      r2_un$r2
    ),

    r2_lo = c(
      r2_ate$r2_lo,
      r2_int$r2_lo,
      r2_pre$r2_lo,
      r2_un$r2_lo
    ),

    r2_hi = c(
      r2_ate$r2_hi,
      r2_int$r2_hi,
      r2_pre$r2_hi,
      r2_un$r2_hi
    ),

    row.names = NULL,
    stringsAsFactors = FALSE
  )

  if (robust == "CR2") {
    small_df <- is.finite(effects$df) & effects$df < 4
    if (any(small_df)) {
      .warn_cr2_small_df(effects$contrast[small_df], effects$df[small_df])
    }
  }

  out <- list(
    model = fit,
    coefficients = tidy,
    effects = effects,
    vcov = vcovM,
    data = df,
    robust = robust,
    family = family,
    cluster = cluster,
    conf_level = conf_level,
    call = match.call()
  )

  class(out) <- "solomon_glm"

  out
}

#' Permutation test for a Solomon contrast
#'
#' Performs a randomization-based test by permuting treatment assignment
#' within pretest strata. This preserves the Solomon four-group design while
#' generating the null distribution for a selected treatment contrast.
#'
#' @details
#' The test statistic is the HC3-studentized contrast. The permutation
#' p-value is a valid test of the sharp null hypothesis that treatment has no
#' effect for any participant; the `+1` correction keeps the Monte Carlo
#' p-value from being zero (Phipson & Smyth, 2010). Studentizing the
#' statistic makes permutation tests asymptotically robust when only an
#' average effect is hypothesized to be zero (DiCiccio & Romano, 2017;
#' Wu & Ding, 2021); for the Pretest x Treatment contrast that robustness
#' should be regarded as approximate.
#'
#' Randomization inference must permute the unit that was randomized.
#' Because this function permutes individual participants, it refuses fits
#' that include a clustering variable. For clustered designs, use the CR2
#' small-sample tests reported by [fit_solomon_glm()].
#'
#' @param object An object returned by \code{fit_solomon_glm()} without a
#'   clustering variable.
#' @param contrast Character string identifying the contrast to test. One of
#'   \code{"ATE (avg over pretest)"}, \code{"Pretest x Treatment"},
#'   \code{"Treatment | pretested"}, or
#'   \code{"Treatment | unpretested"}.
#' @param reps Number of permutations. Default is 5000.
#' @param seed Optional random-number seed for reproducibility. The global
#'   random-number state is restored when the function exits.
#' @param return_dist Logical. If \code{TRUE}, return the permutation
#'   distribution in addition to the observed statistic and p-value.
#'
#' @return A list containing the observed studentized statistic
#'   (\code{z_obs}) and permutation p-value (\code{p_perm}). If
#'   \code{return_dist = TRUE}, the permutation distribution
#'   (\code{z_perm}) is also returned.
#'
#' @references
#' DiCiccio, C. J., & Romano, J. P. (2017). Robust permutation tests for
#' correlation and regression coefficients. *Journal of the American Statistical
#' Association, 112*(519), 1211–1220.
#' https://doi.org/10.1080/01621459.2016.1202117
#'
#' Phipson, B., & Smyth, G. K. (2010). Permutation p-values should never be
#' zero: Calculating exact p-values when permutations are randomly drawn.
#' *Statistical Applications in Genetics and Molecular Biology, 9*(1),
#' Article 39. https://doi.org/10.2202/1544-6115.1585
#'
#' Wu, J., & Ding, P. (2021). Randomization tests for weak null hypotheses in
#' randomized experiments. *Journal of the American Statistical Association,
#' 116*(536), 1898–1913. https://doi.org/10.1080/01621459.2020.1750415
#'
#' @export
perm_solomon <- function(
    object,
    contrast = "ATE (avg over pretest)",
    reps = 5000L,
    seed = NULL,
    return_dist = FALSE
) {

  if (!inherits(object, "solomon_glm")) {
    stop("object must be from fit_solomon_glm().")
  }

  if (!is.null(object$cluster) || identical(object$robust, "CR2")) {
    stop(
      "perm_solomon() permutes individual treatment labels within pretest ",
      "strata, which is not a valid randomization test when treatment was ",
      "assigned to clusters. Use the CR2 small-sample tests reported by ",
      "fit_solomon_glm() for clustered designs; cluster-level randomization ",
      "inference is planned for a future release.",
      call. = FALSE
    )
  }

  valid_contrasts <- c(
    "ATE (avg over pretest)",
    "Pretest x Treatment",
    "Treatment | pretested",
    "Treatment | unpretested"
  )

  if (!contrast %in% valid_contrasts) {
    stop(
      "Unknown contrast. Choose one of: ",
      paste(valid_contrasts, collapse = ", ")
    )
  }

  reps <- as.integer(reps)

  if (length(reps) != 1L || is.na(reps) || reps < 1L) {
    stop("reps must be a positive integer.")
  }

  if (!is.null(seed)) {
    withr::local_seed(seed)
  }

  df <- object$data
  form <- stats::formula(object$model)
  fam <- stats::family(object$model)

  # ------------------------------------------------------------
  # Helper: construct the requested linear contrast
  # ------------------------------------------------------------

  make_contrast <- function(coef_names, contrast) {

    L <- numeric(length(coef_names))
    names(L) <- coef_names

    if (!"treat" %in% coef_names) {
      stop("The fitted model does not contain a treatment coefficient.")
    }

    has_interaction <- "treat:pretested" %in% coef_names

    if (contrast == "ATE (avg over pretest)") {

      # Equal-weighted average of the treatment effects across the
      # pretested and unpretested conditions:
      #
      # beta_treat + 0.5 * beta_treat:pretested

      L["treat"] <- 1

      if (has_interaction) {
        L["treat:pretested"] <- 0.5
      }

    } else if (contrast == "Pretest x Treatment") {

      if (!has_interaction) {
        stop("The fitted model does not contain a treatment-by-pretest interaction.")
      }

      L["treat:pretested"] <- 1

    } else if (contrast == "Treatment | pretested") {

      # beta_treat + beta_treat:pretested

      L["treat"] <- 1

      if (has_interaction) {
        L["treat:pretested"] <- 1
      }

    } else if (contrast == "Treatment | unpretested") {

      # beta_treat

      L["treat"] <- 1
    }

    L
  }

  # ------------------------------------------------------------
  # Helper: studentized statistic for a fitted model
  # ------------------------------------------------------------

  contrast_z <- function(fit, contrast) {

    b <- stats::coef(fit)

    V <- sandwich::vcovHC(
      fit,
      type = "HC3"
    )

    # Force coefficient and covariance-matrix order to match.
    V <- V[names(b), names(b), drop = FALSE]

    L <- make_contrast(
      coef_names = names(b),
      contrast = contrast
    )

    Lm <- matrix(L, nrow = 1)

    estimate <- as.numeric(Lm %*% b)

    variance <- as.numeric(
      Lm %*% V %*% t(Lm)
    )

    if (!is.finite(variance) || variance <= 0) {
      return(NA_real_)
    }

    estimate / sqrt(variance)
  }

  # ------------------------------------------------------------
  # Observed statistic
  # ------------------------------------------------------------

  z_obs <- contrast_z(
    object$model,
    contrast
  )

  if (!is.finite(z_obs)) {
    stop("Could not calculate the observed permutation-test statistic.")
  }

  # ------------------------------------------------------------
  # Permutation distribution
  #
  # Treatment assignment is shuffled WITHIN pretest strata.
  # This preserves the Solomon design and the number assigned
  # to treatment/control within each pretesting condition.
  # ------------------------------------------------------------

  z_perm <- numeric(reps)

  for (i in seq_len(reps)) {

    df_perm <- df

    df_perm$treat <- stats::ave(
      df$treat,
      df$pretested,
      FUN = function(x) sample(x, length(x), replace = FALSE)
    )

    fit_perm <- stats::glm(
      form,
      data = df_perm,
      family = fam,
      na.action = stats::na.exclude
    )

    z_perm[i] <- contrast_z(
      fit_perm,
      contrast
    )
  }

  # In the unlikely event that a permuted model produces an undefined
  # statistic, exclude that replicate explicitly rather than silently
  # allowing NA propagation.

  valid <- is.finite(z_perm)

  if (!any(valid)) {
    stop("No valid permutation statistics were obtained.")
  }

  z_perm_valid <- z_perm[valid]

  # Two-sided randomization p-value with the +1 finite-simulation
  # correction. This prevents an estimated p-value of exactly zero.

  p_perm <- (
    sum(abs(z_perm_valid) >= abs(z_obs)) + 1
  ) / (
    length(z_perm_valid) + 1
  )

  out <- list(
    contrast = contrast,
    z_obs = z_obs,
    p_perm = p_perm,
    reps = reps,
    valid_reps = length(z_perm_valid)
  )

  if (isTRUE(return_dist)) {
    out$z_perm <- z_perm_valid
  }

  class(out) <- "solomon_perm"

  out
}
