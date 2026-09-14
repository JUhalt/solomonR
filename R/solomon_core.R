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
#' Braver & Braver (1988), the Sawilowsky and Markman methodological
#' exchanges, and van Engelenburg (1999).
#'
#' @references
#' Solomon, R. L. (1949). An extension of control group design.
#' *Psychological Bulletin, 46*(2), 137-150.
#'
#' Campbell, D. T., & Stanley, J. C. (1963). *Experimental and
#' quasi-experimental designs for research*. Rand McNally.
#'
#' Huck, S. W., & Sandler, H. M. (1973). A note on the Solomon 4-group
#' design: Appropriate statistical analyses. *The Journal of Experimental
#' Education, 42*(1), 54-55.
#'
#' Braver, M. W., & Braver, S. L. (1988). Statistical treatment of the
#' Solomon four-group design: A meta-analytic approach. *Psychological
#' Bulletin, 104*(1), 150-154.
#'
#' Sawilowsky, S. S., Kelley, D. L., Blair, R. C., & Markman, B. S. (1994).
#' Meta-analysis and the Solomon four-group design. *The Journal of
#' Experimental Education, 62*(4), 361-376.
#'
#' van Engelenburg, G. (1999). *Statistical analysis for the Solomon
#' four-group design* (Research Report 99-06). University of Twente.
#'
#' @keywords Solomon four-group ANCOVA permutation maximum-likelihood SEM
#' @name solomonR
NULL

# ---- helpers ----

#' Convert p-value to Z (one-tailed) for Stouffer's method
#' @param p numeric vector of p-values assumed one-tailed and aligned in the same direction
#' @return numeric Z-scores
#' @references
#' Stouffer, S. A., Suchman, E. A., DeVinney, L. C., Star, S. A., &
#' Williams, R. M., Jr. (1949). *The American soldier: Adjustment during
#' army life* (Vol. 1). Princeton University Press.
#' @export
p_to_z <- function(p) {
  stats::qnorm(1 - p)
}

#' Stouffer's Z combiner (a.k.a. "Test I")
#'
#' Combine one-tailed p-values that test the *same directional* hypothesis into
#' a single Z. This is provided to reproduce the Braver & Braver (1988) option
#' for SFGD. Use cautiously and document assumptions about homogeneity; see
#' the 1988–1990 exchanges for caveats.
#' @param p numeric vector of one-tailed p-values (same direction)
#' @return list with z_meta and p_meta (one-tailed)
#' @references
#' Stouffer, S. A., Suchman, E. A., DeVinney, L. C., Star, S. A., &
#' Williams, R. M., Jr. (1949). *The American soldier: Adjustment during
#' army life* (Vol. 1). Princeton University Press.
#'
#' Braver, M. W., & Braver, S. L. (1988). Statistical treatment of the
#' Solomon four-group design: A meta-analytic approach. *Psychological
#' Bulletin, 104*(1), 150-154.
#'
#' Sawilowsky, S. S., Kelley, D. L., Blair, R. C., & Markman, B. S. (1994).
#' Meta-analysis and the Solomon four-group design. *The Journal of
#' Experimental Education, 62*(4), 361-376.
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
#' - `robust = "none"`: model-based covariance with normal-reference tests.
#' - `robust = "HC3"`: heteroskedasticity-consistent covariance (MacKinnon &
#'   White, 1985), recommended for samples of modest size (Long & Ervin,
#'   2000), with normal-reference tests.
#' - `robust = "CR2"`: bias-reduced cluster-robust covariance (Bell &
#'   McCaffrey, 2002) with Satterthwaite degrees of freedom (Pustejovsky &
#'   Tipton, 2018), computed with the clubSandwich package. Use this when
#'   participants are nested in clusters such as classrooms or sites.
#'
#' The degrees of freedom used for each test are returned in the `df`
#' columns (`Inf` for normal-reference tests). For non-identity links, the
#' contrasts are on the link scale.
#'
#' @param y numeric posttest vector
#' @param treat 0/1 (or logical) treatment indicator (1 = treatment)
#' @param pretested 0/1 (or logical) pretest indicator (1 = group received pretest)
#' @param pretest_score numeric vector for those pretested; NA for others
#' @param covariates optional data.frame of additional covariates
#' @param robust character: "none", "HC3", or "CR2" (cluster-robust; requires `cluster`)
#' @param cluster optional clustering id (e.g., class/site), one value per participant
#' @param family model family (default gaussian())
#' @return An object of class `solomon_glm`: a list with the fitted model,
#'   coefficient and contrast tables (including degrees of freedom), the
#'   covariance matrix, and the settings used.
#' @references
#' Bell, R. M., & McCaffrey, D. F. (2002). Bias reduction in standard errors
#' for linear regression with multi-stage samples. *Survey Methodology,
#' 28*(2), 169-181.
#'
#' Lin, W. (2013). Agnostic notes on regression adjustments to experimental
#' data: Reexamining Freedman's critique. *The Annals of Applied Statistics,
#' 7*(1), 295-318.
#'
#' Long, J. S., & Ervin, L. H. (2000). Using heteroscedasticity consistent
#' standard errors in the linear regression model. *The American
#' Statistician, 54*(3), 217-224.
#'
#' MacKinnon, J. G., & White, H. (1985). Some heteroskedasticity-consistent
#' covariance matrix estimators with improved finite sample properties.
#' *Journal of Econometrics, 29*(3), 305-325.
#'
#' Pustejovsky, J. E., & Tipton, E. (2018). Small-sample methods for
#' cluster-robust variance estimation and hypothesis testing in fixed
#' effects models. *Journal of Business & Economic Statistics, 36*(4),
#' 672-683.
#'
#' Solomon, R. L. (1949). An extension of control group design.
#' *Psychological Bulletin, 46*(2), 137-150.
#' @export
fit_solomon_glm <- function(y, treat, pretested, pretest_score = NULL,
                            covariates = NULL, robust = c("none","HC3","CR2"),
                            cluster = NULL, family = stats::gaussian()) {
  robust <- match.arg(robust)

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

    # clubSandwich cannot use the NA-padded residuals of an na.exclude fit,
    # so CR2 quantities come from a complete-case refit of the same rows
    # (identical coefficients).
    fit_cr <- stats::glm(fml, data = df[used, , drop = FALSE], family = family)
    vcovM <- clubSandwich::vcovCR(fit_cr, cluster = cluster_fit, type = "CR2")
  } else {
    vcovM <- stats::vcov(fit)
  }

  # --- tidy coefficients with chosen vcov ---
  # CR2 tests use Satterthwaite degrees of freedom; the other covariance
  # options use a normal reference distribution (df = Inf).
  tidy <- broom::tidy(fit)
  se_vec <- sqrt(diag(vcovM))
  tidy$std.error <- unname(se_vec[match(tidy$term, names(se_vec))])
  tidy$statistic <- tidy$estimate / tidy$std.error
  if (robust == "CR2") {
    ct <- clubSandwich::coef_test(fit_cr, vcov = vcovM, test = "Satterthwaite")
    tidy$df <- ct$df_Satt[match(tidy$term, ct$Coef)]
  } else {
    tidy$df <- Inf
  }
  tidy$p.value <- 2 * stats::pt(-abs(tidy$statistic), df = tidy$df)

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
      Inf
    }
    p   <- 2 * stats::pt(-abs(z), df = dfc)
    c(estimate = est, std.error = se, statistic = z, p.value = p, df = dfc)
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

  # semi-partial R^2 for each contrast
  r2_ate <- contrast_r2_ci(fit, L_ate, vcovM)
  r2_int <- contrast_r2_ci(fit, L_int, vcovM)
  r2_pre <- contrast_r2_ci(fit, L_pre, vcovM)
  r2_un  <- contrast_r2_ci(fit, L_un,  vcovM)

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

  out <- list(
    model = fit,
    coefficients = tidy,
    effects = effects,
    vcov = vcovM,
    data = df,
    robust = robust,
    family = family,
    cluster = cluster,
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
#' correlation and regression coefficients. *Journal of the American
#' Statistical Association, 112*(519), 1211-1220.
#'
#' Phipson, B., & Smyth, G. K. (2010). Permutation p-values should never be
#' zero: Calculating exact p-values when permutations are randomly drawn.
#' *Statistical Applications in Genetics and Molecular Biology, 9*(1),
#' Article 39.
#'
#' Wu, J., & Ding, P. (2021). Randomization tests for weak null hypotheses
#' in randomized experiments. *Journal of the American Statistical
#' Association, 116*(536), 1898-1913.
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

#' Power simulation for Solomon designs (experimental)
#'
#' **Experimental.** This helper is scheduled to be rebuilt and validated in a
#' later release. Its simulator has not been validated, and its results
#' should not be used for study planning.
#'
#' Simulates normally distributed Solomon four-group data and estimates the
#' rejection rate of several Solomon tests at alpha = .05. Pretest scores are
#' standard normal. The posttest residual standard deviation is `sigma` in
#' every cell, and the pretest-posttest correlation in the pretested cells is
#' `rho`. The treatment effect is `delta` among unpretested participants and
#' `delta + sens` among pretested participants, so the equal-weighted average
#' treatment effect is `delta + sens / 2`.
#'
#' @param n Cell sizes: a single number used for all four cells, or a list
#'   with elements `n1` (pretested treatment), `n2` (pretested control),
#'   `n3` (unpretested treatment), and `n4` (unpretested control).
#' @param delta Treatment effect among unpretested participants, on the
#'   posttest scale.
#' @param rho Pretest-posttest correlation in the pretested cells.
#' @param sens Sensitization: the additional treatment effect among pretested
#'   participants (0 = none).
#' @param sigma Posttest residual standard deviation in all cells.
#' @param sims Number of Monte Carlo replicates.
#' @param stouffer Logical; if `TRUE`, also estimate the rejection rate of the
#'   historical Stouffer Test I, evaluated one-tailed (treatment > control)
#'   as in [fit_solomon_classic()].
#' @return A data frame with the estimated rejection rate for the 2x2 ANOVA
#'   interaction, the GLM average treatment effect, the two simple treatment
#'   effects, and Test I (`NA` when `stouffer = FALSE`).
#' @export
power_solomon <- function(n = list(n1=50,n2=50,n3=50,n4=50),
                          delta = 0.3, rho = 0.5, sens = 0, sigma = 1, sims = 2000,
                          stouffer = TRUE) {
  warning(
    "power_solomon() is experimental and has not been validated; do not use ",
    "it for study planning. A rebuilt planning framework is scheduled for a ",
    "later release.",
    call. = FALSE
  )

  if (length(n) == 1) n <- as.list(rep(n, 4))
  if (is.null(names(n))) names(n) <- paste0("n", 1:4)

  out <- matrix(NA_real_, nrow = sims, ncol = 5)
  colnames(out) <- c("A_interaction","A_ATE","A_pre","A_unpre","A_stouffer")

  for (s in seq_len(sims)) {
    # simulate: pretested cells use a bivariate-normal construction so that
    # cor(pre, post) = rho and the residual SD is sigma in every cell
    pre1 <- stats::rnorm(n$n1); pre2 <- stats::rnorm(n$n2)
    e1 <- stats::rnorm(n$n1); e2 <- stats::rnorm(n$n2)
    post1 <- delta + sens + sigma * (rho*pre1 + sqrt(1-rho^2)*e1)
    post2 <- 0            + sigma * (rho*pre2 + sqrt(1-rho^2)*e2)
    post3 <- delta + stats::rnorm(n$n3, 0, sigma)
    post4 <- 0     + stats::rnorm(n$n4, 0, sigma)

    y_post <- c(post1, post2, post3, post4)
    treat  <- c(rep(1,n$n1), rep(0,n$n2), rep(1,n$n3), rep(0,n$n4))
    preind <- c(rep(1,n$n1 + n$n2), rep(0, n$n3 + n$n4))
    y_pre  <- c(pre1, pre2, rep(NA, n$n3+n$n4))

    # (i) interaction from 2x2 ANOVA
    aov_tab <- broom::tidy(stats::aov(y_post ~ factor(treat)*factor(preind)))
    out[s,"A_interaction"] <- as.numeric(aov_tab$p.value[aov_tab$term=="factor(treat):factor(preind)"] < .05)

    # (ii) unified GLM ATE
    g <- fit_solomon_glm(y_post, treat, preind, y_pre, robust = "HC3")
    ATEp <- g$effects$p.value[match("ATE (avg over pretest)", g$effects$contrast)]
    out[s,"A_ATE"] <- as.numeric(ATEp < .05)

    # (iii) simple effects
    preP <- g$effects$p.value[match("Treatment | pretested", g$effects$contrast)]
    unP  <- g$effects$p.value[match("Treatment | unpretested", g$effects$contrast)]
    out[s,"A_pre"]   <- as.numeric(preP < .05)
    out[s,"A_unpre"] <- as.numeric(unP  < .05)

    # (iv) optional historical Stouffer Test I: one-tailed (treatment >
    # control) ANCOVA in groups 1-2 combined with a t test in groups 3-4
    if (stouffer) {
      pre_rows <- preind == 1
      anc <- stats::lm(y_post[pre_rows] ~ treat[pre_rows] + y_pre[pre_rows])
      t_anc <- summary(anc)$coefficients[2, "t value"]
      p_anc <- stats::pt(t_anc, df = anc$df.residual, lower.tail = FALSE)
      un <- stats::lm(y_post[!pre_rows] ~ treat[!pre_rows])
      t_un <- summary(un)$coefficients[2, "t value"]
      p_un <- stats::pt(t_un, df = un$df.residual, lower.tail = FALSE)
      eps <- .Machine$double.eps
      p_one <- pmin(pmax(c(p_anc, p_un), eps), 1 - eps)
      out[s,"A_stouffer"] <- as.numeric(stouffer_solomon(p_one)$p_meta_one_tailed < .05)
    }
  }

  data.frame(
    metric = colnames(out),
    power  = colMeans(out),
    row.names = NULL
  )
}
