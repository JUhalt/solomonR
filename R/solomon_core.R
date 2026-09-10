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
#' @keywords Solomon four-group ANCOVA permutation maximum-likelihood SEM
#' @name solomonR
NULL

# ---- helpers ----

#' Convert p-value to Z (one-tailed) for Stouffer's method
#' @param p numeric vector of p-values assumed one-tailed and aligned in the same direction
#' @return numeric Z-scores
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
#' @param y numeric posttest vector
#' @param treat 0/1 indicator (1 = treatment)
#' @param pretested 0/1 indicator (1 = group received pretest)
#' @param pretest_score numeric vector for those pretested; NA for others
#' @param covariates optional data.frame of additional covariates
#' @param robust character: "none","HC3","CR2" (cluster-robust via clubSandwich if `cluster` supplied)
#' @param cluster optional clustering id (e.g., class/site)
#' @param family model family (default gaussian())
#' @return a list with model, tidy tables, and predefined contrasts (ATE, interaction, simple effects)
#' @export
fit_solomon_glm <- function(y, treat, pretested, pretest_score = NULL,
                            covariates = NULL, robust = c("none","HC3","CR2"),
                            cluster = NULL, family = stats::gaussian()) {
  robust <- match.arg(robust)

  df <- data.frame(
    y = y,
    treat = as.integer(treat),
    pretested = as.integer(pretested)
  )

  # Safe pretest covariate: equals pretest_score in pretested rows, 0 otherwise
  if (!is.null(pretest_score)) {
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
  if (robust == "HC3") {
    vcovM <- sandwich::vcovHC(fit, type = "HC3")
  } else if (robust == "CR2") {
    if (is.null(cluster)) stop("CR2 requested but 'cluster' is NULL.")
    vcovM <- clubSandwich::vcovCR(fit, cluster = cluster, type = "CR2")
  } else {
    vcovM <- stats::vcov(fit)
  }

  # --- tidy coefficients with chosen vcov ---
  tidy <- broom::tidy(fit)
  se_vec <- sqrt(diag(vcovM))
  tidy$std.error <- unname(se_vec[match(tidy$term, names(se_vec))])
  tidy$statistic <- tidy$estimate / tidy$std.error
  tidy$p.value <- 2 * stats::pnorm(-abs(tidy$statistic))

  # --- linear contrasts (one row each) ---
  cf <- stats::coef(fit)
  cn <- names(cf)
  Z <- function() { v <- numeric(length(cn)); names(v) <- cn; v }

  lin_contrast <- function(L) {
    L <- matrix(L, nrow = 1)
    est <- as.numeric(L %*% cf)
    se  <- sqrt(as.numeric(L %*% vcovM %*% t(L)))
    z   <- est / se
    p   <- 2 * stats::pnorm(-abs(z))
    c(estimate = est, std.error = se, statistic = z, p.value = p)
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
#' @param object An object returned by \code{fit_solomon_glm()}.
#' @param contrast Character string identifying the contrast to test. One of
#'   \code{"ATE (avg over pretest)"}, \code{"Pretest x Treatment"},
#'   \code{"Treatment | pretested"}, or
#'   \code{"Treatment | unpretested"}.
#' @param reps Number of permutations. Default is 5000.
#' @param seed Optional random-number seed for reproducibility.
#' @param return_dist Logical. If \code{TRUE}, return the permutation
#'   distribution in addition to the observed statistic and p-value.
#'
#' @return A list containing the observed studentized statistic
#'   (\code{z_obs}) and permutation p-value (\code{p_perm}). If
#'   \code{return_dist = TRUE}, the permutation distribution
#'   (\code{z_perm}) is also returned.
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
    set.seed(seed)
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

#' Power simulation for Solomon designs
#' @param n list with n1..n4 per cell (or a single n per cell)
#' @param delta average treatment effect (on posttest scale)
#' @param rho correlation(pre, post) in pretested cells
#' @param sens pretest sensitization add-on to treatment in pretested cells (0 = none)
#' @param sigma SD of errors
#' @param sims number of Monte Carlo replicates
#' @return data.frame with estimated power for: interaction, ATE, simple effects, and (optionally) Stouffer Z
#' @export
#' @param stouffer Logical; if `TRUE`, also estimate power for the optional
#'   Stouffer meta-analytic procedure.
power_solomon <- function(n = list(n1=50,n2=50,n3=50,n4=50),
                          delta = 0.3, rho = 0.5, sens = 0, sigma = 1, sims = 2000,
                          stouffer = TRUE) {
  if (length(n) == 1) n <- as.list(rep(n, 4)); names(n) <- paste0("n",1:4)
  out <- matrix(0, nrow = sims, ncol = 5)
  colnames(out) <- c("A_interaction","A_ATE","A_pre","A_unpre","A_stouffer")

  for (s in seq_len(sims)) {
    # simulate
    pre1 <- stats::rnorm(n$n1); pre2 <- stats::rnorm(n$n2)
    # induce pre-post corr via bivariate normal construction
    e1 <- stats::rnorm(n$n1); e2 <- stats::rnorm(n$n2)
    post1 <- delta + sens + rho*pre1 + sqrt(1-rho^2)*e1
    post2 <- 0 + 0     + rho*pre2 + sqrt(1-rho^2)*e2
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

    # (iv) optional Stouffer
    if (stouffer) {
      # ANCOVA on 1&2; t on 3&4
      anc <- stats::lm(y_post[preind==1] ~ treat[preind==1] + y_pre[preind==1])
      pA  <- broom::tidy(anc)$p.value[2]
      tR  <- stats::t.test(y_post[preind==0] ~ treat[preind==0])
      pT  <- tR$p.value
      zM  <- stouffer_solomon(c(pA/2, pT/2))$p_meta_one_tailed*2  # report two-tailed-ish
      out[s,"A_stouffer"] <- as.numeric(zM < .05)
    }
  }

  data.frame(
    metric = colnames(out),
    power  = colMeans(out)
  )
}
