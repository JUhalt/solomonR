# Smallest cell size at which fit_solomon_ml() Wald inference is not flagged
# in printed output. Set from the package's simulation validation (issues #10
# and #22) by the threshold rule posted on #22 before the results were
# examined.
.solomon_ml_small_cell <- 40L

# Short description of the inference of a fit_solomon_ml() fit, for figure
# captions and equivalence tests, with its sources. Fits saved before the
# `inference` argument existed used Wald inference.
.ml_inference_label <- function(fit) {
  if (identical(fit$inference, "satterthwaite")) {
    "maximum likelihood; Satterthwaite inference (Satterthwaite, 1946; Welch, 1947)"
  } else {
    "maximum likelihood; Wald inference (van Engelenburg, 1999)"
  }
}

# Variance and reference degrees of freedom for a linear combination `L` of
# the ML parameters (a, bX, bT, bP, bTP, in that order), using the stored
# `inference_parts` of a fit. Wald inference uses the inverse observed
# information with a normal reference. Satterthwaite inference adds the
# variances from the separate unpretested and pretested regressions, with
# Welch-Satterthwaite degrees of freedom; these reduce to the residual degrees
# of freedom of one regression when the other contributes nothing.
#
# With `center = TRUE`, as for the pretest effects (issue #104), the
# combination is evaluated at the estimated mean pretest of the pretested
# participants, which moves it by L["bP"] * bX per unit, so it gains
# (L["bP"] bX)^2 Var(mean): with Wald inference the ML variance of the mean,
# and with the small-sample option s^2 / n as a third Welch-Satterthwaite
# component with n - 1 degrees of freedom. Under the model the mean is
# independent of the regression estimates.
.ml_combination <- function(L, parts, center = FALSE) {
  L <- as.numeric(L)
  k <- if (center && !is.null(parts$center)) {
    L[match("bP", rownames(parts$V))] * parts$center$slope
  } else {
    0
  }
  if (identical(parts$inference, "wald")) {
    added <- if (k != 0) k^2 * parts$center$variance else 0
    return(c(variance = as.numeric(t(L) %*% parts$V %*% L) + added, df = Inf))
  }
  c_u <- as.numeric(L %*% parts$M_u)
  c_p <- as.numeric(L %*% parts$M_p)
  v_u <- as.numeric(t(c_u) %*% parts$V_u %*% c_u)
  v_p <- as.numeric(t(c_p) %*% parts$V_p %*% c_p)
  v_x <- if (k != 0) k^2 * parts$center$variance else 0
  df_x <- if (k != 0) parts$center$df else Inf
  c(variance = v_u + v_p + v_x,
    df = .welch_df(c(v_u, v_p, v_x), c(parts$df_u, parts$df_p, df_x)))
}

# Estimate, standard error, reference df, and interval for each row of `L`,
# a matrix whose columns name ML parameters (unnamed parameters count as 0),
# computed exactly as fit_solomon_ml() computes its own contrasts.
.ml_linear_combination <- function(fit, L) {
  if (is.null(fit$inference_parts)) {
    stop("This fit was made by an earlier version of solomonR; refit it with ",
         "fit_solomon_ml().", call. = FALSE)
  }
  b <- stats::setNames(fit$coefficients$estimate, fit$coefficients$term)
  full <- matrix(0, nrow = nrow(L), ncol = length(b), dimnames = list(NULL, names(b)))
  full[, colnames(L)] <- L

  rows <- lapply(seq_len(nrow(full)), function(i) {
    comb <- .ml_combination(full[i, ], fit$inference_parts)
    estimate <- sum(full[i, ] * b)
    std.error <- sqrt(comb[["variance"]])
    ci <- .wald_ci(estimate, std.error, comb[["df"]], fit$conf_level)
    data.frame(
      estimate = estimate,
      std.error = std.error,
      df = comb[["df"]],
      conf.low = unname(ci[, "conf.low"]),
      conf.high = unname(ci[, "conf.high"])
    )
  })
  do.call(rbind, rows)
}


#' Full-information ML analysis for a Solomon Four-Group Design
#'
#' `r lifecycle::badge("stable")`
#' Fits the maximum-likelihood regression model described by van Engelenburg
#' (1999). Pretest information is incorporated for the pretested groups while
#' structurally missing pretests in the unpretested groups are handled through
#' a separate residual variance.
#'
#' The model estimates the treatment effect, pretest effect,
#' Pretest x Treatment interaction, pretest-posttest slope, and separate
#' residual standard deviations for pretested and unpretested participants.
#' The pretest enters as a deviation from its mean among pretested
#' participants (returned as `pretest_mean`), so the pretest effect `bP`
#' compares pretested and unpretested controls at that mean. The effects
#' table reports the four Solomon contrasts and then the pretest effects
#' among controls, among treated participants, and averaged over the two, as
#' [fit_solomon_glm()] does (see its section "The pretest effect").
#'
#' @param y_post Numeric posttest scores.
#' @param treat Treatment indicator coded 0 = control and 1 = treatment.
#'   Designs with several treatments are not supported; see
#'   [fit_solomon_glm()].
#' @param pretested Pretest indicator coded 0 = unpretested and 1 = pretested.
#' @param y_pre Numeric pretest scores. These should be missing by design for
#'   participants assigned to the unpretested groups.
#' @param weights Character. How to define the average treatment effect across
#'   pretest conditions. Currently \code{"equal"} gives equal weight to the
#'   pretested and unpretested treatment effects.
#' @param control Optional list passed to \code{stats::optim()}.
#' @param conf_level Confidence level for intervals. Default is 0.95.
#' @param inference How standard errors, tests, and intervals are computed:
#'   `"satterthwaite"` (default) for t tests with residual or
#'   Welch-Satterthwaite degrees of freedom (Satterthwaite, 1946; Welch,
#'   1947), or `"wald"` for van Engelenburg's (1999) large-sample Wald
#'   inference. The point estimates are the same. See the Inference options
#'   section.
#' @param data Optional data frame. When supplied, the other data arguments
#'   are looked up in it first, as bare column names (`y_post = post`) or as
#'   strings (`y_post = "post"`).
#'
#' @section Inference options:
#' The point estimates are maximum-likelihood estimates, which coincide with
#' separate regressions in the pretested and unpretested groups.
#'
#' - `inference = "satterthwaite"` (default): standard errors use unbiased
#'   residual variances within each pretest condition. Contrasts within one
#'   condition use t tests with that condition's residual degrees of
#'   freedom, and contrasts that combine the conditions (the ATE and
#'   Pretest x Treatment) use Welch-Satterthwaite degrees of freedom
#'   (Satterthwaite, 1946; Welch, 1947).
#' - `inference = "wald"` follows van Engelenburg (1999): standard errors
#'   come from the observed information matrix, and tests and intervals use a
#'   normal reference distribution, the usual large-sample basis for
#'   maximum-likelihood inference.
#'
#' The pretest effects are evaluated at the mean pretest of the pretested
#' participants, an estimate whose sampling variance their standard errors
#' include (see the section "The pretest effect" of [fit_solomon_glm()]).
#' Under the model the mean is independent of the regression estimates, so
#' Wald inference adds `bX`^2 times the maximum-likelihood variance of the
#' mean, and the small-sample option adds `bX`^2 s^2 / n, for n pretested
#' participants whose pretests have variance s^2, as a third component with
#' n - 1 degrees of freedom. The coefficient table gives `bP` with the
#' standard error of van Engelenburg (1999), which treats the mean as fixed.
#' The pretest effects were not part of the validation below.
#'
#' In the package's simulation validation (issues #10 and #22; 84 scenarios
#' with 2,000 replications each, reported in the article "Validating
#' fit_solomon_ml()" on the package website), both options recovered the
#' Solomon contrasts without bias. Satterthwaite inference had mean coverage
#' of nominal 95% intervals of 0.949 to 0.950 and Type I error of 0.050 to
#' 0.053 at every cell size studied, from 6 to 100 participants per cell.
#' Wald intervals were too narrow in small samples: mean coverage was 0.893
#' with 6 participants per cell, 0.920 with 10, 0.936 with 20, 0.941 with
#' 30, and 0.948 with 100, and the Pretest x Treatment test rejected a true
#' null hypothesis in 9.9% of samples with 6 per cell and 5.9% with 30.
#'
#' Printed output gives the maximum-likelihood residual standard deviations
#' of the two pretest conditions (the square roots of SSE / n, stored in
#' `sigma`), which Wald standard errors use. For Satterthwaite inference it
#' also gives the square roots of the unbiased residual variances (SSE /
#' residual df, stored in `sigma_unbiased`), which its standard errors use.
#'
#' Satterthwaite inference became the default because of these results
#' (issue #115). In solomonR 0.8.0 and earlier the default was `"wald"`;
#' supply `inference = "wald"` to reproduce results from those versions.
#'
#' When `inference = "wald"` and the smallest cell has fewer than 40
#' participants, printed output notes that Wald intervals were too narrow at
#' such sizes in the validation. The threshold follows a rule set before the
#' validation results were examined: 40 is the smallest cell size at which
#' Wald inference had mean coverage of at least 0.940 and Type I error of at
#' most 0.060, with equal and unequal residual variances, at that size and
#' every larger size studied. Even at and above it, Wald inference was
#' approximately adequate rather than exact.
#'
#' @return An object of class \code{solomon_ml}, with the coefficients, the
#'   `effects` table (the four Solomon contrasts, then the three pretest
#'   effects, in the columns `contrast`, `estimate`, `std.error`,
#'   `statistic`, `df`, `p.value`, `conf.low`, and `conf.high`; `df` is
#'   `Inf` with `inference = "wald"`), the residual standard deviations,
#'   `pretest_mean`, and the settings used. The `coefficients` table has the
#'   same columns, with `term` for `contrast`.
#'
#' @references
#' Satterthwaite, F. E. (1946). An approximate distribution of estimates of
#' variance components. *Biometrics Bulletin, 2*(6), 110–114.
#' https://doi.org/10.2307/3002019
#'
#' van Engelenburg, G. (1999). *Statistical analysis for the Solomon four-group
#' design* (Research Report 99-06). University of Twente. ERIC.
#' https://eric.ed.gov/?id=ED435692
#'
#' Welch, B. L. (1947). The generalization of "Student's" problem when several
#' different population variances are involved. *Biometrika, 34*(1–2), 28–35.
#' https://doi.org/10.1093/biomet/34.1-2.28
#'
#' @examples
#' # Default: Satterthwaite inference.
#' with(solomon_example, fit_solomon_ml(y_post, treat, pretested, y_pre))
#'
#' # van Engelenburg's (1999) large-sample Wald inference. With 30
#' # participants per cell, the printed output notes that these intervals
#' # were too narrow at such sizes in the package's validation.
#' with(solomon_example, fit_solomon_ml(y_post, treat, pretested, y_pre,
#'                                      inference = "wald"))
#' @export
fit_solomon_ml <- function(
    y_post,
    treat,
    pretested,
    y_pre,
    weights = c("equal"),
    control = list(),
    conf_level = 0.95,
    inference = c("satterthwaite", "wald"),
    data = NULL
) {
  .solomon_data_args(
    data, c("y_post", "treat", "pretested", "y_pre"),
    environment(), parent.frame()
  )
  .stop_ngroup_unsupported(treat, "fit_solomon_ml")

  weights <- match.arg(weights)
  inference <- match.arg(inference)
  .check_conf_level(conf_level)

  if (length(y_post) != length(treat) ||
      length(y_post) != length(pretested) ||
      length(y_post) != length(y_pre)) {
    stop("y_post, treat, pretested, and y_pre must have equal lengths.")
  }

  treat <- as.integer(treat)
  pretested <- as.integer(pretested)

  if (any(!stats::na.omit(treat) %in% c(0L, 1L))) {
    stop("treat must be coded 0/1.")
  }

  if (any(!stats::na.omit(pretested) %in% c(0L, 1L))) {
    stop("pretested must be coded 0/1.")
  }

  df <- data.frame(
    y_post = y_post,
    treat = treat,
    pretested = pretested,
    y_pre = y_pre
  )

  if (any(is.na(df$y_post))) {
    stop(
      "fit_solomon_ml() currently requires observed posttest scores. ",
      "Incidental outcome missingness will be addressed separately."
    )
  }

  if (any(df$pretested == 1L & is.na(df$y_pre))) {
    stop(
      "Pretest scores are unexpectedly missing among pretested participants."
    )
  }

  if (any(df$pretested == 0L & !is.na(df$y_pre))) {
    warning(
      "Observed y_pre values were supplied for unpretested participants; ",
      "these values are ignored by the Solomon ML model."
    )
  }

  idx_pre <- df$pretested == 1L
  idx_un <- df$pretested == 0L

  if (!all(c(0L, 1L) %in% unique(df$treat[idx_pre]))) {
    stop("Both treatment conditions must occur among pretested participants.")
  }

  if (!all(c(0L, 1L) %in% unique(df$treat[idx_un]))) {
    stop("Both treatment conditions must occur among unpretested participants.")
  }

  cell_sizes <- table(
    factor(df$pretested, levels = c(1L, 0L)),
    factor(df$treat, levels = c(1L, 0L))
  )
  min_cell_n <- as.integer(min(cell_sizes))
  small_sample <- min_cell_n < .solomon_ml_small_cell

  # Center X over all participants for whom the pretest was administered,
  # following van Engelenburg's deviation-score formulation.
  x_bar <- mean(df$y_pre[idx_pre])

  df$x_c <- NA_real_
  df$x_c[idx_pre] <- df$y_pre[idx_pre] - x_bar

  # ----------------------------------------------------------
  # Starting values (and the separate regressions used by
  # Satterthwaite inference)
  # ----------------------------------------------------------

  un_fit <- stats::lm(
    y_post ~ treat,
    data = df[idx_un, , drop = FALSE]
  )

  pre_fit <- stats::lm(
    y_post ~ treat + x_c,
    data = df[idx_pre, , drop = FALSE]
  )

  bu <- stats::coef(un_fit)
  bp <- stats::coef(pre_fit)

  a_start <- unname(bu["(Intercept)"])
  bT_start <- unname(bu["treat"])

  a_pre_start <- unname(bp["(Intercept)"])
  bT_pre_start <- unname(bp["treat"])
  bX_start <- unname(bp["x_c"])

  bP_start <- a_pre_start - a_start
  bTP_start <- bT_pre_start - bT_start

  r_un <- stats::residuals(un_fit)
  r_pre <- stats::residuals(pre_fit)

  sigma_R_start <- sqrt(mean(r_un^2))
  sigma_E_start <- sqrt(mean(r_pre^2))

  sigma_R_start <- max(sigma_R_start, .Machine$double.eps^0.25)
  sigma_E_start <- max(sigma_E_start, .Machine$double.eps^0.25)

  start <- c(
    a = a_start,
    bX = bX_start,
    bT = bT_start,
    bP = bP_start,
    bTP = bTP_start,
    log_sigma_R = log(sigma_R_start),
    log_sigma_E = log(sigma_E_start)
  )

  # ----------------------------------------------------------
  # Full-information negative log likelihood
  # ----------------------------------------------------------

  nll <- function(par) {

    a <- par["a"]
    bX <- par["bX"]
    bT <- par["bT"]
    bP <- par["bP"]
    bTP <- par["bTP"]

    sigma_R <- exp(par["log_sigma_R"])
    sigma_E <- exp(par["log_sigma_E"])

    mu <- numeric(nrow(df))

    # Unpretested groups: O and T
    mu[idx_un] <-
      a +
      bT * df$treat[idx_un]

    # Pretested groups: P and TP
    mu[idx_pre] <-
      a +
      bX * df$x_c[idx_pre] +
      bT * df$treat[idx_pre] +
      bP +
      bTP * df$treat[idx_pre]

    ll_un <- stats::dnorm(
      df$y_post[idx_un],
      mean = mu[idx_un],
      sd = sigma_R,
      log = TRUE
    )

    ll_pre <- stats::dnorm(
      df$y_post[idx_pre],
      mean = mu[idx_pre],
      sd = sigma_E,
      log = TRUE
    )

    -sum(ll_un) - sum(ll_pre)
  }

  opt <- stats::optim(
    par = start,
    fn = nll,
    method = "BFGS",
    hessian = TRUE,
    control = control
  )

  if (opt$convergence != 0L) {
    warning(
      "Maximum-likelihood optimization did not report successful convergence. ",
      "optim code = ",
      opt$convergence
    )
  }

  est_full <- opt$par

  # First five parameters are on their natural scales.
  b <- est_full[c("a", "bX", "bT", "bP", "bTP")]

  H <- opt$hessian

  V_full <- try(
    solve(H),
    silent = TRUE
  )

  if (inherits(V_full, "try-error")) {
    stop(
      "The observed information matrix could not be inverted; ",
      "standard errors are unavailable."
    )
  }

  dimnames(V_full) <- list(
    names(est_full),
    names(est_full)
  )

  V <- V_full[
    names(b),
    names(b),
    drop = FALSE
  ]

  Z <- function() {
    stats::setNames(
      numeric(length(b)),
      names(b)
    )
  }

  # ----------------------------------------------------------
  # Variance and degrees of freedom for a linear combination
  # of the ML parameters
  # ----------------------------------------------------------

  # Everything needed for the variance and reference df of any linear
  # combination of the ML parameters. It is kept in the fit so that figures
  # use exactly the fit's inference (see .ml_combination()).
  parts <- list(inference = inference, V = V)

  # The mean pretest at which the pretest effects are evaluated is an
  # estimate (issue #104): its ML variance for Wald inference, and s^2 / n
  # with n - 1 degrees of freedom for the small-sample option.
  n_pre <- sum(idx_pre)
  ss_x <- sum(df$x_c[idx_pre]^2)
  parts$center <- list(
    slope = unname(b["bX"]),
    variance = if (inference == "wald") ss_x / n_pre^2 else ss_x / (n_pre * (n_pre - 1)),
    df = if (inference == "wald") Inf else n_pre - 1
  )

  if (inference == "satterthwaite") {

    V_u <- stats::vcov(un_fit)
    V_p <- stats::vcov(pre_fit)

    # Each ML parameter (rows: a, bX, bT, bP, bTP) as a combination of the
    # coefficients of the separate unpretested and pretested regressions.
    parts$M_u <- matrix(
      c(1, 0,
        0, 0,
        0, 1,
        -1, 0,
        0, -1),
      nrow = 5, byrow = TRUE,
      dimnames = list(names(b), colnames(V_u))
    )

    parts$M_p <- matrix(
      c(0, 0, 0,
        0, 0, 1,
        0, 0, 0,
        1, 0, 0,
        0, 1, 0),
      nrow = 5, byrow = TRUE,
      dimnames = list(names(b), colnames(V_p))
    )

    parts$V_u <- V_u
    parts$V_p <- V_p
    parts$df_u <- stats::df.residual(un_fit)
    parts$df_p <- stats::df.residual(pre_fit)
  }

  combination_variance <- function(L, center = FALSE) {
    .ml_combination(L, parts, center)[["variance"]]
  }

  combination_df <- function(L, center = FALSE) {
    .ml_combination(L, parts, center)[["df"]]
  }

  unit <- function(name) {
    L <- Z()
    L[name] <- 1
    L
  }

  coef_se <- vapply(names(b), function(name) sqrt(combination_variance(unit(name))), numeric(1))
  coef_df <- vapply(names(b), function(name) combination_df(unit(name)), numeric(1))
  coef_statistic <- b / coef_se
  coef_p <- 2 * stats::pt(-abs(coef_statistic), df = coef_df)
  coef_ci <- .wald_ci(unname(b), unname(coef_se), unname(coef_df), conf_level)

  coefficients <- data.frame(
    term = names(b),
    estimate = unname(b),
    std.error = unname(coef_se),
    statistic = unname(coef_statistic),
    df = unname(coef_df),
    p.value = unname(coef_p),
    conf.low = unname(coef_ci[, "conf.low"]),
    conf.high = unname(coef_ci[, "conf.high"]),
    row.names = NULL
  )

  # ----------------------------------------------------------
  # Solomon estimands
  # ----------------------------------------------------------

  # The pretest effects are evaluated at the estimated mean pretest, whose
  # variance they also carry (.ml_combination()).
  contrast <- function(L, label) {

    L <- L[names(b)]
    center <- label %in% .solomon_pretest_order

    estimate <- sum(L * b)
    std.error <- sqrt(combination_variance(L, center))
    df <- combination_df(L, center)
    statistic <- estimate / std.error
    p.value <- 2 * stats::pt(-abs(statistic), df = df)
    ci <- .wald_ci(estimate, std.error, df, conf_level)

    data.frame(
      contrast = label,
      estimate = estimate,
      std.error = std.error,
      statistic = statistic,
      df = df,
      p.value = p.value,
      conf.low = unname(ci[, "conf.low"]),
      conf.high = unname(ci[, "conf.high"]),
      row.names = NULL
    )
  }

  # Treatment among unpretested participants
  L_un <- Z()
  L_un["bT"] <- 1

  # Treatment among pretested participants
  L_pre <- Z()
  L_pre["bT"] <- 1
  L_pre["bTP"] <- 1

  # Sensitization
  L_int <- Z()
  L_int["bTP"] <- 1

  # Equal-weight average treatment effect
  L_ate <- Z()
  L_ate["bT"] <- 1
  L_ate["bTP"] <- 0.5

  # The pretest (testing) effects (issue #104): pretested minus unpretested
  # participants with the centered pretest at zero, that is, at the pretested
  # participants' mean pretest. Among controls bP; among treated
  # participants bP + bTP; and their average.
  L_pc <- Z()
  L_pc["bP"] <- 1
  L_pt <- L_pc
  L_pt["bTP"] <- 1
  L_pm <- L_pc
  L_pm["bTP"] <- 0.5

  effects <- .effects_table(rbind(
    contrast(
      L_ate,
      "ATE (avg over pretest)"
    ),
    contrast(
      L_int,
      "Pretest x Treatment"
    ),
    contrast(
      L_pre,
      "Treatment | pretested"
    ),
    contrast(
      L_un,
      "Treatment | unpretested"
    ),
    contrast(
      L_pc,
      "Pretest effect | control"
    ),
    contrast(
      L_pt,
      "Pretest effect | treated"
    ),
    contrast(
      L_pm,
      "Pretest main effect"
    )
  ))

  sigma_R <- exp(est_full["log_sigma_R"])
  sigma_E <- exp(est_full["log_sigma_E"])

  structure(
    list(
      coefficients = coefficients,
      effects = effects,
      sigma = c(
        unpretested = unname(sigma_R),
        pretested = unname(sigma_E)
      ),
      # Square roots of the unbiased residual variances (SSE / residual df),
      # which Satterthwaite standard errors use.
      sigma_unbiased = c(
        unpretested = stats::sigma(un_fit),
        pretested = stats::sigma(pre_fit)
      ),
      pretest_mean = x_bar,
      logLik = -opt$value,
      convergence = opt$convergence,
      optimizer = opt,
      vcov = V,
      inference_parts = parts,
      data = df,
      call = match.call(),
      conf_level = conf_level,
      inference = inference,
      min_cell_n = min_cell_n,
      small_sample = small_sample,
      method = "van Engelenburg (1999) full-information ML"
    ),
    class = "solomon_ml"
  )
}
