# Smallest cell size at which fit_solomon_ml() Wald inference is not flagged.
# Set from the package's simulation validation (issues #10 and #22) by the
# threshold rule posted on #22 before the results were examined.
.solomon_ml_small_cell <- 40L

.warn_ml_small_sample <- function(min_cell) {
  warning(structure(
    class = c("solomonR_small_sample_warning", "warning", "condition"),
    list(
      message = paste0(
        "The smallest Solomon cell has ", min_cell, " participants. In the ",
        "package's simulation validation, maximum-likelihood Wald intervals were ",
        "too narrow with fewer than ", .solomon_ml_small_cell, " participants per ",
        "cell. Consider inference = \"satterthwaite\", or supply ",
        "inference = \"wald\" to keep the default without this warning. ",
        "See ?fit_solomon_ml."
      ),
      call = NULL
    )
  ))
}

# Variance and reference degrees of freedom for a linear combination `L` of
# the ML parameters (a, bX, bT, bP, bTP, in that order), using the stored
# `inference_parts` of a fit. Wald inference uses the inverse observed
# information with a normal reference. The small-sample option adds the
# variances from the separate unpretested and pretested regressions, with
# Welch-Satterthwaite degrees of freedom; these reduce to the residual degrees
# of freedom of one regression when the other contributes nothing.
.ml_combination <- function(L, parts) {
  L <- as.numeric(L)
  if (identical(parts$inference, "wald")) {
    return(c(variance = as.numeric(t(L) %*% parts$V %*% L), df = Inf))
  }
  c_u <- as.numeric(L %*% parts$M_u)
  c_p <- as.numeric(L %*% parts$M_p)
  v_u <- as.numeric(t(c_u) %*% parts$V_u %*% c_u)
  v_p <- as.numeric(t(c_p) %*% parts$V_p %*% c_p)
  denominator <- sum(
    if (v_u > 0) v_u^2 / parts$df_u else 0,
    if (v_p > 0) v_p^2 / parts$df_p else 0
  )
  c(variance = v_u + v_p, df = (v_u + v_p)^2 / denominator)
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
#' Fits the maximum-likelihood regression model described by van Engelenburg
#' (1999). Pretest information is incorporated for the pretested groups while
#' structurally missing pretests in the unpretested groups are handled through
#' a separate residual variance.
#'
#' The model estimates the treatment effect, pretest effect,
#' Treatment x Pretest interaction, pretest-posttest slope, and separate
#' residual standard deviations for pretested and unpretested participants.
#'
#' @param y_post Numeric posttest scores.
#' @param treat Treatment indicator coded 0 = control and 1 = treatment.
#' @param pretested Pretest indicator coded 0 = unpretested and 1 = pretested.
#' @param y_pre Numeric pretest scores. These should be missing by design for
#'   participants assigned to the unpretested groups.
#' @param weights Character. How to define the average treatment effect across
#'   pretest conditions. Currently \code{"equal"} gives equal weight to the
#'   pretested and unpretested treatment effects.
#' @param control Optional list passed to \code{stats::optim()}.
#' @param conf_level Confidence level for intervals. Default is 0.95.
#' @param inference How standard errors, tests, and intervals are computed:
#'   `"wald"` (default) for van Engelenburg's (1999) large-sample Wald
#'   inference, or `"satterthwaite"` for the small-sample option. The point
#'   estimates are the same. See the Inference options section.
#'
#' @section Inference options:
#' The point estimates are maximum-likelihood estimates, which coincide with
#' separate regressions in the pretested and unpretested groups.
#'
#' - `inference = "wald"` (default) follows van Engelenburg (1999): standard
#'   errors come from the observed information matrix, and tests and
#'   intervals use a normal reference distribution, the usual large-sample
#'   basis for maximum-likelihood inference.
#' - `inference = "satterthwaite"` is a small-sample option. Standard errors
#'   use unbiased residual variances within each pretest condition. Contrasts
#'   within one condition use t tests with that condition's residual degrees
#'   of freedom, and contrasts that combine the conditions (the ATE and
#'   Pretest x Treatment) use Welch-Satterthwaite degrees of freedom
#'   (Satterthwaite, 1946; Welch, 1947).
#'
#' In the package's simulation validation (issues #10 and #22; 84 scenarios
#' with 2,000 replications each, reported in the article "Validating
#' fit_solomon_ml()" on the package website), both options recovered the
#' Solomon contrasts without bias. Wald intervals were too narrow in small
#' samples: mean coverage of nominal 95% intervals was 0.893 with 6
#' participants per cell, 0.920 with 10, 0.936 with 20, 0.941 with 30, and
#' 0.948 with 100, and the Pretest x Treatment test rejected a true null
#' hypothesis in 9.9% of samples with 6 per cell and 5.9% with 30. The
#' small-sample option had mean coverage of 0.949 to 0.950 and Type I error of
#' 0.050 to 0.053 at every cell size studied, from 6 to 100 per cell.
#'
#' When the smallest cell has fewer than 40 participants and `inference` is
#' not supplied, `fit_solomon_ml()` issues a warning of class
#' `solomonR_small_sample_warning` that suggests the small-sample option.
#' Supplying `inference = "wald"` explicitly keeps the default without the
#' warning. The threshold follows a rule set before the validation results
#' were examined: 40 is the smallest cell size at which Wald inference had
#' mean coverage of at least 0.940 and Type I error of at most 0.060, with
#' equal and unequal residual variances, at that size and every larger size
#' studied.
#'
#' @return An object of class \code{solomon_ml}.
#'
#' @references
#' Satterthwaite, F. E. (1946). An approximate distribution of estimates of
#' variance components. *Biometrics Bulletin, 2*(6), 110–114.
#' https://doi.org/10.2307/3002019
#'
#' van Engelenburg, G. (1999). *Statistical analysis for the Solomon four-group
#' design* (Research Report 99-06). University of Twente.
#'
#' Welch, B. L. (1947). The generalization of "Student's" problem when several
#' different population variances are involved. *Biometrika, 34*(1–2), 28–35.
#' https://doi.org/10.1093/biomet/34.1-2.28
#'
#' @export
fit_solomon_ml <- function(
    y_post,
    treat,
    pretested,
    y_pre,
    weights = c("equal"),
    control = list(),
    conf_level = 0.95,
    inference = c("wald", "satterthwaite")
) {

  inference_supplied <- !missing(inference)
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
  # Starting values (and the separate regressions used by the
  # small-sample option)
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

  combination_variance <- function(L) {
    .ml_combination(L, parts)[["variance"]]
  }

  combination_df <- function(L) {
    .ml_combination(L, parts)[["df"]]
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
    p.value = unname(coef_p),
    df = unname(coef_df),
    conf.low = unname(coef_ci[, "conf.low"]),
    conf.high = unname(coef_ci[, "conf.high"]),
    row.names = NULL
  )

  # ----------------------------------------------------------
  # Solomon estimands
  # ----------------------------------------------------------

  contrast <- function(L, label) {

    L <- L[names(b)]

    estimate <- sum(L * b)
    std.error <- sqrt(combination_variance(L))
    df <- combination_df(L)
    statistic <- estimate / std.error
    p.value <- 2 * stats::pt(-abs(statistic), df = df)
    ci <- .wald_ci(estimate, std.error, df, conf_level)

    data.frame(
      contrast = label,
      estimate = estimate,
      std.error = std.error,
      statistic = statistic,
      p.value = p.value,
      df = df,
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

  effects <- rbind(
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
    )
  )

  sigma_R <- exp(est_full["log_sigma_R"])
  sigma_E <- exp(est_full["log_sigma_E"])

  out <- structure(
    list(
      coefficients = coefficients,
      effects = effects,
      sigma = c(
        unpretested = unname(sigma_R),
        pretested = unname(sigma_E)
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

  if (inference == "wald" && small_sample && !inference_supplied) {
    .warn_ml_small_sample(min_cell_n)
  }

  out
}
