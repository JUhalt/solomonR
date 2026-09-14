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
#'
#' @return An object of class \code{solomon_ml}.
#'
#' @param conf_level Confidence level for Wald intervals. Default is 0.95.
#'
#' @section Standard errors and intervals:
#' Standard errors come from the observed information matrix. Tests and Wald
#' confidence intervals use a normal reference distribution, the usual
#' large-sample basis for maximum-likelihood inference. The point estimates
#' coincide with separate regressions in the pretested and unpretested
#' groups. In small samples the resulting intervals can be too narrow;
#' simulation validation of interval coverage is in progress.
#'
#' @references
#' van Engelenburg, G. (1999). Statistical analysis for the Solomon
#' four-group design. University of Twente Research Report 99-06.
#'
#' @export
fit_solomon_ml <- function(
    y_post,
    treat,
    pretested,
    y_pre,
    weights = c("equal"),
    control = list(),
    conf_level = 0.95
) {

  weights <- match.arg(weights)
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

  # Center X over all participants for whom the pretest was administered,
  # following van Engelenburg's deviation-score formulation.
  x_bar <- mean(df$y_pre[idx_pre])

  df$x_c <- NA_real_
  df$x_c[idx_pre] <- df$y_pre[idx_pre] - x_bar

  # ----------------------------------------------------------
  # Starting values
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

  se <- sqrt(diag(V))
  z <- b / se
  p <- 2 * stats::pnorm(-abs(z))

  coef_ci <- .wald_ci(unname(b), unname(se), Inf, conf_level)

  coefficients <- data.frame(
    term = names(b),
    estimate = unname(b),
    std.error = unname(se),
    statistic = unname(z),
    p.value = unname(p),
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

    variance <- as.numeric(
      t(L) %*% V %*% L
    )

    std.error <- sqrt(variance)
    statistic <- estimate / std.error
    p.value <- 2 * stats::pnorm(-abs(statistic))

    ci <- .wald_ci(estimate, std.error, Inf, conf_level)

    data.frame(
      contrast = label,
      estimate = estimate,
      std.error = std.error,
      statistic = statistic,
      p.value = p.value,
      conf.low = unname(ci[, "conf.low"]),
      conf.high = unname(ci[, "conf.high"]),
      row.names = NULL
    )
  }

  Z <- function() {
    stats::setNames(
      numeric(length(b)),
      names(b)
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

  structure(
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
      data = df,
      call = match.call(),
      conf_level = conf_level,
      method = "van Engelenburg (1999) full-information ML"
    ),
    class = "solomon_ml"
  )
}
