# Rules of perm_solomon() for supported models, degenerate labelings, and
# the statistic at the exact estimate (issue #134). The conventions are
# those of the reference set posted on #134 (Section 5, decisions D1-D16).

# Classed errors of perm_solomon(). Every one also has the class
# solomonR_perm_error, so one handler can catch them all.
.perm_stop <- function(class, ...) {
  stop(structure(
    class = c(class, "solomonR_perm_error", "error", "condition"),
    list(message = paste0(...), call = NULL)
  ))
}

# The four models the package documents and validates (decision D1): the
# Gaussian (identity), binomial (logit), and Poisson (log) models of glm(),
# and the NB2 model of family = "negative_binomial" (MASS::glm.nb(), theta
# by maximum likelihood). Checked before the participant and cluster paths
# split, so that both refuse the same fits.
.perm_check_family <- function(fit) {
  fam <- stats::family(fit$model)
  ok <- (identical(fam$family, "gaussian") && identical(fam$link, "identity")) ||
    (identical(fam$family, "binomial") && identical(fam$link, "logit")) ||
    (identical(fam$family, "poisson") && identical(fam$link, "log")) ||
    (inherits(fit$model, "negbin") && identical(fam$link, "log"))
  if (ok) return(invisible(NULL))
  hint <- if (fam$family %in% c("quasipoisson", "quasibinomial")) {
    sprintf(" Refit with family = %s(), which gives the same coefficients and HC3 standard errors.",
            sub("^quasi", "", fam$family))
  } else if (.is_negbin(fam)) {
    " For a negative binomial model, fit with family = \"negative_binomial\" (theta estimated by MASS::glm.nb())."
  } else {
    ""
  }
  .perm_stop(
    "solomonR_perm_family_error",
    "perm_solomon() supports four models: Gaussian (identity link), binomial (logit link), ",
    "Poisson (log link), and family = \"negative_binomial\". This fit uses ", fam$family,
    " with the ", fam$link, " link, for which the permutation p-value has not been checked ",
    "(issue #134).", hint
  )
}

# The contrast's weights on the coefficients named `coef_names`, which are
# those of the treat * pretested coding.
.perm_contrast_vector <- function(coef_names, contrast) {
  L <- stats::setNames(numeric(length(coef_names)), coef_names)
  if (!"treat" %in% coef_names) stop("The fitted model does not contain a treatment coefficient.")
  has_interaction <- "treat:pretested" %in% coef_names
  if (contrast == "Pretest x Treatment" && !has_interaction) {
    stop("The fitted model does not contain a treatment-by-pretest interaction.")
  }
  w <- switch(contrast,
              "ATE (avg over pretest)" = c(1, 0.5),
              "Pretest x Treatment" = c(0, 1),
              "Treatment | pretested" = c(1, 1),
              "Treatment | unpretested" = c(1, 0))
  L["treat"] <- w[1]
  if (has_interaction) L["treat:pretested"] <- w[2]
  L
}

# The outcome as glm() fits it: a factor counts its first level as 0 and the
# others as 1, and a logical outcome is 0 or 1 (binomial()$initialize), so
# that a binary outcome meets the same rules however it is coded.
.perm_outcome <- function(y) {
  if (is.factor(y)) return(as.numeric(y != levels(y)[1L]))
  as.numeric(y)
}

# No variation in a pretest condition the contrast uses (the acceptance
# criterion of #134 and decision D5). `y` (as .perm_outcome() gives it),
# `pretested`, `shared` (the model-matrix columns of the covariates the two
# conditions share, NULL or none without them), and `offset` (NULL without
# an exposure) are those of the rows in the model.
#
# D5: an error when the contrast is 0 or undefined under every labeling
# (ERR). With the outcome constant in a condition the contrast uses, that
# holds when
# - the constant is at a bound (0 or 1 for a binary outcome, 0 for a count)
#   in a condition the contrast uses: the contrast is undefined;
# - every condition the contrast uses is fitted exactly whatever the labels,
#   with a treatment effect of 0: its outcome and offset are constant, and no
#   shared covariate links it to the other condition (there is none, or each
#   is constant within the condition, where its own intercept absorbs it);
# - the outcome and the offset are constant in both conditions: the model
#   without treatment fits every row, so every treatment coefficient is 0.
# Otherwise the contrast still varies, through the other condition, a shared
# covariate, or the exposure (ERR-S): the p-value, with a classed warning
# (D5's proposal).
.perm_check_variation <- function(y, pretested, contrast, family, shared = NULL, offset = NULL) {
  uses <- .perm_contrast_weights(contrast) != 0
  rows <- list(pretested = pretested == 1, unpretested = pretested == 0)
  one_value <- function(v) length(unique(v)) == 1L
  constant <- vapply(rows, function(r) one_value(y[r]), logical(1))
  if (!any(constant[uses])) return(invisible(NULL))
  flat <- constant & vapply(rows, function(r) is.null(offset) || one_value(offset[r]), logical(1))
  apart <- vapply(rows, function(r) {
    is.null(shared) || ncol(shared) == 0L ||
      all(vapply(seq_len(ncol(shared)), function(j) one_value(shared[r, j]), logical(1)))
  }, logical(1))
  at_bound <- constant & vapply(rows, function(r) {
    v <- y[r][1L]
    switch(family$family, gaussian = FALSE, binomial = v %in% c(0, 1), isTRUE(v == 0))
  }, logical(1))
  and_list <- function(x) {
    if (length(x) <= 2L) return(paste(x, collapse = " and "))
    paste0(paste(x[-length(x)], collapse = ", "), ", and ", x[length(x)])
  }
  named <- function(keep) {
    x <- names(rows)[keep & uses]
    paste0(and_list(x), " condition", if (length(x) > 1L) "s" else "")
  }
  undefined <- any(at_bound & uses)
  if (undefined || all((flat & apart)[uses]) || all(flat)) {
    .perm_stop(
      "solomonR_perm_no_variation_error",
      if (all(constant) && !undefined) {
        paste0("The outcome does not vary in either pretest condition, so the contrast \"",
               contrast, "\" is 0 under every labeling and the test has no information.")
      } else {
        paste0("The outcome does not vary in the ", named(if (undefined) at_bound else constant),
               ", which the contrast \"", contrast, "\" uses, so the contrast is ",
               if (undefined) "undefined" else "0",
               " under every labeling and the test has no information.",
               if (!all(constant)) " Choose a contrast on a condition whose outcome varies.")
      }
    )
  }
  # The contrast still varies (decision D5's proposal): a p-value, with a
  # classed warning that names the constant condition and what the contrast
  # varies through.
  through <- c(
    if (any(!constant[uses])) "the other condition",
    if (any((constant & !apart)[uses])) "the covariates the conditions share",
    if (any((constant & !flat)[uses])) "the exposure"
  )
  warning(structure(
    class = c("solomonR_perm_no_variation_warning", "warning", "condition"),
    list(message = paste0(
      "The outcome does not vary in the ", named(constant), "; the contrast \"", contrast,
      "\" varies only through ", and_list(through), "."), call = NULL)
  ))
  invisible(NULL)
}

# The refit frame (issue #134, items 6 and 9). The outcome of an
# identity-link fit, and each column of a numeric matrix covariate, are
# moved by one of their own values (the lower median of the rows in the
# model). That changes only the intercept, and it is exact when the values
# lie within a factor of two of it (Sterbenz), so scores far from zero give
# the same refits as scores near it.
.perm_center_frame <- function(df, used, family) {
  lower_median <- function(x) {
    x <- sort(x[is.finite(x)])
    if (length(x)) x[ceiling(length(x) / 2)] else 0
  }
  if (identical(family$family, "gaussian") && identical(family$link, "identity")) {
    df$y <- df$y - lower_median(df$y[used])
  }
  for (name in setdiff(names(df), c("y", "treat", "pretested"))) {
    x <- df[[name]]
    if (is.matrix(x) && is.numeric(x)) {
      for (j in seq_len(ncol(x))) x[, j] <- x[, j] - lower_median(x[used, j])
      df[[name]] <- x
    }
  }
  df
}

# Whether the contrast L is estimable from the rows of a design whose
# pivoted QR decomposition is `q` (L orthogonal to the null space).
.perm_estimable <- function(L, q, tol = 1e-7) {
  p <- length(L)
  r <- q$rank
  if (r == p) return(TRUE)
  if (r == 0) return(FALSE)
  piv1 <- q$pivot[seq_len(r)]
  piv2 <- q$pivot[-seq_len(r)]
  R <- qr.R(q)
  t1 <- backsolve(R[seq_len(r), seq_len(r), drop = FALSE], R[seq_len(r), r + seq_along(piv2), drop = FALSE])
  gap <- L[piv2] - drop(crossprod(t1, L[piv1]))
  max(abs(gap)) / max(1, abs(L)) <= tol
}

# The value of a statistic whose HC3 variance `v` may be zero in exact
# arithmetic: constant arms or an exact count fit give an infinite statistic
# when the contrast is not 0 (counted, the acceptance criterion), and 0 when
# it is (0/0, decision D11). `scale` is the model-based variance of the
# contrast, against which v is judged zero.
.perm_finish <- function(est, v, scale) {
  if (!is.finite(scale) || scale <= 0) scale <- 1
  if (!is.finite(v) || v <= 1e-20 * scale) {
    return(if (abs(est) / sqrt(scale) <= 1e-10) 0 else Inf)
  }
  est / sqrt(v)
}

# Newton's method for a model fitted by iteration, from glm()'s stopping
# point (decision D8): the binomial and Poisson models have canonical links,
# so this is IRLS; for the negative binomial with theta fixed it uses the
# observed information, (y + theta) theta mu / (mu + theta)^2, so that it
# converges quadratically too. Returns the coefficients and the change in
# the linear predictor at the last step, which marks rows whose fitted value
# runs to the boundary (they move by about 1 per step; others by rounding).
.perm_newton <- function(X, y, off, pw, fam, beta, theta, steps = 6L) {
  deta <- rep(0, nrow(X))
  for (s in seq_len(steps)) {
    eta <- drop(X %*% beta) + off
    mu <- fam$linkinv(eta)
    me <- fam$mu.eta(eta)
    va <- fam$variance(mu)
    g <- pw * (y - mu) * me / va
    w <- if (is.null(theta)) pw * me^2 / va else pw * (y + theta) * theta * mu / (mu + theta)^2
    w <- pmax(w, 1e-300)
    sw <- sqrt(w)
    step <- qr.coef(qr(sw * X, tol = 1e-11), g / sw)
    step[is.na(step)] <- 0
    beta <- beta + step
    deta <- drop(X %*% step)
    # From glm()'s stopping point one step reaches the estimate to rounding
    # error; a second confirms it. Rows still moving run to the boundary.
    if (s >= 2L && max(abs(deta)) < 1e-6) break
  }
  list(beta = beta, deta = deta)
}

# The statistic of a fit by iteration at coefficients `beta`: the contrast,
# or the contrast divided by its HC3 standard error (MacKinnon & White,
# 1985) with the Fisher working weights, as sandwich computes it for a glm
# (Zeileis, 2006), here at the estimate itself.
.perm_glm_stat <- function(X, y, off, pw, fam, beta, L, statistic) {
  est <- sum(L * beta)
  if (statistic == "difference") return(est)
  eta <- drop(X %*% beta) + off
  mu <- fam$linkinv(eta)
  me <- fam$mu.eta(eta)
  w <- pw * me^2 / fam$variance(mu)
  res <- (y - mu) / me
  q <- qr(sqrt(w) * X, tol = 1e-11)
  r <- seq_len(q$rank)
  Q <- qr.Q(q)[, r, drop = FALSE]
  hat <- rowSums(Q^2)
  B <- backsolve(qr.R(q)[r, r, drop = FALSE], L[q$pivot[r]], transpose = TRUE)
  infl <- drop(Q %*% B) * sqrt(w) * res
  term <- ifelse(hat > 1 - 1e-8, 0, (infl / (1 - hat))^2)
  .perm_finish(est, sum(term), sum(B^2))
}

# One labeling's statistic, with its code. The codes follow the reference
# set: "ALIAS" (a covariate collinear with the cells leaves the contrast not
# estimable; left out, decision D16), "H1" (studentized: a row with leverage
# 1 whose HC3 term enters the contrast's variance; an error, decision D3),
# "INF" (an emptied arm, separation, or constant arms: +Inf, counted;
# decisions D4 and D10), "LIM" (a degenerate cell the contrast does not use:
# the value of the limiting model, decision D9), "FAIL" (glm() failed), and
# "ok". `nonconv` marks a fit that did not converge. `X` is the design of
# the refit (the model matrix with the labeling's treatment columns), `y` the
# outcome, `off` the offset, `pw` the prior weights, `L` the contrast.
.perm_statistic <- function(X, y, off, pw, fam, theta, L, statistic) {
  out <- function(z, code, nonconv = FALSE) list(z = as.numeric(z), code = code, nonconv = nonconv)
  gaussian <- identical(fam$family, "gaussian")
  # The models fitted by iteration go through glm.fit(), the fitter of glm(),
  # on the design built once; the linear model needs only the QR below.
  fit <- NULL
  if (!gaussian) {
    fit <- tryCatch(withCallingHandlers(
      stats::glm.fit(X, y, weights = pw, offset = off, family = fam,
                     control = stats::glm.control(maxit = 200)),
      warning = function(w) if (startsWith(conditionMessage(w), "glm.fit:")) invokeRestart("muffleWarning")
    ), error = function(e) NULL)
    if (is.null(fit)) return(out(NA_real_, "FAIL"))
  }
  qx <- qr(X, tol = 1e-7)
  if (qx$rank < ncol(X) && !.perm_estimable(L, qx)) return(out(NA_real_, "ALIAS"))
  r <- seq_len(qx$rank)
  Q <- NULL
  if (statistic == "studentized") {
    Q <- qr.Q(qx)[, r, drop = FALSE]
    h <- rowSums(Q^2)
    lev1 <- which(h > 1 - 1e-8)
    if (length(lev1)) {
      a <- drop(Q %*% backsolve(qr.R(qx)[r, r, drop = FALSE], L[qx$pivot[r]], transpose = TRUE))
      if (any(abs(a[lev1]) > 1e-8 * max(abs(a)))) return(out(NA_real_, "H1"))
    }
  }
  if (gaussian) {
    # Least squares by the QR decomposition glm() uses (LINPACK, tolerance
    # 1e-11); with full rank it is the decomposition above.
    q <- if (qx$rank == ncol(X)) qx else qr(X, tol = 1e-11)
    b <- qr.coef(q, y - off)
    b[is.na(b)] <- 0
    est <- sum(L * b)
    if (statistic == "difference") return(out(est, "ok"))
    rq <- seq_len(q$rank)
    Qq <- if (identical(q, qx)) Q else qr.Q(q)[, rq, drop = FALSE]
    hat <- rowSums(Qq^2)
    B <- backsolve(qr.R(q)[rq, rq, drop = FALSE], L[q$pivot[rq]], transpose = TRUE)
    e <- y - off - drop(X %*% b)
    infl <- drop(Qq %*% B) * e
    term <- ifelse(hat > 1 - 1e-8, 0, (infl / (1 - hat))^2)
    scale <- sum(backsolve(qr.R(qx)[r, r, drop = FALSE], L[qx$pivot[r]], transpose = TRUE)^2) *
      stats::var(y)
    z <- .perm_finish(est, sum(term), scale)
    return(out(z, if (is.infinite(z)) "INF" else "ok"))
  }
  b <- fit$coefficients
  keep <- !is.na(b)
  # Rows whose fitted values run to the boundary. An arm whose outcomes are
  # all 0 (or all 1 for a binary outcome) is found from its counts, exactly
  # (Albert & Anderson, 1984); separation by a covariate is found from rows
  # that keep moving under Newton's method once the others have converged.
  # The limiting model is refitted without them (Geyer, 2009): a contrast it
  # does not estimate diverges or has no limit, and is +Inf (decisions D4 and
  # D10); one it estimates takes its value (D9).
  cell <- X[, "treat"] + 2 * X[, "pretested"]
  K <- rep(TRUE, length(y))
  for (k in unique(cell)) {
    in_k <- cell == k
    if (all(y[in_k] <= 0) || (identical(fam$family, "binomial") && all(y[in_k] >= 1))) K[in_k] <- FALSE
  }
  beta <- b[keep]
  cols <- which(keep)
  for (round in 1:5) {
    if (!all(K)) {
      qK <- qr(X[K, , drop = FALSE], tol = 1e-7)
      if (!.perm_estimable(L, qK)) return(out(Inf, "INF"))
      cols <- sort(qK$pivot[seq_len(qK$rank)])
      lim <- tryCatch(suppressWarnings(stats::glm.fit(X[K, cols, drop = FALSE], y[K], weights = pw[K],
                                                      offset = off[K], family = fam,
                                                      control = stats::glm.control(maxit = 200))),
                      error = function(e) NULL)
      if (is.null(lim)) return(out(NA_real_, "FAIL"))
      beta <- lim$coefficients
      beta[is.na(beta)] <- 0
    }
    XK <- X[K, cols, drop = FALSE]
    nt <- .perm_newton(XK, y[K], off[K], pw[K], fam, beta, theta)
    moving <- abs(nt$deta) > 1e-3
    if (!any(moving)) break
    yK <- y[K]
    side <- if (identical(fam$family, "binomial")) (yK <= 0 & nt$deta < 0) | (yK >= 1 & nt$deta > 0) else
      (yK <= 0 & nt$deta < 0)
    if (!any(moving & side)) break
    K[which(K)[moving & side]] <- FALSE
  }
  z <- .perm_glm_stat(XK, y[K], off[K], pw[K], fam, nt$beta, L[cols], statistic)
  code <- if (is.infinite(z)) "INF" else if (all(K)) "ok" else "LIM"
  out(z, code, nonconv = any(moving) || (all(K) && !fit$converged))
}
