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

# The participants marginal_solomon() standardizes over with
# `pretest_center = "observed"`, with their design rows (issue #135). The
# unpretested cells are standardized over the unpretested participants in
# the model, as with the default. The pretested cells are standardized over
# every pretested participant with a pretest and the model's other
# covariates, in the model or not: with or without an outcome, and with or
# without a recorded condition, because each participant's risk is
# predicted under treatment and under control. Their mean pretest is the
# "observed" center of fit_solomon_glm() when every pretested participant
# with a pretest has the model's covariates; one without a covariate is in
# that center but not in this set, whose risks need the covariates. (The
# center itself plays no part in the standardization.) Without a pretest
# covariate, and
# when every such participant is in the model, the rows are the model's,
# the participants the default standardizes over.
#
# Returns `X`, one row per participant in the order of the fit's data;
# `pretested` and `treat` (NA where the condition was not recorded);
# `in_model`, which flags the rows the model used; `index`, their rows in
# the fit's data; and `n`, the number of pretested participants the
# default standardizes over (`analyzed`, those in the model), the number
# this set standardizes over (`observed`), and the number of those outside
# the model who have no outcome (`no_outcome`).
.marginal_rows <- function(fit) {
  model <- fit$model
  X_model <- stats::model.matrix(model)
  d <- fit$data
  used <- .glm_rows_used(fit)
  extra <- rep(FALSE, nrow(d))
  if ("pre_obs" %in% colnames(X_model)) {
    extra <- .in_pretest_center(d) & !used
    covariates <- setdiff(names(d), c("y", "treat", "pretested", "pre_obs", "log_exposure"))
    if (length(covariates)) {
      extra <- extra & stats::complete.cases(d[, covariates, drop = FALSE])
    }
  }
  analyzed <- sum(X_model[, "pretested"] == 1)
  n <- c(analyzed = analyzed, observed = analyzed + sum(extra),
         no_outcome = sum(extra & is.na(d$y)))
  storage.mode(n) <- "integer"
  if (!any(extra)) {
    return(list(
      X = X_model, pretested = as.integer(X_model[, "pretested"]),
      treat = as.integer(X_model[, "treat"]), in_model = rep(TRUE, nrow(X_model)),
      index = which(used), n = n
    ))
  }

  # The design rows of the participants outside the model. The treatment
  # terms are removed from the linear predictor before a cell is predicted
  # (.marginal_risks()), so any value serves where the condition was not
  # recorded; the offset does not enter the design matrix.
  frame <- d[extra, , drop = FALSE]
  frame$treat[is.na(frame$treat)] <- 0L
  if ("log_exposure" %in% names(frame)) frame$log_exposure[is.na(frame$log_exposure)] <- 0
  tt <- stats::delete.response(stats::terms(model))
  mf <- stats::model.frame(tt, frame, na.action = stats::na.pass, xlev = model$xlevels)
  X_extra <- stats::model.matrix(tt, mf, contrasts.arg = model$contrasts)

  rows <- used | extra
  in_model <- used[rows]
  X <- matrix(NA_real_, sum(rows), ncol(X_model), dimnames = list(NULL, colnames(X_model)))
  X[in_model, ] <- X_model
  # The columns are those of the model, in its order, and are taken by
  # position: two can share a name (.scaled_columns()).
  X[!in_model, ] <- if (identical(colnames(X_extra), colnames(X_model))) X_extra else
    X_extra[, colnames(X_model), drop = FALSE]
  list(
    X = X, pretested = as.integer(d$pretested[rows]), treat = as.integer(d$treat[rows]),
    in_model = in_model, index = which(rows), n = n
  )
}

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

# The four Solomon contrasts and the three pretest effects (issue #104) on
# one scale from the four marginal risks. The average treatment effect
# compares risks averaged over equal numbers of pretested and unpretested
# participants, and the pretest main effect risks averaged over equal
# numbers of treated and control participants. Ratio scales are undefined
# when a risk is 0 or 1 (within 1e-8, as under separation), or when a count
# outcome's rate is 0.
# (A function, because R/solomon_plot_effects.R, which defines the names,
# is collated after this file.)
.marginal_contrast_names <- function() c(.solomon_contrast_order, .solomon_pretest_order)

.marginal_contrasts <- function(r, scale, count = FALSE) {
  g <- switch(scale, difference = identity, ratio = log, odds_ratio = stats::qlogis)
  degenerate <- if (count) any(r < 1e-8) else .degenerate_risk(r)
  if (scale != "difference" && degenerate) {
    return(stats::setNames(rep(NA_real_, length(.marginal_contrast_names())),
                           .marginal_contrast_names()))
  }
  pre <- g(r[["t1p1"]]) - g(r[["t0p1"]])
  un <- g(r[["t1p0"]]) - g(r[["t0p0"]])
  ate <- g((r[["t1p1"]] + r[["t1p0"]]) / 2) - g((r[["t0p1"]] + r[["t0p0"]]) / 2)
  pretest_control <- g(r[["t0p1"]]) - g(r[["t0p0"]])
  pretest_treated <- g(r[["t1p1"]]) - g(r[["t1p0"]])
  pretest_main <- g((r[["t1p1"]] + r[["t0p1"]]) / 2) - g((r[["t1p0"]] + r[["t0p0"]]) / 2)
  stats::setNames(c(ate, pre - un, pre, un, pretest_control, pretest_treated, pretest_main),
                  .marginal_contrast_names())
}

.marginal_all <- function(b, X, pretested, scales, count = FALSE) {
  r <- .marginal_risks(b, X, pretested, linkinv = if (count) exp else stats::plogis)
  unlist(lapply(scales, function(s) .marginal_contrasts(r, s, count)), use.names = FALSE)
}

# The gradient of .marginal_all() with respect to the coefficients `b`: one
# row per contrast, in its order, and one column per coefficient. A risk
# r(t, p) is the mean over the participants of pretest condition p of
# linkinv(eta_i), with eta_i their linear predictor with the treatment
# terms set to t (.marginal_risks()), so its derivative is the mean of
# mu.eta(eta_i) x_i(t), with x_i(t) their design row with those terms set
# likewise. Each contrast is a difference of g() of risks or of their
# means, with g the identity, the log, or the logit, and takes the
# derivatives of the risks times g'. Contrasts that .marginal_contrasts()
# leaves undefined have rows of NA. The delta method took numerical
# derivatives instead, with steps of 1e-6 times each coefficient: with
# nearly collinear covariates, whose standardized coefficients are large,
# the steps were too large, and so were the standard errors (issue #137).
.marginal_gradient <- function(b, X, pretested, scales, count = FALSE) {
  b[is.na(b)] <- 0
  linkinv <- if (count) exp else stats::plogis
  slope_link <- if (count) exp else stats::dlogis
  interaction <- "treat:pretested" %in% names(b)
  b_t <- b[["treat"]]
  b_tp <- if (interaction) b[["treat:pretested"]] else 0
  base <- drop(X %*% b) - X[, "treat"] * b_t -
    if (b_tp != 0) X[, "treat:pretested"] * b_tp else 0
  pre <- pretested == 1L
  cell <- function(t, rows, p) {
    Xc <- X[rows, , drop = FALSE]
    Xc[, "treat"] <- t
    if (interaction) Xc[, "treat:pretested"] <- t * p
    eta <- base[rows] + t * b_t + t * p * b_tp
    list(risk = mean(linkinv(eta)), slope = colMeans(slope_link(eta) * Xc))
  }
  cells <- list(t1p1 = cell(1, pre, 1), t0p1 = cell(0, pre, 1), t1p0 = cell(1, !pre, 0),
                t0p0 = cell(0, !pre, 0))
  r <- vapply(cells, `[[`, numeric(1), "risk")
  D <- do.call(rbind, lapply(cells, `[[`, "slope"))
  rows <- lapply(scales, function(s) {
    degenerate <- if (count) any(r < 1e-8) else .degenerate_risk(r)
    if (s != "difference" && degenerate) {
      return(matrix(NA_real_, length(.marginal_contrast_names()), ncol(D)))
    }
    g <- function(x) switch(s, difference = 1, ratio = 1 / x, odds_ratio = 1 / (x * (1 - x)))
    treated <- (r[["t1p1"]] + r[["t1p0"]]) / 2
    control <- (r[["t0p1"]] + r[["t0p0"]]) / 2
    tested <- (r[["t1p1"]] + r[["t0p1"]]) / 2
    untested <- (r[["t1p0"]] + r[["t0p0"]]) / 2
    J <- rbind(
      c(g(treated), -g(control), g(treated), -g(control)) / 2,
      c(g(r[["t1p1"]]), -g(r[["t0p1"]]), -g(r[["t1p0"]]), g(r[["t0p0"]])),
      c(g(r[["t1p1"]]), -g(r[["t0p1"]]), 0, 0),
      c(0, 0, g(r[["t1p0"]]), -g(r[["t0p0"]])),
      c(0, g(r[["t0p1"]]), 0, -g(r[["t0p0"]])),
      c(g(r[["t1p1"]]), 0, -g(r[["t1p0"]]), 0),
      c(g(tested), g(tested), -g(untested), -g(untested)) / 2
    )
    J %*% D
  })
  out <- do.call(rbind, rows)
  dimnames(out) <- NULL
  out
}

# The sampling variance of the standardization, for the pretest effects of
# marginal_solomon()'s delta method (issue #104). `rows` index the stacked
# contrasts (each scale with the contrasts of .marginal_contrast_names()) and
# must be pretest effects; `grad` is the delta method's gradient with respect
# to the coefficients `b`. Each marginal risk is a mean of predicted risks
# over one sample, the pretested or the unpretested participants, so to
# first order a contrast moves by the sum over the two samples of the mean of
# a_i: its derivatives with respect to the four risks applied to participant
# i's predicted risks in the cells of their pretest condition.
# .center_terms() gives the variance of that sum (`v`, with degrees of
# freedom `df`) and its covariance with the coefficients, returned as the
# term 2 grad'cross (`cross`). `model` is the fit the covariance was computed
# on.
#
# `X` and `pretested` run over the participants standardized over. With
# `pretest_center = "observed"` (.marginal_rows()) some of them may be
# outside the model (issue #135): `in_model` flags those the model used, and
# `cluster` gives the cluster of each for a CR2 fit. The variance of each
# standardization is then taken over every participant in it, and its
# covariance with the coefficients over those of them in the model, as
# .center_terms() describes. By default, every participant is in the model.
.marginal_standardization_variance <- function(fit, model, b, X, pretested, scale, count,
                                               rows, grad, in_model = NULL, cluster = NULL) {
  linkinv <- if (count) exp else stats::plogis
  b0 <- b
  b0[is.na(b0)] <- 0
  b_t <- b0[["treat"]]
  b_tp <- if ("treat:pretested" %in% names(b0)) b0[["treat:pretested"]] else 0
  base <- drop(X %*% b0) - X[, "treat"] * b_t -
    if (b_tp != 0) X[, "treat:pretested"] * b_tp else 0
  pre <- pretested == 1L
  f <- cbind(
    t1p1 = ifelse(pre, linkinv(base + b_t + b_tp), 0),
    t0p1 = ifelse(pre, linkinv(base), 0),
    t1p0 = ifelse(pre, 0, linkinv(base + b_t)),
    t0p0 = ifelse(pre, 0, linkinv(base))
  )
  r <- .marginal_risks(b, X, pretested, linkinv)
  slope <- function(s, x) {
    switch(s, difference = 1, ratio = 1 / x, odds_ratio = 1 / (x * (1 - x)))
  }
  n_contrasts <- length(.marginal_contrast_names())
  set <- ifelse(pre, "pretested", "unpretested")
  if (is.null(cluster) && identical(fit$robust, "CR2")) cluster <- fit$cluster[.glm_rows_used(fit)]
  keep <- !is.na(b)

  parts <- lapply(rows, function(row) {
    s <- scale[(row - 1L) %/% n_contrasts + 1L]
    contrast <- .marginal_contrast_names()[(row - 1L) %% n_contrasts + 1L]
    J <- c(t1p1 = 0, t0p1 = 0, t1p0 = 0, t0p0 = 0)
    if (contrast == .solomon_pretest_order[1]) {
      J[c("t0p1", "t0p0")] <- c(slope(s, r[["t0p1"]]), -slope(s, r[["t0p0"]]))
    } else if (contrast == .solomon_pretest_order[2]) {
      J[c("t1p1", "t1p0")] <- c(slope(s, r[["t1p1"]]), -slope(s, r[["t1p0"]]))
    } else {
      J[c("t1p1", "t0p1")] <- slope(s, (r[["t1p1"]] + r[["t0p1"]]) / 2) / 2
      J[c("t1p0", "t0p0")] <- -slope(s, (r[["t1p0"]] + r[["t0p0"]]) / 2) / 2
    }
    terms <- .center_terms(model, fit$robust, a = drop(f %*% J), set = set, cluster = cluster,
                           used = in_model)
    c(v = terms$v, cross = 2 * sum(grad[row, keep] * terms$cross), df = terms$df)
  })
  parts <- do.call(rbind, parts)
  list(v = parts[, "v"], cross = parts[, "cross"], df = parts[, "df"])
}

# A Solomon cell with no counts at all has no rate ratio (the log-link
# estimate diverges).
.cell_empty <- function(y, treat, pretested) {
  any(tapply(y, list(treat, pretested), function(v) all(v == 0)))
}

# Classed warning for link-scale Solomon contrasts from a noncollapsible link
# with the pretest as a covariate (#43).
#
# marginal_solomon() takes the fit of a four-group design, so the fit of a
# design with several treatments (`ngroup = TRUE`) is told to subset first.
.warn_noncollapsible <- function(link, ngroup = FALSE) {
  advice <- if (ngroup) {
    paste0(
      "marginal_solomon() compares the effects on a common scale, and takes ",
      "the fit of a four-group design: subset the data to one treatment and ",
      "the control, fit the subset with fit_solomon_glm(), and give that fit ",
      "to marginal_solomon()."
    )
  } else {
    "Use marginal_solomon() to compare the effects on a common scale."
  }
  warning(structure(
    class = c("solomonR_noncollapsible_warning", "warning", "condition"),
    list(
      message = paste0(
        "With the ", link, " link and a pretest covariate, the Pretest x Treatment ",
        "contrast compares a treatment effect conditional on the pretest ",
        "(pretested participants) with a marginal one (unpretested participants). ",
        "These differ whenever the pretest predicts the outcome, even without ",
        "sensitization (Daniel et al., 2021). ", advice
      ),
      call = NULL
    )
  ))
}

# Classed warning when a link-scale contrast that involves pretesting is
# tested on a fit with a pretest covariate (issue #114): the Pretest x
# Treatment contrast on a noncollapsible link (any but the identity and the
# log; Daniel et al., 2021), and the pretest effects (issue #104) on any
# link but the identity, because a fitted mean at the mean pretest is not
# the mean over the pretests. Used by equivalence_solomon(), perm_solomon(),
# and plot_solomon_effects(bounds =).
.warn_link_scale <- function(fit, contrast) {
  if (!inherits(fit, c("solomon_glm", "solomon_ngroup")) ||
      !"pre_obs" %in% names(fit$data)) {
    return(invisible(NULL))
  }
  family <- stats::family(fit$model)
  link <- family$link
  pretest <- contrast %in% .solomon_pretest_order
  affected <- if (pretest) {
    !identical(link, "identity")
  } else {
    identical(contrast, "Pretest x Treatment") && !link %in% c("identity", "log")
  }
  if (!affected) return(invisible(NULL))

  # The pretest effects' gap is that between a fitted mean at the mean
  # pretest and the mean of the fitted means over the pretests, which a
  # nonlinear link opens; it is stated without a source. The sensitization
  # contrast's is the noncollapsibility Daniel et al. (2021) describe.
  why <- if (pretest) {
    paste0(
      "compares the pretested participants' fitted mean at the mean pretest with ",
      "the unpretested participants' mean over their unmeasured pretests. Under a ",
      "nonlinear link, a mean at the average pretest is not the average of the means ",
      "over the pretests, so these differ even when the pretest has no effect."
    )
  } else {
    paste0(
      "compares a treatment effect conditional on the pretest (pretested ",
      "participants) with a marginal one (unpretested participants). On this scale ",
      "they differ whenever the pretest predicts the outcome, even without sensitization ",
      "(Daniel et al., 2021)."
    )
  }
  advice <- if (inherits(fit, "solomon_ngroup")) {
    paste0(
      "marginal_solomon() estimates the contrast on a common scale, from standardized ",
      "risks or rates; it takes the fit of a four-group design, so subset the data to ",
      "one treatment and the control first."
    )
  } else {
    "marginal_solomon() estimates the contrast on a common scale, from standardized risks or rates."
  }
  warning(structure(
    class = c("solomonR_link_scale_warning", "warning", "condition"),
    list(
      message = paste0(
        "The ", contrast, " contrast is on the ", .contrast_scale(family), " scale (",
        link, " link). With the pretest as a covariate, it ", why, " ", advice
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
# `epsilon` within `maxit` iterations; ?glm.control). A lean version for
# bootstrap refits, warm-started at `start`. Each step solves its weighted
# least-squares problem by QR, with the tolerance of glm.fit(); solving the
# normal equations, which square the design's condition number, failed in
# every resample with nearly collinear covariates, such as a covariate and
# its copy in single precision (issue #137).
# Each step is a full step: glm.fit() halves a step that makes the deviance
# infinite (Marschner, 2011), while this refit stops. So it can diverge
# from `start` on a resample whose maximum-likelihood estimate exists.
# A fit that fails gives its `reason` (issue #154). "risk": fitted risks
# reach 0 or 1 in floating point, which makes the deviance, the working
# response, or the coefficients non-finite; or the weighted design loses
# rank while `X` has full rank and a fitted risk is within 1e-8 of 0 or 1,
# because the weights mu (1 - mu) of those rows vanish. That happens under
# separation and when the refit diverges. "rank": `X` itself is not of
# full rank at that tolerance, as when a resample lacks a level of a
# factor or covariates are collinear to within it, or the weighted design
# loses rank with no fitted risk that close to 0 or 1. `X` is checked
# only once the weighted design has lost rank. "maxit": the rule
# is not met in `maxit` iterations; `kappa` is the 2-norm condition number
# of the last weighted design (?kappa, exact = TRUE). The arithmetic is
# that of the version before #154, so a fit that converged then converges
# now with the same coefficients, and a fit that failed then fails now.
.logistic_irls <- function(X, y, start, epsilon = 1e-8, maxit = 25L) {
  failed <- function(reason) list(converged = FALSE, reason = reason)
  events <- y == 1
  deviance <- function(mu) -2 * (sum(log(mu[events])) + sum(log1p(-mu[!events])))
  b <- start
  eta <- drop(X %*% b)
  mu <- stats::plogis(eta)
  dev <- deviance(mu)
  if (!is.finite(dev)) return(failed("risk"))
  for (iter in seq_len(maxit)) {
    w <- mu * (1 - mu)
    z <- eta + (y - mu) / w
    if (!all(is.finite(z))) return(failed("risk"))
    root <- sqrt(w)
    q <- tryCatch(qr(root * X, tol = 1e-11), error = function(e) NULL)
    if (is.null(q)) return(failed("risk"))
    if (q$rank < ncol(X)) {
      vanished <- .degenerate_risk(mu) && qr(X, tol = 1e-11)$rank == ncol(X)
      return(failed(if (vanished) "risk" else "rank"))
    }
    b <- tryCatch(qr.coef(q, root * z), error = function(e) NULL)
    if (is.null(b) || !all(is.finite(b))) return(failed("risk"))
    eta <- drop(X %*% b)
    mu <- stats::plogis(eta)
    dev_new <- deviance(mu)
    if (!is.finite(dev_new)) return(failed("risk"))
    if (abs(dev_new - dev) / (abs(dev_new) + 0.1) < epsilon) {
      return(list(coefficients = b, fitted.values = mu, converged = TRUE, iter = iter))
    }
    dev <- dev_new
  }
  list(coefficients = b, fitted.values = mu, converged = FALSE, reason = "maxit",
       kappa = kappa(q, exact = TRUE))
}

# The bootstrap's design in a well-conditioned basis (issue #154), for the
# refits that do not converge on the model's design. The columns other
# than the treatment terms, X_N, are replaced by the orthonormal columns Q
# of their Householder QR decomposition, X_N = Q R (Golub, 1965), which
# span the same space; the treatment columns, whose coefficients
# .marginal_risks() reads, are kept. Each IRLS step is a weighted
# least-squares fit, whose fitted values depend only on the space the
# columns span (Green, 1984), so in exact arithmetic every iterate, fitted
# risk, and contrast is the same in either basis, from the start `b`
# mapped to it (R b_N). In floating point the fitted values of a
# least-squares step move by about the unit roundoff times the design's
# condition number (Higham, 2002, Chapter 20). With nearly collinear
# covariates (a condition number near 1e14), between iterations 5 and 25
# the deviance changed by 4e-7 to 7e-4 of itself, against the 1e-8 that
# the rule allows, and the fitted risks by up to 0.06, so almost no
# resample met the rule. Q is orthonormal and the QR is backward stable
# (Higham, 2002, Theorem 19.4), so the condition number left for the IRLS
# steps is that of the treatment terms against the other columns, and
# the refits converge. That is a statement about convergence, not
# accuracy: Q is computed from X_N in floating point, so its columns span
# those of X_N only to within about the unit roundoff times the condition
# number of X_N (Higham, 2002, Theorem 19.4 and Chapter 20), an angle of
# 0.02 to 0.08 radians on the data of #154. A refit in the basis is
# accurate to about that, as the model's own fit is, and a resample whose
# fitted risks come near the 1e-8 limit can be counted differently than
# under an exact reparameterization.
#
# Rank. A refit that fails as "rank" is not fitted again, and the
# resample counts as failed: its design is not of full rank at qr()'s
# tolerance, 1e-11, on the rows of the resample, as when the resample
# lacks a level of a factor or makes nearly collinear covariates
# collinear. In the basis, a direction that the resample makes collinear
# is carried only by rounding error, which a rank check cannot tell from
# a direction in the data, and a refit there can converge on it. In the
# data of a test of this in test-covariate-location.R (seed 4344), z2
# equals z1 except in two participants with an outcome, and in one
# pretested participant without, by 1e-6 to 2e-6, so a resample without
# the first two has z2 == z1 on the rows it is fitted to. Fitted again in
# the basis, as an earlier version of this fix did, the 15 such resamples
# of 199 converged, with risk differences up to 0.03 from those of the
# refit without z2, which is the same model on those rows. Under
# pretest_center = "observed", the 31 such resamples converged too,
# though the participant without an outcome, who is standardized over,
# makes their contrasts depend on the coefficient of z2, which they do
# not identify. The exact reparameterization (z1, (z2 - z1) / 1e-6) fails
# in the same resamples.
#
# A refit that fails as "maxit" passed that rank check at every
# iteration, but that does not show that its resample identifies the
# coefficients. When a covariate is a combination of others with large
# coefficients, the rounding error of the stored and standardized values,
# multiplied by those coefficients, can keep a design that is not of full
# rank in the data above 1e-11 of its norm. In the example of
# .representation_error(), a time of testing that is the time of
# enrollment plus a delay, with the times in seconds since 1970, the
# testing times of two participants were corrected, so that the relation
# fails only for them (seed 6101 in test-covariate-location.R). On a
# resample without the two, tested - enrolled == delay on every row. Its
# refit on the model's design passed the rank check at every iteration,
# with condition numbers of 1e15 to 1e16, and did not converge. Fitted
# again in the basis, as an earlier version of this fix did, each of the
# 43 such resamples of 199 converged, with risk differences up to 0.088
# from those of the refit without the delay (0.18 under pretest_center =
# "observed"), and the intervals were reported where the exact
# reparameterization (enrolled, tested - enrolled, delay) withholds them.
# Before a refit that failed as "maxit" is fitted again, its resample is
# therefore checked by the rule that fit_solomon_glm() applies to the
# model (.collinear_terms()), on the rows of the resample, in the
# model's working weights and with the rounding error the fit computed
# for their entries (.bootstrap_rank_check(), .resample_collinear()). If
# a term is collinear with the others by that rule, the refit counts as
# failed with the reason "rank" and is not fitted again. The check is
# the fit's rule applied to the resample's rows, not a fit to them:
# fit_solomon_glm() given those rows would weight and standardize them
# anew. What is left of that term on those rows after the others is then
# within the rounding error the rule allows for (4 B_k of
# .collinear_terms()), or within glm()'s tolerance, 1e-11 of its norm.
# In the first case, as on each of the 43 timestamp resamples, a refit
# in the basis could converge on rounding error. Near the edge of the
# fit's check, as in the next paragraph, the rule stops refits by the
# 1e-11 tolerance alone; what is left there is the data's own direction,
# and those resamples are the cost described there. On the timestamp
# data fit_solomon_glm() stops on the rows of each of the 43 resamples,
# at the delay; near the edge of its check it need not. With the data of
# the next paragraph at g = 1.1e-11 to 1.3e-11, over bootstrap seeds 1
# to 20, it went ahead on 47 of the 136 resamples that the rule stopped,
# on all of which the exact reparameterization's refit converged. On the
# timestamp data the bootstrap fails in the same resamples as the exact
# reparameterization, 43 of 199 under "analyzed" and 33 under
# "observed", and withholds the intervals as it does. Only refits that
# did not converge are checked. Applied to every resample, the check
# found no term collinear on the resamples whose refit converged on
# those data, on the data of the next paragraph at g = 1.2e-11 and above
# (seed = 2024), and on the nested covariates of #154; at g = 1.1e-11 to
# 1.15e-11 it found one on 4 or 5 of the 62 to 87 resamples whose refit
# converged, which are accepted as before.
#
# The cost is a band just inside the collinearity check of
# fit_solomon_glm(), where rounding error can make nearly collinear
# covariates collinear, to within qr()'s tolerance or by the fit's rule,
# on a resample whose design has full rank in exact arithmetic. With
# z2 = z1 + g e1 (seed 5151, as in another test there; the fit stops at
# g = 1.09e-11 and below on Windows) and seed = 2024, 3, 24, 61, 86, and
# 106 of 199 resamples failed so at g = 1.3e-11, 1.2e-11, 1.15e-11,
# 1.12e-11, and 1.1e-11, where the exact reparameterization
# (z1, (z2 - z1) / g) has no failures. At the last two, 3 of them did not
# converge and were collinear by the fit's rule. With other seeds the
# rule also stops a few refits that did not converge from g = 1.15e-11
# up to at least 1.4e-11, about 30% above where the fit stops: over seeds
# 1 to 20, 2, 19, and 25 in all at 1.3e-11, 1.2e-11, and 1.15e-11, and
# over seeds 21 to 80, 2 at 1.4e-11 and none at 1.45e-11. That is at
# most 1.25 per seed on average, small next to the resamples that fail
# qr()'s rank check there, 3 and 24 with seed = 2024 at 1.3e-11 and
# 1.2e-11. With z3 = z1 + z2 + g e1 (seed
# 7501), whose fit stops at g = 1.75e-11 and goes ahead at 1.8e-11, the
# rule stopped 1 refit at g = 2e-11 with seed = 2024, and 4 over seeds 1
# to 20. A refit can also fail as "rank" far from the edge of the
# fit's check, when one participant alone carries a direction of the
# design: the refit moves that participant's fitted risk toward the
# observed outcome, and the participant's weight in the weighted design
# vanishes, as when separation drives a coefficient toward infinity
# (Albert & Anderson, 1984), until the weighted design falls short of
# full rank at 1e-11 while no fitted risk is yet within 1e-8 of 0 or 1.
# In the seed 4343 data of that test, with z2 = z1 except in the two
# participants, by delta and 1.5 delta, 30 to 105 of 199 resamples that
# hold one of the two failed so at delta = 5e-8 to 1e-9, where the
# reparameterization's refits converge, before this fix as well.
#
# Q is formed from every column of X_N, with no rank decision of its own:
# at tol = 0, qr() moves no column to the end. (qr.Q() applies only as
# many reflections as the rank, so a column moved to the end would be
# replaced by a direction outside the span of X_N.) fit_solomon_glm()
# stops when a term is collinear with the others, in the model as given
# and standardized (.collinear_terms()), so X_N has full rank on the
# model's rows as glm() and that check judge it, in the working weights
# of the fit, and so on any set of rows that contains them. Without
# those weights X_N can fall short of full rank at 1e-11: on the data
# above at g = 1.1e-11, qr(X_N, tol = 1e-11) has rank 4 of 5. A basis
# formed at that tolerance did not exist, so the 34 refits that did not
# converge failed, with the 103 not of full rank; formed at tol = 0, it
# lets the 31 of them that are not collinear by the fit's rule converge.
# Where qr(X_N, tol = 1e-11) has full rank, it moves no column, and Q and
# R are the same, bit for bit, at either tolerance. Returns the design
# `X` and the start `start` in the basis, or NULL when X_N is not finite
# or a column's residual is exactly zero, as for a column of zeros. An
# exact relation computed in floating point leaves a residual of rounding
# error, so this does not replace the check of fit_solomon_glm().
.bootstrap_basis <- function(X, b) {
  treatment <- match(c("treat", "treat:pretested"), colnames(X))
  treatment <- treatment[!is.na(treatment)]
  others <- setdiff(seq_len(ncol(X)), treatment)
  if (!all(is.finite(X[, others]))) return(NULL)
  q <- qr(X[, others, drop = FALSE], tol = 0)
  R <- qr.R(q)
  if (any(diag(R) == 0)) return(NULL)
  X[, others] <- qr.Q(q)
  b[others] <- drop(R %*% b[others])
  list(X = X, start = b)
}

# The collinearity check of fit_solomon_glm() (.collinear_terms()),
# prepared for the bootstrap resamples of marginal_solomon() whose refit
# does not converge on the model's design (issue #154;
# .resample_collinear()). `X` is the design that the bootstrap refits,
# one row per participant it resamples; `model` is the model whose design
# it is, standardized when the fit standardized anything; and `in_model`
# flags the rows of `X` that `model` used, which are its rows, in its
# order. Returns `X`; `noise`, the rounding error of each of its entries,
# computed by .representation_error() from the model as given and the
# model standardized, as fit_solomon_glm() computed it for its own check;
# `w`, the weight of each row; `last`, which flags the covariates'
# columns, checked last; and `intercept`, the position of the intercept's
# column.
#
# The weights are the working weights of the model, those in which glm()
# and the check of fit_solomon_glm() decided that its terms were not
# collinear. A resample's row takes the weight of the participant it
# draws, so a participant drawn twice enters twice with that weight, as
# a row given twice to glm() would. A refit starts from the model's
# coefficients, so its first step weights the rows nearly so. A relation
# that holds exactly in the data holds row by row, and so in any positive
# weights, and the bound it is compared with is computed in the same
# weights. The weights of the refit's last step are not used: after 25
# steps without converging, its coefficients can have moved along the
# direction that the resample makes collinear, and the check would then
# depend on how far they moved. The rows outside the model, under
# `pretest_center = "observed"`, are never refitted; they get a weight of
# zero, and rows with a weight of zero are left out of the check, as
# .collinear_terms() leaves them out.
.bootstrap_rank_check <- function(fit, model, X, in_model) {
  noise <- matrix(0, nrow(X), ncol(X))
  noise[in_model, ] <- .representation_error(model, fit$model, fit$covariate_scaling,
                                             .null_na(fit$pretest_mean))
  w <- numeric(nrow(X))
  w[in_model] <- model$weights
  pretest <- is.null(fit$pretest_mean) || is.finite(fit$pretest_mean)
  covariates <- setdiff(names(fit$data), c("y", "treat", "pretested", "log_exposure",
                                           if (pretest) "pre_obs"))
  assign <- attr(stats::model.matrix(model), "assign")
  labels <- attr(stats::terms(model), "term.labels")
  list(X = X, noise = noise, w = w,
       last = assign %in% match(.formula_names(covariates), labels),
       intercept = match(0L, assign))
}

# Whether the rows `rows` of the bootstrap's design, with repeats, as a
# resample draws them, make a term collinear with the others by the check
# of fit_solomon_glm(), from `check` (.bootstrap_rank_check()): each term
# is regressed on those before it, the covariates last, in the weights
# of these rows and with each column but the intercept centered at its
# weighted mean over them, and is taken to be collinear when its residual
# is within 4 times the rounding error that its regression carries into
# it, or within 1e-11 of its norm (.collinear_columns()).
.resample_collinear <- function(check, rows) {
  length(.collinear_columns(check$X[rows, , drop = FALSE], check$w[rows],
                            check$noise[rows, , drop = FALSE], check$last,
                            check$intercept)) > 0L
}

# The warning of a bootstrap in which more than 10% of the resamples
# failed, with the failures by reason (issue #154): `reasons` holds one
# of "separated", "risk", "rank", or "maxit" per failed resample, and
# `kappa` the condition numbers of those that did not converge (NA for
# the others). Each cause is named only for the reasons that occurred.
# The condition number is reported, not judged: no threshold is applied,
# and the sentence does not say that rounding caused the failures. In
# double precision a least-squares step carries a relative error of about
# the unit roundoff (1.1e-16) times the condition number (Higham, 2002,
# Chapter 20), about 1% near 1e14, where the refits of #154 did not
# converge.
.bootstrap_failure_message <- function(reasons, R, kappa) {
  n <- vapply(c("separated", "risk", "rank", "maxit"), function(r) sum(reasons == r), integer(1))
  parts <- c(
    separated = "%d had a Solomon cell with only events or only non-events",
    risk = "%d had a fitted risk within 1e-8 of 0 or 1",
    rank = paste("%d had a design not of full rank, as when a resample lacks a level of a factor",
                 "or makes nearly collinear covariates collinear"),
    maxit = "%d did not converge in 25 iterations"
  )
  parts <- sprintf(parts, n)[n > 0]
  if (length(parts) > 1L) {
    parts[length(parts)] <- paste("and", parts[length(parts)])
  }
  reasons_text <- paste0(paste(parts, collapse = if (length(parts) > 2L) ", " else " "), ".")
  sparse_text <- if (n[["separated"]] > 0) {
    " Cells with only events or only non-events come from sparse cells."
  } else {
    ""
  }
  risk_text <- if (n[["risk"]] > 0) {
    paste0(" Fitted risks near 0 or 1 come from separation or near separation, from ",
           "extreme covariate values, or from a refit that diverged: each refit takes ",
           "full steps from the model's coefficients, without step halving.")
  } else {
    ""
  }
  maxit_text <- if (n[["maxit"]] > 0) {
    sprintf(paste0(" For the refits that did not converge, the largest condition number ",
                   "of the weighted design was %.2g. Rounding error in each step grows ",
                   "with it, and near 1e14, as with nearly collinear covariates, it can ",
                   "keep a refit from converging."),
            max(kappa, na.rm = TRUE))
  } else {
    ""
  }
  sprintf(
    paste0("The logistic fit failed in %d of %d bootstrap resamples (more than 10%%), ",
           "so intervals are not reported. Of these, %s%s%s%s method = \"delta\" does not ",
           "resample; in the package's simulation study (issue #43), its ",
           "risk-difference intervals met the registered tolerances in most designs ",
           "in which the bootstrap failed this often."),
    length(reasons), R, reasons_text, sparse_text, risk_text, maxit_text
  )
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
#'   risks;
#' - `Pretest effect | control` and `Pretest effect | treated`: the pretested
#'   cell compared with the unpretested cell, among control and among
#'   treated participants (issue #104);
#' - `Pretest main effect`: the pretest effect in a population with equal
#'   numbers of treated and control participants, computed from the averaged
#'   risks.
#'
#' The pretest effects compare two different sets of participants, because
#' the unpretested participants have no pretest from which to predict their
#' risk had they been pretested. Each set is standardized over its own
#' participants, and with random assignment both estimate the same
#' population: the expected pretest of the unpretested participants equals
#' that of the pretested ones (Solomon, 1949, p. 141). Both
#' risks are then marginal, so the comparison is on a common scale, unlike
#' the link-scale pretest effects of [fit_solomon_glm()].
#'
#' **Missing outcomes** (issue #135). When outcomes are missing, the
#' pretested cells can be standardized over two sets of pretested
#' participants, which `pretest_center` chooses, as it chooses the center of
#' [fit_solomon_glm()] (see "The pretest effect when posttests are missing"
#' there, for what each assumes): by default over the pretested participants
#' in the model (`"analyzed"`), or over every pretested participant with a
#' pretest and the model's covariates (`"observed"`; a participant without
#' a covariate is in the observed center of [fit_solomon_glm()] but not in
#' this set). On the risk or rate scale the choice moves
#' Treatment | pretested, the average treatment effect, and Pretest x
#' Treatment as well as the pretest effects, because a marginal risk depends
#' on the participants it is averaged over; only Treatment | unpretested does
#' not move. Whenever the two sets differ, the print gives the moved
#' contrasts under both, with their difference, and `pretest_sensitivity`
#' returns every contrast under both: with full inference for the other set
#' under the delta method, and estimates only under the bootstrap, so that
#' the resamples are those of the set used. No delta line is given on these
#' scales, because no verified source gives one there.
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
#' risk differences, log risk ratios, or log odds ratios.
#'
#' Bootstrap resamples in which the logistic fit fails are excluded and
#' counted: those with a Solomon cell of only events or only non-events, a
#' fitted probability within 1e-8 of 0 or 1, a design not of full rank (as
#' when a resample lacks a level of a factor or makes nearly collinear
#' covariates collinear), or no convergence in 25 iterations by the rule of
#' [stats::glm.control()]. When more than 10% fail, intervals are not
#' reported, and the warning gives the failures by reason.
#' Each resample is refitted by iteratively reweighted least squares from
#' the model's coefficients, with full steps. Where a step takes fitted
#' probabilities to 0 or 1, the refit stops and the resample fails;
#' [stats::glm.fit()] would instead halve a step that makes the deviance
#' infinite (Marschner, 2011). So a fitted probability near 0 or 1 comes
#' from complete or quasicomplete separation (Albert & Anderson, 1984) or
#' near separation, from an extreme covariate value, or from a refit that
#' diverged.
#' A refit that does not converge is fitted again with the intercept, the
#' pretesting indicator, the pretest, and the covariates replaced by the
#' orthonormal columns of their QR decomposition (Golub, 1965). These span
#' the same space, and the treatment terms are kept, so in exact arithmetic
#' it is the same fit (Green, 1984). With nearly collinear covariates,
#' rounding error can keep a refit from converging, and the new columns
#' let it converge. The new columns are formed whenever [fit_solomon_glm()]
#' fitted the model, since it stops at a term collinear with the others.
#' A refit whose design is not of full rank at the tolerance of [qr()]
#' (1e-11) is not fitted again, and the resample counts as failed. In the
#' new columns, a direction that a resample makes collinear is carried
#' only by rounding error, and a refit could converge on it. Rounding
#' error can also keep a design that is not of full rank above that
#' tolerance, as when a covariate is computed from others, such as a time
#' of testing that is the time of enrollment plus a delay, with the times
#' in seconds. So before a refit that did not converge is fitted again,
#' its resample is checked by the collinearity check of
#' [fit_solomon_glm()] (see "Covariates" there), applied to the rows of the
#' resample with the weights and rounding error of the model's own fit: a
#' resample on whose rows a term is collinear with the others by that
#' check counts as failed, as a design not of full rank, and is not fitted
#' again. A refit that converged is not checked. Just inside that check,
#' rounding error can also make nearly collinear covariates collinear on
#' some resamples, at the tolerance of [qr()] or by the check, and those
#' resamples fail too, though an exact reparameterization of the
#' covariates would fit them. A
#' resample can also fail as a design not of full rank when only one of its
#' participants departs from a relation that holds among the covariates of
#' all the others: the refit moves that participant's fitted probability
#' toward the observed outcome, and the participant's weight in the refit
#' vanishes before the probability comes within 1e-8 of 0 or 1. The new
#' columns cannot help with a treatment term nearly collinear with the
#' other terms, because the contrasts are read from the treatment
#' coefficients.
#' They are computed in floating point, so they span the old ones only to
#' within about the unit roundoff times the condition number of the
#' columns they replace (Higham, 2002, Theorem 19.4 and Chapter 20). A
#' refit in them is accurate to about that, as the model's own fit is, and
#' a resample with a fitted probability near the 1e-8 limit can be counted
#' as failed, or not, differently than in exact arithmetic. For the refits
#' that still do not converge, the warning gives the largest condition
#' number of their weighted design: the rounding error of each step grows
#' with it (Higham, 2002, Chapter 20). A bootstrap whose refits all
#' converge is computed as before.
#'
#' For the pretest effects, the pretested and the unpretested risks are
#' standardized over different participants, so the sampling variance of
#' each standardization does not cancel, as it does for the treatment
#' contrasts, which standardize both arms over the same participants. The
#' delta method therefore adds it for the pretest effects, with its
#' covariance with the coefficients under HC3 or CR2 covariance, by stacking
#' the estimating equations of the standardized means with those of the model
#' (Stefanski & Boos, 2002), as [fit_solomon_glm()] does for the mean
#' pretest; for clustered fits, the degrees of freedom are combined by the
#' Welch-Satterthwaite formula (Satterthwaite, 1946; Welch, 1947). The
#' bootstrap resamples participants and so includes this variance.
#'
#' Under HC3, a delta-method contrast that depends on a participant whose
#' leverage is 1 has an `NA` standard error, interval, and p-value, as in
#' the fit (see "HC3 at leverage 1" in [fit_solomon_glm()]).
#'
#' Count outcomes: for a fit with `family = poisson()` or
#' `family = "negative_binomial"`, rates per unit of
#' exposure are standardized in the same way and compared as rate differences
#' or rate ratios. Log-link rate ratios are collapsible (Daniel et al., 2021),
#' so the treatment contrasts' rate ratios agree with the fitted model's
#' coefficients when there are no other covariates; rate differences depend on
#' the covariate distribution. The pretest effects' rate ratios do not agree
#' with the model's pretest effects when the pretest is a covariate, because
#' the model compares a rate at the mean pretest with a marginal rate. For
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
#' degrees of freedom. It gives the four treatment contrasts only, the
#' comparisons of treated and control clusters that the study validated. With pretesting assigned within clusters, each cluster
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
#' @param pretest_center For a fit with a pretest covariate, the pretested
#'   participants the pretested cells are standardized over (issue #135):
#'   `"analyzed"`, the pretested participants in the model; or `"observed"`,
#'   every pretested participant with a pretest and the model's other
#'   covariates, with or without an outcome or a recorded condition.
#'   Defaults to the fit's `pretest_center` (`"analyzed"` for fits made
#'   before it existed). The centering constant itself plays no role.
#'
#' @return An object of class `solomon_marginal`, a list with:
#'   - `effects`: one row for each scale and contrast, in the columns
#'     `scale`, `contrast`, `estimate`, `std.error`, `statistic`, `df`,
#'     `p.value`, `conf.low`, and `conf.high`. On each scale the rows are
#'     the four treatment contrasts, `ATE (avg over pretest)`,
#'     `Pretest x Treatment`, `Treatment | pretested`, and
#'     `Treatment | unpretested`, followed by the three pretest effects,
#'     `Pretest effect | control`, `Pretest effect | treated`, and
#'     `Pretest main effect`; `method = "cluster_summary"` gives the four
#'     treatment contrasts. The estimate and the interval are on the
#'     reporting scale (a difference or a ratio). The standard error and the
#'     test statistic are on the analysis scale, the difference or the log
#'     of the ratio, so for a ratio `statistic` is
#'     `log(estimate) / std.error`. `df` is `Inf` for a normal reference
#'     distribution, which every fit without CR2 covariance uses.
#'   - `conf_level`: the confidence level of the intervals.
#'   - `risks` (binary outcomes) or `rates` (counts, per unit of exposure):
#'     the value in each of the four cells (`cell`, `treat`, `pretested`,
#'     and `risk` or `rate`). The other of the two is `NULL`.
#'   - `outcome`: `"binary"` or `"count"`.
#'   - `method`: the method of inference, and `vcov`, a description of the
#'     covariance the delta method used (`NA` for the other methods).
#'   - `R` and `failures`: with `method = "bootstrap"`, the number of
#'     resamples and the number that failed.
#'   - `dispersion` and `theta`: as in the fit.
#'   - `pretest_adjusted`: whether the fit adjusted for the pretest.
#'   - `pretest_center`: the set the pretested cells were standardized over,
#'     `"analyzed"` or `"observed"`. Absent with `method = "cluster_summary"`,
#'     which uses neither set, so that `pretest_center` is ignored there.
#'   - `standardized`: the numbers of pretested participants each set holds
#'     (`analyzed` and `observed`), and of those outside the model, how many
#'     have no outcome (`no_outcome`). Absent with
#'     `method = "cluster_summary"`.
#'   - `attrition`: each Solomon group's `n`, outcomes missing
#'     (`posttest_missing`), and `share`, from the fit, for the attrition
#'     line of the print and the report; `NULL` for fits made before 1.0.0.
#'   - `pretest_sensitivity`: for a fit with a pretest covariate, every
#'     contrast under both sets, a list with `effects` (keyed by `scale` and
#'     `center`; the other set's standard errors and intervals are `NA`
#'     under the bootstrap), `difference` (the estimates under the
#'     `analyzed` and `observed` sets and their `difference`), `differ`, and
#'     `note`; `NULL` without a pretest covariate and with
#'     `method = "cluster_summary"`. Outside the stable interface (see
#'     [solomon_output]).
#'   - `design` and `clusters` (the number of clusters in each cell), with
#'     `method = "cluster_summary"`.
#'
#'   The `effects` table, `conf_level`, and [`tidy()`][solomon_output],
#'   which returns the table, are the stable interface of the result; see
#'   [solomon_output].
#'
#' @references
#' Albert, A., & Anderson, J. A. (1984). On the existence of maximum
#' likelihood estimates in logistic regression models. *Biometrika, 71*(1),
#' 1–10. https://doi.org/10.1093/biomet/71.1.1
#'
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
#' Golub, G. (1965). Numerical methods for solving linear least squares
#' problems. *Numerische Mathematik, 7*(3), 206–216.
#' https://doi.org/10.1007/BF01436075
#'
#' Green, P. J. (1984). Iteratively reweighted least squares for maximum
#' likelihood estimation, and some robust and resistant alternatives.
#' *Journal of the Royal Statistical Society: Series B (Methodological),
#' 46*(2), 149–170. https://doi.org/10.1111/j.2517-6161.1984.tb01288.x
#'
#' Hayes, R. J., & Moulton, L. H. (2017). *Cluster randomised trials* (2nd
#' ed.). Chapman and Hall/CRC. https://doi.org/10.4324/9781315370286
#'
#' Higham, N. J. (2002). *Accuracy and stability of numerical algorithms* (2nd
#' ed.). Society for Industrial and Applied Mathematics.
#' https://doi.org/10.1137/1.9780898718027
#'
#' Localio, A. R., Margolis, D. J., & Berlin, J. A. (2007). Relative risks and
#' confidence intervals were easily computed indirectly from multivariable
#' logistic regression. *Journal of Clinical Epidemiology, 60*(9), 874–882.
#' https://doi.org/10.1016/j.jclinepi.2006.12.001
#'
#' Marschner, I. C. (2011). glm2: Fitting generalized linear models with
#' convergence problems. *The R Journal, 3*(2), 12–15.
#' https://doi.org/10.32614/RJ-2011-012
#'
#' Pustejovsky, J. E., & Tipton, E. (2018). Small-sample methods for
#' cluster-robust variance estimation and hypothesis testing in fixed effects
#' models. *Journal of Business & Economic Statistics, 36*(4), 672–683.
#' https://doi.org/10.1080/07350015.2016.1247004
#'
#' Satterthwaite, F. E. (1946). An approximate distribution of estimates of
#' variance components. *Biometrics Bulletin, 2*(6), 110–114.
#' https://doi.org/10.2307/3002019
#'
#' Solomon, R. L. (1949). An extension of control group design. *Psychological
#' Bulletin, 46*(2), 137–150. https://doi.org/10.1037/h0062958
#'
#' Stefanski, L. A., & Boos, D. D. (2002). The calculus of M-estimation. *The
#' American Statistician, 56*(1), 29–38.
#' https://doi.org/10.1198/000313002753631330
#'
#' Tipton, E. (2015). Small sample adjustments for robust variance estimation
#' with meta-regression. *Psychological Methods, 20*(3), 375–393.
#' https://doi.org/10.1037/met0000011
#'
#' Welch, B. L. (1947). The generalization of "Student's" problem when several
#' different population variances are involved. *Biometrika, 34*(1–2), 28–35.
#' https://doi.org/10.1093/biomet/34.1-2.28
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
                             seed = NULL, conf_level = fit$conf_level,
                             pretest_center = fit$pretest_center) {

  .stop_ngroup_fit(fit, "marginal_solomon")
  if (!inherits(fit, "solomon_glm")) {
    stop("`fit` must come from fit_solomon_glm().", call. = FALSE)
  }
  # Fits made before issue #135 have no `pretest_center`; theirs is the
  # default.
  if (is.null(pretest_center)) pretest_center <- "analyzed"
  pretest_center <- match.arg(pretest_center, c("analyzed", "observed"))
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
    out <- .marginal_cluster_summary(fit, conf_level)
    out$attrition <- .fit_attrition(fit)
    return(out)
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

  # With a pretest or numeric covariates the fit's inference was computed
  # with the covariates standardized and the pretest scaled (issue #137;
  # .covariate_scaling()), and so is this: the model, its coefficients and
  # covariance, and the design rows are those of the standardized model, in
  # which the delta method (.marginal_gradient(), .combination_variance())
  # and the refits of the bootstrap are accurate however far a covariate
  # lies from zero, and whatever the scale of a covariate or the pretest.
  # The risks and the contrasts are the same. `V` is the covariance of `b`.
  # The refitted model codes each factor as the fit did, so that its design
  # rows go with the fit's coefficients (issue #168; .standardized_model()).
  scaling <- fit$covariate_scaling
  model <- fit$model
  b <- stats::coef(model)
  V <- fit$vcov
  if (!is.null(scaling)) {
    model <- .standardized_model(fit)
    b <- scaling$coefficients
    V <- scaling$vcov
  }
  y <- model$y
  n_contrasts <- length(.marginal_contrast_names())
  linkinv <- if (count) exp else stats::plogis

  # The participants standardized over (issue #135). By default, those in
  # the model: the pretested cells are standardized over the pretested
  # participants the model used. With `pretest_center = "observed"`, the
  # pretested cells are standardized over every pretested participant with
  # a pretest and the model's other covariates, with or without an outcome
  # (.marginal_rows()). `X`, `pretested`, and `treat` run over them;
  # `in_model` flags the rows the model used, to which `y` belongs.
  model_set <- function() {
    X <- stats::model.matrix(model)
    list(X = X, pretested = as.integer(X[, "pretested"]), treat = as.integer(X[, "treat"]),
         in_model = NULL, index = NULL)
  }
  observed_set <- .marginal_rows(fit)
  if (!is.null(scaling)) {
    observed_set$X <- .standardized_columns(observed_set$X, scaling, fit$model)
  }
  standardized <- observed_set$n
  set <- if (pretest_center == "analyzed") model_set() else observed_set
  X <- set$X
  pretested <- set$pretested
  treat <- set$treat
  in_model <- if (is.null(set$in_model)) rep(TRUE, nrow(X)) else set$in_model

  # The estimates over a set of participants, with the ratio scales set to
  # NA when a cell is separated (binary) or empty (counts).
  estimate_at <- function(set) {
    keep <- if (is.null(set$in_model)) rep(TRUE, nrow(set$X)) else set$in_model
    estimate <- .marginal_all(b, set$X, set$pretested, scale, count)
    risks <- .marginal_risks(b, set$X, set$pretested, linkinv = linkinv)
    separated <- if (count) {
      .cell_empty(y, set$treat[keep], set$pretested[keep]) || any(risks < 1e-8)
    } else {
      .cell_separated(y, set$treat[keep], set$pretested[keep]) || .degenerate_risk(risks)
    }
    if (separated) {
      estimate[rep(scale, each = n_contrasts) != "difference"] <- NA_real_
    }
    list(estimate = estimate, risks = risks, separated = separated)
  }
  at <- estimate_at(set)
  estimate <- at$estimate
  risks <- at$risks
  separated <- at$separated
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

  # The delta method over a set of participants: the standard errors on the
  # analysis scale, the degrees of freedom, and the intervals. `cluster` is
  # the cluster of each participant of the set, for CR2 fits whose set has
  # participants outside the model. The gradient is exact
  # (.marginal_gradient()). Its variance is computed through the QR
  # decomposition of the model (.combination_variance()), or under CR2 with
  # a numeric covariate in its orthonormal basis (.cr2_basis()): with nearly
  # collinear covariates the gradient weights their standardized
  # coefficients nearly equally, and grad' V grad canceled the large
  # entries of their covariance, so that the standard errors were NaN or up
  # to 78 times too large (issue #137).
  delta_at <- function(set, estimate, cluster = NULL) {
    X <- set$X
    pretested <- set$pretested
    df <- rep(Inf, length(estimate))
    grad <- .marginal_gradient(b, X, pretested, scale, count)
    colnames(grad) <- names(b)
    basis <- NULL
    if (cr2) {
      # The complete-case refit that the CR2 covariance was computed on,
      # with each factor coded as in the fit (issue #168).
      fit_cr <- stats::glm(stats::formula(model), data = stats::model.frame(model),
                           family = fit$family, contrasts = fit$model$contrasts)
      if (.numeric_covariates(scaling)) {
        basis <- .cr2_basis(fit_cr, fit$cluster[.glm_rows_used(fit)])
      }
    }
    std.error <- rep(NA_real_, length(estimate))
    finite <- is.finite(estimate) & rowSums(!is.finite(grad)) == 0
    if (any(finite)) {
      G <- grad[finite, , drop = FALSE]
      std.error[finite] <- sqrt(if (!cr2) {
        .combination_variance(model, fit$robust, G)
      } else if (is.null(basis)) {
        rowSums((G %*% as.matrix(V)) * G)
      } else {
        U <- backsolve(basis$R, t(G), transpose = TRUE)
        colSums(U * (.plain_matrix(basis$vcov) %*% U))
      })
    }
    if (cr2) {
      # Satterthwaite degrees of freedom for the delta method's linear
      # approximation (Pustejovsky & Tipton, 2018).
      ok <- is.finite(estimate) & is.finite(std.error) & std.error > 0
      df[ok] <- if (is.null(basis)) {
        vapply(which(ok), function(i) {
          as.data.frame(clubSandwich::linear_contrast(
            fit_cr, vcov = V, contrasts = grad[i, , drop = FALSE], test = "Satterthwaite"
          ))$df
        }, numeric(1))
      } else {
        .cr2_df(grad[ok, , drop = FALSE], fit_cr, V, basis)
      }
    }
    # The pretest effects compare risks standardized over two different
    # samples, the pretested and the unpretested participants, so the
    # sampling variance of each sample's standardization does not cancel as
    # it does for the treatment contrasts, which standardize both arms over
    # the same participants. It is added for them (issue #104), over the
    # participants each sample holds, in the model or not (issue #135).
    pretest_rows <- which(rep(.marginal_contrast_names(), length(scale)) %in%
                            .solomon_pretest_order & is.finite(estimate) & is.finite(std.error))
    if (length(pretest_rows)) {
      added <- .marginal_standardization_variance(
        fit, if (cr2) fit_cr else model, b, X, pretested, scale, count, pretest_rows, grad,
        in_model = set$in_model, cluster = cluster
      )
      v_delta <- std.error[pretest_rows]^2 + added$cross
      if (cr2) {
        df[pretest_rows] <- vapply(seq_along(pretest_rows), function(j) {
          .welch_df(c(v_delta[j], added$v[j]), c(df[pretest_rows[j]], added$df[j]))
        }, numeric(1))
      }
      std.error[pretest_rows] <- sqrt(v_delta + added$v)
    }
    q <- ifelse(is.finite(df), stats::qt(1 - alpha / 2, pmax(df, 1e-8)), stats::qnorm(1 - alpha / 2))
    ci <- cbind(estimate - q * std.error, estimate + q * std.error)
    list(std.error = std.error, df = df, ci = ci)
  }

  # A CR2 fit needs the cluster of every participant standardized over; a
  # participant outside the model may lack it (issue #135).
  set_clusters <- function(set) {
    if (!cr2 || is.null(set$index)) return(list(cluster = NULL, note = NULL))
    cluster <- fit$cluster[set$index]
    note <- .check_center_clusters(cluster, !set$in_model, stop = FALSE)
    list(cluster = cluster, note = note)
  }

  # Bootstrap refits (issue #154). Each resample is refitted on the
  # model's design, from its coefficients, as before. A refit that does
  # not converge in 25 iterations is fitted again on the same columns in
  # an orthonormal basis (.bootstrap_basis()). In exact arithmetic that is
  # the same fit, so it can succeed only where rounding error kept the
  # first from converging; a bootstrap whose refits all converge is
  # unchanged. A resample that fails again is counted by the reason of the
  # refit in the basis. A refit whose design is not of full rank at the
  # tolerance of qr() is not fitted again, and the resample counts as
  # failed. Nor is a refit that did not converge on a resample whose rows
  # make a term collinear by the check of fit_solomon_glm(), applied in
  # the model's working weights (.resample_collinear(), prepared once by
  # .bootstrap_rank_check()); it counts as failed with the reason "rank"
  # (see .bootstrap_basis() for why). `status` records how each resample
  # ended, "ok" or the reason it failed, and `condition_number` that of a
  # refit that did not converge.
  basis <- if (method == "bootstrap") .bootstrap_basis(X, b)
  rank_check <- NULL
  status <- character(if (method == "bootstrap") R else 0L)
  condition_number <- rep(NA_real_, length(status))
  refit_resample <- function(i, rows, keep, y_fit) {
    Xi <- X[rows, , drop = FALSE]
    refit <- .logistic_irls(Xi[keep, , drop = FALSE], y_fit, start = b)
    if (identical(refit$reason, "maxit") && !is.null(basis)) {
      if (is.null(rank_check)) rank_check <<- .bootstrap_rank_check(fit, model, X, in_model)
      if (.resample_collinear(rank_check, rows[keep])) {
        refit <- list(converged = FALSE, reason = "rank")
      } else {
        Xi <- basis$X[rows, , drop = FALSE]
        refit <- .logistic_irls(Xi[keep, , drop = FALSE], y_fit, start = basis$start)
      }
    }
    if (!is.null(refit$fitted.values) && .degenerate_risk(refit$fitted.values)) {
      status[i] <<- "risk"
    } else if (!isTRUE(refit$converged)) {
      status[i] <<- refit$reason
      if (!is.null(refit$kappa)) condition_number[i] <<- refit$kappa
    } else {
      status[i] <<- "ok"
      return(.marginal_all(refit$coefficients, Xi, pretested[rows], scale))
    }
    rep(NA_real_, length(estimate))
  }

  if (method == "bootstrap" && pretest_center == "analyzed") {
    cells <- split(seq_along(y), interaction(treat, pretested, drop = TRUE))
    sizes <- lengths(cells)
    ends <- cumsum(sizes)
    run <- function() {
      t(vapply(seq_len(R), function(i) {
        idx <- unlist(lapply(cells, function(rows) rows[sample.int(length(rows), replace = TRUE)]),
                      use.names = FALSE)
        # Resampling is within cells, so the cells occupy consecutive blocks.
        events <- diff(c(0, cumsum(y[idx])[ends]))
        if (any(events == 0 | events == sizes)) {
          status[i] <<- "separated"
          return(rep(NA_real_, length(estimate)))
        }
        refit_resample(i, idx, TRUE, y[idx])
      }, numeric(length(estimate))))
    }
  } else if (method == "bootstrap") {
    # Participants are resampled within the four Solomon cells. A pretested
    # participant who is standardized over without being in the model (no
    # outcome) is resampled with the cell of their condition, and those with
    # no recorded condition as a stratum of their own, after the cells. Each
    # resample is refitted on its members who are in the model and
    # standardized over all of them, so the number with an outcome varies
    # between resamples as it does between samples (issue #135).
    y_all <- numeric(nrow(X))
    y_all[in_model] <- y
    cells <- split(seq_len(nrow(X)), interaction(treat, pretested, drop = TRUE))
    solomon_cells <- seq_along(cells)
    if (anyNA(treat)) cells <- c(cells, list(which(is.na(treat))))
    ends <- cumsum(lengths(cells))
    run <- function() {
      t(vapply(seq_len(R), function(i) {
        idx <- unlist(lapply(cells, function(members) {
          members[sample.int(length(members), replace = TRUE)]
        }), use.names = FALSE)
        # Resampling is within cells, so the cells occupy consecutive blocks.
        refitted <- in_model[idx]
        events <- diff(c(0, cumsum(y_all[idx])[ends]))[solomon_cells]
        outcomes <- diff(c(0, cumsum(refitted)[ends]))[solomon_cells]
        if (any(events == 0 | events == outcomes)) {
          status[i] <<- "separated"
          return(rep(NA_real_, length(estimate)))
        }
        refit_resample(i, idx, refitted, y_all[idx][refitted])
      }, numeric(length(estimate))))
    }
  }

  if (method == "bootstrap") {
    boot <- if (is.null(seed)) run() else withr::with_seed(seed, run())
    # A refit whose every contrast is undefined counts as failed, as before.
    status[status == "ok" & apply(boot, 1, function(row) all(is.na(row)))] <- "risk"
    failed <- status != "ok"
    failures <- sum(failed)
    boot <- boot[!failed, , drop = FALSE]
    std.error <- apply(boot, 2, stats::sd, na.rm = TRUE)
    ci <- t(apply(boot, 2, stats::quantile, probs = c(alpha / 2, 1 - alpha / 2),
                  na.rm = TRUE, names = FALSE))
    if (failures > 0.1 * R) {
      ci[] <- NA_real_
      warning(structure(
        class = c("solomonR_bootstrap_warning", "warning", "condition"),
        list(message = .bootstrap_failure_message(status[failed], R, condition_number[failed]),
             call = NULL)
      ))
    }
  } else {
    clusters <- set_clusters(set)
    if (!is.null(clusters$note)) stop(clusters$note, call. = FALSE)
    inference <- delta_at(set, estimate, clusters$cluster)
    std.error <- inference$std.error
    df <- inference$df
    ci <- inference$ci
  }

  # The test of each contrast is on the analysis scale: the difference, or
  # the log of the ratio.
  scale_col <- rep(scale, each = n_contrasts)
  ratio <- scale_col != "difference"
  report <- function(x) ifelse(ratio, exp(x), x)
  table_at <- function(estimate, std.error, df, ci) {
    statistic <- estimate / std.error
    p.value <- 2 * stats::pt(-abs(statistic), df)
    .effects_table(data.frame(
      scale = unname((if (count) .marginal_rate_scales else .marginal_scales)[scale_col]),
      contrast = rep(.marginal_contrast_names(), length(scale)),
      estimate = report(estimate),
      std.error = std.error,
      statistic = statistic,
      df = df,
      p.value = p.value,
      conf.low = report(ci[, 1]),
      conf.high = report(ci[, 2]),
      stringsAsFactors = FALSE
    ), keys = "scale")
  }
  effects <- table_at(estimate, std.error, df, ci)

  if (cr2) {
    small_df <- is.finite(df) & df < 4
    if (any(small_df)) {
      .warn_cr2_small_df(paste(effects$scale[small_df], effects$contrast[small_df], sep = ": "),
                         df[small_df])
    }
  }

  # Every contrast standardized over both sets of participants (issue #135):
  # on the risk or rate scale the pretested cells move with the set, and so
  # do Treatment | pretested, the average treatment effect, and Pretest x
  # Treatment, as well as the pretest effects. With the delta method the
  # other set has full inference; with the bootstrap only its estimates, so
  # that the resamples are those of the set used. When the two sets are the
  # same participants, the rows are copied.
  sensitivity <- NULL
  pretest_adjusted <- "pre_obs" %in% colnames(X)
  if (pretest_adjusted) {
    other_center <- .other_center(pretest_center)
    other_effects <- effects
    note <- NULL
    if (standardized[["observed"]] > standardized[["analyzed"]]) {
      other_set <- if (other_center == "analyzed") model_set() else observed_set
      other_estimate <- estimate_at(other_set)$estimate
      missing <- rep(NA_real_, length(other_estimate))
      inference <- list(std.error = missing, df = missing, ci = cbind(missing, missing))
      if (method == "delta") {
        clusters <- set_clusters(other_set)
        note <- clusters$note
        if (is.null(note)) {
          inference <- suppressWarnings(delta_at(other_set, other_estimate, clusters$cluster))
        } else {
          note <- paste(note, "The standard errors of that standardization are not given.")
        }
      }
      other_effects <- table_at(other_estimate, inference$std.error, inference$df, inference$ci)
    }
    effects$center <- pretest_center
    other_effects$center <- other_center
    both <- if (pretest_center == "analyzed") rbind(effects, other_effects) else
      rbind(other_effects, effects)
    effects$center <- NULL
    both <- .effects_table(both, keys = c("scale", "center"))
    analyzed <- both[both$center == "analyzed", , drop = FALSE]
    observed <- both[both$center == "observed", , drop = FALSE]
    difference <- data.frame(
      scale = analyzed$scale, contrast = analyzed$contrast,
      analyzed = analyzed$estimate, observed = observed$estimate,
      stringsAsFactors = FALSE
    )
    difference$difference <- difference$analyzed - difference$observed
    sensitivity <- list(
      effects = both, difference = difference,
      differ = standardized[["observed"]] > standardized[["analyzed"]], note = note
    )
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
      pretest_adjusted = pretest_adjusted,
      pretest_center = pretest_center,
      standardized = standardized,
      pretest_sensitivity = sensitivity,
      attrition = .fit_attrition(fit)
    ),
    class = "solomon_marginal"
  )
}


#' @export
print.solomon_marginal <- function(x, digits = 3, ...) {
  count <- identical(x$outcome, "count")
  what <- if (count) "Rates" else "Risks"
  cat("Marginal Solomon contrasts for a", if (count) "count" else "binary", "outcome\n")
  cat(if (x$pretest_adjusted && identical(x$pretest_center, "observed")) {
    paste(what, "standardized over the pretest among every pretested participant with a",
          "pretest (pretest_center = \"observed\").\n")
  } else if (x$pretest_adjusted) {
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
  # The two standardizations of the pretested cells when outcomes are
  # missing, or else the attrition line whenever an outcome is missing
  # (issue #135).
  block <- .marginal_block(x, digits)
  if (length(block)) {
    cat("\n", paste(block, collapse = "\n"), "\n", sep = "")
  } else {
    line <- .attrition_line(x$attrition, "Outcomes missing")
    if (length(line)) cat("\n", .wrap_lines(line), "\n", sep = "")
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
#' @return An object of class `solomon_fisher`, a list with:
#'   - `tests`: one row for each pretest condition and one for both
#'     combined, in the columns `condition`, `events_treatment`,
#'     `n_treatment`, `events_control`, `n_control`, `risk_treatment`,
#'     `risk_control`, `chisq` (the uncorrected Pearson chi-square),
#'     `chisq_p`, and `fisher_p` (Fisher's exact p-value).
#'   - `sensitization`: the verdict of the historical rule, `TRUE` or
#'     `FALSE`.
#'   - `alpha`: the significance level of the rule.
#'
#'   The tests are of two-by-two tables, not of the Solomon contrasts, so
#'   the result has no effects table; see [solomon_output].
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
