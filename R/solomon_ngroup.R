# Solomon N-group designs (issue #45): k treatments and a control, each with
# and without a pretest, giving 2(k + 1) groups (Steyn, 2009). The joint
# model of fit_solomon_glm() is fitted with one indicator per treatment, and
# the four Solomon contrasts are estimated for each comparison, with
# omnibus Wald tests and a multiplicity adjustment within each family of
# comparisons (Holm, 1979, by default).

# Classify the treatment argument of a fitting function.
#
# A 0/1 (or logical) indicator, or a factor or character vector with two
# conditions, describes a four-group design: `k = 1` and `treat` is returned
# as a 0/1 integer with NA kept. A factor or character vector with three or
# more conditions describes an N-group design: `k` treatments and the
# control named by `control`, with `condition` a factor whose first level
# is the control.
.solomon_conditions <- function(treat, control = NULL) {

  if (is.factor(treat) || is.character(treat)) {
    f <- droplevels(as.factor(treat))
    lev <- levels(f)

    if (length(lev) < 2L) {
      stop(
        "`treat` must contain at least two conditions; it contains ",
        if (length(lev)) paste0("only ", dQuote(lev, FALSE)) else "none",
        ".",
        call. = FALSE
      )
    }

    if (is.null(control)) {
      stop(
        "`treat` has the conditions ", paste(dQuote(lev, FALSE), collapse = ", "),
        ". Name the control condition with `control`.",
        call. = FALSE
      )
    }

    if (!is.atomic(control) || length(control) != 1L || is.na(control) ||
        !as.character(control) %in% lev) {
      stop(
        "`control` must be one of the conditions of `treat`: ",
        paste(dQuote(lev, FALSE), collapse = ", "), ".",
        call. = FALSE
      )
    }

    control <- as.character(control)
    treatments <- setdiff(lev, control)

    if (length(treatments) == 1L) {
      return(list(
        k = 1L,
        treat = as.integer(f != control),
        control = control,
        treatments = treatments
      ))
    }

    return(list(
      k = length(treatments),
      condition = factor(as.character(f), levels = c(control, treatments)),
      control = control,
      treatments = treatments
    ))
  }

  if (!is.null(control) &&
      !(is.atomic(control) && length(control) == 1L && !is.na(control) &&
        isTRUE(suppressWarnings(as.numeric(control)) == 0))) {
    stop(
      "With a 0/1 `treat`, the control condition is 0; leave `control` unset.",
      call. = FALSE
    )
  }

  if (is.numeric(treat) && length(unique(stats::na.omit(treat))) > 2L) {
    stop(
      "`treat` has more than two values. For a design with several ",
      "treatments, give `treat` as a factor or character vector and name ",
      "the control condition with `control`.",
      call. = FALSE
    )
  }

  list(k = 1L, treat = .solomon_indicator(treat, "treat"))
}


# Stop with a clear message when a function that analyzes four-group
# designs receives a treatment with more than two conditions.
.stop_ngroup_unsupported <- function(treat, fun) {
  if ((is.factor(treat) || is.character(treat)) &&
      length(unique(stats::na.omit(as.character(treat)))) > 2L) {
    stop(structure(
      class = c("solomonR_ngroup_unsupported", "error", "condition"),
      list(
        message = paste0(
          "`", fun, "()` analyzes a Solomon four-group design: one treatment ",
          "and a control, each with and without a pretest. `treat` has ",
          length(unique(stats::na.omit(as.character(treat)))), " conditions. ",
          "Use `fit_solomon_glm(..., control = )` for designs with several ",
          "treatments, or analyze one treatment and the control at a time."
        ),
        call = NULL
      )
    ))
  }
  invisible(NULL)
}


# Comparison weights over the conditions (control first).
#
# `contrasts` is "control" (each treatment against the control), "pairwise"
# (every pair of conditions), or a named list of weight vectors named by
# condition. Each set of weights must sum to zero. Returns a matrix with one
# row per comparison and one column per condition.
.ngroup_weights <- function(contrasts, control, treatments) {

  lev <- c(control, treatments)
  unit <- function(plus, minus) {
    w <- stats::setNames(numeric(length(lev)), lev)
    w[plus] <- 1
    w[minus] <- -1
    w
  }

  if (is.character(contrasts)) {
    contrasts <- match.arg(contrasts, c("control", "pairwise"))
    rows <- lapply(treatments, function(t) unit(t, control))
    names(rows) <- paste(treatments, "vs", control)
    if (contrasts == "pairwise") {
      pairs <- utils::combn(treatments, 2L, simplify = FALSE)
      extra <- lapply(pairs, function(p) unit(p[1], p[2]))
      names(extra) <- vapply(pairs, function(p) paste(p[1], "vs", p[2]), character(1))
      rows <- c(rows, extra)
    }
    return(do.call(rbind, rows))
  }

  if (!is.list(contrasts) || !length(contrasts) || is.null(names(contrasts)) ||
      any(!nzchar(names(contrasts))) || anyDuplicated(names(contrasts))) {
    stop(
      "`contrasts` must be \"control\", \"pairwise\", or a list of weight ",
      "vectors with unique names, one per comparison.",
      call. = FALSE
    )
  }

  rows <- lapply(names(contrasts), function(nm) {
    w <- contrasts[[nm]]
    if (!is.numeric(w) || is.null(names(w)) || any(!is.finite(w))) {
      stop(
        "Contrast ", dQuote(nm, FALSE), " must be a numeric vector named by ",
        "condition, for example c(", deparse(treatments[1]), " = 1, ",
        deparse(control), " = -1).",
        call. = FALSE
      )
    }
    unknown <- setdiff(names(w), lev)
    if (length(unknown)) {
      stop(
        "Contrast ", dQuote(nm, FALSE), " names conditions not in `treat`: ",
        paste(dQuote(unknown, FALSE), collapse = ", "), ".",
        call. = FALSE
      )
    }
    if (anyDuplicated(names(w))) {
      stop("Contrast ", dQuote(nm, FALSE), " names a condition twice.", call. = FALSE)
    }
    full <- stats::setNames(numeric(length(lev)), lev)
    full[names(w)] <- w
    if (all(full == 0)) {
      stop("Contrast ", dQuote(nm, FALSE), " has no nonzero weights.", call. = FALSE)
    }
    if (abs(sum(full)) > 1e-8) {
      stop(
        "The weights of contrast ", dQuote(nm, FALSE), " sum to ",
        format(sum(full)), "; the weights of a comparison must sum to zero.",
        call. = FALSE
      )
    }
    full
  })
  names(rows) <- names(contrasts)
  do.call(rbind, rows)
}


# The names of the Solomon groups flagged by `sel`, as "Pretested, RP":
# `counts` is the table of analyzable participants, with one row per
# condition and the columns "1" (pretested) and "0" (unpretested).
.group_names <- function(counts, sel) {
  at <- which(sel, arr.ind = TRUE)
  sprintf(
    "%s, %s", c("1" = "Pretested", "0" = "Unpretested")[colnames(counts)[at[, 2]]],
    rownames(counts)[at[, 1]]
  )
}

# Stop when a coefficient of the fitted model is inestimable: an empty
# group, or a term collinear with the others. `empty` names the groups
# without an analyzable participant (.group_names()). Used by the four-group
# and the several-treatment paths of fit_solomon_glm(). With `noise`, the
# rounding error of each entry of the model matrix (.representation_error()),
# a term that is collinear with the others but for that error is collinear
# too (.collinear_terms()); `covariates` names the covariates' terms, which
# are checked last.
.stop_inestimable <- function(fit, empty, noise = NULL, covariates = character()) {
  aliased <- names(stats::coef(fit))[is.na(stats::coef(fit))]
  if (!length(aliased) && !is.null(noise)) aliased <- .collinear_terms(fit, noise, covariates)
  if (!length(aliased)) return(invisible(NULL))
  if (length(empty)) {
    stop(
      "The model cannot be estimated because these groups have no analyzable ",
      "participants: ", paste(empty, collapse = "; "), ".",
      call. = FALSE
    )
  }
  aliased[aliased == "pre_obs"] <- "pre_obs (the pretest, `y_pre`)"
  stop(
    "The model cannot be estimated because these terms are collinear with ",
    "other terms of the model: ", paste(aliased, collapse = ", "),
    ". A covariate or pretest that is constant, or that repeats another ",
    "term, causes this.",
    call. = FALSE
  )
}

# The rounding error that each entry of the model matrix of the
# standardized model `fit` carries from the values the model was given, for
# .collinear_terms() (issue #137). `given` is the model as given, `scaling`
# comes from .covariate_scaling(), and `center` is the mean the pretest was
# centered at, NA without a pretest. A value x is stored with an error of up
# to |x| eps / 2, with eps = .Machine$double.eps, and a covariate computed
# from others, such as a time of testing that is the time of enrolment plus
# a delay, carries the errors of that computation, of the same order.
# Standardizing a covariate divides that error by the covariate's standard
# deviation and adds an error of its own, of up to the standardized value
# times eps / 2, and dividing the pretest by a power of two is exact. The
# error of entry (i, j) is therefore taken to be
#   eps (|x_ij| / s_j + |z_ij|),
# with x_ij the value as given, s_j the number column j was divided by (1
# for a column that is not standardized), and z_ij the value in `fit`. The
# pretest enters centered, but it was stored before centering, so its x_ij
# is the pretest as given, the centered value plus the center, in the
# pretested rows; the 0 of the unpretested rows is exact. The intercept,
# the treatment and pretest indicators, their interactions, and the
# indicators of factors hold 0s and 1s, exactly, and get eps (|x| + |z|)
# too. Returns a matrix of the dimensions of the model matrix.
.representation_error <- function(fit, given, scaling, center = NA_real_) {
  Z <- stats::model.matrix(fit)
  X <- stats::model.matrix(given)
  s <- rep(1, ncol(Z))
  if (length(scaling$means)) {
    columns <- .scaled_columns(fit, scaling)
    s[columns] <- scaling$scales
    if ("pre_obs" %in% names(columns) && is.finite(center)) {
      j <- columns[["pre_obs"]]
      X[, j] <- X[, j] + center * X[, match("pretested", colnames(X))]
    }
  }
  .Machine$double.eps * (sweep(abs(X), 2L, s, "/") + abs(Z))
}

# The terms of the standardized model `fit` that are collinear with the
# other terms but for the rounding error of the values the model was given
# (issue #137), from `noise`, the error of each entry of its model matrix
# (.representation_error()).
#
# Terms that are exactly collinear in the data as given, sum_j c_j x_j = 0
# with the intercept among the x_j, satisfy it once stored only up to their
# rounding errors: sum_j c_j x_j = r, with |r_i| <= sum_j |c_j| e_ij for the
# errors e_ij of the entries. Standardizing changes the units of the terms,
# not the relation. Regressed on the other terms, as glm()'s QR
# decomposition does in the metric of the fit's weights (the working
# weights, times the prior weights), a term k of the relation then has a
# residual of at most || r || / |c_k|, and its coefficients are those of
# the relation, b_kj = -c_j / c_k. So its residual is at most
#   B_k = || e_k + sum_j |b_kj| e_j ||,
# the rounding error the regression carries into term k from the terms it
# is regressed on: in the units of the terms as given, with coefficients
# beta_kj, entry i is eps (|x_ik| + sum_j |beta_kj| |x_ij|) / s_k, which
# grows with the coefficients and so with the spreads of the far terms
# relative to term k's own. A tolerance relative to each term's own spread,
# the same for every term, missed a time of testing that is the time of
# enrolment plus a delay, with the times as dates in seconds: the error of
# the times reached the delay amplified by their spreads relative to its
# spread, 1,200 times its own, and the fit went on to meaningless standard
# errors. A single tolerance large enough to cover that amplification
# would instead stop nearly collinear covariates near zero whenever another
# covariate lies far from zero, though glm() and master fit them.
#
# Each term is therefore compared with its own bound: it is taken to be
# collinear with the terms before it when its residual, in the metric of
# the weights and with each column but the intercept centered at its
# weighted mean, is at most 4 B_k, or at most 1e-11 of its norm, the
# tolerance of glm(). The bound holds with the coefficients of an exact
# relation, and the regression's coefficients equal them to within the
# relation's own rounding error; the margin of 4 covers errors several
# times eps / 2 per entry, as computed covariates carry, and the rounding
# of the decomposition. A term that is not collinear is stopped only if its
# residual is within 4 B_k: near zero B_k is of the order of eps, so that
# glm()'s tolerance decides, as in master, and a covariate far from zero
# that is unrelated to the term enters its regression with a coefficient
# near zero and adds nothing to its bound.
#
# The terms are checked in turn, each regressed on those before it, with
# the covariates (the terms named in `covariates`) after the treatment and
# pretest terms, so that a relation is found at its last term in that order
# whatever the order of the covariates; a term found collinear is left out
# and the rest are checked again. Returns the names of the collinear terms'
# columns, in the order of the model, or an empty vector.
.collinear_terms <- function(fit, noise, covariates = character(), margin = 4, floor = 1e-11) {
  X <- stats::model.matrix(fit)
  w <- fit$weights
  if (is.null(w) || length(w) != nrow(X)) w <- rep(1, nrow(X))
  labels <- attr(stats::terms(fit), "term.labels")
  last <- attr(X, "assign") %in% match(covariates, labels)
  order <- c(which(!last), which(last))
  keep <- w > 0
  A <- X[keep, order, drop = FALSE]
  E <- noise[keep, order, drop = FALSE]
  w <- w[keep]
  j <- order != match(0L, attr(X, "assign"))
  A[, j] <- sweep(A[, j, drop = FALSE], 2L, colSums(w * A[, j, drop = FALSE]) / sum(w))
  A <- sqrt(w) * A
  E <- sqrt(w) * E
  # The first term, in this order, collinear with the terms before it; NA
  # for none. qr() with a tolerance of 0 keeps the columns in order.
  first <- function(A, E) {
    R <- qr.R(qr(A, tol = 0))
    for (k in seq_len(ncol(A))[-1L]) {
      before <- seq_len(k - 1L)
      b <- backsolve(R[before, before, drop = FALSE], R[before, k])
      bound <- E[, k] + drop(E[, before, drop = FALSE] %*% abs(b))
      residual <- abs(R[k, k])
      if (!is.finite(residual) || !all(is.finite(bound)) ||
          residual <= max(margin * sqrt(sum(bound^2)), floor * sqrt(sum(A[, k]^2)))) {
        return(k)
      }
    }
    NA_integer_
  }
  active <- seq_len(ncol(A))
  found <- integer()
  repeat {
    k <- first(A[, active, drop = FALSE], E[, active, drop = FALSE])
    if (is.na(k)) break
    found <- c(found, active[k])
    active <- active[-k]
  }
  colnames(X)[sort(order[found])]
}

# Fit the joint model of an N-group design. Called by fit_solomon_glm()
# after its arguments have been resolved; `cond` comes from
# .solomon_conditions(). `pretest_center` is the center of the pretest
# (issue #135), as in the four-group model.
.fit_solomon_ngroup <- function(y_post, cond, pretested, y_pre, covariates,
                                robust, cluster, family, negbin, conf_level,
                                exposure, contrasts, adjust, call,
                                pretest_center = "analyzed") {

  pretested <- .solomon_indicator(pretested, "pretested")

  .solomon_check_lengths(
    y_post = y_post,
    treat = cond$condition,
    pretested = pretested,
    y_pre = y_pre,
    covariates = covariates,
    cluster = cluster,
    exposure = exposure
  )

  k <- cond$k
  lev <- c(cond$control, cond$treatments)
  terms <- make.names(paste0("treat_", cond$treatments), unique = TRUE)
  weights <- .ngroup_weights(contrasts, cond$control, cond$treatments)

  df <- data.frame(
    y = y_post,
    condition = cond$condition,
    pretested = pretested
  )
  for (j in seq_len(k)) {
    df[[terms[j]]] <- as.integer(cond$condition == cond$treatments[j])
  }

  if (!is.null(y_pre)) {

    n_incidental <- sum(pretested == 1L & is.na(y_pre), na.rm = TRUE)

    if (n_incidental > 0L) {
      warning(
        n_incidental, " pretested participant(s) have missing pretest scores ",
        "and are excluded from the model. Incidental missingness is not imputed.",
        call. = FALSE
      )
    }

    if (any(pretested == 0L & !is.na(y_pre), na.rm = TRUE)) {
      warning(
        "Observed `y_pre` values were supplied for unpretested ",
        "participants; these values are ignored.",
        call. = FALSE
      )
    }

    df$pre_obs <- ifelse(df$pretested == 1, y_pre, 0)
  }

  if (!is.null(covariates)) {
    # `log_exposure` is reserved with or without `exposure`: report_solomon()
    # and the clustered permutation test read a column of that name in the
    # fit's data as the exposure offset.
    clash <- intersect(names(covariates), c(names(df), "log_exposure"))
    if (length(clash)) {
      stop(
        "Rename these covariates, whose names the model uses: ",
        paste(clash, collapse = ", "), ".",
        call. = FALSE
      )
    }
    df <- cbind(df, covariates)
  }

  rhs <- c(
    sprintf("(%s) * pretested", paste(terms, collapse = " + ")),
    if (!is.null(y_pre)) "pre_obs" else NULL,
    if (!is.null(covariates)) .formula_names(names(covariates)) else NULL
  )
  if (!is.null(exposure)) {
    if (!is.numeric(exposure) || any(exposure <= 0, na.rm = TRUE)) {
      stop("`exposure` must contain positive numbers.", call. = FALSE)
    }
    if (!identical(family$link, "log")) {
      stop("`exposure` requires a log-link family, such as poisson().", call. = FALSE)
    }
    df$log_exposure <- log(exposure)
    rhs <- c(rhs, "offset(log_exposure)")
  }
  fml <- stats::as.formula(paste("y ~", paste(rhs, collapse = " + ")))

  # Center the pretest at its mean among the pretested participants in the
  # model, as in the four-group model (issue #104), or at the mean of every
  # observed pretest among the pretested participants (issue #135). Either
  # center pools the pretested participants of every condition.
  pretest_mean <- NA_real_
  centers <- c(analyzed = NA_real_, observed = NA_real_)
  if (!is.null(y_pre)) {
    centers <- .pretest_centers(df)
    pretest_mean <- centers[[pretest_center]]
    df$pre_obs <- ifelse(df$pretested == 1, df$pre_obs - pretest_mean, 0)
  }

  # The model as given, and the model refitted with the numeric covariates
  # standardized and the pretest scaled, which the inference comes from and
  # which gives the warnings, as in the four-group model (issue #137).
  standardize <- .standardizes(df, names(covariates))
  fit <- .fit_model(fml, df, family, negbin, quiet = standardize)
  used <- rep(TRUE, nrow(df))
  if (!is.null(fit$na.action)) used[fit$na.action] <- FALSE
  fit_given <- fit
  scaling <- .covariate_scaling(df, names(covariates), used)
  if (standardize) fit <- .fit_model(fml, scaling$data, family, negbin)
  if (negbin) family <- stats::family(fit)

  if (!is.null(y_pre) && !stats::family(fit)$link %in% c("identity", "log")) {
    .warn_noncollapsible(stats::family(fit)$link, ngroup = TRUE)
  }

  # Analyzable participants (the rows the model kept) in each of the
  # 2(k + 1) groups.
  counts <- table(
    factor(df$condition[used], levels = lev),
    factor(df$pretested[used], levels = c(1L, 0L))
  )

  # An empty group, or a term collinear with the others, leaves a
  # coefficient inestimable, in the model as given or standardized, the
  # latter checked as well for terms collinear but for the rounding error of
  # the values given, as in the four-group model (.collinear_terms()).
  .stop_inestimable(fit_given, .group_names(counts, counts == 0))
  .stop_inestimable(fit, .group_names(counts, counts == 0),
                    noise = .representation_error(fit, fit_given, scaling, pretest_mean),
                    covariates = .formula_names(names(covariates)))

  # A group with one participant is fitted exactly: its hat value is 1, and
  # HC3 divides the squared residual by (1 - h)^2.
  if (robust == "HC3" && any(counts == 1)) {
    stop(
      "HC3 standard errors are undefined because these groups have a single ",
      "analyzable participant: ", paste(.group_names(counts, counts == 1), collapse = "; "), ". ",
      "Use `robust = \"none\"` or leave the condition out; see validate_solomon().",
      call. = FALSE
    )
  }

  # Robust covariance, as in the four-group model.
  fit_cr <- fit
  n_clusters <- NULL
  if (robust == "HC3") {
    vcovM <- sandwich::vcovHC(fit, type = "HC3")
  } else if (robust == "CR2") {
    if (is.null(cluster)) stop("CR2 requested but 'cluster' is NULL.", call. = FALSE)

    cluster_fit <- cluster[used]
    if (anyNA(cluster_fit)) {
      stop("`cluster` is missing for participants included in the model.", call. = FALSE)
    }
    # The clusters the CR2 covariance is computed from; `cluster` keeps one
    # value per input row.
    n_clusters <- length(unique(cluster_fit))
    single <- .cluster_structure(
      df$condition[used], df$pretested[used], cluster_fit,
      cells = .solomon_cells(lev)
    )$single_cluster
    if (length(single)) {
      .stop_confounded_clusters(single)
    }

    fit_cr <- stats::glm(fml, data = scaling$data[used, , drop = FALSE], family = family)
    vcovM <- clubSandwich::vcovCR(fit_cr, cluster = cluster_fit, type = "CR2")
  } else {
    vcovM <- stats::vcov(fit)
  }

  # The sampling variance of the mean pretest at which the pretest effects
  # are evaluated, as in the four-group model (issues #104 and #135), and the
  # other center's, when the two differ, for `pretest_sensitivity`, in the
  # units of the scaled pretest of `fit`.
  center <- NULL
  center_other <- NULL
  center_note <- NULL
  parts <- NULL
  differ <- FALSE
  if (!is.null(y_pre)) {
    both <- .center_terms_both(scaling$data, y_pre, used, if (robust == "CR2") fit_cr else fit,
                               robust, cluster, if (robust == "CR2") cluster_fit, pretest_center)
    center <- both$center
    center_other <- both$other
    center_note <- both$note
    parts <- both$parts
    differ <- both$differ
  }

  dispersion_fixed <- .fixed_dispersion(stats::family(fit))
  df_model <- if (dispersion_fixed) Inf else stats::df.residual(fit)

  coefficients <- .coefficient_table(fit, fit_cr, vcovM, scaling, robust, df_model, conf_level)

  cf <- stats::coef(fit)
  cn <- names(cf)

  # The pretest effects also carry the variance of the mean pretest at which
  # they are evaluated (.center_adjusted()); `terms` is that mean's
  # .center_terms(), of the center used or of the other center.
  lin_contrast <- function(L, terms = center) {
    Lm <- matrix(L, nrow = 1)
    est <- as.numeric(Lm %*% cf)
    dfc <- if (robust == "CR2") {
      as.data.frame(
        clubSandwich::linear_contrast(
          fit_cr,
          vcov = vcovM,
          contrasts = Lm,
          test = "Satterthwaite"
        )
      )$df
    } else {
      df_model
    }
    adj <- .center_adjusted(L, cf, vcovM, terms, robust, dfc)
    se  <- sqrt(adj$variance)
    dfc <- adj$df
    z   <- est / se
    p   <- 2 * stats::pt(-abs(z), df = dfc)
    ci  <- .wald_ci(est, se, dfc, conf_level)
    c(estimate = est, std.error = se, statistic = z, p.value = p, df = dfc,
      conf.low = unname(ci[, "conf.low"]), conf.high = unname(ci[, "conf.high"]),
      added = adj$added)
  }

  # The four Solomon contrasts of a comparison with weights w over the
  # conditions. The control's coefficients are zero, so only the treatment
  # weights enter:
  #   unpretested: sum_j w_j b_j
  #   sensitization: sum_j w_j b_{j:pretested}
  #   pretested: the two together; average: unpretested + sensitization / 2.
  contrast_rows <- function(w) {
    main <- stats::setNames(numeric(length(cn)), cn)
    inter <- main
    main[terms] <- w[cond$treatments]
    inter[paste0(terms, ":pretested")] <- w[cond$treatments]
    list(
      "ATE (avg over pretest)" = main + 0.5 * inter,
      "Pretest x Treatment" = inter,
      "Treatment | pretested" = main + inter,
      "Treatment | unpretested" = main
    )
  }
  types <- c("ATE (avg over pretest)", "Pretest x Treatment",
             "Treatment | pretested", "Treatment | unpretested")

  conventional <- robust == "none"
  L_by_comparison <- lapply(seq_len(nrow(weights)), function(r) contrast_rows(weights[r, ]))

  effect_row <- function(L, comparison, contrast) {
    est <- lin_contrast(L)
    r2 <- contrast_r2_ci(fit, L, vcovM, conf_level, conventional,
                         added_variance = est[["added"]])
    data.frame(
      comparison = comparison,
      contrast = contrast,
      estimate = unname(est["estimate"]),
      std.error = unname(est["std.error"]),
      statistic = unname(est["statistic"]),
      df = unname(est["df"]),
      p.value = unname(est["p.value"]),
      conf.low = unname(est["conf.low"]),
      conf.high = unname(est["conf.high"]),
      r2 = r2$r2,
      r2_lo = r2$r2_lo,
      r2_hi = r2$r2_hi,
      stringsAsFactors = FALSE
    )
  }

  effects <- do.call(rbind, lapply(types, function(type) {
    rows <- lapply(seq_len(nrow(weights)), function(r) {
      effect_row(L_by_comparison[[r]][[type]], rownames(weights)[r], type)
    })
    block <- do.call(rbind, rows)
    block$p.adjusted <- stats::p.adjust(block$p.value, method = adjust)
    block
  }))

  # The pretest effect in each condition (issue #104): pretested minus
  # unpretested participants with the centered pretest at zero. The
  # control's is the pretesting coefficient, and a treatment's adds that
  # treatment's interaction with pretesting. The main effect weights the
  # k + 1 conditions equally. The treatments' pretest effects are adjusted
  # across the treatments, like each contrast across the comparisons.
  base <- stats::setNames(numeric(length(cn)), cn)
  base["pretested"] <- 1
  L_pretest <- c(
    list(base),
    lapply(terms, function(term) {
      L <- base
      L[paste0(term, ":pretested")] <- 1
      L
    }),
    list(replace(base, paste0(terms, ":pretested"), 1 / (k + 1)))
  )
  pretest_block <- do.call(rbind, Map(
    effect_row, L_pretest,
    c(cond$control, cond$treatments, "All conditions"),
    c(.solomon_pretest_order[1], rep(.solomon_pretest_order[2], k), .solomon_pretest_order[3])
  ))
  treated <- pretest_block$contrast == .solomon_pretest_order[2]
  pretest_block$p.adjusted <- pretest_block$p.value
  pretest_block$p.adjusted[treated] <- stats::p.adjust(pretest_block$p.value[treated],
                                                       method = adjust)
  effects <- rbind(effects, pretest_block)
  # The adjusted p-values follow the columns that every effects table has
  # (issue #110).
  effects <- effects[, c("comparison", .solomon_effect_columns, "p.adjusted",
                         "r2", "r2_lo", "r2_hi")]
  effects <- .effects_table(effects, keys = "comparison")

  # The pretest effects at both centers (issue #135), computed for the
  # other center by shifting the contrasts, as in the four-group model: the
  # shift in the units of the scaled pretest, the slope in those of the
  # pretest as given.
  sensitivity <- NULL
  if (!is.null(y_pre)) {
    pretest_labels <- c(cond$control, cond$treatments, "All conditions")
    pretest_types <- c(.solomon_pretest_order[1], rep(.solomon_pretest_order[2], k),
                       .solomon_pretest_order[3])
    columns <- c("comparison", .solomon_effect_columns, "p.adjusted")
    used_rows <- effects[effects$contrast %in% .solomon_pretest_order, columns, drop = FALSE]
    other_rows <- NULL
    if (differ) {
      other_rows <- .other_center_rows(
        L_pretest, pretest_types,
        (centers[[.other_center(pretest_center)]] - pretest_mean) / .pretest_scale(scaling),
        function(L) lin_contrast(L, center_other), comparison = pretest_labels,
        no_se = !is.null(center_note)
      )
      other_treated <- other_rows$contrast == .solomon_pretest_order[2]
      other_rows$p.adjusted <- other_rows$p.value
      other_rows$p.adjusted[other_treated] <- stats::p.adjust(other_rows$p.value[other_treated],
                                                              method = adjust)
      other_rows <- other_rows[, columns]
    }
    identity <- identical(stats::family(fit)$link, "identity")
    slope <- cf[["pre_obs"]] / .pretest_scale(scaling)
    sensitivity <- .pretest_sensitivity(
      used_rows, other_rows, pretest_center, "comparison", parts, centers, slope, identity,
      .pretest_attrition(as.character(df$condition), lev, df$pretested, df$y, used, slope,
                         centers, identity),
      note = center_note
    )
  }

  if (robust == "CR2") {
    small_df <- is.finite(effects$df) & effects$df < 4
    if (any(small_df)) {
      .warn_cr2_small_df(
        paste0(effects$comparison[small_df], ": ", effects$contrast[small_df]),
        effects$df[small_df]
      )
    }
  }

  # Omnibus Wald tests that all conditions are equal on each contrast type:
  # k constraints, one per treatment against the control.
  basis <- lapply(cond$treatments, function(t) {
    w <- stats::setNames(numeric(length(lev)), lev)
    w[t] <- 1
    w[cond$control] <- -1
    contrast_rows(w)
  })
  omnibus_labels <- c(
    "ATE (avg over pretest)" = "Condition (avg over pretest)",
    "Pretest x Treatment" = "Pretest x Condition",
    "Treatment | pretested" = "Condition | pretested",
    "Treatment | unpretested" = "Condition | unpretested"
  )
  omnibus <- do.call(rbind, lapply(types, function(type) {
    L <- do.call(rbind, lapply(basis, `[[`, type))
    b <- as.numeric(L %*% cf)
    W <- as.numeric(t(b) %*% solve(L %*% vcovM %*% t(L), b))
    if (robust == "CR2") {
      wt <- clubSandwich::Wald_test(fit_cr, constraints = L, vcov = vcovM, test = "HTZ")
      data.frame(test = omnibus_labels[[type]], statistic = wt$Fstat,
                 df1 = wt$df_num, df2 = wt$df_denom, p.value = wt$p_val,
                 reference = "F", stringsAsFactors = FALSE)
    } else if (is.finite(df_model)) {
      data.frame(test = omnibus_labels[[type]], statistic = W / k,
                 df1 = k, df2 = df_model,
                 p.value = stats::pf(W / k, k, df_model, lower.tail = FALSE),
                 reference = "F", stringsAsFactors = FALSE)
    } else {
      data.frame(test = omnibus_labels[[type]], statistic = W,
                 df1 = k, df2 = Inf,
                 p.value = stats::pchisq(W, k, lower.tail = FALSE),
                 reference = "chisq", stringsAsFactors = FALSE)
    }
  }))
  rownames(omnibus) <- NULL

  # The dispersion and theta come from the fit the inference came from, as
  # in the four-group model.
  dispersion <- if (dispersion_fixed) {
    sum(stats::residuals(fit, type = "pearson")^2, na.rm = TRUE) / stats::df.residual(fit)
  } else {
    NA_real_
  }

  out <- list(
    model = fit_given,
    coefficients = coefficients$table,
    effects = effects,
    omnibus = omnibus,
    conditions = data.frame(
      condition = lev,
      role = c("control", rep("treatment", k)),
      term = c(NA_character_, terms),
      stringsAsFactors = FALSE
    ),
    weights = weights,
    adjust = adjust,
    dispersion = dispersion,
    theta = if (negbin) c(theta = fit$theta, std.error = fit$SE.theta, alpha = 1 / fit$theta),
    vcov = coefficients$vcov,
    covariate_scaling = .covariate_scaling_element(scaling, fit, vcovM),
    data = df,
    pretest_mean = pretest_mean,
    pretest_center = pretest_center,
    pretest_mean_analyzed = centers[["analyzed"]],
    pretest_mean_observed = centers[["observed"]],
    pretest_n = parts$n,
    pretest_sensitivity = sensitivity,
    robust = robust,
    family = family,
    cluster = cluster,
    n_clusters = n_clusters,
    conf_level = conf_level,
    call = call
  )

  class(out) <- "solomon_ngroup"

  out
}


# Name of a multiplicity adjustment for printing and reporting.
.adjust_label <- function(adjust) {
  switch(
    adjust,
    holm = "Holm's (1979) procedure",
    bonferroni = "the Bonferroni procedure",
    none = "none"
  )
}


# Print a data frame of formatted columns with a header row.
.print_columns <- function(heads, cols, left = 1L) {
  widths <- mapply(function(h, v) max(nchar(h), nchar(v)), heads, cols)
  is_left <- seq_along(heads) <= left
  print_row <- function(values) {
    cells <- ifelse(
      is_left,
      sprintf("%-*s", widths, values),
      sprintf("%*s", widths, values)
    )
    cat(paste(cells, collapse = "  "), "\n", sep = "")
  }
  print_row(heads)
  for (i in seq_along(cols[[1]])) {
    print_row(vapply(cols, `[`, character(1), i))
  }
}


#' @export
print.solomon_ngroup <- function(x, digits = 3, ...) {

  cond <- x$conditions
  k <- sum(cond$role == "treatment")
  treatments <- cond$condition[cond$role == "treatment"]
  control <- cond$condition[cond$role == "control"]

  cat("Solomon GLM (unified model), N-group design\n")
  cat(sprintf(
    "Conditions: %s; control: %s. With and without a pretest: %d groups.\n",
    paste(treatments, collapse = ", "), control, 2L * (k + 1L)
  ))
  f <- tryCatch(stats::formula(x$model), error = function(e) NULL)
  if (!is.null(f)) cat("Formula: ", .formula_line(f), "\n", sep = "")
  cat("Covariance: ", .solomon_vcov_label(x), "\n", sep = "")
  if (!is.null(x$theta)) {
    cat(sprintf("Negative binomial (NB2): theta = %.3g (SE %.3g); alpha = 1/theta = %.3g\n",
                x$theta[["theta"]], x$theta[["std.error"]], x$theta[["alpha"]]))
  }

  level <- if (is.null(x$conf_level)) 0.95 else x$conf_level
  ci_label <- sprintf("%s%% CI", format(100 * level))
  show_df <- any(is.finite(x$effects$df))

  cat("\nOmnibus tests\n")
  om <- x$omnibus
  stat <- ifelse(om$reference == "F",
                 sprintf("F(%s, %s) = %.2f", .df_fmt(om$df1), .df_fmt(om$df2), om$statistic),
                 sprintf("chi2(%s) = %.2f", .df_fmt(om$df1), om$statistic))
  .print_columns(c("Test", "Statistic", "p"), list(om$test, stat, p_fmt(om$p.value)))

  adj <- .adjust_label(x$adjust)
  cat("\nContrasts\n")
  e <- x$effects
  # The treatment contrasts are adjusted across the comparisons, so with one
  # comparison they are not adjusted; the treatments' pretest effects are
  # adjusted across the k >= 2 treatments. Older fits have no pretest rows.
  n_comp <- length(unique(e$comparison[e$contrast %in% .solomon_contrast_order]))
  pretest_rows <- any(e$contrast %in% .solomon_pretest_order)
  show_adj <- x$adjust != "none" && (n_comp > 1L || pretest_rows)
  heads <- c("Comparison", "Contrast", "Est (SE)", if (show_df) "t" else "z")
  cols <- list(e$comparison, e$contrast, estse_str(e$estimate, e$std.error, digits),
               sprintf("%.2f", e$statistic))
  if (show_df) {
    heads <- c(heads, "df")
    cols <- c(cols, list(.df_fmt(e$df)))
  }
  heads <- c(heads, "p", if (show_adj) "p adj.", ci_label)
  cols <- c(
    cols,
    list(p_fmt(e$p.value)),
    if (show_adj) list(p_fmt(e$p.adjusted)),
    list(sprintf("[%.*f, %.*f]", digits, e$conf.low, digits, e$conf.high))
  )
  .print_columns(heads, cols, left = 2L)

  if (!show_adj) {
    cat(if (n_comp == 1L) {
      "\nWith one comparison, the p-values need no adjustment for multiple comparisons.\n"
    } else {
      "\nNo adjustment for multiple comparisons.\n"
    })
  } else {
    cat(if (n_comp == 1L) {
      "\nWith one comparison, the treatment contrasts need no adjustment for multiple comparisons.\n"
    } else {
      sprintf("\np adj.: adjusted by %s within each contrast, across the %d comparisons.\n",
              adj, n_comp)
    })
    if (pretest_rows) {
      cat(sprintf("%s of the treatments' pretest effects: adjusted by %s across the %d treatments.\n",
                  "p adj.", adj, k))
    }
    cat("Confidence intervals are not adjusted.\n")
  }
  if (pretest_rows) {
    cat(if (!is.finite(.null_na(x$pretest_mean))) {
      "Pretest effects: pretested minus unpretested participants in each condition.\n"
    } else if (identical(x$pretest_center, "observed")) {
      paste0(.wrap_lines(sprintf(paste(
        "Pretest effects: pretested minus unpretested participants in each condition,",
        "at the mean of every observed pretest among the pretested participants (%.*f),",
        "with standard errors that include the sampling variance of that mean."
      ), digits, x$pretest_mean)), "\n")
    } else {
      sprintf(paste0(
        "Pretest effects: pretested minus unpretested participants in each condition,\n",
        "at the pretested participants' mean pretest (%.*f), with standard errors\n",
        "that include the sampling variance of that mean.\n"
      ), digits, x$pretest_mean)
    })
    .link_pretest_note(x)
  }
  # The two pretest centers when posttests are missing (issue #135).
  .center_output(x, digits)
  cat(.ngroup_lifecycle_note(), sep = "\n")
  invisible(x)
}

# The lifecycle label of the analysis of designs with several treatments.
# The protocol of its simulation study (issue #45) required the printed
# output to say that the analysis is experimental and to name the scenarios
# that failed the rule for error control.
.ngroup_lifecycle_note <- function() {
  c(
    "",
    "Experimental: in the package's simulation study (issue #45), the omnibus",
    "tests of Condition | pretested and Condition | unpretested rejected in up",
    "to 6.9% of replications at the .05 level with three treatments and 10",
    "participants per group. No other test, and no family of adjusted",
    "comparisons, failed the study's rule. See ?fit_solomon_glm."
  )
}


#' @export
summary.solomon_ngroup <- function(object, ...) {
  object
}
