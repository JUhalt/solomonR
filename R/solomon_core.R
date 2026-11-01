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
#' @section Lifecycle:
#' Each exported function's help page, and the reference index of the
#' package website, shows the function's lifecycle stage, in the stages of
#' the lifecycle package (Henry & Wickham, 2026):
#' - **Stable.** The interface is settled. Any change goes through
#'   deprecation: the old name keeps working, with a warning, through v1.x,
#'   when the arguments that follow a former argument name are also named.
#' - **Experimental.** The function is tested, but its interface or defaults
#'   may change. [fit_solomon_sem()] has no simulation study yet at Solomon
#'   sample sizes. The issue #55 study found no measurement-invariance
#'   criterion that holds its false-rejection rate in Solomon-sized groups,
#'   which affects [fit_solomon_sem_latent()] and [invariance_solomon()].
#'   The analysis of designs with several treatments in [fit_solomon_glm()] is
#'   experimental, because the rule for error control set before its
#'   simulation study (issue #45) was not met for two omnibus tests with
#'   three treatments and 10 participants per group. The pre-specified
#'   rules of their simulation studies were not met for [fit_solomon_mmrm()]
#'   (issue #57) and for [fit_solomon_mi()], [tipping_point_solomon()], and
#'   [plot_tipping_point()] (issue #82). [analysis_plan_solomon()] stays
#'   experimental until researchers have used it (issue #81).
#' - **Deprecated.** A former name that still works, with a warning. Its
#'   help page names the replacement.
#'
#' solomonR 0.8.0 settled the interface (issue #83). The outcome arguments
#' are `y_post` and `y_pre` everywhere, a fitted model is passed as `fit`,
#' and the planning functions share one argument order: `n` (or `power`),
#' `delta`, `sens`, `rho`, `sigma`, `alpha`. Each function that takes data
#' vectors also takes an optional `data` data frame.
#'
#' solomonR 1.0.0 settled the output (issue #110). The fits and analyses
#' listed in [solomon_output] return the Solomon contrasts in a table named
#' `effects` (`results` in [compare_solomon_methods()] and
#' [tipping_point_solomon()], whose rows come from several fits), with one
#' order of columns and one set of contrast labels, beside the confidence
#' level `conf_level`, and [`tidy()`][solomon_output] returns that table.
#' [solomon_output] also names the results that have no such table, and
#' states which parts of a result are stable and which are not.
#'
#' @references
#' Campbell, D. T., & Stanley, J. C. (1966). *Experimental and
#' quasi-experimental designs for research*. Rand McNally. (Original work
#' published 1963)
#'
#' Henry, L., & Wickham, H. (2026). *lifecycle: Manage the life cycle of
#' your package functions* (Version 1.0.5) \[R package\].
#' https://doi.org/10.32614/CRAN.package.lifecycle
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
#' Steyn, R. (2009). Re-designing the Solomon four-group: Can we improve on
#' this exemplary model? *Design Principles and Practices: An International
#' Journal—Annual Review, 3*(1), 383–394. https://doi.org/10.18848/1833-1874/CGP/v03i01/37588
#'
#' van Engelenburg, G. (1999). *Statistical analysis for the Solomon four-group
#' design* (Research Report 99-06). University of Twente. ERIC.
#' https://eric.ed.gov/?id=ED435692
#'
#' Walton Braver, M. C., & Braver, S. L. (1988). Statistical treatment of the
#' Solomon four-group design: A meta-analytic approach. *Psychological Bulletin,
#' 104*(1), 150–154. https://doi.org/10.1037/0033-2909.104.1.150
#'
#' @keywords Solomon four-group ANCOVA permutation maximum-likelihood SEM
#' @name solomonR
NULL

# ---- helpers ----

# The mean pretest among the pretested participants a model uses: the
# complete cases of `df`, which holds the model's variables, including the
# uncentered `pre_obs`. fit_solomon_glm(), its several-treatment path, and
# fit_solomon_mmrm() center the pretest at this mean, as fit_solomon_ml()
# does (van Engelenburg, 1999), so that the pretesting coefficient compares
# pretested and unpretested participants at the pretested participants'
# mean pretest (issue #104). Centering changes neither the fitted values nor
# the treatment contrasts. Without any analyzable pretested participant
# there is nothing to center, and the model's own checks report the problem.
.pretest_center <- function(df) {
  used <- stats::complete.cases(df) & df$pretested %in% 1
  if (!any(used)) return(0)
  mean(df$pre_obs[used])
}

# The participants whose pretests define the other center,
# `pretest_center = "observed"` (issue #135): every pretested participant
# with an observed pretest, whether or not the posttest, the condition, a
# covariate, or the exposure was observed. `df` holds the model's variables,
# with `pre_obs` the pretest (centered or not) in the pretested rows; one
# logical value per row. The participants of .pretest_center() are among
# them, so the two centers are the same number when no pretested
# participant with a pretest is outside the model.
.in_pretest_center <- function(df) {
  df$pretested %in% 1 & !is.na(df$pre_obs)
}

# Both centers, from `df` with the uncentered `pre_obs`: `analyzed`, the
# default, by .pretest_center() itself, and `observed`, the mean over
# .in_pretest_center(), with the same fallback to 0 when there is no such
# participant (the model's own checks then report the problem).
.pretest_centers <- function(df) {
  observed <- .in_pretest_center(df)
  c(analyzed = .pretest_center(df),
    observed = if (any(observed)) mean(df$pre_obs[observed]) else 0)
}

# The participants behind the two centers (issue #135). `x` is the
# uncentered pretest; `in_center` flags the pretested participants with a
# pretest (the "observed" set), `in_model` those in the model, and
# `no_posttest` those without a posttest. Returns the mean pretest of those
# of them outside the model (`mean_outside`; NA when there are none) and
# `n`: how many are in the model (`analyzed`), how many there are in all
# (`observed`), and how many of those outside the model have no posttest
# (`no_posttest`; the rest lack the condition, a covariate, or the
# exposure).
.pretest_center_parts <- function(x, in_center, in_model, no_posttest) {
  inside <- in_center & in_model
  outside <- in_center & !in_model
  list(
    mean_outside = if (any(outside)) mean(x[outside]) else NA_real_,
    n = c(analyzed = sum(inside), observed = sum(in_center),
          no_posttest = sum(outside & no_posttest))
  )
}

# The other center of a fit (issue #135): "observed" for the default, and
# "analyzed" for pretest_center = "observed".
.other_center <- function(center) {
  if (identical(center, "observed")) "analyzed" else "observed"
}

# The attrition of each Solomon group and, on the identity link, the
# offsets at which the default center is right (issue #135), for
# `pretest_sensitivity`. `condition` names each participant's condition (NA
# when it was not recorded), `levels` the conditions in order (the control
# first), `pretested` and `y` the indicator and the posttest, and `in_model`
# the rows the model used. `slope` is the pretest slope and `centers` the
# two centers.
#
# - `attrition`: for each group, the number of participants, the number
#   without a posttest, and their share.
# - `delta`: for each unpretested group, its number of participants, its
#   share outside the model (m0), and the difference between its dropouts'
#   and its completers' mean posttest at which the pretest effect at the
#   "observed" center equals the default's, b (c_B - c_A) / m0 (National
#   Research Council, 2010, pp. 88-89; Mavridis et al., 2015, sec. 2): NA
#   where m0 is zero. The delta parameterization is on the scale of the
#   posttest, so `delta` is NULL off the identity link.
.pretest_attrition <- function(condition, levels, pretested, y, in_model, slope, centers,
                               identity) {
  group_condition <- rep(levels, 2L)
  group_pretested <- rep(c(1L, 0L), each = length(levels))
  lost <- is.na(y)
  counts <- vapply(seq_along(group_condition), function(i) {
    member <- pretested %in% group_pretested[i] & condition %in% group_condition[i]
    c(sum(member), sum(member & lost))
  }, numeric(2))
  n <- as.integer(counts[1, ])
  missing <- as.integer(counts[2, ])
  attrition <- list2DF(list(
    group = paste0(ifelse(group_pretested == 1L, "Pretested", "Unpretested"), ", ",
                   group_condition),
    condition = group_condition, pretested = group_pretested,
    n = n, posttest_missing = missing, share = ifelse(n > 0L, missing / n, NA_real_)
  ))
  delta <- NULL
  if (identity) {
    outside <- vapply(levels, function(level) {
      member <- pretested %in% 0L & condition %in% level
      c(sum(member), sum(member & !in_model))
    }, numeric(2))
    n0 <- as.integer(outside[1, ])
    m0 <- ifelse(n0 > 0L, unname(outside[2, ]) / n0, NA_real_)
    gap <- slope * (centers[["observed"]] - centers[["analyzed"]])
    delta <- list2DF(list(
      condition = levels, n = n0, m0 = m0,
      delta = ifelse(is.finite(m0) & m0 > 0, gap / m0, NA_real_)
    ))
  }
  list(attrition = attrition, delta = delta)
}

# `pretest_sensitivity` of a fit with a pretest covariate (issue #135): the
# pretest effects at both centers and their difference, with what the
# printed block and check_solomon_missing() report. `rows` holds the pretest
# rows of the effects table at the center used (`center`) and at the other
# center (`other`; NULL when the two centers are the same number, in which
# case the rows are copied). `keys` names the columns that index the rows
# before `contrast` (`comparison` for several treatments, `occasion` for
# repeated posttests). `parts` comes from .pretest_center_parts(), `centers`
# from .pretest_centers(), and `attrition` from .pretest_attrition();
# `slope` is the pretest slope (one per occasion for repeated posttests,
# named by occasion), `identity` whether the link is the identity, and
# `note` the reason the other center's standard errors are missing, if they
# are.
.pretest_sensitivity <- function(used, other, center, keys, parts, centers, slope, identity,
                                 attrition, note = NULL) {
  differ <- parts$n[["observed"]] > parts$n[["analyzed"]]
  if (!differ || is.null(other)) other <- used
  analyzed <- if (center == "analyzed") used else other
  observed <- if (center == "analyzed") other else used
  columns <- names(used)
  both <- lapply(stats::setNames(columns, columns), function(column) {
    c(analyzed[[column]], observed[[column]])
  })
  both$center <- rep(c("analyzed", "observed"), each = nrow(used))
  effects <- .effects_table(list2DF(both), keys = c(keys, "center"))
  difference <- list2DF(c(
    lapply(stats::setNames(c(keys, "contrast"), c(keys, "contrast")), function(column) {
      analyzed[[column]]
    }),
    list(analyzed = analyzed$estimate, observed = observed$estimate,
         difference = analyzed$estimate - observed$estimate)
  ))
  m1 <- if (parts$n[["observed"]] > 0L) {
    (parts$n[["observed"]] - parts$n[["analyzed"]]) / parts$n[["observed"]]
  } else {
    NA_real_
  }
  list(
    effects = effects,
    difference = difference,
    slope = slope,
    m1 = m1,
    pretest_means = c(analyzed = centers[["analyzed"]], not_analyzed = parts$mean_outside),
    delta = attrition$delta,
    predicted = if (identity && length(slope) == 1L) {
      slope * (parts$mean_outside - centers[["analyzed"]])
    } else {
      NA_real_
    },
    attrition = attrition$attrition,
    differ = differ,
    note = note
  )
}

# Each analyzed participant's influence on the coefficients of a GLM fit,
# (X'WX)^{-1} x_i w_i e_i, with w_i the working weight and e_i the working
# residual, so that the sum of the cross-products is the HC0 covariance; and
# the hat values h_i, by which HC3 divides each influence by 1 - h_i
# (MacKinnon & White, 1985). One row per row the model used.
#
# Both come from the QR decomposition W^(1/2) X = QR that the fit computed:
# the influence is R^{-1} q_i w_i^(1/2) e_i and the hat value is ||q_i||^2,
# with q_i the ith row of Q. Forming X'WX instead squares the condition
# number of the design, so that a covariate far from zero relative to its
# spread made it numerically singular (issue #137). The decomposition holds
# the rows whose working weight is positive: a row of weight zero, which
# glm.fit() leaves out, has no influence and a hat value of zero. Should the
# fit hold no decomposition of exactly those rows, the weighted design is
# decomposed again, with the fit's tolerance.
#
# With `map`, a matrix whose columns are the coefficients, the influence is
# that on the combinations map %*% b, one column for each row of `map`:
# q_i' R^{-T} map' w_i^(1/2) e_i, with R^{-T} map' from backsolve(). That
# keeps a combination accurate whose terms are large and cancel, such as the
# intercept of covariates as given computed from their standardized
# coefficients (.given_vcov()), where the influences on the terms, combined
# afterwards, would carry their rounding error into it.
#
# The columns of the decomposition are the coefficients q$pivot, taken by
# position: two coefficients can share a name (.scaled_columns()).
.glm_influence <- function(model, map = NULL) {
  cf <- stats::coef(model)
  w <- model$weights
  e <- model$residuals
  good <- w > 0
  q <- model$qr
  if (is.null(q) || nrow(q$qr) != sum(good)) {
    X <- stats::model.matrix(model)[good, , drop = FALSE]
    q <- qr(sqrt(w[good]) * X, tol = if (is.null(q$tol)) 1e-11 else q$tol)
  }
  rank <- seq_len(q$rank)
  Q <- qr.Q(q)[, rank, drop = FALSE]
  hat <- numeric(length(w))
  hat[good] <- rowSums(Q^2)
  pivot <- q$pivot[rank]
  if (!is.null(map)) {
    B <- backsolve(qr.R(q)[rank, rank, drop = FALSE], t(map[, pivot, drop = FALSE]),
                   transpose = TRUE)
    influence <- matrix(0, length(w), nrow(map), dimnames = list(NULL, rownames(map)))
    influence[good, ] <- (Q * (sqrt(w[good]) * e[good])) %*% B
    return(list(influence = influence, hat = hat))
  }
  influence <- matrix(0, length(w), length(cf), dimnames = list(NULL, names(cf)))
  influence[good, pivot] <- t(backsolve(qr.R(q)[rank, rank, drop = FALSE],
                                        t(Q * (sqrt(w[good]) * e[good]))))
  list(influence = influence[, !is.na(cf), drop = FALSE], hat = hat)
}

# The variance that estimated sample means add to a contrast evaluated at
# them (issue #104). A pretest effect compares the pretested participants at
# their mean pretest with the unpretested participants, whose expected
# pretest that mean estimates. The mean is itself an estimate: treating it
# as a fixed number leaves out b^2 Var(mean), with b the pretest slope.
#
# Here a contrast depends, to first order, on the sum over the samples in
# `set` of each sample's mean of `a`. `a` and `set` run over participants;
# `set` is NA for those in no sample. Stacking the estimating equations of
# the means with the model's (Stefanski & Boos, 2002) gives the variance of
# that sum, `v`, and its covariance with the coefficients, `cross` (a vector
# over the coefficients), in the covariance type of the fit:
# - "none" (model-based): s^2 / n in each sample, and no covariance with the
#   coefficients, which is zero when the model is correctly specified;
# - "HC3": the jackknife form, sum (a_i - mean)^2 / (n - 1)^2 in each
#   sample, and the covariance from the HC3-scaled influences (which equal
#   the jackknife changes in the coefficients of a linear model);
# - "CR2": the CR2 covariance of the sample means (Bell & McCaffrey, 2002)
#   with Satterthwaite degrees of freedom (Pustejovsky & Tipton, 2018), and
#   the covariance from cluster sums with the small-sample factor G / (G - 1).
# `df` is the degrees of freedom of `v`: n - 1 for one sample, combined by
# the Welch-Satterthwaite formula over several, or Satterthwaite's for CR2.
#
# A sample may hold participants who are not in the model: the mean pretest
# of `pretest_center = "observed"` is taken over every pretested participant
# with a pretest, with or without a posttest (issue #135). `used` then flags,
# among the participants `a` and `set` run over, those the model used, in
# the model's row order; by default all of them, which is the computation
# for the default center. Each sum runs over the participants who
# contribute to it:
# - `v` over every member of a sample, in the model or not, with n their
#   number;
# - `cross` over the members who are in the model, because a participant
#   outside the model has no influence on the coefficients (and leaving one
#   out changes the mean and not the coefficients). With CR2, a cluster's
#   sum of the deviations takes every member of the sample in the cluster,
#   and G counts the clusters with a participant in the model or in a
#   sample.
# `cluster` runs over the same participants as `a`.
.center_terms <- function(model, robust, a, set, cluster = NULL, used = NULL) {
  if (is.null(used)) used <- rep(TRUE, length(a))
  in_set <- !is.na(set)
  set <- factor(set)
  n_s <- as.numeric(table(set)[as.character(set)])
  dev <- numeric(length(a))
  dev[in_set] <- a[in_set] - stats::ave(a[in_set], set[in_set])
  n_coef <- sum(!is.na(stats::coef(model)))
  cross <- stats::setNames(numeric(n_coef), names(stats::coef(model))[!is.na(stats::coef(model))])
  # Means of constants, such as standardized risks without covariates, add
  # nothing.
  if (!any(in_set) || max(abs(dev)) <= 1e-12 * max(1, abs(a[in_set]))) {
    return(list(v = 0, cross = cross, df = Inf))
  }

  # Welch-Satterthwaite degrees of freedom over the samples.
  sample_df <- function(v_by_set) {
    n <- as.numeric(table(set))
    keep <- v_by_set > 0
    if (!any(keep)) return(Inf)
    sum(v_by_set[keep])^2 / sum(v_by_set[keep]^2 / (n[keep] - 1))
  }

  if (robust == "CR2") {
    # A sample whose values are all equal adds no variance, and would
    # distort the Satterthwaite degrees of freedom, which depend on the
    # design and not on the residuals.
    varies <- tapply(abs(dev[in_set]), set[in_set], max) > 1e-12 * max(1, abs(a[in_set]))
    in_set <- in_set & set %in% names(varies)[varies]
    rows <- data.frame(a = a[in_set], set = droplevels(set[in_set]))
    one <- nlevels(rows$set) == 1L
    m <- if (one) stats::lm(a ~ 1, data = rows) else stats::lm(a ~ 0 + set, data = rows)
    Vm <- clubSandwich::vcovCR(m, cluster = cluster[in_set], type = "CR2")
    lc <- as.data.frame(clubSandwich::linear_contrast(
      m, vcov = Vm, contrasts = matrix(1, nrow = 1, ncol = length(stats::coef(m))),
      test = "Satterthwaite"
    ))
    # Each participant's influence on the coefficients: zero outside the
    # model.
    infl <- matrix(0, length(a), n_coef)
    infl[used, ] <- .glm_influence(model)$influence
    u <- numeric(length(a))
    u[in_set] <- dev[in_set] / n_s[in_set]
    # The participants of the stacked equations: in the model or in a sample.
    stacked <- used | in_set
    by_cluster <- cluster[stacked]
    G <- length(unique(by_cluster))
    cross[] <- G / (G - 1) * colSums(
      rowsum(infl[stacked, , drop = FALSE], by_cluster) * rowsum(u[stacked], by_cluster)[, 1]
    )
    return(list(v = lc$SE^2, cross = cross, df = lc$df))
  }

  if (robust == "HC3") {
    u <- numeric(length(a))
    u[in_set] <- dev[in_set] / pmax(n_s[in_set] - 1, 1)
    infl <- .glm_influence(model)
    cross[] <- colSums(infl$influence / (1 - infl$hat) * u[used])
  } else {
    u <- numeric(length(a))
    u[in_set] <- dev[in_set] / sqrt(n_s[in_set] * pmax(n_s[in_set] - 1, 1))
  }
  v_by_set <- as.numeric(tapply(u[in_set]^2, set[in_set], sum))
  list(v = sum(v_by_set), cross = cross, df = sample_df(v_by_set))
}

# The clusters a clustered (CR2) fit needs for the "observed" center (issue
# #135): those of every pretested participant whose pretest enters it,
# including those who are not in the model (`outside` flags them), because
# the cluster-robust variance of the mean pretest, and its covariance with
# the coefficients, are built from cluster sums. A participant is outside
# the model for want of a posttest, a recorded condition, a covariate, or
# the exposure. With `pretest_center = "observed"` a missing cluster stops
# the fit; with the default it only leaves the other center's standard
# errors missing (`stop = FALSE` returns the message instead).
.check_center_clusters <- function(cluster, outside, stop = TRUE) {
  missing <- sum(is.na(cluster[outside]))
  if (!missing) return(NULL)
  message <- paste0(
    "`cluster` is missing for ", missing, " pretested participant(s) who have a ",
    "pretest and are not in the model (no posttest, condition, covariate, or exposure). ",
    "With pretest_center = \"observed\" their pretests enter the analysis, whose ",
    "cluster-robust variance needs their clusters."
  )
  if (stop) stop(message, call. = FALSE)
  message
}

# The terms of the center used and of the other center (issue #135), for
# the four-group and the several-treatment paths of fit_solomon_glm(). `df`
# holds the model's variables with the centered `pre_obs`, `x` is the
# uncentered pretest, `in_model` flags the rows the model used, `model` is
# the fit the covariance comes from (the complete-case refit for CR2),
# `cluster` the cluster of every row, and `cluster_fit` those of the rows
# in the model (CR2 only). The analyzed center's terms are computed as
# before #135: over the model's rows, with their clusters. The other
# center's are computed only when the two centers differ; with the default
# and CR2, a missing cluster of a participant outside the model leaves them
# out, with `note` saying so, while with pretest_center = "observed" it
# stops the fit. Returns `center`, `other`, `note`, `parts`
# (.pretest_center_parts()), and `differ`.
.center_terms_both <- function(df, x, in_model, model, robust, cluster, cluster_fit,
                               pretest_center) {
  in_observed <- .in_pretest_center(df)
  parts <- .pretest_center_parts(x, in_observed, in_model, is.na(df$y))
  differ <- parts$n[["observed"]] > parts$n[["analyzed"]]
  terms_analyzed <- function() {
    .center_terms(
      model, robust,
      a = df$pre_obs[in_model],
      set = ifelse(df$pretested[in_model] == 1L, "pretested", NA_character_),
      cluster = if (robust == "CR2") cluster_fit
    )
  }
  terms_observed <- function() {
    .center_terms(
      model, robust,
      a = df$pre_obs,
      set = ifelse(in_observed, "pretested", NA_character_),
      cluster = if (robust == "CR2") cluster,
      used = in_model
    )
  }
  center <- NULL
  other <- NULL
  note <- NULL
  if (pretest_center == "analyzed") {
    center <- terms_analyzed()
    if (differ) {
      if (robust == "CR2") {
        note <- .check_center_clusters(cluster, in_observed & !in_model, stop = FALSE)
        if (!is.null(note)) {
          note <- paste(note, "The standard errors at that center are not given.")
        }
      }
      if (is.null(note)) other <- terms_observed()
    }
  } else {
    if (robust == "CR2") .check_center_clusters(cluster, in_observed & !in_model)
    center <- terms_observed()
    if (differ) other <- terms_analyzed()
  }
  list(center = center, other = other, note = note, parts = parts, differ = differ)
}

# The pretest effects at the other center (issue #135), without a refit:
# moving the center by `shift` moves a pretest effect by its weight on the
# pretesting coefficient times the pretest slope times `shift`, so in the
# fitted model it is the contrast L + shift L["pretested"] e_pre_obs.
# `Ls` holds the contrasts of the pretest effects, `contrast` their labels,
# `comparison` (for several treatments) their comparisons, and `estimate` a
# function that estimates a contrast with the other center's terms, as in
# the effects table. With `no_se`, the other center has no standard errors
# (a CR2 fit without a cluster; see .center_terms_both()).
.other_center_rows <- function(Ls, contrast, shift, estimate, comparison = NULL, no_se = FALSE) {
  rows <- do.call(rbind, Map(function(L, label) {
    L[["pre_obs"]] <- L[["pre_obs"]] + shift * L[["pretested"]]
    est <- estimate(L)
    data.frame(contrast = label, estimate = est[["estimate"]],
               std.error = est[["std.error"]], statistic = est[["statistic"]],
               df = est[["df"]], p.value = est[["p.value"]],
               conf.low = est[["conf.low"]], conf.high = est[["conf.high"]],
               stringsAsFactors = FALSE)
  }, unname(Ls), contrast))
  if (!is.null(comparison)) rows <- cbind(comparison = comparison, rows, stringsAsFactors = FALSE)
  if (no_se) rows[setdiff(.solomon_effect_columns, c("contrast", "estimate"))] <- NA_real_
  rownames(rows) <- NULL
  rows
}

# Welch-Satterthwaite degrees of freedom of a sum of independent variance
# components `v` with degrees of freedom `df` (Satterthwaite, 1946; Welch,
# 1947).
.welch_df <- function(v, df) {
  v <- pmax(v, 0)
  keep <- v > 0
  if (!any(keep)) return(min(df))
  sum(v[keep])^2 / sum(v[keep]^2 / df[keep])
}

# The variance a GLM fit's contrast `L` (a vector over the coefficients)
# adds when it is evaluated at the estimated mean pretest: the contrast
# moves by L["pretested"] * b_x for each unit the mean moves, so it gains
# (L["pretested"] b_x)^2 v + 2 L["pretested"] b_x L'cross, from
# .center_terms() (`center`). Returns the variance of L'b, `variance`, the
# added part, `added`, and, with `df_contrast`, the degrees of freedom: those
# of the contrast for HC3 and model-based fits, whose reference
# distribution has the residual degrees of freedom of the model, and for CR2
# their Welch-Satterthwaite combination with those of the mean.
.center_adjusted <- function(L, cf, V, center, robust, df_contrast) {
  variance <- as.numeric(t(L) %*% V %*% L)
  k <- if (is.null(center) || !"pretested" %in% names(L)) 0 else L[["pretested"]] * cf[["pre_obs"]]
  if (!is.finite(k) || k == 0) return(list(variance = variance, added = 0, df = df_contrast))
  # `cross` runs over the coefficients in their order, as `L` does; by
  # position, since two coefficients can share a name (.scaled_columns()).
  weights <- if (length(center$cross) == length(L)) L else L[names(center$cross)]
  v_contrast <- variance + 2 * k * sum(weights * center$cross)
  v_center <- k^2 * center$v
  df <- if (robust == "CR2") {
    .welch_df(c(v_contrast, v_center), c(df_contrast, center$df))
  } else {
    df_contrast
  }
  list(variance = v_contrast + v_center, added = v_contrast + v_center - variance, df = df)
}

# NA for an element that objects made by earlier versions lack.
.null_na <- function(x) if (is.null(x)) NA_real_ else x

#' Convert p-value to Z (one-tailed) for Stouffer's method
#'
#' `r lifecycle::badge("stable")`
#'
#' @param p numeric vector of p-values assumed one-tailed and aligned in the same direction
#' @return A numeric vector of z scores, one for each p-value: the standard
#'   normal quantile of `1 - p`.
#' @references
#' Stouffer, S. A., Suchman, E. A., DeVinney, L. C., Star, S. A., & Williams, R.
#' M., Jr. (1949). *The American soldier: Adjustment during army life* (Vol. 1).
#' Princeton University Press.
#' @examples
#' p_to_z(c(0.05, 0.025, 0.5))
#' @export
p_to_z <- function(p) {
  stats::qnorm(1 - p)
}

#' Stouffer's Z combiner (Walton Braver & Braver, 1988, Test I)
#'
#' `r lifecycle::badge("stable")`
#' Combine one-tailed p-values that test the *same directional* hypothesis into
#' a single Z. This is provided to reproduce the Walton Braver & Braver (1988)
#' meta-analytic option (Test I) for the Solomon four-group design. Use cautiously and document assumptions about homogeneity; see
#' the 1988–1990 exchanges for caveats.
#'
#' Walton Braver and Braver (1988, p. 152) convert the one-tailed p-value of
#' each test to z and refer the combined z to a normal table. Their worked
#' example reports z = 2.05 with p = .040 (p. 153), the two-tailed value, so
#' [fit_solomon_classic()] judges Test I by `p_meta_two_tailed`. The one-tailed
#' value is also returned.
#' @param p numeric vector of one-tailed p-values (same direction)
#' @return A list with:
#'   - `z_meta`: the combined z, the sum of the z scores divided by the
#'     square root of their number.
#'   - `p_meta_two_tailed`: its two-tailed p-value.
#'   - `p_meta_one_tailed`: its one-tailed p-value, in the direction of the
#'     p-values supplied.
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
#'
#' # Walton Braver and Braver's (1988, p. 153) example: the ANCOVA (p = .0993)
#' # and the t test (p = .2127), halved to one-tailed p-values in the
#' # direction of the effect, give z = 2.05 and p = .040.
#' stouffer_solomon(c(0.0993, 0.2127) / 2)
#' @export
stouffer_solomon <- function(p) {
  stopifnot(all(is.finite(p)), all(p > 0 & p < 1))
  z <- p_to_z(p)
  z_meta <- sum(z) / sqrt(length(z))
  list(z_meta = z_meta, p_meta_two_tailed = 2 * stats::pnorm(-abs(z_meta)),
       p_meta_one_tailed = 1 - stats::pnorm(z_meta))
}

#' Fit the unified GLM for a Solomon Four-Group design
#'
#' `r lifecycle::badge("stable")`
#' Fits one generalized linear model to all four Solomon groups and reports
#' four Solomon contrasts: the equal-weighted average treatment effect, the
#' Pretest x Treatment (sensitization) contrast, and the treatment effect
#' within each pretesting condition. It also reports the pretest (testing)
#' effect among controls, among treated participants, and averaged over the
#' two (see "The pretest effect").
#'
#' When `y_pre` is supplied, it enters the model as `pre_obs`: in the
#' pretested groups, the pretest score minus the mean pretest of the
#' pretested participants in the model (returned as `pretest_mean`), and in
#' the unpretested groups 0, so the structurally absent pretests do not
#' remove Groups 3 and 4.
#' This is the model Newman et al. (1990) proposed for the design, which they
#' call a pseudo-analysis of covariance (pp. 94, 98), written with treatment
#' and pretest indicators in place of their group indicators (equation 7,
#' p. 98); the two forms give the same fitted values. Each group keeps its
#' own mean, so the pretest slope is estimated within the pretested groups,
#' and the unpretested groups' means are not adjusted (p. 98). They tested
#' the interaction as the restriction b1 = b3 + b2 - b4 on the group
#' coefficients (p. 98), the difference of differences that the Pretest x
#' Treatment contrast estimates. With `robust = "none"`, the fit reproduces
#' their pretest slope (.55264, p. 98), their within-groups sum of squares,
#' and their tests of the treatment (F = 21.01) and the interaction
#' (F = 0.21, printed .22) (Table 3, p. 100). Centering changes neither the
#' fitted values nor the treatment contrasts; it makes the pretesting
#' coefficient the pretest effect among controls at the pretested
#' participants' mean pretest, as in [fit_solomon_ml()], the mean at which
#' Newman et al. took their adjusted means (p. 101). Without centering, the
#' coefficient would compare the two control groups at a pretest score of 0,
#' which treats the unpretested participants as if they had scored 0, the
#' adjustment Newman et al. cautioned against (p. 98).
#' solomonR adds the named contrasts with confidence intervals, the pretest
#' effects, robust and cluster-robust inference, other families, and designs
#' with several treatments. Regression adjustment for baseline covariates in
#' randomized experiments, and the case for pairing it with
#' heteroskedasticity-robust standard errors, is discussed by Lin (2013).
#' Pretested participants with a missing pretest score are excluded with a
#' warning; incidental missingness is never imputed. Their exclusion biases
#' the pretest effects when it relates to baseline (issue #142).
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
#'   participants are nested in clusters such as classrooms or sites. In the
#'   package's simulation of cluster-randomized Solomon designs (issue #19),
#'   CR2 tests exceeded the nominal level with four clusters per arm (Type I
#'   error up to 0.068 with no effect in any cluster); with few clusters, the
#'   cluster-level randomization test of [perm_solomon()] is preferred.
#'
#' Tests and confidence intervals use the same reference distribution. For
#' `"HC3"` and `"none"` it is the t distribution with residual degrees of
#' freedom when the dispersion is estimated, as in Gaussian models, which is
#' the conventional choice when t approximations are used with robust
#' standard errors (Imbens & Kolesár, 2016; Rajh-Weber et al., 2025). Families
#' with a fixed dispersion (binomial, Poisson) use the normal distribution.
#' With a noncollapsible link (such as the logit) and `y_pre`, the
#' link-scale contrasts compare a treatment effect conditional on the pretest
#' among pretested participants with a marginal effect among unpretested
#' participants, who have no pretest. When the pretest predicts the outcome,
#' the Pretest x Treatment contrast is then nonzero even without
#' sensitization (Daniel et al., 2021): in the package's simulation study
#' (issue #43) it averaged 0.08 to 0.21 on the log-odds scale with no
#' sensitization present. Such fits give the classed warning
#' `solomonR_noncollapsible_warning`; for binary outcomes, estimate the
#' Solomon contrasts on a common scale with [marginal_solomon()]. It takes
#' the fit of a four-group design: for a design with several treatments,
#' subset the data to one treatment and the control and fit the subset.
#' Identity and log links are collapsible and are not affected.
#'
#' Count outcomes: with `family = poisson()`, the default HC3 covariance gives
#' the robust (quasi-likelihood) inference that Cameron and Trivedi (2013)
#' describe. The Poisson estimator stays consistent when the counts are not
#' Poisson, provided the mean is correctly specified (p. 72), whereas
#' model-based standard errors should not be used under overdispersion. The
#' Pearson dispersion statistic is returned as `dispersion` for description;
#' values well above 1 indicate overdispersion. Log-link rate ratios are
#' collapsible, so the noncollapsibility caution above does not apply.
#' Use `exposure` for counts observed over different times or exposures, and
#' [marginal_solomon()] for rate differences. No published Solomon study with
#' a count outcome has been identified, so applying these count-data methods
#' to the Solomon contrasts is a solomonR extension; its simulation validation
#' is on issue #44.
#'
#' `family = "negative_binomial"` fits the NB2 model, with variance
#' \eqn{\mu + \alpha\mu^2}, by maximum likelihood using `MASS::glm.nb()`
#' (Venables & Ripley, 2002). Its coefficient estimates stay consistent when the
#' counts are not negative binomial, provided the mean is correctly specified,
#' but its model-based standard errors do not, so robust standard errors are
#' advised (Cameron & Trivedi, 2013, pp. 84--85); the default HC3 covariance
#' provides them. The estimated \eqn{\theta = 1/\alpha} is returned as
#' `theta`. When the counts show little overdispersion, \eqn{\theta} does not
#' converge and a classed warning (`solomonR_theta_boundary_warning`) is
#' given. Clustered (CR2) negative-binomial fits are not supported.
#'
#' In the package's pre-registered simulation (issue #62), which reused the
#' datasets of issue #44, the NB2 fit with HC3 met the coverage and Type I
#' tolerances in 85% of the overdispersed contrasts with 50 or 100
#' participants per cell, against 82% for robust Poisson. That fell short of
#' the 90% set in advance for recommending it, although its estimates were
#' on average 3.5% more precise with strong overdispersion. Robust Poisson
#' therefore remains the recommendation, and Poisson is preferred when there
#' is no evidence of overdispersion. Model-based NB2 standard errors fell
#' outside the coverage tolerance in 162 of 384 contrasts and should not be
#' used.
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
#' @section The pretest effect:
#' Solomon (1949) added the unpretested groups to separate the effect of
#' taking the pretest from the effect of the treatment, and Campbell and
#' Stanley (1963/1966, p. 25) list the main effect of testing among the
#' quantities the design estimates. The effects table reports it after the
#' four treatment contrasts, as pretested minus unpretested participants:
#' - `Pretest effect | control`: among control participants;
#' - `Pretest effect | treated`: among treated participants;
#' - `Pretest main effect`: the equal-weighted average of the two.
#'
#' The two pretest effects differ by the Pretest x Treatment contrast: an
#' interaction can be read as a treatment effect that depends on pretesting
#' or as a pretest effect that depends on treatment.
#'
#' With `y_pre`, the pretested groups are compared with the unpretested
#' groups at the pretested participants' mean pretest, where the centered
#' pretest is zero. Unpretested participants were never measured, but with
#' random assignment their expected pretest equals that of the pretested
#' participants (Solomon, 1949, p. 141), whose combined mean is its best
#' estimate (Solomon & Lessac, 1968, p. 147). Each pretest effect is therefore
#' the pretest-adjusted mean of a pretested group minus the mean of the
#' unpretested group in the same treatment condition: the differences between
#' the adjusted means that [plot_sensitization()] draws. Without `y_pre`, the
#' pretest effects are differences between the fitted cell means. When
#' posttests are missing, the mean pretest can be taken over two sets of
#' pretested participants; see "The pretest effect when posttests are
#' missing".
#'
#' That mean pretest is an estimate of the expected pretest, not a fixed
#' value, and the pretest effects shift by b for each point it shifts, where
#' b is the pretest slope. Their standard errors therefore include its sampling
#' variance, about b^2 s^2 / n for n pretested participants whose pretests
#' have variance s^2. The estimating equation of the mean is stacked with
#' those of the model (Stefanski & Boos, 2002): each pretest effect gains b^2
#' times the variance of the mean and 2b times the covariance of the mean
#' with the contrast, estimated in the fit's covariance type:
#' - with model-based covariance, the variance of the mean is s^2 / n, and
#'   the covariance is zero, as it is when the model is correctly specified;
#' - with HC3, the variance is in its jackknife form, and the covariance is
#'   estimated from the HC3-scaled influence of each participant on the
#'   coefficients;
#' - with CR2, the variance is the CR2 variance of the mean, and the
#'   covariance is estimated from cluster sums with the factor G / (G - 1)
#'   for G clusters. The degrees of freedom combine the Satterthwaite degrees
#'   of freedom of the contrast and of the mean by the Welch-Satterthwaite
#'   formula (Satterthwaite, 1946; Welch, 1947).
#'
#' With HC3 and model-based covariance, the reference distribution keeps the
#' residual degrees of freedom of the model. The treatment contrasts do not
#' depend on the mean, so their standard errors are unchanged. The coefficient
#' table reports the pretesting coefficient with the model's standard error,
#' which treats the mean as fixed, and the Wald R-squared of a pretest effect
#' uses its full standard error and has no interval.
#'
#' Treating the mean as fixed leaves its variance out at every sample size.
#' In a simulation check of this correction (not a pre-registered study;
#' 1,000 to 4,000 replications in each of 10 scenarios), with a
#' pretest-posttest correlation of .8 the 95% intervals of the pretest main
#' effect that treat the mean as fixed covered 0.90 of the time, with 30 and
#' with 100 participants per cell. With its variance included, the intervals
#' of the pretest effects covered 0.937 to 0.958 of the time, within the
#' Monte Carlo tolerance in all 82 cells: HC3, model-based, and CR2 fits of
#' four-group designs, including one whose pretest slope differed between
#' the treatment conditions, a design with two treatments, [fit_solomon_ml()]
#' with either inference, [fit_solomon_mmrm()], and the delta method of
#' [marginal_solomon()].
#'
#' On a link other than the identity, with `y_pre`, the pretest effects
#' compare the pretested participants' fitted mean at the mean pretest with
#' the unpretested participants' mean over their unmeasured pretests. With a
#' nonlinear link, a mean at the average pretest is not the average of the
#' means over the pretests, so these contrasts differ from the marginal pretest
#' effect even when the pretest has no effect. This holds for the log link
#' too, although its treatment rate ratios are collapsible.
#' [marginal_solomon()] estimates the pretest effects on a common scale, from
#' standardized risks or rates (Daniel et al., 2021).
#'
#' @section The pretest effect when posttests are missing:
#' Each pretest effect compares a pretested group, at a mean pretest (the
#' center), with the unpretested group in the same treatment condition. The
#' center stands in for the unpretested groups' expected pretest, which
#' random assignment makes equal to the pretested groups' (Solomon, 1949, p.
#' 141). With every posttest observed it is the mean pretest of both
#' pretested groups (Newman et al., 1990, p. 101). When posttests are
#' missing there are two candidates, and `pretest_center` chooses between
#' them:
#' - `"analyzed"` (the default): the mean pretest of the pretested
#'   participants in the model. It assumes that, within each treatment
#'   condition, the pretested and the unpretested groups lost the same kinds
#'   of participants. Random assignment makes that the natural case, unless
#'   taking the pretest itself changed who dropped out.
#' - `"observed"`: the mean of every observed pretest among the pretested
#'   participants, with or without a posttest. It assumes that the pretested
#'   groups' losses depend at most on the observed pretest, and that the
#'   unpretested groups' losses are unrelated to baseline as well as to the
#'   posttest. In the pretested groups this is missingness at random given
#'   the pretest (Seaman et al., 2013, Definition 3; the mapping is the
#'   package's inference), the assumption under which [fit_solomon_mi()]
#'   imputes; in one check on issue #135, normal-theory full-information
#'   maximum likelihood with the pretest gave the same pretest effects to
#'   three decimals. It needs the rows of the pretested participants without
#'   a posttest, with `y_post` missing; without them it gives the default.
#'
#' The two centers differ by m1 (x_in - x_out), and the pretest effects at
#' the two centers by b m1 (x_in - x_out), where b is the pretest slope, m1
#' the share of the pretested participants with a pretest who are not in
#' the model, and x_in and x_out the mean pretests of those in and out of
#' the model. The data cannot tell which center is right, because the
#' unpretested groups have no pretest (Campbell, 1957, p. 305), and no
#' published source decides between them. The default is a judgment, set out
#' with its reasons in the decision record posted on issue #135:
#' - Under the design's null hypothesis, that taking the pretest changes
#'   neither the posttest nor who stays, random assignment makes all groups
#'   lose the same kinds of participants, whatever drives the loss, and the
#'   default's test of a zero pretest effect keeps its level. In the record's
#'   check (500 data sets of 200 participants per group, about a third of
#'   the posttests lost by baseline in all four groups), the test of the
#'   pretest main effect at the .05 level rejected in 5% of the data sets at
#'   the default and in 54% at the observed center.
#' - The default fails only when taking the pretest changed who stayed,
#'   which is itself an effect of pretesting; it then reports that effect as
#'   an effect on scores, and can be wrong by up to the whole gap, as when
#'   taking the pretest itself caused the selective loss. When only the
#'   pretested groups lose posttests, by their pretest, the observed center
#'   is unbiased and the default is not. When pretesting only adds losses at
#'   random, the default is exact and the observed center is off by the
#'   whole gap. When the pretested groups lose by their pretest and the
#'   unpretested groups lose at random at the same rate, the loss rates do
#'   not differ and no diagnostic warns: the printed observed center is
#'   then the only check.
#' - The default keeps the design's own pretest effect. For the main effect
#'   it is, when every pretested participant has a pretest (see issue #142)
#'   and up to a term that vanishes when the two pretested groups keep equal
#'   numbers, the difference between the posttest means of the participants
#'   analyzed, as the design defines the effect (Solomon & Lessac, 1968, p.
#'   147). It does not depend on whether the data keep rows for participants
#'   without a posttest.
#'
#' Whenever the two centers differ, the printed fit gives the pretest
#' effects at both, their difference, and the shares of posttests missing,
#' and `pretest_sensitivity` returns them, with each group's attrition in
#' `pretest_sensitivity$attrition`; [report_solomon()] names the center,
#' its assumption, and the other center. On the identity link the print also
#' gives the delta line: for each unpretested group, the difference between
#' its dropouts' and its completers' mean posttest that the default assumes,
#' which is the difference the pretested dropouts' pretests predict for
#' them, scaled by the ratio of the loss rates; `"observed"` assumes 0. This
#' is the delta parameterization of the National Research Council (2010, pp.
#' 88--89) and Mavridis et al. (2015, sec. 2), and it assumes that the
#' pretested losses are at random given the pretest. [check_solomon_missing()]
#' gives the attrition checks.
#'
#' **Choosing.** Choose the center before the outcomes are seen, and report
#' the other as a sensitivity analysis (National Research Council, 2010,
#' Recommendation 9, p. 4, and Recommendation 15, p. 106; White et al.,
#' 2011). Specify `"observed"` in advance when there are grounds to believe
#' that the unpretested groups' losses are unrelated to baseline and to the
#' posttest. The clear case is
#' unpretested groups that lost no posttests, with a row for everyone
#' randomized to them. A burdensome pretest that adds losses makes the
#' default doubtful, but it does not make `"observed"` right. Loss rates
#' that differ between the pretest conditions suggest that pretesting
#' changed who dropped out, in which case neither center is free of
#' assumptions.
#'
#' **Cautions.**
#' - Moving the center changes what the coefficient of `pretested`
#'   measures, and nothing else in the model. Only the pretest effects, and
#'   that coefficient with its standard error and covariances, depend on the
#'   center; the treatment contrasts, Pretest x Treatment, and the fitted
#'   values do not, on the model's scale. On the risk or rate
#'   scale of [marginal_solomon()], Treatment | pretested, the average
#'   treatment effect, and Pretest x Treatment depend on it as well as the
#'   pretest effects.
#' - The gap above assumes one pretest slope for all conditions, as the
#'   model has.
#' - If the treatment conditions lost different kinds of participants,
#'   neither center gives both conditional pretest effects: in the record's
#'   case S3 (losses by baseline in the control arms and at random in the
#'   treated arms, alike in both pretest conditions) the default misstated
#'   both, although not their average when
#'   the pretested groups kept equal numbers, and the observed center
#'   misstated one of them and their average. Pretest x Treatment is then
#'   biased at every center (issue #143). [check_solomon_missing()] tests
#'   whether the pretested conditions' selection differs and then gives the
#'   raw comparison of posttest means in each condition.
#' - Pretested participants with a posttest and no pretest are dropped at
#'   either center. When that loss relates to baseline, both centers are
#'   biased alike, and a fit without `y_pre` is not (issue #142).
#' - When loss depends on the posttest itself, both centers understate a
#'   real pretest effect: in a derivation reported on issue #135, a true
#'   effect of 0.3 gave 0.234 at the default and 0.174 at the observed
#'   center. Neither center protects a test of equivalence; see
#'   [equivalence_solomon()].
#' - The default's target is the pretest effect among the participants who
#'   would be analyzed in either pretest condition. It equals the effect in
#'   the randomized population only when the pretest effect does not vary
#'   with what drives the loss.
#'
#' @section Covariates:
#' `covariates` enter the model as main effects, beside the pretest. The
#' inference is computed with each numeric covariate standardized, centered
#' at its mean and divided by its standard deviation over the participants
#' in the model, and with the pretest, centered already, divided by the
#' power of two nearest its standard deviation over the pretested
#' participants in the model (issue #137). That is a reparameterization of
#' the same model, so it changes no contrast, standard error, degrees of
#' freedom, or test, and keeps them accurate when a covariate, such as a
#' calendar year, lies far from zero relative to its spread, or when a
#' covariate or the pretest is on a scale far from 1, where the covariances
#' of the model as given lose accuracy or cannot be computed. `model`,
#' `coefficients`, and `vcov` are those of the covariates and the pretest as
#' given, the intercept included. With a pretest or a numeric covariate, a
#' CR2 `vcov` is therefore a matrix rather than the clubSandwich object,
#' whose small-sample adjustments belong to the standardized model; that
#' model's coefficients and covariance are in `covariate_scaling`. The
#' warnings, such as one that a negative-binomial theta did not converge,
#' come from the standardized model too, as do `theta` and `dispersion`,
#' which do not depend on the parameterization. Far from zero, a fit to a
#' covariate as given can warn where the same model standardized does not.
#'
#' Covariates collinear with the other terms stop the fit wherever they
#' lie. Far from zero, the rounding error of the stored values can make
#' [stats::glm()] take collinear terms for distinct ones, such as `z + 1e6`
#' and `2 * z + 3`, or a time of testing that is the time of enrollment plus
#' a delay, with the times as dates, so the standardized model is checked
#' too. Each term is regressed on the terms before it, the covariates after
#' the treatment and pretest terms, weighted as [stats::glm()] weighs them.
#' It is taken to be collinear with them when its residual is within 4
#' times the rounding error that the regression carries into it, or within
#' 1e-11 of its norm, the tolerance of [stats::glm()]. That rounding error
#' is the term's own plus each other term's times the magnitude of its
#' coefficient, relative to the term's standard deviation, with the error
#' of each value taken to be `.Machine$double.eps` times its magnitude as
#' given (the pretest's before it is centered), and that of standardizing,
#' centering, and weighting it. The residual is computed from the stored
#' values, accumulated in twice the working precision, rather than taken
#' from the QR decomposition, whose own rounding the coefficients would
#' multiply too.
#'
#' The check guarantees two things.
#'
#' - An exact relation among the terms (the covariates, the pretest, the
#'   treatment and pretest indicators, their interaction, and the
#'   indicators of factors) stops the fit, however far from zero its terms
#'   lie, however their spreads differ, and in any order of the covariates,
#'   as long as its values carry no more rounding error than the check
#'   allows for. Values stored from exact data, or computed from one another
#'   in a few operations, do not exceed it. The message names a term of the
#'   relation: the last in that order, or the one that [stats::glm()] takes
#'   to be collinear in the model as given, which it does when the spread of
#'   a term falls below about 1e-11 of its distance from zero.
#' - A term that is not exactly collinear stops the fit only when its
#'   residual is within that bound or below the tolerance of
#'   [stats::glm()]. Near zero the bound grows with the coefficients: it is
#'   about 2.7 * 10^-15 of the term's norm times 1 plus the sum of the
#'   magnitudes of its coefficients, each times the norm of its term over
#'   the term's own. A covariate that nearly repeats another, or a
#'   combination of others with coefficients of order 1, such as a
#'   covariate and its copy in single precision, or covariates 10^-9
#'   standard deviations apart, has a bound of about 5 * 10^-15, far below
#'   that tolerance, which then decides as it does for the model as given;
#'   such covariates are fitted as [stats::glm()] fits them, whatever other
#'   covariates lie far from zero. The bound exceeds that tolerance once the
#'   sum exceeds about 4,000, as in nested near-collinearity. With
#'   `z2 = z1 + 1e-6 * e1` and `z3 = e1 + delta * e2`, all near zero, the
#'   sum is about 2 * 10^6 and the bound about 5 * 10^-9 of the norm of
#'   `z3`, so that the fit goes ahead at `delta = 1e-8` and stops at
#'   `delta = 3e-9` and below, where [stats::glm()] still fits the model: a
#'   residual within the bound cannot be told from that of an exact
#'   relation between rounded values.
#'
#' Its limits:
#'
#' - A relation computed with more rounding error, such as one among values
#'   rounded to fewer digits, is not exact in double precision, and its
#'   terms are fitted as nearly collinear ones.
#' - Nearly collinear covariates far enough from zero that their difference
#'   falls below about 1e-11 of their distance from zero, such as
#'   covariates 10^-6 standard deviations apart at 10^6 standard deviations
#'   from zero, are taken to be collinear by [stats::glm()] as given, which
#'   stops the fit too.
#' - The check runs once the model as given is fitted. When that fit fails
#'   itself, as [MASS::glm.nb()] did with "NA/NaN/Inf in 'x'" for an exact
#'   relation with a covariate 3 * 10^10 standard deviations from zero, the
#'   fit stops with that error, as in earlier versions, before the check
#'   can name a term.
#' - Fitting nearly collinear covariates is not the same as estimating
#'   their standard errors well. With covariates that differ by about 10^-7
#'   of their standard deviation or less, as a covariate near zero and its
#'   copy in single precision do, the HC3 and CR2 standard errors of the
#'   contrasts lose accuracy, as [sandwich::vcovHC()] and
#'   [clubSandwich::vcovCR()] form the covariance from the inverse of `X'WX`.
#'   In one check, under HC3 and CR2, they were off by up to 3% and 4% for
#'   such a copy (2% and 3% in earlier versions), 0.3% and 1.5% at 10^-7,
#'   and 76% and 28% at 10^-8. At 10^-9, HC3 standard errors were off by up
#'   to a factor of 11 or `NaN`, with a warning, and CR2 fits gave `NaN` or
#'   stopped with "the leading minor ... is not positive", as some did at
#'   10^-8. Earlier versions stopped with "system is computationally
#'   singular" from 10^-8 near zero, and from 10^-7 for covariates 5
#'   standard deviations from zero. The model-based standard errors, and
#'   those of the coefficients under every covariance, keep their accuracy.
#'
#' A numeric covariate is standardized whatever its class, such as that of a
#' variable with value labels read from an SPSS file, or of one wrapped in
#' `I()`, and so are dates and times, as the numbers they enter the model
#' as; factors, logical columns, and character columns enter as given.
#'
#' @section Designs with several treatments:
#' A Solomon N-group design crosses k treatments and a control with
#' pretesting, giving 2(k + 1) groups: six for two treatments and eight for
#' three (Steyn, 2009). Give `treat` as a factor or character vector of
#' conditions and name the control with `control`. One model is then fitted
#' to all the groups, with an indicator for each treatment, and the four
#' Solomon contrasts are estimated for each comparison:
#' - each treatment against the control (`contrasts = "control"`, the
#'   default);
#' - every pair of conditions (`contrasts = "pairwise"`);
#' - or comparisons given as weights over the conditions, such as the main
#'   effects of two treatments crossed factorially:
#'   `contrasts = list(Lecture = c(Lecture = 0.5, Both = 0.5, Service = -0.5,
#'   None = -0.5))`. The weights of each comparison must sum to zero.
#'
#' Omnibus Wald tests ask whether the conditions differ on each contrast. The
#' Pretest x Condition test asks whether pretesting changes the effect of any
#' treatment.
#'
#' The p-values of each contrast are adjusted across the comparisons, by
#' Holm's (1979) procedure by default. It controls the familywise error rate
#' "for any combination of true hypotheses" (p. 65). The confidence
#' intervals are not adjusted. The result has class `solomon_ngroup`.
#'
#' The effects table ends with the pretest effect in each condition (see
#' "The pretest effect"), with `comparison` naming the condition:
#' `Pretest effect | control` for the control and `Pretest effect | treated`
#' for each treatment, whose p-values are adjusted across the treatments.
#' The `Pretest main effect`, with `comparison` `"All conditions"`, is their
#' equal-weighted average over the k + 1 conditions. Steyn (2009) tests the
#' pretest main effect separately for each intervention, in the two-way
#' analysis of variance of that intervention's groups and the control groups
#' that [fit_solomon_steyn()] carries out; this row averages over all the
#' conditions and is adjusted for the pretest score.
#'
#' Published studies with several treatments analyzed them as overlapping
#' four-group designs: one for each treatment against the control (McCarthy &
#' Tucker, 2002), or one for each pair of conditions (Mai et al., 2020). Those
#' analyses reuse the same groups, so their tests are dependent, and each
#' extra analysis adds to the chance of a false finding. The joint model asks
#' each question once. In Mai et al.'s (2020) six-group study, the Pretest x
#' Condition test of the posttests alone (without the pretest as a
#' covariate), with conventional covariance, gives F(2, 127) = 1.86,
#' p = .161. [fit_solomon_steyn()] carries out the sequence of tests Steyn
#' (2009) proposed for these designs.
#'
#' `r lifecycle::badge("experimental")` The analysis of designs with several
#' treatments is experimental. In the package's simulation study (issue #45;
#' 112 scenarios with two or three treatments and 10 to 50 participants per
#' group, 5,000 replications each):
#' - the treatment contrasts were unbiased, and coverage of their 95%
#'   intervals was 0.939 to 0.967;
#' - the familywise error rates of the Holm-adjusted comparisons were at most
#'   0.059;
#' - the Pretest x Condition test and the test of the conditions averaged
#'   over pretest rejected a true null hypothesis in 0.036 to 0.055 of
#'   replications;
#' - with three treatments and 10 participants per group, the omnibus tests
#'   of Condition | pretested and Condition | unpretested rejected in 0.055
#'   to 0.069 of replications at the .05 level. The rule for error control
#'   set before the study was therefore not met, which is why the analysis
#'   is experimental. With groups that small, judge those two questions by
#'   the adjusted comparisons.
#'
#' The study did not cover the pretest effects, which were added later
#' (issue #104), binary or count outcomes, clustered designs, or comparisons
#' given as weights. It is reported in the article "Designs With Several
#' Treatments: Validating the Joint Model".
#'
#' @param y_post numeric posttest vector
#' @param treat 0/1 (or logical) treatment indicator (1 = treatment); or, for
#'   a design with several treatments, a factor or character vector of
#'   conditions, with the control named by `control`
#' @param pretested 0/1 (or logical) pretest indicator (1 = group received pretest)
#' @param y_pre numeric pretest vector: the pretest score for pretested
#'   participants and `NA` for the others
#' @param covariates optional data frame of additional covariates, or, with
#'   `data`, the names of its columns to use
#' @param robust character: "HC3" (default), "none", or "CR2" (cluster-robust; requires `cluster`)
#' @param cluster optional clustering id (e.g., class/site), one value per
#'   participant. CR2 fits refuse designs in which a Solomon cell contains a
#'   single cluster, because cluster and condition are then confounded; see
#'   [validate_solomon()]. A classed warning (`solomonR_small_df_warning`)
#'   flags Solomon contrasts whose Satterthwaite degrees of freedom are below
#'   4, where Tipton (2015) advises that p-values not be trusted.
#' @param family model family (default gaussian()), given as a family object,
#'   a family function, or its name; or `"negative_binomial"` for the NB2
#'   model for counts.
#' @param conf_level confidence level for intervals (default 0.95)
#' @param exposure optional positive exposure (for example, observation time)
#'   for each participant, entered as a log offset; requires a log-link family
#'   such as `poisson()` or `"negative_binomial"`. Contrasts are then log rate
#'   ratios per unit of exposure.
#' @param control the control condition, when `treat` is a factor or
#'   character vector. With two conditions the design is a four-group design
#'   and the result is the same as with a 0/1 `treat`; with three or more it
#'   is an N-group design (see "Designs with several treatments").
#' @param contrasts for designs with several treatments: `"control"` (each
#'   treatment against the control), `"pairwise"` (every pair of
#'   conditions), or a named list of weight vectors named by condition, one
#'   per comparison, each summing to zero.
#' @param adjust for designs with several treatments: the adjustment of the
#'   p-values of each contrast across the comparisons, `"holm"` (Holm, 1979;
#'   the default), `"bonferroni"`, or `"none"`.
#' @param data optional data frame. When supplied, the other data arguments
#'   are looked up in it first: give them as bare column names
#'   (`y_post = post`) or as strings (`y_post = "post"`).
#' @param pretest_center with `y_pre`, the mean pretest at which the
#'   pretested groups are compared with the unpretested groups (issue #135).
#'   `"analyzed"` (the default): the mean pretest of the pretested
#'   participants in the model. `"observed"`: the mean of every observed
#'   pretest among the pretested participants, with or without a posttest,
#'   condition, covariate, or exposure; their rows must be in the data, with
#'   `y_post` missing, or this gives the default. Only the pretest effects
#'   (and the coefficient of `pretested`, with its standard error and
#'   covariances) depend on it; the treatment contrasts, Pretest x
#'   Treatment, and the fitted values do not, on the model's scale. Ignored
#'   without `y_pre`.
#' @param y,pretest_score `r lifecycle::badge("deprecated")` Use `y_post`
#'   and `y_pre`.
#' @return An object of class `solomon_glm`, a list with:
#'   - `effects`: the four treatment contrasts, `ATE (avg over pretest)`,
#'     `Pretest x Treatment`, `Treatment | pretested`, and
#'     `Treatment | unpretested`, followed by the three pretest effects,
#'     `Pretest effect | control`, `Pretest effect | treated`, and
#'     `Pretest main effect`, in the columns `contrast`, `estimate`,
#'     `std.error`, `statistic`, `df` (`Inf` for a normal reference
#'     distribution), `p.value`, `conf.low`, and `conf.high`, followed by
#'     the Wald partial R-squared and its interval (`r2`, `r2_lo`, and
#'     `r2_hi`).
#'   - `conf_level`: the confidence level of the intervals.
#'   - `coefficients`: the coefficients of the model, in the same columns
#'     through `conf.high`, with `term` for `contrast`.
#'   - `model`: the fitted model, a `glm` object.
#'   - `vcov`: the covariance matrix of the coefficients, as `robust` sets
#'     it. With CR2 it is the clubSandwich object, or, with a pretest or a
#'     numeric covariate, a matrix (see "Covariates").
#'   - `covariate_scaling`: with a pretest or numeric covariates, the
#'     parameterization in which the inference was computed (see
#'     "Covariates"): the `means` and `scales` of the pretest, `pre_obs`
#'     (a mean of 0 and the power of two it was divided by), and of the
#'     numeric covariates (their means and standard deviations over the
#'     participants in the model), and the `coefficients` and the covariance
#'     `vcov` of the model with them standardized; otherwise `NULL`. Outside
#'     the stable interface (see [solomon_output]).
#'   - `dispersion`: the Pearson dispersion statistic for binomial, Poisson,
#'     and negative-binomial fits, and `NA` for the others.
#'   - `theta`: for negative-binomial fits, its estimate, standard error,
#'     and \eqn{\alpha = 1/\theta}; otherwise `NULL`.
#'   - `pretest_mean`: the mean pretest at which the pretest is centered, the
#'     center used; `NA` without `y_pre`.
#'   - `pretest_center`: the center used, `"analyzed"` or `"observed"`.
#'   - `pretest_mean_analyzed` and `pretest_mean_observed`: the two centers,
#'     the mean pretest of the pretested participants in the model and the
#'     mean of every observed pretest among the pretested participants; `NA`
#'     without `y_pre`.
#'   - `pretest_n`: the numbers of pretested participants with a pretest in
#'     the model (`analyzed`) and in all (`observed`), and, of those outside
#'     the model, how many have no posttest (`no_posttest`; the others lack
#'     the condition, a covariate, or the exposure). `NULL` without `y_pre`.
#'   - `pretest_sensitivity`: the pretest effects at both centers (see "The
#'     pretest effect when posttests are missing"), a list with `effects`
#'     (the pretest rows of the effects table at each `center`, the center
#'     used giving the rows of `effects`); `difference` (each pretest
#'     effect's estimates at the `analyzed` and `observed` centers and their
#'     `difference`); `slope`, the pretest slope; `m1`, the share of the
#'     pretested participants with a pretest who are not in the model;
#'     `pretest_means`, the mean pretest of those `analyzed` and
#'     `not_analyzed`; `delta`, on the identity link only, each unpretested
#'     group's `condition`, `n`, share not analyzed `m0`, and the offset
#'     `delta` the default assumes (`NA` where `m0` is 0), and otherwise
#'     `NULL`; `predicted`, what the pretested dropouts' pretests predict for
#'     them; `attrition`, each group's `n`, `posttest_missing`, and `share`;
#'     `differ`, whether the two centers differ; and `note`, the reason the
#'     other center's standard errors are missing, if they are. The other
#'     center's effects are computed without a refit, by moving the
#'     contrasts, with the variance of its mean stacked as for the center
#'     used. `NULL` without `y_pre`; outside the stable interface (see
#'     [solomon_output]).
#'   - `robust`, `family`, `cluster` (the cluster of each participant), and
#'     `n_clusters`: the settings used. The last two are `NULL` without
#'     `cluster`.
#'   - `call`, and `data` (the data the model was fitted to), which
#'     solomonR's own functions use.
#'
#'   For a design with several treatments, an object of class
#'   `solomon_ngroup`, with the same elements and these changes:
#'   - `effects` has a `comparison` column before `contrast` and the
#'     adjusted p-values `p.adjusted` after `conf.high`. It has the four
#'     treatment contrasts for each comparison, then the pretest effect in
#'     each condition and the pretest main effect.
#'   - `omnibus`: the omnibus tests, in the columns `test`, `statistic`,
#'     `df1`, `df2`, `p.value`, and `reference` (`"F"` or `"chisq"`).
#'   - `conditions`: the control and the treatments, with their model
#'     terms (`condition`, `role`, and `term`).
#'   - `weights`: the weights of each comparison over the conditions, one
#'     row for each comparison.
#'   - `adjust`: the adjustment of the p-values.
#'   - The tables of `pretest_sensitivity` have `comparison` first, as
#'     `effects` does, and its `effects` has `p.adjusted`; `delta` has a row
#'     for each unpretested group.
#'
#'   The `effects` table, `conf_level`, and [`tidy()`][solomon_output],
#'   which returns the table, are the stable interface of the result; see
#'   [solomon_output].
#' @references
#' Bell, R. M., & McCaffrey, D. F. (2002). Bias reduction in standard errors for
#' linear regression with multi-stage samples. *Survey Methodology, 28*(2),
#' 169–181.
#'
#' Cameron, A. C., & Trivedi, P. K. (2013). *Regression analysis of count data*
#' (2nd ed.). Cambridge University Press.
#' https://doi.org/10.1017/CBO9781139013567
#'
#' Campbell, D. T. (1957). Factors relevant to the validity of experiments in
#' social settings. *Psychological Bulletin, 54*(4), 297–312.
#' https://doi.org/10.1037/h0040950
#'
#' Campbell, D. T., & Stanley, J. C. (1966). *Experimental and
#' quasi-experimental designs for research*. Rand McNally. (Original work
#' published 1963)
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
#' Holm, S. (1979). A simple sequentially rejective multiple test procedure.
#' *Scandinavian Journal of Statistics, 6*(2), 65–70.
#' https://www.jstor.org/stable/4615733
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
#' Mai, N. N., Takahashi, Y., & Oo, M. M. (2020). Testing the effectiveness of
#' transfer interventions using Solomon four-group designs. *Education
#' Sciences, 10*(4), Article 92. https://doi.org/10.3390/educsci10040092
#'
#' Mavridis, D., White, I. R., Higgins, J. P. T., Cipriani, A., & Salanti, G.
#' (2015). Allowing for uncertainty due to missing continuous outcome data in
#' pairwise and network meta-analysis. *Statistics in Medicine, 34*(5),
#' 721–741. https://doi.org/10.1002/sim.6365
#'
#' McCarthy, A. M., & Tucker, M. L. (2002). Encouraging community service
#' through service learning. *Journal of Management Education, 26*(6), 629–647.
#' https://doi.org/10.1177/1052562902238322
#'
#' National Research Council. (2010). *The prevention and treatment of missing
#' data in clinical trials*. The National Academies Press.
#' https://doi.org/10.17226/12955
#'
#' Newman, I., Benz, C., & Williams, J. D. (1990). Alternatives in analyzing the
#' Solomon four group design. *Multiple Linear Regression Viewpoints, 17*(2),
#' 91–103. https://ojs.lib.ua.edu/glmj/article/view/125
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
#' Satterthwaite, F. E. (1946). An approximate distribution of estimates of
#' variance components. *Biometrics Bulletin, 2*(6), 110–114.
#' https://doi.org/10.2307/3002019
#'
#' Seaman, S., Galati, J., Jackson, D., & Carlin, J. (2013). What is meant by
#' "missing at random"? *Statistical Science, 28*(2), 257–268.
#' https://doi.org/10.1214/13-STS415
#'
#' Solomon, R. L. (1949). An extension of control group design. *Psychological
#' Bulletin, 46*(2), 137–150. https://doi.org/10.1037/h0062958
#'
#' Solomon, R. L., & Lessac, M. S. (1968). A control group design for
#' experimental studies of developmental processes. *Psychological Bulletin,
#' 70*(3, Pt. 1), 145–150. https://doi.org/10.1037/h0026147
#'
#' Stefanski, L. A., & Boos, D. D. (2002). The calculus of M-estimation. *The
#' American Statistician, 56*(1), 29–38.
#' https://doi.org/10.1198/000313002753631330
#'
#' Steiger, J. H. (2004). Beyond the F test: Effect size confidence intervals
#' and tests of close fit in the analysis of variance and contrast analysis.
#' *Psychological Methods, 9*(2), 164–182.
#' https://doi.org/10.1037/1082-989X.9.2.164
#'
#' Steyn, R. (2009). Re-designing the Solomon four-group: Can we improve on
#' this exemplary model? *Design Principles and Practices: An International
#' Journal—Annual Review, 3*(1), 383–394. https://doi.org/10.18848/1833-1874/CGP/v03i01/37588
#'
#' Tipton, E. (2015). Small sample adjustments for robust variance estimation
#' with meta-regression. *Psychological Methods, 20*(3), 375–393.
#' https://doi.org/10.1037/met0000011
#'
#' Venables, W. N., & Ripley, B. D. (2002). *Modern applied statistics with S*
#' (4th ed.). Springer. https://doi.org/10.1007/978-0-387-21706-2
#'
#' Welch, B. L. (1947). The generalization of "Student's" problem when several
#' different population variances are involved. *Biometrika, 34*(1–2), 28–35.
#' https://doi.org/10.1093/biomet/34.1-2.28
#'
#' White, I. R., Horton, N. J., Carpenter, J., & Pocock, S. J. (2011). Strategy
#' for intention to treat analysis in randomised trials with missing outcome
#' data. *BMJ, 342*, Article d40. https://doi.org/10.1136/bmj.d40
#' @examples
#' fit <- fit_solomon_glm(y_post, treat, pretested, y_pre, data = solomon_example)
#' fit
#'
#' # The four treatment contrasts and the three pretest effects, as a data
#' # frame; tidy(fit) returns the same table.
#' fit$effects
#'
#' # A six-group design: two treatments and a control (Mai et al., 2020).
#' fit6 <- fit_solomon_glm(post_behavior, condition, pretested, pre_behavior,
#'                         control = "Control", data = mai2020)
#' fit6
#'
#' # Posttests are missing in mai2020, so the print gives the pretest effects
#' # at both centers; pretest_sensitivity returns them.
#' fit6$pretest_sensitivity$difference
#' @export
fit_solomon_glm <- function(y_post, treat, pretested, y_pre = NULL,
                            covariates = NULL, robust = c("HC3", "none", "CR2"),
                            cluster = NULL, family = stats::gaussian(),
                            conf_level = 0.95, exposure = NULL,
                            control = NULL, contrasts = "control",
                            adjust = c("holm", "bonferroni", "none"),
                            data = NULL,
                            pretest_center = c("analyzed", "observed"),
                            y = deprecated(),
                            pretest_score = deprecated()) {
  .solomon_data_args(
    data,
    c("y_post", "treat", "pretested", "y_pre", "covariates", "cluster",
      "exposure", "y", "pretest_score"),
    environment(), parent.frame(), as_frame = "covariates"
  )
  if (lifecycle::is_present(y)) {
    .renamed_arg(!missing(y_post), "y", "y_post", "fit_solomon_glm")
    y_post <- y
  }
  if (lifecycle::is_present(pretest_score)) {
    .renamed_arg(!is.null(y_pre), "pretest_score", "y_pre", "fit_solomon_glm")
    y_pre <- pretest_score
  }
  robust <- match.arg(robust)
  pretest_center <- match.arg(pretest_center)
  .check_conf_level(conf_level)
  negbin <- identical(family, "negative_binomial")
  if (negbin) {
    if (robust == "CR2") {
      stop("Clustered (CR2) negative-binomial fits are not supported.", call. = FALSE)
    }
    # NB2 uses the log link; the Poisson family stands in for the link checks
    # below and is replaced by the fitted family afterwards.
    family <- stats::poisson()
  }
  if (is.character(family)) family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) family <- family()
  adjust <- match.arg(adjust)

  conditions <- .solomon_conditions(treat, control)
  if (conditions$k > 1L) {
    return(.fit_solomon_ngroup(
      y_post, conditions, pretested, y_pre, covariates, robust, cluster,
      family, negbin, conf_level, exposure, contrasts, adjust, match.call(),
      pretest_center = pretest_center
    ))
  }
  if (!(is.character(contrasts) && length(contrasts) == 1L &&
        contrasts %in% c("control", "pairwise"))) {
    stop(
      "`contrasts` applies to designs with more than one treatment; this ",
      "design has one treatment and a control.",
      call. = FALSE
    )
  }

  treat <- conditions$treat
  pretested <- .solomon_indicator(pretested, "pretested")

  .solomon_check_lengths(
    y_post = y_post,
    treat = treat,
    pretested = pretested,
    y_pre = y_pre,
    covariates = covariates,
    cluster = cluster,
    exposure = exposure
  )

  df <- data.frame(
    y = y_post,
    treat = treat,
    pretested = pretested
  )

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

    # Safe pretest covariate: equals y_pre in pretested rows, 0 otherwise
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

  # One model: treatment + pretest indicator + their interaction + (optional) pre_obs
  rhs <- c("treat*pretested",
           if (!is.null(y_pre)) "pre_obs" else NULL,
           if (!is.null(covariates)) .formula_names(names(covariates)) else NULL)
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
  # model, so that the pretesting coefficient compares pretested and
  # unpretested participants at that score (issue #104); or, with
  # `pretest_center = "observed"`, at the mean of every observed pretest
  # among the pretested participants, in the model or not (issue #135).
  pretest_mean <- NA_real_
  centers <- c(analyzed = NA_real_, observed = NA_real_)
  if (!is.null(y_pre)) {
    centers <- .pretest_centers(df)
    pretest_mean <- centers[[pretest_center]]
    df$pre_obs <- ifelse(df$pretested == 1, df$pre_obs - pretest_mean, 0)
  }

  # The model as given, which the result returns as `model`. With a pretest
  # or a numeric covariate the inference comes from the model refitted with
  # the covariates standardized and the pretest scaled, which is accurate
  # however far a covariate lies from zero, and whatever the scale of a
  # covariate or the pretest (issue #137; .covariate_scaling()). That fit
  # gives the warnings, and this one is fitted silently: far from zero, a
  # covariate can make the fit as given warn, for instance that theta did
  # not converge, where the same model standardized does not.
  standardize <- .standardizes(df, names(covariates))
  fit <- .fit_model(fml, df, family, negbin, quiet = standardize)
  in_model <- rep(TRUE, nrow(df))
  if (!is.null(fit$na.action)) in_model[fit$na.action] <- FALSE
  fit_given <- fit
  scaling <- .covariate_scaling(df, names(covariates), in_model)
  if (standardize) fit <- .fit_model(fml, scaling$data, family, negbin)
  if (negbin) family <- stats::family(fit)

  # An empty group, or a term collinear with the others, leaves a
  # coefficient inestimable, in the model as given or standardized: stop as
  # the several-treatment path does. The groups are named by the conditions
  # when `control` named them. With `pretest_center = "observed"` the center
  # exists even when a pretested group has no posttest, so the model's own
  # check is needed. The standardized model is checked as well, for terms
  # collinear but for the rounding error of the values given, which glm()
  # can take for distinct ones far from zero (.collinear_terms()).
  arms <- if (is.null(conditions$control)) {
    c("treatment", "control")
  } else {
    c(conditions$treatments, conditions$control)
  }
  group_counts <- table(
    factor(df$treat[in_model], levels = c(1L, 0L), labels = arms),
    factor(df$pretested[in_model], levels = c(1L, 0L))
  )
  .stop_inestimable(fit_given, .group_names(group_counts, group_counts == 0))
  .stop_inestimable(fit, .group_names(group_counts, group_counts == 0),
                    noise = .representation_error(fit, fit_given, scaling, pretest_mean),
                    covariates = .formula_names(names(covariates)))

  # A group with one participant is fitted exactly: its hat value is 1, and
  # HC3 divides the squared residual by (1 - h)^2, so its standard errors
  # are undefined. Stop as the several-treatment path does.
  if (robust == "HC3" && any(group_counts == 1)) {
    stop(
      "HC3 standard errors are undefined because these groups have a single ",
      "analyzable participant: ", paste(.group_names(group_counts, group_counts == 1),
                                         collapse = "; "), ". ",
      "Use `robust = \"none\"`; see validate_solomon().",
      call. = FALSE
    )
  }

  if (!is.null(y_pre) && !stats::family(fit)$link %in% c("identity", "log")) {
    .warn_noncollapsible(stats::family(fit)$link)
  }

  # Robust VCOV
  fit_cr <- fit
  n_clusters <- NULL
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
    # The clusters the CR2 covariance is computed from; `cluster` keeps one
    # value per input row.
    n_clusters <- length(unique(cluster_fit))
    single <- .cluster_structure(df$treat[used], df$pretested[used], cluster_fit)$single_cluster
    if (length(single)) {
      .stop_confounded_clusters(single)
    }

    # clubSandwich cannot use the NA-padded residuals of an na.exclude fit,
    # so CR2 quantities come from a complete-case refit of the same rows
    # (identical coefficients), with the covariates standardized as `fit`.
    fit_cr <- stats::glm(fml, data = scaling$data[used, , drop = FALSE], family = family)
    vcovM <- clubSandwich::vcovCR(fit_cr, cluster = cluster_fit, type = "CR2")
  } else {
    vcovM <- stats::vcov(fit)
  }

  # The sampling variance of the mean pretest at which the pretest effects
  # are evaluated, and its covariance with the coefficients (issue #104).
  # For the "observed" center the variance is taken over every pretested
  # participant with a pretest, and the covariance over those of them in the
  # model (issue #135). The other center's terms are computed too when the
  # two centers differ, for `pretest_sensitivity`. They are in the units of
  # the scaled pretest of `fit` (`scaling$data`), as are its coefficients.
  center <- NULL
  center_other <- NULL
  center_note <- NULL
  parts <- NULL
  differ <- FALSE
  if (!is.null(y_pre)) {
    rows <- rep(TRUE, nrow(df))
    if (!is.null(fit$na.action)) rows[fit$na.action] <- FALSE
    both <- .center_terms_both(scaling$data, y_pre, rows, if (robust == "CR2") fit_cr else fit,
                               robust, cluster, if (robust == "CR2") cluster_fit, pretest_center)
    center <- both$center
    center_other <- both$other
    center_note <- both$note
    parts <- both$parts
    differ <- both$differ
  }

  # --- reference distribution ---
  # CR2 tests use Satterthwaite degrees of freedom. Otherwise, as in
  # summary.glm(), tests use t with residual df when the dispersion is
  # estimated (e.g., Gaussian models) and the normal distribution when it is
  # fixed (binomial, Poisson).
  dispersion_fixed <- .fixed_dispersion(stats::family(fit))
  df_model <- if (dispersion_fixed) Inf else stats::df.residual(fit)

  # --- tidy coefficients with chosen vcov, for the covariates as given ---
  coefficients <- .coefficient_table(fit, fit_cr, vcovM, scaling, robust, df_model, conf_level)

  # --- linear contrasts (one row each) ---
  cf <- stats::coef(fit)
  cn <- names(cf)
  Z <- function() { v <- numeric(length(cn)); names(v) <- cn; v }

  # The pretest effects also carry the variance of the mean pretest at which
  # they are evaluated (.center_adjusted()); `added` is that part, which the
  # Wald R-squared uses too. `terms` is that mean's .center_terms(): the
  # center used, or the other center for `pretest_sensitivity`.
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

  # The pretest (testing) effects, pretested minus unpretested participants
  # with the centered pretest at zero, that is, at the pretested
  # participants' mean pretest (issue #104):
  # among controls, beta_pretested; among treated participants,
  # beta_pretested + beta_treat:pretested; and their equal-weighted average.
  L_pc <- Z()
  L_pc["pretested"] <- 1
  L_pt <- L_pc
  L_pm <- L_pc
  if ("treat:pretested" %in% cn) {
    L_pt["treat:pretested"] <- 1
    L_pm["treat:pretested"] <- 0.5
  }

  Ls <- stats::setNames(
    list(L_ate, L_int, L_pre, L_un, L_pc, L_pt, L_pm),
    c(.solomon_contrast_order, .solomon_pretest_order)
  )

  # Wald-based partial R^2 for each contrast (intervals only for
  # conventional covariance)
  conventional <- robust == "none"
  effects <- do.call(rbind, lapply(names(Ls), function(label) {
    est <- lin_contrast(Ls[[label]])
    r2 <- contrast_r2_ci(fit, Ls[[label]], vcovM, conf_level, conventional,
                         added_variance = est[["added"]])
    data.frame(
      contrast = label,
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
  }))
  effects <- .effects_table(effects)

  if (robust == "CR2") {
    small_df <- is.finite(effects$df) & effects$df < 4
    if (any(small_df)) {
      .warn_cr2_small_df(effects$contrast[small_df], effects$df[small_df])
    }
  }

  # The pretest effects at both centers (issue #135). The other center's are
  # computed without a refit (.other_center_rows()), with the stacked
  # variance of the other mean. This is exact, because the fitted values,
  # the HC3, CR2, and model-based covariances, and the Satterthwaite degrees
  # of freedom of a linear function of the estimates do not depend on where
  # the pretest is centered. When the two centers are the same number, the
  # rows are copied. The shift between the centers is in the units of the
  # scaled pretest, and the slope in those of the pretest as given.
  sensitivity <- NULL
  if (!is.null(y_pre)) {
    used_rows <- effects[effects$contrast %in% .solomon_pretest_order,
                         .solomon_effect_columns, drop = FALSE]
    other_rows <- NULL
    if (differ) {
      other_rows <- .other_center_rows(
        Ls[.solomon_pretest_order], .solomon_pretest_order,
        (centers[[.other_center(pretest_center)]] - pretest_mean) / .pretest_scale(scaling),
        function(L) lin_contrast(L, center_other), no_se = !is.null(center_note)
      )
    }
    identity <- identical(stats::family(fit)$link, "identity")
    slope <- cf[["pre_obs"]] / .pretest_scale(scaling)
    sensitivity <- .pretest_sensitivity(
      used_rows, other_rows, pretest_center, character(), parts, centers, slope, identity,
      .pretest_attrition(ifelse(df$treat == 1L, arms[1], arms[2]), rev(arms), df$pretested,
                         df$y, rows, slope, centers, identity),
      note = center_note
    )
  }

  # The dispersion and theta do not depend on the parameterization; they
  # come from the fit the inference came from.
  dispersion <- if (.fixed_dispersion(stats::family(fit))) {
    sum(stats::residuals(fit, type = "pearson")^2, na.rm = TRUE) / stats::df.residual(fit)
  } else {
    NA_real_
  }

  out <- list(
    model = fit_given,
    coefficients = coefficients$table,
    effects = effects,
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
    call = match.call()
  )

  class(out) <- "solomon_glm"

  out
}

#' Permutation test for a Solomon contrast
#'
#' `r lifecycle::badge("stable")`
#' Performs a randomization-based test by permuting treatment assignment
#' within pretest strata. This preserves the Solomon four-group design while
#' generating the null distribution for a selected treatment contrast.
#' Participants are permuted in unclustered designs and whole clusters in
#' clustered ones.
#'
#' @details
#' **Unclustered designs.** Treatment labels of participants are permuted
#' within pretest strata and the model is refitted for each permutation. The
#' default statistic is the HC3-studentized contrast. The permutation
#' p-value is a valid test of the sharp null hypothesis that treatment has no
#' effect for any participant; the `+1` correction keeps the Monte Carlo
#' p-value from being zero (Phipson & Smyth, 2010). That result is for the
#' count of permuted statistics at least as extreme as the observed one, so
#' a permutation whose statistic equals the observed one is counted. Such
#' ties occur when a permutation repeats the observed assignment or swaps
#' arms of equal size, and often when scores are tied, as with binary
#' outcomes, counts, and ratings. So that rounding error does not decide,
#' statistics that differ by less than a small relative tolerance are
#' treated as equal: 1e-10 for a linear model and 1e-6 for models fitted by
#' iteration. The tolerance is relative to the observed statistic or, when
#' that is smaller, to 1 for the studentized statistic and to the largest
#' permuted difference for the difference statistic. In a model fitted by
#' iteration, equal fits usually agree to about 1e-7; with a covariate they
#' can differ by more, most with a link other than the canonical one, and
#' such a tie can still be missed (issue #134). Studentizing the
#' statistic makes permutation tests asymptotically robust when only an
#' average effect is hypothesized to be zero (DiCiccio & Romano, 2017;
#' Wu & Ding, 2021); for the Pretest x Treatment contrast that robustness
#' should be regarded as approximate.
#'
#' **The difference statistic.** `statistic = "difference"` tests only the
#' sharp null hypothesis. A permutation test of a difference in means is
#' exact when the two arms' outcomes have the same distribution, but when
#' only the average effect is zero it keeps its level, even asymptotically,
#' only if the arms are equal in size or in variance (Romano, 1990). In a
#' check for issue #113 (1,000 replications of 199 permutations; the average
#' treatment effect; 8 treated and 24 control participants in each pretest
#' condition; a treated standard deviation twice the control one; no
#' effect), the Type I error at .05 was 0.144 (Monte Carlo standard error
#' 0.011) for the difference statistic and 0.064 (0.008) for the studentized
#' statistic. A classed warning (`solomonR_unbalanced_arms_warning`) is
#' therefore given when the difference statistic is chosen and treated and
#' control participants differ in number in a pretest condition that the
#' contrast uses. The studentized default is recommended.
#'
#' **Clustered designs.** Randomization inference must permute the unit that
#' was randomized, so for fits with a `cluster` variable the treatment labels
#' of whole clusters are permuted. Two assignment mechanisms are supported:
#' - whole clusters assigned to the four Solomon conditions, as in Kvalem et
#'   al. (1996), where treatment labels are permuted among clusters within
#'   each pretest condition; and
#' - treatment assigned to clusters and pretesting to participants within
#'   clusters, where treatment labels are permuted among all clusters and
#'   every cluster must contain pretested and unpretested participants.
#'
#' Designs that assign treatment to participants within clusters are refused,
#' and stratified or restricted randomization of clusters is not supported.
#'
#' The statistic is built from cluster-level summaries (Gail et al., 1996;
#' Hayes & Moulton, 2017, ch. 10). A Stage 1 model with every term of the fit
#' except treatment (the pretest indicator, the pretest score, covariates and
#' any exposure offset) gives each cluster a covariate-adjusted difference
#' residual: observed minus expected, divided by the number of participants
#' or, for counts, by the total exposure (Hayes & Moulton, 2017, pp.
#' 221--224, following Bennett et al., 2002). Because Stage 1 ignores
#' treatment, the test remains exact under the sharp null hypothesis. When
#' pretesting is assigned within clusters, each cluster has a pretested and
#' an unpretested residual, and each contrast compares treated and control
#' clusters on one combination of the two.
#'
#' The contrast is estimated from unweighted means of the cluster residuals,
#' so each cluster counts once (Hayes & Moulton, 2017, pp. 202--205). It
#' therefore estimates the average effect across clusters, which can differ
#' from the participant-weighted contrast of [fit_solomon_glm()] when cluster
#' sizes vary and effects depend on cluster size. The studentized statistic
#' divides the contrast by its separate-variances standard error (Hayes &
#' Moulton, 2017, p. 212), which is the studentization Wu and Ding (2021)
#' use for weak null hypotheses, applied here with clusters as the units.
#' Gail et al. (1996, p. 1079) showed that the unstudentized difference can
#' exceed the nominal level under the weak null hypothesis when the arms have
#' unequal numbers of clusters and unequal variances.
#'
#' When the number of possible allocations is at most `reps`, all of them are
#' enumerated and the p-value is exact; the smallest attainable p-value is
#' then reported. Hayes and Moulton (2017, p. 239) note that at least four
#' clusters per arm are needed for a two-sided p below .05.
#'
#' **Simulation evidence.** In the package's pre-registered simulation study
#' (issue #19; 96 scenarios, 2,000 replications each):
#' - Under the sharp null hypothesis, both statistics had Type I errors of at
#'   most 0.063, consistent with the exactness of randomization tests. With
#'   four clusters per arm, few allocations exist and the test is
#'   conservative: the attainable level at .05 is 2/70, about 0.03.
#' - When treatment made treated clusters four times as variable as control
#'   clusters, so that only the average effect was zero, the studentized
#'   statistic's Type I error reached 0.065 with equal numbers of treated and
#'   control clusters, 0.0685 with 15 treated and 47 control clusters per
#'   pretest condition, and 0.0805 with 4 and 8. The difference statistic
#'   reached 0.158 with 15 and 47, as Gail et al. (1996) found for unbalanced
#'   designs. A classed warning (`solomonR_unbalanced_clusters_warning`) is
#'   therefore given, for either statistic, whenever treated and control
#'   clusters differ in number. It also has the class
#'   `solomonR_unbalanced_arms_warning` of the participant-level warning, so
#'   one handler can catch both.
#' - The CR2 tests of [fit_solomon_glm()] reached 0.068 with four clusters per
#'   arm even under the sharp null, where this test is exact, so the
#'   permutation test is preferred for designs with few clusters.
#'
#' @param fit An object returned by \code{fit_solomon_glm()}. Designs with
#'   several treatments are not supported; see [fit_solomon_glm()].
#' @param contrast Character string identifying the contrast to test. One of
#'   \code{"ATE (avg over pretest)"}, \code{"Pretest x Treatment"},
#'   \code{"Treatment | pretested"}, or
#'   \code{"Treatment | unpretested"}. The pretest effects of the fit are not
#'   tested: permuting treatment labels says nothing about them, and with a
#'   pretest covariate the pretest labels cannot be permuted. For a binary
#'   fit with a pretest covariate on a noncollapsible link such as the logit,
#'   \code{"Pretest x Treatment"} gives a classed warning
#'   (`solomonR_link_scale_warning`): on that scale the contrast is nonzero
#'   whenever the pretest predicts the outcome, even without sensitization
#'   (Daniel et al., 2021), so a rejection need not reflect sensitization;
#'   see [marginal_solomon()].
#' @param reps Number of permutations. Default is 5000. In clustered designs
#'   with at most `reps` possible allocations, every allocation is used.
#' @param seed Optional random-number seed for reproducibility. The global
#'   random-number state is restored when the function exits.
#' @param return_dist Logical. If \code{TRUE}, return the permutation
#'   distribution in addition to the observed statistic and p-value.
#' @param statistic `"studentized"` (the default, and recommended) divides
#'   the contrast by its standard error; `"difference"` uses the contrast
#'   itself, which tests only the sharp null hypothesis (see "The difference
#'   statistic").
#' @param object `r lifecycle::badge("deprecated")` Use `fit`.
#'
#' @return An object of class `solomon_perm`, a list with:
#'   - `contrast`: the label of the contrast tested.
#'   - `statistic`: the statistic permuted, `"studentized"` or
#'     `"difference"`.
#'   - `level`: the level permuted, `"participant"` or `"cluster"`.
#'   - `estimate`: the estimated contrast.
#'   - `z_obs`: the observed statistic (the contrast itself when
#'     `statistic = "difference"`).
#'   - `p_perm`: the permutation p-value.
#'   - `reps` and `valid_reps`: the number of permutations, and the number
#'     that gave a usable statistic.
#'   - `exact`: whether every possible allocation was used.
#'   - `design`, `n_allocations` (the number of possible allocations),
#'     `min_p` (the smallest attainable p-value when exact), and `clusters`
#'     (the numbers of treated and control clusters in each stratum), for
#'     clustered fits.
#'   - `z_perm`: the permutation distribution, with `return_dist = TRUE`.
#'
#'   The result tests one contrast and gives no interval, so it has no
#'   effects table; see [solomon_output].
#'
#' @references
#' Bennett, S., Parpia, T., Hayes, R., & Cousens, S. (2002). Methods for the
#' analysis of incidence rates in cluster randomized trials. *International
#' Journal of Epidemiology, 31*(4), 839–846.
#' https://doi.org/10.1093/ije/31.4.839
#'
#' Daniel, R., Zhang, J., & Farewell, D. (2021). Making apples from oranges:
#' Comparing noncollapsible effect estimators and their standard errors after
#' adjustment for different covariate sets. *Biometrical Journal, 63*(3),
#' 528–557. https://doi.org/10.1002/bimj.201900297
#'
#' DiCiccio, C. J., & Romano, J. P. (2017). Robust permutation tests for
#' correlation and regression coefficients. *Journal of the American Statistical
#' Association, 112*(519), 1211–1220.
#' https://doi.org/10.1080/01621459.2016.1202117
#'
#' Gail, M. H., Mark, S. D., Carroll, R. J., Green, S. B., & Pee, D. (1996).
#' On design considerations and randomization-based inference for community
#' intervention trials. *Statistics in Medicine, 15*(11), 1069–1092.
#' https://doi.org/10.1002/(SICI)1097-0258(19960615)15:11%3C1069::AID-SIM220%3E3.0.CO;2-Q
#'
#' Hayes, R. J., & Moulton, L. H. (2017). *Cluster randomised trials* (2nd
#' ed.). Chapman and Hall/CRC. https://doi.org/10.4324/9781315370286
#'
#' Kvalem, I. L., Sundet, J. M., Rivø, K. I., Eilertsen, D. E., & Bakketeig,
#' L. S. (1996). The effect of sex education on adolescents' use of condoms:
#' Applying the Solomon four-group design. *Health Education Quarterly,
#' 23*(1), 34–47. https://doi.org/10.1177/109019819602300103
#'
#' Phipson, B., & Smyth, G. K. (2010). Permutation p-values should never be
#' zero: Calculating exact p-values when permutations are randomly drawn.
#' *Statistical Applications in Genetics and Molecular Biology, 9*(1),
#' Article 39. https://doi.org/10.2202/1544-6115.1585
#'
#' Romano, J. P. (1990). On the behavior of randomization tests without a
#' group invariance assumption. *Journal of the American Statistical
#' Association, 85*(411), 686–692.
#' https://doi.org/10.1080/01621459.1990.10474928
#'
#' Wu, J., & Ding, P. (2021). Randomization tests for weak null hypotheses in
#' randomized experiments. *Journal of the American Statistical Association,
#' 116*(536), 1898–1913. https://doi.org/10.1080/01621459.2020.1750415
#'
#' @examples
#' fit <- with(solomon_example, fit_solomon_glm(y_post, treat, pretested, y_pre))
#' # Few permutations keep the example fast; use the default for analyses.
#' perm_solomon(fit, reps = 199, seed = 1)
#' @export
perm_solomon <- function(
    fit,
    contrast = "ATE (avg over pretest)",
    reps = 5000L,
    seed = NULL,
    return_dist = FALSE,
    statistic = c("studentized", "difference"),
    object = deprecated()
) {
  if (lifecycle::is_present(object)) {
    .renamed_arg(!missing(fit), "object", "fit", "perm_solomon")
    fit <- object
  }

  .stop_ngroup_fit(fit, "perm_solomon")

  if (!inherits(fit, "solomon_glm")) {
    stop("`fit` must be from fit_solomon_glm().", call. = FALSE)
  }

  statistic <- match.arg(statistic)

  valid_contrasts <- c(
    "ATE (avg over pretest)",
    "Pretest x Treatment",
    "Treatment | pretested",
    "Treatment | unpretested"
  )

  if (length(contrast) == 1L && contrast %in% .solomon_pretest_order) {
    stop(
      "perm_solomon() permutes treatment labels within pretest conditions, so it ",
      "tests the treatment contrasts, not the pretest effects, which compare the ",
      "pretest conditions. With a pretest covariate, pretest labels cannot be ",
      "permuted, because unpretested participants have no pretest score. The fit's ",
      "own tests of the pretest effects are in `fit$effects`.",
      call. = FALSE
    )
  }

  if (!contrast %in% valid_contrasts) {
    stop(
      "Unknown contrast. Choose one of: ",
      paste(valid_contrasts, collapse = ", ")
    )
  }

  # The Pretest x Treatment contrast on a noncollapsible link (issue #114).
  .warn_link_scale(fit, contrast)

  reps <- as.integer(reps)

  if (length(reps) != 1L || is.na(reps) || reps < 1L) {
    stop("reps must be a positive integer.")
  }

  if (!is.null(seed)) {
    withr::local_seed(seed)
  }

  if (!is.null(fit$cluster)) {
    return(.perm_solomon_cluster(fit, contrast, reps, return_dist, statistic))
  }

  # The model is refitted with the numeric covariates standardized and the
  # pretest scaled as the fit's inference was (issue #137;
  # .covariate_scaling()). The statistic is the same, and accurate however
  # far a covariate lies from zero.
  df <- fit$data
  if (!is.null(fit$covariate_scaling)) df <- .standardized_frame(df, fit$covariate_scaling)
  form <- stats::formula(fit$model)
  fam <- stats::family(fit$model)

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

    L <- make_contrast(
      coef_names = names(b),
      contrast = contrast
    )

    Lm <- matrix(L, nrow = 1)

    estimate <- as.numeric(Lm %*% b)

    if (statistic == "difference") {
      return(estimate)
    }

    V <- sandwich::vcovHC(
      fit,
      type = "HC3"
    )

    # Force coefficient and covariance-matrix order to match.
    V <- V[names(b), names(b), drop = FALSE]

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

  estimate <- as.numeric(
    make_contrast(names(stats::coef(fit$model)), contrast) %*%
      stats::coef(fit$model)
  )

  # One routine refits the model for the observed and for the permuted
  # labels, so that labelings that give the same fit in exact arithmetic
  # are computed the same way. Convergence stays at glm()'s default, and
  # the tie tolerance allows for it (see .perm_tie_tolerance(); issue #131).
  refit <- function(d) {
    stats::glm(
      form,
      data = d,
      family = fam,
      na.action = stats::na.exclude,
      control = stats::glm.control(maxit = 200)
    )
  }
  tie_tol <- .perm_tie_tolerance(fam)

  z_obs <- contrast_z(
    refit(df),
    contrast
  )

  if (!is.finite(z_obs)) {
    stop("Could not calculate the observed permutation-test statistic.")
  }

  # The difference statistic tests only the sharp null hypothesis; with
  # arms of unequal size in a pretest condition the contrast uses, it can
  # reject too often when only the average effect is zero (Romano, 1990;
  # issue #113). The counts are of the participants in the model.
  if (statistic == "difference") {
    used <- rep(TRUE, nrow(df))
    if (!is.null(fit$model$na.action)) used[fit$model$na.action] <- FALSE
    strata <- c(pretested = 1L, unpretested = 0L)[.perm_contrast_weights(contrast) != 0]
    arm_count <- function(arm) {
      vapply(strata, function(p) sum(df$treat[used] == arm & df$pretested[used] == p), numeric(1))
    }
    n_treated <- arm_count(1L)
    n_control <- arm_count(0L)
    if (any(n_treated != n_control)) .warn_unbalanced_arms(n_treated, n_control)
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

    df_perm$treat <- .perm_treat(df$treat, df$pretested)

    z_perm[i] <- contrast_z(
      refit(df_perm),
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
  # Ties with the observed statistic are counted (see .perm_at_least()).

  p_perm <- (
    sum(.perm_at_least(z_perm_valid, z_obs, .perm_unit(z_perm_valid, statistic), tie_tol)) + 1
  ) / (
    length(z_perm_valid) + 1
  )

  out <- list(
    contrast = contrast,
    statistic = statistic,
    level = "participant",
    estimate = estimate,
    z_obs = z_obs,
    p_perm = p_perm,
    reps = reps,
    valid_reps = length(z_perm_valid),
    exact = FALSE
  )

  if (isTRUE(return_dist)) {
    out$z_perm <- z_perm_valid
  }

  class(out) <- "solomon_perm"

  out
}
