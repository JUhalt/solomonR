# Pretested participants with a posttest and no pretest (issue #142).
#
# fit_solomon_glm() handles them in one of two ways, which `missing_pretest`
# chooses.
# - "drop": they are left out of the model, as in every earlier version.
#   The pretest effects then compare the pretested participants who have a
#   pretest with all the unpretested participants.
# - "indicator": they are kept by the missing-indicator method (Groenwold et
#   al., 2012; Carpenter & Kenward, 2007, pp. 37-38). Their centered pretest
#   is set to 0, and a 0/1 column, .missing_pretest_column, marks them. The
#   indicator is common to all the conditions. Both pretest centers are
#   computed from observed pretests only. Each pretest effect is the
#   `pretested` contrast plus the indicator's mean among the pretested
#   participants in the model times the indicator's coefficient, which
#   standardizes it over those participants. That step is the package's
#   own, by analogy with Morris et al. (2022, p. 14).
# The decision record posted on issue #142 sets the method out (its sections
# 3 to 5), with the identities this file's comments name (i) and (iii).
#
# Notation, as in the record: b is the pretest slope, g the indicator's
# coefficient, m the share of the pretested participants in the model who
# have no pretest (the indicator's mean), and p = 1 - m.

# The default of `missing_pretest`, for the functions that follow the
# default of fit_solomon_glm(): its first value.
.missing_pretest_default <- function() {
  eval(formals(fit_solomon_glm)$missing_pretest)[[1L]]
}

# Whether a fit holds the indicator: it was fitted with missing_pretest =
# "indicator" and at least one pretested participant in the model has no
# pretest. With none, the fit is the "drop" fit, with no indicator term. A
# fit made by an earlier version can have a covariate named as the
# indicator's column; it has no `missing_pretest`, and is not taken for one.
.has_missing_indicator <- function(fit) {
  identical(fit$missing_pretest, "indicator") &&
    isTRUE(sum(fit$n_missing_pretest) > 0L) &&
    .missing_pretest_column %in% names(fit$data)
}

# The rows of a fit's data that the indicator marks: the pretested
# participants in the model who have no pretest. NULL for a fit without the
# indicator.
.kept_rows <- function(fit) {
  if (!.has_missing_indicator(fit)) return(NULL)
  fit$data[[.missing_pretest_column]] %in% 1L
}

# The share of the pretested participants in the model whom the indicator
# marks, for a fit with the indicator.
.kept_share <- function(fit) fit$missing_pretest_share

# The names of the columns of a fit's data that are not covariates: the
# outcome, the treatment and pretest terms, the pretest, the exposure, and,
# in a fit with the indicator, the indicator.
.design_names <- function(fit) {
  c("y", "treat", "pretested", "pre_obs", "log_exposure",
    if (.has_missing_indicator(fit)) .missing_pretest_column)
}

# The pretested participants whom "indicator" keeps: those without a pretest
# who have every other variable of the model, so that only the pretest keeps
# them out of the "drop" fit. `df` holds the model's variables, with
# `pre_obs` still missing for them; one logical value per row. All FALSE
# without a pretest.
.missing_pretest_rows <- function(df) {
  if (!"pre_obs" %in% names(df)) return(rep(FALSE, nrow(df)))
  others <- df[setdiff(names(df), "pre_obs")]
  stats::complete.cases(others) & df$pretested %in% 1L & is.na(df$pre_obs)
}

# The model data `df` with the indicator's column after the pretest: 1 for
# the participants `kept` flags and 0 for everyone else.
.add_missing_indicator <- function(df, kept) {
  column <- .missing_pretest_column
  df[[column]] <- as.integer(kept)
  at <- match("pre_obs", names(df))
  df[c(names(df)[seq_len(at)], column, setdiff(names(df)[-seq_len(at)], column))]
}

# The classed warning on pretested participants without a pretest, which
# names the handling. `n` counts them, and `with_posttest` those of them
# who have a posttest, whom "indicator" keeps.
.warn_missing_pretest <- function(n, with_posttest, handling) {
  text <- if (handling == "drop") {
    paste0(
      n, " pretested participant(s) have missing pretest scores and are excluded ",
      "from the model (missing_pretest = \"drop\"). Incidental missingness is not ",
      "imputed. Excluding them assumes that, in each pretested group, those without ",
      "a pretest would have had the same mean posttest as those with one; ",
      "missing_pretest = \"indicator\" keeps those who have a posttest. See ",
      "\"Pretested participants without a pretest\" in ?fit_solomon_glm."
    )
  } else if (with_posttest > 0L) {
    paste0(
      n, " pretested participant(s) have missing pretest scores. ",
      if (with_posttest == n) "They" else paste("The", with_posttest, "with a posttest"),
      " are kept by the missing-indicator method (missing_pretest = \"indicator\"; ",
      "Groenwold et al., 2012): their centered pretest is set to 0, and the column ",
      .missing_pretest_column, " marks them. This assumes that losing the pretest did ",
      "not depend on the treatment condition. See \"Pretested participants without a ",
      "pretest\" in ?fit_solomon_glm."
    )
  } else {
    paste0(
      n, " pretested participant(s) have missing pretest scores. None of them has a ",
      "posttest, so none is in the model, and missing_pretest = \"indicator\" has ",
      "no one to keep."
    )
  }
  warning(structure(
    class = c("solomonR_missing_pretest_warning", "warning", "condition"),
    list(message = text, call = NULL)
  ))
}

# Whether the indicator separates the outcome: every participant it marks
# has the same outcome, at a bound of the family's mean. The indicator's
# coefficient then has no finite maximum-likelihood estimate (Albert &
# Anderson, 1984): it runs to minus infinity when all their outcomes are 0
# (a binary outcome, or a count under the log link), and to plus infinity
# when all are 1 (a binary outcome). Returns -1, 1, or 0 for no separation.
# `y` holds the outcomes of the marked participants as glm() fits them.
.indicator_separation <- function(family, y) {
  if (!length(y) || anyNA(y)) return(0L)
  name <- family$family
  binary <- name %in% c("binomial", "quasibinomial")
  count <- (name %in% c("poisson", "quasipoisson") || .is_negbin(family)) &&
    identical(family$link, "log")
  if ((binary || count) && all(y == 0)) return(-1L)
  if (binary && !identical(family$link, "log") && all(y == 1)) return(1L)
  0L
}

# The classed warning of a fit whose indicator separates the outcome
# (.indicator_separation()): the estimates that weight the indicator's
# coefficient, the pretest effects, are infinite on the model's scale and
# have no standard error.
.warn_indicator_separation <- function(k, direction, family) {
  outcome <- if (family$family %in% c("binomial", "quasibinomial")) {
    if (direction < 0) "the outcome 0" else "the outcome 1"
  } else {
    "a count of 0"
  }
  warning(structure(
    class = c("solomonR_missing_pretest_separation_warning", "warning", "condition"),
    list(message = paste0(
      if (k == 1L) {
        paste0("The pretested participant kept without a pretest has ", outcome)
      } else {
        paste0("All ", k, " pretested participants kept without a pretest have ", outcome)
      },
      ", so the coefficient of ", .missing_pretest_column, " has no finite estimate on ",
      "the ", .separation_scale(family), " scale (separation; Albert & Anderson, 1984). ",
      "The pretest effects with them kept are ", if (direction < 0) "minus " else "",
      "infinity on that scale and have no standard error, test, or interval; the ",
      "treatment contrasts keep theirs. marginal_solomon() estimates the pretest ",
      "effects from standardized ", if (.is_count_family(family)) "rates" else "risks",
      ", which are defined."
    ), call = NULL)
  ))
}

.is_count_family <- function(family) {
  family$family %in% c("poisson", "quasipoisson") || .is_negbin(family)
}

# The scale of a model's coefficients, as in "log-odds", for the messages
# on separation.
.separation_scale <- function(family) {
  switch(family$link, logit = "log-odds", log = "log", paste(family$link, "link"))
}

# The coefficient table and the covariance of a fit whose indicator
# separates the outcome (`coefficients`, from .coefficient_table()): the
# indicator's coefficient is plus or minus infinity (`direction`), the limit
# of the maximum-likelihood fit, with no standard error, test, or interval,
# and its entries of the covariance are NA. The model itself, and the
# standardized coefficients and covariance the inference came from, keep
# the values at which the fit stopped, which the functions that take the
# fit compute from.
.separated_coefficients <- function(coefficients, direction) {
  row <- match(.missing_pretest_column, coefficients$table$term)
  coefficients$table$estimate[row] <- direction * Inf
  for (column in c("std.error", "statistic", "p.value", "conf.low", "conf.high")) {
    coefficients$table[[column]][row] <- NA_real_
  }
  at <- match(.missing_pretest_column, rownames(coefficients$vcov))
  coefficients$vcov[at, ] <- NA_real_
  coefficients$vcov[, at] <- NA_real_
  coefficients
}

# `pretest_sensitivity` of a fit whose indicator separates the outcome: the
# pretest effects are infinite at both centers, and their difference, which
# is that of two infinities, is set to its value in the fitted model, the
# slope times the share with a pretest (`slope_gap`) times the difference
# between the centers.
.separated_sensitivity <- function(sensitivity, slope_gap, centers) {
  sensitivity$difference$difference <- rep(
    slope_gap * (centers[["analyzed"]] - centers[["observed"]]),
    nrow(sensitivity$difference)
  )
  sensitivity
}

# The variance of a sum of sample means, and its covariance with the
# coefficients of a GLM fit, when the samples can share participants.
#
# .center_terms() gives them for samples that share no participant. With the
# indicator, a pretest effect depends, to first order, on two means over
# overlapping samples: the mean pretest, over the pretested participants
# with a pretest, and the indicator's mean, over all the pretested
# participants in the model. Their estimating equations are stacked with
# the model's (Stefanski & Boos, 2002), as .center_terms() stacks one
# mean's, in the covariance type of the fit.
#
# `samples` is a list, one element for each mean, each with `a` (one value
# per participant, already multiplied by the mean's weight in the contrast)
# and `member` (which participants the sample holds). `used` flags the
# participants the model used, in the model's row order, and `cluster`
# gives the cluster of each participant for CR2. Returns `v`, the variance
# of the sum of the means; `cross`, its covariance with the coefficients;
# and, for CR2, `df`, the Satterthwaite degrees of freedom of `v`.
# - "none" (model-based): each participant contributes the sum, over the
#   samples that hold them, of their deviation from the sample's mean
#   divided by sqrt(n (n - 1)), and `v` is the sum of the squares. For one
#   sample that is s^2 / n. No covariance with the coefficients.
# - "HC3": the jackknife form, with the deviations divided by n - 1, and
#   the covariance from the HC3-scaled influences.
# - "CR2": the CR2 variance (Bell & McCaffrey, 2002) of the sum of the
#   means of a regression of the stacked values on their sample, clustered
#   by the participants' clusters, with Satterthwaite degrees of freedom
#   (Pustejovsky & Tipton, 2018); the covariance from cluster sums with the
#   factor G / (G - 1).
# A sample whose values are all equal adds nothing and is left out.
.stacked_mean_terms <- function(model, robust, samples, cluster = NULL, used = NULL) {
  n_all <- length(samples[[1L]]$a)
  if (is.null(used)) used <- rep(TRUE, n_all)
  estimated <- !is.na(stats::coef(model))
  cross <- stats::setNames(numeric(sum(estimated)), names(stats::coef(model))[estimated])
  parts <- lapply(samples, function(s) {
    member <- s$member %in% TRUE
    dev <- numeric(n_all)
    if (any(member)) dev[member] <- s$a[member] - mean(s$a[member])
    varies <- any(member) &&
      max(abs(dev)) > 1e-12 * max(1, abs(s$a[member]))
    list(a = s$a, member = member, dev = dev, n = sum(member), varies = varies)
  })
  parts <- parts[vapply(parts, `[[`, logical(1), "varies")]
  if (!length(parts)) return(list(v = 0, cross = cross, df = Inf))

  if (robust == "CR2") {
    long <- do.call(rbind, lapply(seq_along(parts), function(j) {
      p <- parts[[j]]
      data.frame(a = p$a[p$member], set = j, row = which(p$member))
    }))
    long$set <- factor(long$set)
    m <- if (nlevels(long$set) == 1L) {
      stats::lm(a ~ 1, data = long)
    } else {
      stats::lm(a ~ 0 + set, data = long)
    }
    Vm <- clubSandwich::vcovCR(m, cluster = cluster[long$row], type = "CR2")
    lc <- as.data.frame(clubSandwich::linear_contrast(
      m, vcov = Vm, contrasts = matrix(1, nrow = 1, ncol = length(stats::coef(m))),
      test = "Satterthwaite"
    ))
    infl <- matrix(0, n_all, sum(estimated))
    infl[used, ] <- .glm_influence(model)$influence
    u <- Reduce(`+`, lapply(parts, function(p) p$dev / p$n))
    stacked <- used | Reduce(`|`, lapply(parts, `[[`, "member"))
    by_cluster <- cluster[stacked]
    G <- length(unique(by_cluster))
    cross[] <- G / (G - 1) * colSums(
      rowsum(infl[stacked, , drop = FALSE], by_cluster) * rowsum(u[stacked], by_cluster)[, 1]
    )
    return(list(v = lc$SE^2, cross = cross, df = lc$df))
  }

  u <- Reduce(`+`, lapply(parts, function(p) {
    p$dev / if (robust == "HC3") max(p$n - 1, 1) else sqrt(p$n * max(p$n - 1, 1))
  }))
  # A participant at leverage 1 adds its limit, 0, to the covariance, as in
  # .center_terms() (issue #173).
  if (robust == "HC3") cross[] <- colSums(.hc3_influence(model)$G * u[used])
  list(v = sum(u^2), cross = cross, df = NA_real_)
}

# What a fit with the indicator needs for its pretest effects, from the
# model `model` its covariance comes from (the complete-case refit for CR2),
# whose data `df` hold the centered `pre_obs`, 0 for the participants the
# indicator marks (`kept`), in the units the model's coefficients use.
# `in_model` flags the rows the model used, `cluster` gives the cluster of
# every row and `cluster_fit` those of the model's rows (CR2 only), `x` is
# the pretest as given, `centers` the two centers, and `scale` the number
# the pretest was divided by.
#
# A pretest effect is L'beta + m g, with beta the model's coefficients.
# - Moving the center by D moves only the participants with a pretest, so
#   it moves the effect by b p D, not b D: `pre_obs` at the new center is
#   `pre_obs` less D (`pretested` less the indicator), a combination of the
#   model's own columns, so the fit is the same, the coefficient of
#   `pretested` gains b D, and g loses b D.
# - m is a sample mean, as the center is, and the effect moves by g for
#   each unit m moves.
# So at the center used the effect gains the variance of b p (mean pretest)
# + g m and twice its covariance with the contrast (.stacked_mean_terms()).
# At the other center, which lies `shift` away in the model's units, the
# effect is the contrast L + m e_g + p shift e_b, and it moves by
# g - b shift for each unit m moves.
#
# Returns `m` and `p`; `shift`; `parts` and `differ`, as
# .center_terms_both() returns them; `center`, the terms at the center
# used; `other`, those at the other center (NULL when the centers are the
# same number, or when `note` says why they are missing, as
# .center_terms_both() does); and `note`.
.indicator_terms <- function(df, x, kept, in_model, model, robust, cluster, cluster_fit,
                             pretest_center, centers, scale) {
  cf <- stats::coef(model)
  b <- cf[["pre_obs"]]
  g <- cf[[.missing_pretest_column]]
  pretested <- df$pretested %in% 1L
  in_m <- in_model & pretested
  m <- sum(kept & in_model) / sum(in_m)
  p <- 1 - m
  in_observed <- .in_pretest_center(df, kept)
  in_analyzed <- in_observed & in_model
  parts <- .pretest_center_parts(x, in_observed, in_model, is.na(df$y))
  differ <- parts$n[["observed"]] > parts$n[["analyzed"]]
  other_center <- .other_center(pretest_center)
  shift <- (centers[[other_center]] - centers[[pretest_center]]) / scale
  terms <- function(center, g_at) {
    .stacked_mean_terms(
      model, robust,
      samples = list(
        list(a = b * p * df$pre_obs, member = if (center == "analyzed") in_analyzed else in_observed),
        list(a = g_at * df[[.missing_pretest_column]], member = in_m)
      ),
      cluster = if (robust == "CR2") cluster,
      used = in_model
    )
  }
  note <- NULL
  other <- NULL
  outside <- in_observed & !in_model
  if (pretest_center == "analyzed") {
    center <- terms("analyzed", g)
    if (differ) {
      if (robust == "CR2") {
        note <- .check_center_clusters(cluster, outside, stop = FALSE)
        if (!is.null(note)) note <- paste(note, "The standard errors at that center are not given.")
      }
      if (is.null(note)) other <- terms("observed", g - b * shift)
    }
  } else {
    if (robust == "CR2") .check_center_clusters(cluster, outside)
    center <- terms("observed", g)
    if (differ) other <- terms("analyzed", g - b * shift)
  }
  list(m = m, p = p, shift = shift, parts = parts, differ = differ, center = center,
       other = other, note = note)
}

# The variance of a contrast `L` of a fit with the indicator, with what the
# estimated center and the indicator's estimated mean add to a pretest
# effect (`terms`, from .indicator_terms(), which are for a weight of 1 on
# `pretested`, the weight of every pretest effect). `L` already holds the
# indicator's mean on the indicator's coefficient. A contrast without
# weight on `pretested` gains nothing. Returns what .center_adjusted()
# returns: `variance`, the part `added`, and `df`.
.indicator_adjusted <- function(L, V, terms, robust, df_contrast) {
  variance <- as.numeric(t(L) %*% V %*% L)
  w <- if ("pretested" %in% names(L)) L[["pretested"]] else 0
  if (is.null(terms) || !is.finite(w) || w == 0) {
    return(list(variance = variance, added = 0, df = df_contrast))
  }
  weights <- if (length(terms$cross) == length(L)) L else L[names(terms$cross)]
  v_contrast <- variance + 2 * w * sum(weights * terms$cross)
  v_means <- w^2 * terms$v
  df <- if (robust == "CR2") {
    .welch_df(c(v_contrast, v_means), c(df_contrast, terms$df))
  } else {
    df_contrast
  }
  list(variance = v_contrast + v_means, added = v_contrast + v_means - variance, df = df)
}

# The estimates of the pretest effects under the handling a fit did not
# use, for `missing_pretest_sensitivity`: a fit of the same model to the
# same data with the other handling, from which only the coefficients are
# taken.
# - `to = "indicator"` (the fit dropped the participants `kept`): they are
#   added at a centered pretest of 0 with the indicator, and each pretest
#   effect is its contrast plus the indicator's mean times its coefficient;
#   plus or minus infinity when the indicator separates the outcome
#   (.indicator_separation()).
# - `to = "drop"` (the fit kept them): they are left out, with the
#   indicator.
# `fml` is the formula of the other handling's model, `data` the data the
# inference came from, and `Ls` the contrasts of the pretest effects over
# the fit's coefficients. The centers are those of the fit: both handlings
# take them from observed pretests only. Returns `estimates`, one per
# contrast, NA where the other fit cannot be made or leaves a coefficient
# inestimable, and `direction`, the separation of the indicator in a fit
# that kept them (0 otherwise).
.other_handling_estimates <- function(fml, data, family, negbin, kept, to, Ls) {
  column <- .missing_pretest_column
  none <- list(estimates = rep(NA_real_, length(Ls)), direction = 0L)
  if (to == "indicator") {
    data$pre_obs[kept] <- 0
    data[[column]] <- as.integer(kept)
  } else {
    data$pre_obs[kept] <- NA
    data[[column]] <- NULL
  }
  fit <- tryCatch(.fit_model(fml, data, family, negbin, quiet = TRUE), error = function(e) NULL)
  if (is.null(fit)) return(none)
  cf <- stats::coef(fit)
  used <- rep(TRUE, nrow(data))
  if (!is.null(fit$na.action)) used[fit$na.action] <- FALSE
  m <- 0
  direction <- 0L
  if (to == "indicator") {
    if (!column %in% names(cf) || is.na(cf[[column]])) return(none)
    m <- sum(kept & used) / sum(used & data$pretested %in% 1L)
    direction <- .indicator_separation(stats::family(fit), fit$y[kept[used]])
  }
  estimates <- vapply(Ls, function(L) {
    terms <- names(L)[L != 0]
    if (!all(terms %in% names(cf)) || anyNA(cf[terms])) return(NA_real_)
    if (direction != 0L) return(direction * Inf)
    estimate <- sum(L[terms] * cf[terms])
    if (to == "indicator") estimate <- estimate + L[["pretested"]] * m * cf[[column]]
    estimate
  }, numeric(1), USE.NAMES = FALSE)
  list(estimates = estimates, direction = direction)
}

# Hedges' g of two samples with means `m1` and `m0`, standard deviations
# `s1` and `s0`, and sizes `n1` and `n0`: the difference of the means over
# the pooled standard deviation, times the small-sample correction
# 1 - 3 / (4 (n1 + n0 - 2) - 1) (What Works Clearinghouse, 2022, eqs. E.1 to
# E.4, pp. 162-163).
.hedges_g <- function(m1, m0, s1, s0, n1, n0) {
  pooled <- sqrt(((n1 - 1) * s1^2 + (n0 - 1) * s0^2) / (n1 + n0 - 2))
  (1 - 3 / (4 * (n1 + n0 - 2) - 1)) * (m1 - m0) / pooled
}

# The arm check (the decision record on issue #142, section 5, item 7d),
# for each treatment condition against the control, among the pretested
# participants analyzed: those in a "indicator" fit, with a pretest or
# without. `x` is their pretest (NA without one) and `condition` their
# condition; `levels` names the conditions, the control first.
# - The detector: the mean observed pretest by arm, their difference
#   (treatment minus control) as Hedges' g (.hedges_g()), and Welch's
#   (1947) t test. Regressing the covariate on the arm among complete cases
#   is the check Kayembe et al. (2020, sec. 5) give for loss that depends
#   on the arm and the covariate together.
# - A description, not a detector: the number and share without a pretest
#   by arm, with Fisher's exact test, which is the package's own choice.
#   Different rates alone do not bias the indicator method when the loss is
#   not selective (Kayembe et al., 2020, sec. 3.3).
# One row for each treatment. A test that cannot be computed is NA: Welch's
# needs two pretests in each arm, with some spread, and Fisher's needs both
# arms and both kinds of participant.
.missing_pretest_arm_check <- function(x, condition, levels) {
  control <- levels[1L]
  rows <- lapply(levels[-1L], function(treatment) {
    in_c <- condition %in% control
    in_t <- condition %in% treatment
    x_c <- x[in_c & !is.na(x)]
    x_t <- x[in_t & !is.na(x)]
    n <- c(length(x_c), length(x_t))
    means <- c(if (n[1]) mean(x_c) else NA_real_, if (n[2]) mean(x_t) else NA_real_)
    sds <- c(if (n[1] > 1L) stats::sd(x_c) else NA_real_,
             if (n[2] > 1L) stats::sd(x_t) else NA_real_)
    welch <- c(statistic = NA_real_, df = NA_real_, p.value = NA_real_)
    g <- NA_real_
    if (all(n > 1L) && all(is.finite(sds)) && sum(sds) > 0) {
      test <- stats::t.test(x_t, x_c, var.equal = FALSE)
      welch <- c(statistic = unname(test$statistic), df = unname(test$parameter),
                 p.value = test$p.value)
      g <- .hedges_g(means[2], means[1], sds[2], sds[1], n[2], n[1])
    }
    without <- c(sum(in_c & is.na(x)), sum(in_t & is.na(x)))
    total <- c(sum(in_c), sum(in_t))
    tab <- rbind(c(without[1], n[1]), c(without[2], n[2]))
    fisher <- if (all(rowSums(tab) > 0) && all(colSums(tab) > 0)) {
      stats::fisher.test(tab)$p.value
    } else {
      NA_real_
    }
    data.frame(
      comparison = paste(treatment, "vs", control),
      n_control = n[1], mean_control = means[1], sd_control = sds[1],
      n_treatment = n[2], mean_treatment = means[2], sd_treatment = sds[2],
      difference = means[2] - means[1], g = g,
      statistic = welch[["statistic"]], df = welch[["df"]], p.value = welch[["p.value"]],
      missing_control = without[1], share_control = if (total[1]) without[1] / total[1] else NA_real_,
      missing_treatment = without[2],
      share_treatment = if (total[2]) without[2] / total[2] else NA_real_,
      fisher_p = fisher,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# The posttests of the pretested participants analyzed, by condition: the
# number and the mean posttest of those without a pretest and of those with
# one. `y` is the posttest as the model fits it, `has_pretest` whether the
# participant has a pretest, and `condition` the condition of each; `levels`
# names the conditions, the control first.
.missing_pretest_posttests <- function(y, has_pretest, condition, levels) {
  rows <- lapply(levels, function(l) {
    member <- condition %in% l
    without <- member & !has_pretest
    with <- member & has_pretest
    data.frame(
      condition = l,
      n_without = sum(without), mean_without = if (any(without)) mean(y[without]) else NA_real_,
      n_with = sum(with), mean_with = if (any(with)) mean(y[with]) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# What a fit returns on the pretested participants without a pretest: the
# elements `missing_pretest` (the handling), `n_missing_pretest` (how many
# the model keeps, or would keep, by condition), `missing_pretest_share`
# (their share of the pretested participants analyzed with them kept),
# `missing_pretest_sensitivity`, and `missing_pretest_arm_check`. The last
# two are NULL when there is no such participant, and all but the first are
# NULL or NA without a pretest (`kept` NULL).
#
# `handling` is the value of `missing_pretest`. `kept` flags the
# participants concerned (.missing_pretest_rows()) and `analyzed` the
# pretested participants in the model with them kept. `condition` names
# each participant's condition and `levels` the conditions, the control
# first. `y` is the posttest and `x` the pretest as given. `labels` holds
# the columns that index the pretest effects (`contrast`, after
# `comparison` for several treatments), `kept_estimates` and
# `dropped_estimates` their estimates under each handling, and `direction`
# the separation of the indicator (.indicator_separation()).
.missing_pretest_elements <- function(handling, kept = NULL, analyzed = NULL, condition = NULL,
                                      levels = NULL, y = NULL, x = NULL, labels = NULL,
                                      kept_estimates = NULL, dropped_estimates = NULL,
                                      direction = 0L) {
  if (is.null(kept)) {
    return(list(missing_pretest = handling, n_missing_pretest = NULL,
                missing_pretest_share = NA_real_, missing_pretest_sensitivity = NULL,
                missing_pretest_arm_check = NULL))
  }
  n <- vapply(levels, function(l) sum(kept & condition %in% l), integer(1))
  out <- list(missing_pretest = handling, n_missing_pretest = n,
              missing_pretest_share = if (any(analyzed)) sum(kept) / sum(analyzed) else NA_real_,
              missing_pretest_sensitivity = NULL, missing_pretest_arm_check = NULL)
  if (!any(kept)) return(out)
  effects <- labels
  effects$kept <- kept_estimates
  effects$dropped <- dropped_estimates
  effects$difference <- kept_estimates - dropped_estimates
  rownames(effects) <- NULL
  out$missing_pretest_sensitivity <- list(
    effects = effects,
    posttests = .missing_pretest_posttests(.perm_outcome(y)[analyzed], !kept[analyzed],
                                           condition[analyzed], levels),
    separation = direction,
    one = sum(kept) == 1L
  )
  out$missing_pretest_arm_check <- .missing_pretest_arm_check(
    ifelse(kept, NA_real_, x)[analyzed], condition[analyzed], levels
  )
  out
}
