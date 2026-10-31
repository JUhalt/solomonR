# The diagnostics of check_solomon_missing() for the pretest effects when
# posttests are missing (issue #135). They follow the decision record's
# section 5, items 1 to 8:
#   1. posttests missing by group, overall and differential attrition
#      between the pretest conditions and between the conditions (What Works
#      Clearinghouse, 2022, p. 40);
#   2. the pretest scores of those lost, for each pretested group (Campbell,
#      1957, pp. 301-302; Campbell & Stanley, 1963/1966, pp. 15-16);
#   3. each pretested condition's selection shift, with its standard error
#      and a Wald test of equal shifts, and the raw comparison of posttest
#      means when that test gives p < .05 (issue #143);
#   4. the gap between the two centers, b m1 (completers - dropouts), the
#      adjustment of What Works Clearinghouse (2022, eq. H.2, p. 207);
#   5. for each unpretested group, its share lost and the offset the default
#      assumes (the delta line);
#   6. covariates recorded in all groups: loss by covariate within each
#      pretest condition, and their distribution among the participants
#      analyzed;
#   7. pretested participants with a posttest and no pretest (issue #142);
#   8. what cannot be checked: the unpretested groups' losses against
#      baseline (Campbell, 1957, p. 305).
# Items 4 and 5 come from an identity-link fit_solomon_glm() of the data, and
# items 3 and 7 from one without the pretest. Each fit is wrapped so that
# check_solomon_missing() never fails because the data do not support a
# model. validate_solomon(), which uses only its classification of the
# missing values, skips these diagnostics and their fits.

# A fit for the diagnostics, or NULL when the data do not support one.
.diagnostic_fit <- function(...) {
  tryCatch(suppressWarnings(fit_solomon_glm(...)), error = function(e) NULL)
}

# Item 1: the attrition of each group and the rates by pretest condition and
# by condition. `condition` is each participant's condition (NA when not
# recorded), `levels` the conditions with the control first. Fisher's exact
# test of loss by pretest condition gives the rule for the statement that the
# loss rates differ (p < .05): a test of rates only, which a choice of the
# package, not of the record, fixes (issue #135).
.missing_attrition <- function(y_post, condition, levels, pretested) {
  lost <- is.na(y_post)
  known <- !is.na(condition) & !is.na(pretested)
  groups <- .pretest_attrition(condition, levels, pretested, y_post, rep(TRUE, length(lost)), 0,
                               c(analyzed = 0, observed = 0), identity = FALSE)$attrition
  share <- function(rows) if (any(rows)) mean(lost[rows]) else NA_real_
  by_pretest <- c(pretested = share(known & pretested %in% 1L),
                  unpretested = share(known & pretested %in% 0L))
  by_condition <- list2DF(list(
    condition = levels,
    n = vapply(levels, function(l) sum(known & condition %in% l), integer(1), USE.NAMES = FALSE),
    posttest_missing = vapply(levels, function(l) sum(known & condition %in% l & lost), integer(1),
                              USE.NAMES = FALSE)
  ))
  by_condition$share <- ifelse(by_condition$n > 0L, by_condition$posttest_missing / by_condition$n,
                               NA_real_)
  tab <- table(factor(pretested[known], levels = c(1L, 0L)),
               factor(lost[known], levels = c(FALSE, TRUE)))
  p <- if (all(dim(tab) == 2L) && all(rowSums(tab) > 0) && all(colSums(tab) > 0)) {
    stats::fisher.test(tab)$p.value
  } else {
    NA_real_
  }
  list(
    groups = groups,
    overall = share(known),
    by_pretest = by_pretest,
    pretest_difference = by_pretest[["pretested"]] - by_pretest[["unpretested"]],
    by_condition = by_condition,
    condition_difference = if (any(is.finite(by_condition$share))) {
      diff(range(by_condition$share, na.rm = TRUE))
    } else {
      NA_real_
    },
    p.value = p,
    differ = isTRUE(p < .05)
  )
}

# Items 2 and 3: in each pretested condition, the mean pretest of the
# participants with and without a posttest, and the selection shift: the
# mean pretest of those with a posttest minus that of every participant with
# a pretest. Its standard error comes from its influence function,
# R_i (x_i - mean_with) N / n - (x_i - mean_all) over the N participants
# with a pretest, n of them with a posttest (R_i = 1), as for any smooth
# function of means (Stefanski & Boos, 2002): the sum of its squares over
# N (N - 1). A condition that lost no posttest has a shift of exactly 0,
# with no sampling error.
#
# The shifts "differ" when a Wald test of equal shifts gives p < .05: the
# differences from the first condition's shift, with their covariance (the
# conditions hold different participants), in a chi-square statistic whose
# degrees of freedom are the rank of that covariance, so that conditions
# without error enter as the known value 0. The rule is the package's (the
# decision record on issue #135 leaves it open); its rate when every
# condition loses by the same rule is checked in
# tools/pretest-center-check.R (part 8). Being a ratio of a shift to its
# standard error, it does not depend on the units of the pretest.
.missing_pretest_scores <- function(y_post, y_pre, condition, levels, pretested) {
  lost <- is.na(y_post)
  rows <- lapply(levels, function(l) {
    member <- pretested %in% 1L & condition %in% l & !is.na(y_pre)
    with_post <- member & !lost
    without <- member & lost
    m <- function(v) if (any(v)) mean(y_pre[v]) else NA_real_
    N <- sum(member)
    n <- sum(with_post)
    se <- NA_real_
    if (n > 0L && N > 1L) {
      x <- y_pre[member]
      r <- !lost[member]
      influence <- r * (x - mean(x[r])) * N / n - (x - mean(x))
      se <- sqrt(sum(influence^2) / (N * (N - 1)))
    }
    list(condition = l, n = N, mean_all = m(member),
         n_with = n, mean_with = m(with_post),
         n_without = sum(without), mean_without = m(without), shift_se = se)
  })
  scores <- list2DF(lapply(stats::setNames(names(rows[[1]]), names(rows[[1]])), function(nm) {
    unlist(lapply(rows, `[[`, nm), use.names = FALSE)
  }))
  scores$difference <- scores$mean_without - scores$mean_with
  scores$shift <- ifelse(scores$n_with > 0L, scores$mean_with - scores$mean_all, NA_real_)
  scores <- scores[c(setdiff(names(scores), "shift_se"), "shift_se")]
  test <- .shift_test(scores$shift, scores$shift_se)
  list(scores = scores, test = test, differ = isTRUE(test[["p.value"]] < .05))
}

# The Wald test of equal selection shifts `shift` with standard errors `se`
# (independent estimates): `statistic`, `df`, and `p.value`, NA when fewer
# than two conditions have a shift or none has sampling error.
.shift_test <- function(shift, se) {
  keep <- is.finite(shift) & is.finite(se)
  none <- c(statistic = NA_real_, df = NA_real_, p.value = NA_real_)
  if (sum(keep) < 2L) return(none)
  shift <- shift[keep]
  v <- se[keep]^2
  d <- shift[-1] - shift[1]
  V <- diag(v[-1], nrow = length(d)) + v[1]
  e <- eigen(V, symmetric = TRUE)
  rank <- e$values > 1e-10 * max(e$values)
  if (!any(rank)) return(none)
  z <- crossprod(e$vectors[, rank, drop = FALSE], d)
  statistic <- sum(z^2 / e$values[rank])
  df <- sum(rank)
  c(statistic = statistic, df = df,
    p.value = stats::pchisq(statistic, df, lower.tail = FALSE))
}

# Fisher's test of loss by level, for a table of levels by posttest
# observed and missing: the exact test for two levels, and for more a
# p-value simulated from 10,000 tables with the same margins, under a fixed
# seed that leaves the user's random-number stream as it was, because the
# exact test of a larger table can exhaust its workspace at ordinary sample
# sizes. The choice is the package's; no source fixes one. A failure gives
# NA. Returns the p-value with the test used as its "test" attribute.
.loss_by_level_test <- function(tab) {
  exact <- nrow(tab) == 2L
  p <- tryCatch(
    if (exact) {
      stats::fisher.test(tab)$p.value
    } else {
      withr::with_seed(135, stats::fisher.test(tab, simulate.p.value = TRUE, B = 10000)$p.value)
    },
    error = function(e) NA_real_
  )
  structure(p, test = if (exact) "exact" else "simulated")
}

# Item 6: for each covariate, loss by covariate within each pretest
# condition, and its distribution among the participants with a posttest.
# A covariate with more than two distinct numbers is summarized by its mean
# among those with and those without a posttest; any other is taken by its
# levels, with Fisher's test of loss by level in each pretest condition
# (.loss_by_level_test()).
.missing_covariates <- function(covariates, y_post, condition, pretested) {
  if (is.null(covariates) || !ncol(covariates)) return(NULL)
  lost <- is.na(y_post)
  known <- !is.na(condition) & !is.na(pretested)
  out <- lapply(names(covariates), function(name) {
    v <- covariates[[name]]
    numeric <- is.numeric(v) && length(unique(v[!is.na(v)])) > 2L
    if (numeric) {
      return(list2DF(lapply(list(
        pretested = c(1L, 0L),
        n_with = vapply(c(1L, 0L), function(p) sum(known & pretested == p & !lost & !is.na(v)), 1L),
        mean_with = vapply(c(1L, 0L), function(p) {
          r <- known & pretested == p & !lost & !is.na(v)
          if (any(r)) mean(v[r]) else NA_real_
        }, 1),
        n_without = vapply(c(1L, 0L), function(p) sum(known & pretested == p & lost & !is.na(v)), 1L),
        mean_without = vapply(c(1L, 0L), function(p) {
          r <- known & pretested == p & lost & !is.na(v)
          if (any(r)) mean(v[r]) else NA_real_
        }, 1)
      ), identity)))
    }
    v <- as.character(v)
    levels <- sort(unique(v[!is.na(v)]))
    rows <- expand.grid(level = levels, pretested = c(1L, 0L), stringsAsFactors = FALSE)
    rows <- rows[, c("pretested", "level")]
    rows$n <- mapply(function(p, l) sum(known & pretested == p & v %in% l), rows$pretested, rows$level)
    rows$posttest_missing <- mapply(function(p, l) sum(known & pretested == p & v %in% l & lost),
                                    rows$pretested, rows$level)
    rows$share <- ifelse(rows$n > 0L, rows$posttest_missing / rows$n, NA_real_)
    rows$analyzed_share <- mapply(function(p, l) {
      analyzed <- known & pretested == p & !lost & !is.na(v)
      if (any(analyzed)) mean(v[analyzed] == l) else NA_real_
    }, rows$pretested, rows$level)
    tests <- lapply(c(1L, 0L), function(p) {
      r <- known & pretested == p & !is.na(v)
      tab <- table(factor(v[r], levels = levels), factor(lost[r], levels = c(FALSE, TRUE)))
      tab <- tab[rowSums(tab) > 0, , drop = FALSE]
      if (nrow(tab) >= 2L && all(colSums(tab) > 0)) {
        .loss_by_level_test(tab)
      } else {
        structure(NA_real_, test = NA_character_)
      }
    })
    at <- match(rows$pretested, c(1L, 0L))
    rows$p.value <- vapply(tests, as.numeric, 1)[at]
    rows$test <- vapply(tests, attr, "", "test")[at]
    rownames(rows) <- NULL
    rows
  })
  stats::setNames(out, names(covariates))
}

# Item 7: pretested participants with a posttest and no pretest, whom every
# fit with a pretest covariate drops (issue #142): their number and mean
# posttest in each pretested condition, beside those of the pretested
# participants with both scores.
.missing_pretests <- function(y_post, y_pre, condition, levels, pretested) {
  has_post <- !is.na(y_post)
  rows <- lapply(levels, function(l) {
    member <- pretested %in% 1L & condition %in% l & has_post
    none <- member & is.na(y_pre)
    both <- member & !is.na(y_pre)
    list(condition = l, n_without_pretest = sum(none),
         mean_without_pretest = if (any(none)) mean(y_post[none]) else NA_real_,
         n_with_pretest = sum(both),
         mean_with_pretest = if (any(both)) mean(y_post[both]) else NA_real_)
  })
  list2DF(lapply(stats::setNames(names(rows[[1]]), names(rows[[1]])), function(nm) {
    unlist(lapply(rows, `[[`, nm), use.names = FALSE)
  }))
}

# All the diagnostics, as elements of check_solomon_missing()'s result.
# `treat` and `control` are as the user gave them, for the internal fits;
# `condition` names each participant's condition and `levels` the
# conditions, the control first.
.missing_diagnostics <- function(y_post, treat, control, condition, levels, pretested, y_pre,
                                 covariates) {
  attrition <- .missing_attrition(y_post, condition, levels, pretested)
  out <- list(
    attrition = attrition,
    pretest_scores = NULL,
    selection_test = NULL,
    selection_differ = NULL,
    gap = NULL,
    unpretested = NULL,
    covariates = .missing_covariates(covariates, y_post, condition, pretested),
    missing_pretest = NULL,
    raw_pretest_effects = NULL,
    unverifiable = attrition$groups$posttest_missing[attrition$groups$pretested == 0L] > 0L
  )
  out$unverifiable <- any(out$unverifiable)
  # The checks on the pretest need numeric scores; the classification of
  # the missing values above does not.
  if (is.null(y_pre) || !is.numeric(y_pre)) return(out)

  scores <- .missing_pretest_scores(y_post, y_pre, condition, levels, pretested)
  out$pretest_scores <- scores$scores
  out$selection_test <- scores$test
  out$selection_differ <- scores$differ
  out$missing_pretest <- .missing_pretests(y_post, y_pre, condition, levels, pretested)

  # Items 4 and 5, on the identity link, with the covariates.
  fit <- .diagnostic_fit(y_post, treat, pretested, y_pre, covariates = covariates,
                         robust = "none", control = control)
  s <- if (!is.null(fit)) fit$pretest_sensitivity
  if (!is.null(s)) {
    out$gap <- c(slope = s$slope, m1 = s$m1, analyzed = s$pretest_means[["analyzed"]],
                 not_analyzed = s$pretest_means[["not_analyzed"]],
                 gap = s$difference$difference[1], predicted = s$predicted)
    unpretested <- attrition$groups[attrition$groups$pretested == 0L,
                                    c("condition", "n", "posttest_missing"), drop = FALSE]
    # The fit's conditions come in the same order, the control first, but a
    # four-group fit names them by `control` when it was given, so they are
    # matched by position.
    position <- if (nrow(s$delta) == nrow(unpretested)) seq_len(nrow(unpretested)) else
      match(unpretested$condition, s$delta$condition)
    unpretested$m0 <- s$delta$m0[position]
    unpretested$delta <- s$delta$delta[position]
    rownames(unpretested) <- NULL
    out$unpretested <- unpretested
  }

  # Items 3 and 7: the raw comparison of posttest means in each condition,
  # from a fit without the pretest, when the selection differs between the
  # pretested conditions or pretested participants lack a pretest.
  if (isTRUE(scores$differ) || sum(out$missing_pretest$n_without_pretest) > 0L) {
    raw <- .diagnostic_fit(y_post, treat, pretested, control = control)
    if (!is.null(raw)) {
      e <- raw$effects
      keys <- intersect("comparison", names(e))
      out$raw_pretest_effects <- e[e$contrast %in% .solomon_pretest_order,
                                   c(keys, .solomon_effect_columns), drop = FALSE]
      rownames(out$raw_pretest_effects) <- NULL
    }
  }
  out
}

# The lines print.solomon_missing() adds for the diagnostics, indented by
# two spaces. Nothing when no posttest is missing and no pretested
# participant lacks a pretest.
.missing_diagnostic_lines <- function(x, digits = 2) {
  a <- x$attrition
  if (is.null(a)) return(character(0))
  lost <- sum(a$groups$posttest_missing)
  no_pretest <- if (is.null(x$missing_pretest)) 0L else sum(x$missing_pretest$n_without_pretest)
  if (!lost && !no_pretest) return(character(0))
  num <- function(v) formatC(v, format = "f", digits = digits)
  # A factor of a printed product or difference: in parentheses when
  # negative, as in "0.43 - (-0.61)".
  factor_num <- function(v) ifelse(v < 0, paste0("(", num(v), ")"), num(v))
  wrap <- function(text, indent = 2, exdent = indent) {
    strwrap(text, width = 80, indent = indent, exdent = exdent)
  }
  lines <- "Missing posttests and the pretest effects (see ?fit_solomon_glm)"
  if (lost) {
    by_condition <- a$by_condition
    lines <- c(lines, wrap(sprintf(paste(
      "Posttests missing: %s overall; pretested %s, unpretested %s (a difference of %s",
      "points%s); by condition, %s (the largest difference %s points)."
    ), .pct(a$overall), .pct(a$by_pretest[["pretested"]]), .pct(a$by_pretest[["unpretested"]]),
    formatC(100 * a$pretest_difference, format = "f", digits = 1),
    if (is.finite(a$p.value)) paste0("; Fisher's exact test, ", .apa_p(a$p.value, md = FALSE)) else "",
    .series_and(sprintf("%s %s", by_condition$condition, .pct(by_condition$share))),
    formatC(100 * a$condition_difference, format = "f", digits = 1))))
    lines <- c(lines, wrap(if (a$differ) {
      paste("Loss rates differ between the pretest conditions, which suggests that taking the",
            "pretest changed attrition. Neither center is assured; see the delta line.")
    } else {
      "Equal loss rates do not show that the same kinds of participants were lost."
    }))
  }
  sc <- x$pretest_scores
  if (!is.null(sc) && sum(a$groups$posttest_missing[a$groups$pretested == 1L]) > 0L) {
    shown <- sc[sc$n_without > 0L, , drop = FALSE]
    if (nrow(shown)) {
      lines <- c(lines, wrap(sprintf(
        "Mean pretest of the pretested participants without and with a posttest: %s.",
        paste(sprintf("%s %s (%d) and %s (%d)", shown$condition, num(shown$mean_without),
                      shown$n_without, num(shown$mean_with), shown$n_with), collapse = "; ")
      )))
    }
    test <- x$selection_test
    test_text <- if (!is.null(test) && is.finite(test[["p.value"]])) {
      sprintf("; Wald test of equal shifts, chi-square(%d) = %s, %s", as.integer(test[["df"]]),
              formatC(test[["statistic"]], format = "f", digits = 2),
              .apa_p(test[["p.value"]], md = FALSE))
    } else {
      ""
    }
    se_text <- ifelse(is.finite(sc$shift_se), paste0(" (SE ", num(sc$shift_se), ")"), "")
    lines <- c(lines, wrap(sprintf(paste(
      "Selection shift in each pretested condition, the mean pretest with a posttest minus",
      "every pretest's: %s%s."
    ), .series_and(sprintf("%s %s%s", sc$condition, num(sc$shift), se_text)), test_text)))
    if (isTRUE(x$selection_differ)) {
      lines <- c(lines, wrap(paste(
        "The shifts differ between the pretested conditions (p < .05, the package's rule).",
        "If the conditions lost different kinds of participants, neither center gives both",
        "conditional pretest effects, and Pretest x Treatment is biased at every center",
        "(issue #143). The raw comparison of posttest means in each condition, from a fit",
        "without the pretest, is given below."
      )))
    }
  }
  g <- x$gap
  if (!is.null(g) && lost && is.finite(g[["not_analyzed"]])) {
    lines <- c(lines, wrap(sprintf(paste(
      "Gap between the centers: b m1 (completers - dropouts) = %s x %s x (%s - %s) = %s, the",
      "pretest effects at the analyzed center minus those at the observed center (identity",
      "link%s)."
    ), num(g[["slope"]]), num(g[["m1"]]), num(g[["analyzed"]]), factor_num(g[["not_analyzed"]]),
    formatC(g[["gap"]], format = "f", digits = digits + 1L),
    if (is.null(x$covariates)) "" else ", with the covariates")))
  }
  u <- x$unpretested
  if (!is.null(u) && lost && !is.null(g) && is.finite(g[["not_analyzed"]])) {
    lost_some <- is.finite(u$m0) & u$m0 > 0
    if (any(lost_some)) {
      lines <- c(lines, wrap(sprintf(paste(
        "Unpretested groups, the share without a posttest and the offset the default assumes",
        "for their dropouts: %s; the pretested dropouts' pretests predict %s for them."
      ), paste(sprintf("%s %s and %s", u$condition[lost_some], .pct(u$m0[lost_some]),
                       num(u$delta[lost_some])), collapse = "; "),
      num(g[["predicted"]]))))
    }
    if (any(!lost_some)) {
      several <- sum(!lost_some) > 1L
      lines <- c(lines, wrap(sprintf(paste(
        "The unpretested %s %s lost no posttests: if the data hold everyone randomized to %s,",
        "the default's pretest %s there %s right only if the pretested groups' losses,",
        "taken together, were unrelated to baseline, because the default's center pools",
        "the pretested groups."
      ), .series_and(u$condition[!lost_some]), if (several) "groups" else "group",
      if (several) "them" else "it", if (several) "effects" else "effect",
      if (several) "are" else "is")))
    }
  }
  for (name in names(x$covariates)) {
    cv <- x$covariates[[name]]
    if ("level" %in% names(cv)) {
      condition_text <- vapply(c(1L, 0L), function(p) {
        r <- cv[cv$pretested == p, , drop = FALSE]
        test <- if (!"test" %in% names(r) || is.na(r$p.value[1])) {
          ""
        } else if (identical(r$test[1], "simulated")) {
          paste0(" (Fisher's test, p-value simulated from 10,000 tables, ",
                 .apa_p(r$p.value[1], md = FALSE), ")")
        } else {
          paste0(" (Fisher's exact test, ", .apa_p(r$p.value[1], md = FALSE), ")")
        }
        sprintf("%s %s%s", if (p == 1L) "pretested" else "unpretested",
                .series_and(sprintf("%s %s", .pct(r$share, 0L), r$level)), test)
      }, "")
      analyzed_text <- vapply(c(1L, 0L), function(p) {
        r <- cv[cv$pretested == p, , drop = FALSE]
        sprintf("%s %s", if (p == 1L) "pretested" else "unpretested",
                .series_and(sprintf("%s %s", .pct(r$analyzed_share, 0L), r$level)))
      }, "")
      lines <- c(lines, wrap(sprintf(
        "Recorded in all groups, %s: posttests missing by level, %s. Among those with a posttest: %s.",
        name, paste(condition_text, collapse = "; "), paste(analyzed_text, collapse = "; ")
      )))
    } else {
      lines <- c(lines, wrap(sprintf(
        "Recorded in all groups, %s: mean without and with a posttest, pretested %s and %s; unpretested %s and %s.",
        name, num(cv$mean_without[1]), num(cv$mean_with[1]), num(cv$mean_without[2]),
        num(cv$mean_with[2])
      )))
    }
  }
  mp <- x$missing_pretest
  if (!is.null(mp) && no_pretest) {
    lines <- c(lines, wrap(sprintf(paste(
      "Pretested participants with a posttest and no pretest: %d. Their mean posttest was %s,",
      "against %s for the pretested participants with both scores. Every fit with the pretest",
      "drops them, at either center (issue #142); the raw comparison below keeps them."
    ), no_pretest,
    num(sum(mp$mean_without_pretest * mp$n_without_pretest, na.rm = TRUE) / no_pretest),
    num(sum(mp$mean_with_pretest * mp$n_with_pretest, na.rm = TRUE) / sum(mp$n_with_pretest)))))
  }
  raw <- x$raw_pretest_effects
  if (!is.null(raw)) {
    labels <- .center_labels(raw)
    column <- function(v) formatC(num(v), width = max(nchar(num(c(raw$estimate, raw$conf.low,
                                                                    raw$conf.high)))))
    lines <- c(lines, "  Raw comparison (a fit without the pretest, HC3; 95% CI):",
               paste0("    ", formatC(labels, width = -max(nchar(labels))), "  ",
                      column(raw$estimate), "  [", column(raw$conf.low), ", ",
                      column(raw$conf.high), "]"))
  }
  if (x$unverifiable) {
    lines <- c(lines, wrap(paste0(
      "The unpretested groups' losses cannot be checked against baseline (Campbell, 1957, ",
      "p. 305)", if (sum(a$groups$posttest_missing[a$groups$pretested == 1L]) == 0L) {
        paste("; with only those groups losing posttests, the two centers coincide, but neither",
              "is protected. Only the delta line reaches that case: if an unpretested group's",
              "dropouts would have scored Delta above its completers, the pretest effect in its",
              "condition is m0 x Delta lower, m0 being that group's share without a posttest",
              "(see ?fit_solomon_glm).")
      } else {
        "."
      }
    )))
  }
  lines
}
