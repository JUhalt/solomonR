# The output of the two pretest centers (issue #135): the block that the
# print methods add when the two centers differ, the attrition line they add
# whenever a posttest is missing, and the sentences report_solomon() adds.
# The decision record on issue #135 sets out what is printed (its sections 4
# and 5). Fits made by versions before #135 have no `pretest_center` and
# print and report as before, and so does a fit with every posttest
# observed.

# A percentage with `digits` decimals.
.pct <- function(x, digits = 1L) paste0(formatC(100 * x, format = "f", digits = digits), "%")

# A number with `digits` decimals, as the printed tables show it.
.center_num <- function(x, digits) formatC(x, format = "f", digits = digits)

# The differences between the pretest effects at the two centers, with
# those that are zero but for rounding error set to zero: a difference that
# is zero in exact arithmetic, as for a rate ratio, which does not depend on
# the center, would otherwise print with the sign of its rounding error
# (issue #137). A small real difference keeps its sign. `d` is
# `pretest_sensitivity$difference`. When the effects are infinite, as under
# separation by the indicator of a missing pretest (issue #142), the
# difference is the finite value the fit gives it and is not compared with
# them.
.center_differences <- function(d) {
  size <- pmax(abs(d$analyzed), abs(d$observed))
  size[is.infinite(size)] <- 0
  zero <- abs(d$difference) <= 1e-12 * size
  ifelse(is.finite(d$difference) & zero, 0, d$difference)
}

# The attrition of a fit's Solomon groups: a data frame with `group`,
# `condition`, `pretested`, `n`, `posttest_missing`, and `share`, one row
# for each group. It is `pretest_sensitivity$attrition` when the fit has a
# pretest covariate, and is otherwise counted from the data the model was
# fitted to. NULL for fits made before #135 and for fits whose data do not
# hold the posttests (fit_solomon_mmrm(), which .mmrm_attrition() covers).
# The groups are those of the conditions in a fit with several treatments
# (class `solomon_ngroup`), and of `treat` otherwise: a four-group fit made
# before #174 can have a covariate named condition (issue #174).
.fit_attrition <- function(x) {
  if (is.null(x$pretest_center)) return(NULL)
  # A summary carries the table it was given.
  if (!is.null(x$attrition)) return(x$attrition)
  if (!is.null(x$pretest_sensitivity$attrition)) return(x$pretest_sensitivity$attrition)
  d <- x$data
  if (!is.data.frame(d) || !all(c("y", "pretested") %in% names(d))) return(NULL)
  condition <- as.character(if (inherits(x, "solomon_ngroup")) d$condition else d$treat)
  levels <- sort(unique(condition[!is.na(condition)]))
  .pretest_attrition(condition, levels, d$pretested, d$y, rep(TRUE, nrow(d)), 0,
                     c(analyzed = 0, observed = 0), identity = FALSE)$attrition
}

# The attrition of a fit_solomon_mmrm() fit at its last occasion, in the
# columns of .fit_attrition(): each group's participants (`x$data`) and its
# observed posttests at that occasion (`x$observed`). NULL for fits made
# before #135.
.mmrm_attrition <- function(x) {
  if (is.null(x$pretest_center) || is.null(x$observed) || is.null(x$data)) return(NULL)
  groups <- .solomon_group_labels
  members <- factor(.solomon_group_labels[.solomon_group(x$data$treat, x$data$pretested)],
                    levels = groups)
  n <- as.integer(table(members))
  observed <- as.integer(x$observed[groups, ncol(x$observed)])
  observed[is.na(observed)] <- 0L
  list2DF(list(
    group = groups,
    condition = c("treatment", "control", "treatment", "control"),
    pretested = c(1L, 1L, 0L, 0L),
    n = n,
    posttest_missing = n - observed,
    share = ifelse(n > 0L, (n - observed) / n, NA_real_)
  ))
}

# The shares of missing posttests among the pretested and the unpretested
# participants, their range over the groups, and the numbers missing, from
# the attrition table `a`. NULL when `a` is NULL or empty.
.attrition_shares <- function(a) {
  if (is.null(a)) return(NULL)
  a <- a[a$n > 0L, , drop = FALSE]
  if (!nrow(a)) return(NULL)
  pre <- a$pretested == 1L
  share <- function(rows) {
    if (any(rows)) sum(a$posttest_missing[rows]) / sum(a$n[rows]) else NA_real_
  }
  list(
    pretested = share(pre),
    unpretested = share(!pre),
    range = range(a$posttest_missing / a$n),
    missing = c(pretested = sum(a$posttest_missing[pre]),
                unpretested = sum(a$posttest_missing[!pre]))
  )
}

# The attrition line, as in "Posttests missing: pretested 37.3%, unpretested
# 36.6% (29% to 46% by group).", for an attrition table `a`; `what` names
# what is missing. character(0) when no posttest is missing.
.attrition_line <- function(a, what = "Posttests missing") {
  s <- .attrition_shares(a)
  if (is.null(s) || !sum(s$missing)) return(character(0))
  sprintf("%s: pretested %s, unpretested %s (%s).", what, .pct(s$pretested),
          .pct(s$unpretested), .attrition_spread(s))
}

# The range of the share missing over the groups, as in "29% to 46% by
# group" or "20% in every group", from .attrition_shares() (`s`).
.attrition_spread <- function(s) {
  by_group <- .pct(s$range, 0L)
  if (by_group[1] == by_group[2]) {
    paste(by_group[1], "in every group")
  } else {
    paste(by_group[1], "to", by_group[2], "by group")
  }
}

# Whether only the unpretested groups lost posttests: the two centers then
# coincide, but neither is protected, because those losses cannot be checked
# against the pretest (Campbell, 1957, p. 305; the decision record on issue
# #135, section 5, item 8).
.unpretested_only <- function(a) {
  s <- .attrition_shares(a)
  !is.null(s) && s$missing[["pretested"]] == 0L && s$missing[["unpretested"]] > 0L
}

.unpretested_only_line <- paste(
  "Only the unpretested groups lost posttests, so the two pretest centers coincide;",
  "those losses cannot be checked against the pretest, and neither center is",
  "protected against them: the pretest effects assume that they were unrelated to",
  "baseline and to the posttest. See ?check_solomon_missing."
)

# The pretested participants outside the model, as in "44 without a
# posttest" or "5 not analyzed, 3 of them without a posttest", from
# `pretest_n` (`n`). `posttest` words the lack of a posttest; `of_them`
# words the phrase for the report, which follows a count of participants,
# as in "44 of them without a posttest".
.center_outside_phrase <- function(n, posttest = "without a posttest", of_them = FALSE) {
  outside <- n[["observed"]] - n[["analyzed"]]
  none <- n[["no_posttest"]]
  lead <- if (of_them) "%d of them" else "%d"
  if (none == outside) {
    sprintf(paste(lead, "%s"), outside, posttest)
  } else if (none == 0L) {
    sprintf(paste(lead, "not analyzed"), outside)
  } else if (of_them) {
    sprintf("%d of them not analyzed, %d %s", outside, none, posttest)
  } else {
    sprintf("%d not analyzed, %d of them %s", outside, none, posttest)
  }
}

# The label of each pretest effect in the printed table: the contrast, with
# the condition for designs with several treatments ("Pretest effect | RP")
# and the occasion for repeated posttests.
.center_labels <- function(d) {
  label <- d$contrast
  if ("comparison" %in% names(d)) {
    by_condition <- d$contrast %in% .solomon_pretest_order[1:2]
    label[by_condition] <- paste("Pretest effect |", d$comparison[by_condition])
  }
  if ("occasion" %in% names(d)) label <- sprintf("%s, occasion %s", label, d$occasion)
  label
}

# The table of the pretest effects at both centers and their difference,
# from `pretest_sensitivity$difference` (`d`), as lines indented by two
# spaces (`indent`).
.center_table <- function(d, digits, indent = "  ") {
  cols <- list(.center_labels(d), .center_num(d$analyzed, digits),
               .center_num(d$observed, digits), .center_num(.center_differences(d), digits))
  heads <- c("", "analyzed", "observed", "difference")
  widths <- mapply(function(h, v) max(nchar(h), nchar(v)), heads, cols)
  row <- function(values) {
    paste0(indent, paste(c(sprintf("%-*s", widths[1], values[1]),
                           sprintf("%*s", widths[-1], values[-1])), collapse = "   "))
  }
  c(row(heads), vapply(seq_along(cols[[1]]), function(i) row(vapply(cols, `[`, "", i)), ""))
}

# The scale of a fit's pretest effects for the block, as in "log-odds".
.center_scale <- function(family) {
  link <- if (is.null(family)) "identity" else family$link
  switch(link, logit = "log-odds", log = "log", identity = "outcome", paste(link, "link"))
}

# The sentence on the delta line (the decision record on issue #135,
# section 3, point 4): on the identity link, the offset each unpretested
# group's dropouts would need for the "observed" center to give the
# default's pretest effect, and what the pretested dropouts' pretests
# predict for them. Groups that lost no posttests are named in a sentence of
# their own.
.center_delta_text <- function(s, digits) {
  d <- s$delta
  lost <- is.finite(d$m0) & d$m0 > 0
  text <- character(0)
  if (any(lost)) {
    values <- d$delta[lost]
    below <- all(values <= 0) && any(values < 0)
    listed <- .series_and(sprintf("%s (%s)", .center_num(abs(values), digits), d$condition[lost]))
    predicted <- .center_num(abs(s$predicted), digits)
    text <- sprintf(paste(
      "The default assumes that each unpretested group's dropouts would have scored",
      "%s %s its completers, as the pretested dropouts' pretests predict for them (%s);",
      "\"observed\" assumes 0."
    ), listed, if (below) "below" else "above", predicted)
  }
  c(text, .no_loss_text(lost, d$condition))
}

# The sentence on the unpretested groups that lost no posttests, the one
# case that points to the "observed" center (the decision record on issue
# #135, section 3, point 2): `lost` flags, for each condition in
# `condition`, whether its unpretested group lost any. `lost_none` words
# the loss for repeated posttests. character(0) when every unpretested
# group lost some.
.no_loss_text <- function(lost, condition, lost_none = "lost no posttests",
                          lost_one = "lost a posttest") {
  if (all(lost)) return(character(0))
  if (!any(lost)) {
    return(sprintf(paste(
      "No unpretested group %s. If the data hold everyone randomized to those groups,",
      "this is the case for \"observed\": the default's assumption then holds only if",
      "the pretested groups' losses were unrelated to baseline."
    ), lost_one))
  }
  none <- condition[!lost]
  several <- length(none) > 1L
  sprintf(paste(
    "The unpretested %s %s %s: if the data hold everyone randomized",
    "to %s, the default's pretest %s there %s right only if the pretested",
    "groups' losses, taken together, were unrelated to baseline, because the",
    "default's center pools the pretested groups."
  ), .series_and(none), if (several) "groups" else "group", lost_none,
  if (several) "them" else "it", if (several) "effects" else "effect",
  if (several) "are" else "is")
}

# How a contrast at the two centers (or standardizations) disagrees, row by
# row of `a` and `o`, which hold `estimate`, `conf.low`, and `conf.high`
# for the same contrasts: the record's rule for giving no verdict on a
# pretest effect (issue #135, section 6). `null` is the value of no effect
# for each row, 0 for a difference and 1 for a ratio. One value per row:
# - "interval": the interval excludes the null value at one center only,
#   so whether the contrast differs from it depends on the center;
# - "sign_none": the direction (the sign of the estimate minus the null
#   value) differs and neither interval excludes the null value, so only
#   the direction depends on the center;
# - "sign_both": the direction differs and both intervals exclude it;
# - "sign": the direction differs and an interval is missing (a CR2 fit
#   without a cluster, or the bootstrap at the other standardization);
# - NA: no disagreement.
.disagreement_kind <- function(a, o, null = 0) {
  excludes <- function(r) r$conf.low > null | r$conf.high < null
  ea <- excludes(a)
  eo <- excludes(o)
  known <- !is.na(ea) & !is.na(eo)
  sign_differs <- (sign(a$estimate - null) != sign(o$estimate - null)) %in% TRUE
  ifelse(known & ea != eo, "interval",
         ifelse(!sign_differs, NA_character_,
                ifelse(!known, "sign", ifelse(ea, "sign_both", "sign_none"))))
}

# The kind of disagreement of each pretest effect of a fit (see
# .disagreement_kind()), named by its label. `s` is `pretest_sensitivity`.
.center_disagreement <- function(s) {
  e <- s$effects
  a <- e[e$center == "analyzed", , drop = FALSE]
  o <- e[e$center == "observed", , drop = FALSE]
  stats::setNames(.disagreement_kind(a, o), .center_labels(a))
}

# The sentences on a disagreement, one for each kind present, naming the
# contrasts `labels` of each kind (`kind`, from .disagreement_kind()).
# `between` names what differs ("centers" or "standardizations"), `level`
# is the confidence level, and `null` the value of no effect (0, or 1 for a
# ratio).
.disagreement_sentences <- function(labels, kind, level, between = "centers", null = 0) {
  one <- if (between == "centers") "center" else "standardization"
  pct <- format(100 * level)
  side <- if (null == 0) "the sign differs" else "the direction from 1 differs"
  sentence <- function(k) {
    l <- labels[kind %in% k]
    if (!length(l)) return(NULL)
    several <- length(l) > 1L
    switch(
      k,
      interval = sprintf(paste("For %s, only one %s's %s%% interval excludes %s: whether %s",
                               "from %s depends on the %s."),
                         .series_and(l), one, pct, null,
                         if (several) "they differ" else "it differs",
                         if (null == 0) "zero" else null, one),
      sign_none = sprintf(paste("For %s, %s between the %s, so the direction",
                                "depends on the %s; neither %s's %s%% interval excludes %s."),
                          .series_and(l), side, between, one, one, pct, null),
      sign_both = sprintf(paste("For %s, %s between the %s, so the direction",
                                "depends on the %s, although each %s's %s%% interval excludes %s."),
                          .series_and(l), side, between, one, one, pct, null),
      sign = sprintf("For %s, %s between the %s, so the direction depends on the %s.",
                     .series_and(l), side, between, one)
    )
  }
  unlist(lapply(c("interval", "sign_none", "sign_both", "sign"), sentence))
}

# The lines on the disagreement of a fit's pretest effects, if there is any.
.center_disagreement_text <- function(s, level) {
  kind <- .center_disagreement(s)
  .disagreement_sentences(names(kind), kind, level)
}

# The printed block of a fit whose two centers differ (the decision record
# on issue #135, section 4), as lines: the center used and the other, the
# pretest effects at both and their difference, the attrition, the
# pretested participants' mean pretest in and out of the model, the delta
# line (identity link only), and the cautions. character(0) when the
# centers do not differ or the fit predates #135. `mmrm` words the lack of a
# posttest for repeated posttests, where the delta line is not derived.
.pretest_block <- function(x, digits = 3, mmrm = FALSE) {
  s <- x$pretest_sensitivity
  if (is.null(x$pretest_center) || is.null(s) || !isTRUE(s$differ)) return(character(0))
  n <- x$pretest_n
  m <- max(digits - 1L, 1L)
  posttest <- if (mmrm) "without any posttest" else "without a posttest"
  outside <- .center_outside_phrase(n, posttest)
  analyzed <- .center_num(x$pretest_mean_analyzed, m)
  observed <- .center_num(x$pretest_mean_observed, m)
  centers <- if (identical(x$pretest_center, "observed")) {
    c(sprintf(paste("Centered at %s, the mean of every observed pretest (%d participants, %s;",
                    "pretest_center = \"observed\")."), observed, n[["observed"]], outside),
      sprintf("Pretested participants analyzed%s: %s (%d participants).", .analyzed_with_pretest(x),
              analyzed, n[["analyzed"]]))
  } else {
    c(sprintf("Centered at %s, the mean pretest of the %d pretested participants analyzed%s (default).",
              analyzed, n[["analyzed"]], .analyzed_with_pretest(x)),
      sprintf("Every observed pretest: %s (%d participants, %s).", observed, n[["observed"]],
              outside))
  }
  attrition <- if (mmrm) {
    c(.attrition_line(s$attrition, "Participants without any posttest"),
      .attrition_line(.mmrm_attrition(x),
                      sprintf("Posttests missing at occasion %s", x$occasions[length(x$occasions)])))
  } else {
    .attrition_line(s$attrition)
  }
  all_lost <- n[["no_posttest"]] == n[["observed"]] - n[["analyzed"]]
  means <- if (all_lost) {
    sprintf("Pretested participants' mean pretest: %s %s, %s with one.",
            .center_num(s$pretest_means[["not_analyzed"]], m), posttest,
            .center_num(s$pretest_means[["analyzed"]], m))
  } else {
    sprintf("Pretested participants' mean pretest: %s among those not analyzed, %s among those analyzed.",
            .center_num(s$pretest_means[["not_analyzed"]], m),
            .center_num(s$pretest_means[["analyzed"]], m))
  }
  delta <- if (mmrm) {
    # The case for "observed": an unpretested group in which everyone has a
    # posttest at some occasion, and so is in the model.
    un <- s$attrition[s$attrition$pretested == 0L & s$attrition$n > 0L, , drop = FALSE]
    c("No delta line is derived for repeated posttests.",
      .no_loss_text(un$posttest_missing > 0L, un$condition,
                    "has no participant without any posttest",
                    "has a participant without any posttest"))
  } else if (!is.null(s$delta)) {
    .center_delta_text(s, m)
  } else {
    sprintf("The pretest effects are on the %s scale, where no delta line is given.",
            .center_scale(x$family))
  }
  level <- if (is.null(x$conf_level)) 0.95 else x$conf_level
  wrap <- function(text) strwrap(text, width = 80, indent = 2, exdent = 2)
  c(
    if (n[["no_posttest"]] > 0L) {
      "Missing posttests and the pretest effects"
    } else {
      "Pretested participants not analyzed and the pretest effects"
    },
    unlist(lapply(centers, wrap)),
    .center_table(s$difference, digits),
    unlist(lapply(c(attrition, means, delta, .center_disagreement_text(s, level), s$note,
                    "The data cannot tell which center is right. See ?check_solomon_missing."),
                  wrap))
  )
}

# What a print method adds for the pretest centers: the block when the two
# centers differ; otherwise the attrition line whenever a posttest is
# missing, with the sentence for losses in the unpretested groups only.
# Nothing for complete data or for fits made before #135. Each part is
# preceded by an empty line.
.center_output <- function(x, digits = 3, mmrm = FALSE) {
  block <- .pretest_block(x, digits, mmrm)
  if (length(block)) {
    cat("\n", paste(block, collapse = "\n"), "\n", sep = "")
    return(invisible(NULL))
  }
  a <- if (mmrm) .mmrm_attrition(x) else .fit_attrition(x)
  what <- if (mmrm) {
    sprintf("Posttests missing at occasion %s", x$occasions[length(x$occasions)])
  } else {
    "Posttests missing"
  }
  line <- .attrition_line(a, what)
  if (!length(line)) return(invisible(NULL))
  pretest <- is.finite(.null_na(x$pretest_mean))
  lines <- c(line, if (pretest && .unpretested_only(a)) .unpretested_only_line)
  cat("\n", paste(vapply(lines, .wrap_lines, ""), collapse = "\n"), "\n", sep = "")
  invisible(NULL)
}

# The first line of plot_sensitization()'s caption: where the pretested
# groups' adjusted means are evaluated (`pretest_mean`, NA without a
# pretest covariate). When the two centers differ it names the center used
# and gives the other (issue #135); otherwise it is as before.
.sensitization_adjustment <- function(fit, pretest_mean) {
  if (is.na(pretest_mean)) return("Means from the fitted model")
  num <- function(v) formatC(v, format = "f", digits = 2)
  s <- fit$pretest_sensitivity
  # In a fit with the missing-indicator method (issue #142), the pretested
  # groups' means are also at the indicator's mean, which standardizes them
  # over the pretested participants analyzed, with a pretest or without.
  kept <- if (.has_missing_indicator(fit)) {
    sprintf(paste0(".\nPretested participants kept without a pretest: %d (missing-indicator method;\n",
                   "Groenwold et al., 2012)"), sum(fit$n_missing_pretest))
  } else {
    ""
  }
  if (is.null(fit$pretest_center) || is.null(s) || !isTRUE(s$differ)) {
    return(paste0(sprintf("Adjusted means: pretested groups at the mean pretest (%s)",
                          num(pretest_mean)), kept))
  }
  if (identical(fit$pretest_center, "observed")) {
    paste0(sprintf(paste0("Adjusted means: pretested groups at the mean of every observed pretest\n",
                          "(%s; pretested participants analyzed%s: %s)"),
                   num(fit$pretest_mean_observed), .analyzed_with_pretest(fit), num(fit$pretest_mean_analyzed)), kept)
  } else {
    paste0(sprintf(paste0("Adjusted means: pretested groups at the mean pretest of the %d pretested\n",
                          "participants analyzed%s (%s; every observed pretest: %s)"),
                   fit$pretest_n[["analyzed"]], .analyzed_with_pretest(fit), num(fit$pretest_mean_analyzed),
                   num(fit$pretest_mean_observed)), kept)
  }
}

# The participants each standardization of marginal_solomon() runs over, as
# in "the 74 pretested participants analyzed" or "every pretested
# participant with a pretest (118, 44 without an outcome)", from
# `standardized` (`n`). When the "observed" standardization leaves some of
# them out (issue #171), it is "97 of the 100 pretested participants with
# a pretest (7 without an outcome)".
.marginal_set_phrase <- function(n, center) {
  if (center == "analyzed") {
    return(sprintf("the %d pretested participants analyzed", n[["analyzed"]]))
  }
  outside <- .center_outside_phrase(
    c(analyzed = n[["analyzed"]], observed = n[["observed"]], no_posttest = n[["no_outcome"]]),
    "without an outcome"
  )
  left_out <- .marginal_left_out_n(n)
  if (left_out > 0L) {
    return(sprintf("%d of the %d pretested participants with a pretest (%s)", n[["observed"]],
                   n[["observed"]] + left_out, outside))
  }
  sprintf("every pretested participant with a pretest (%d, %s)", n[["observed"]], outside)
}

# The number of pretested participants with a pretest whom the "observed"
# standardization of marginal_solomon() leaves out (issue #171), from
# `standardized` (`n`); 0 for results made before it was counted.
.marginal_left_out_n <- function(n) {
  sum(n[intersect(c("missing_covariate", "unseen_level"), names(n))])
}

# The rows of marginal_solomon() that move with the standardization: every
# contrast but Treatment | unpretested, whose cells are standardized over
# the unpretested participants either way.
.marginal_moved <- function(d) d[d$contrast != "Treatment | unpretested", , drop = FALSE]

# The contrasts of a marginal_solomon() result whose direction differs
# between the two standardizations, or whose interval excludes the null
# value under one only (0 for a difference, 1 for a ratio): a data frame
# with `scale`, `contrast`, `null`, and `kind` (.disagreement_kind()), one
# row for each. The bootstrap gives no interval for the other
# standardization, so only the direction is compared there.
.marginal_disagreement <- function(s) {
  e <- .marginal_moved(s$effects)
  a <- e[e$center == "analyzed", , drop = FALSE]
  o <- e[e$center == "observed", , drop = FALSE]
  o <- o[match(paste(a$scale, a$contrast), paste(o$scale, o$contrast)), , drop = FALSE]
  null <- ifelse(grepl("difference", a$scale, fixed = TRUE), 0, 1)
  kind <- .disagreement_kind(a, o, null)
  keep <- !is.na(kind)
  data.frame(scale = a$scale[keep], contrast = a$contrast[keep], null = null[keep],
             kind = kind[keep], stringsAsFactors = FALSE)
}

# The marginal sentences of .disagreement_sentences(), scale by scale.
.marginal_disagreement_text <- function(dis, level) {
  unlist(lapply(unique(dis$scale), function(sc) {
    r <- dis[dis$scale == sc, , drop = FALSE]
    .disagreement_sentences(sprintf("%s (%s)", r$contrast, tolower(sc)), r$kind, level,
                            "standardizations", r$null[1])
  }))
}

# The lines print.solomon_marginal() adds when the two standardizations
# differ (issue #135): the participants of each, the moved contrasts under
# both with their difference, and the cautions. On the risk or rate scale
# the standardization moves Treatment | pretested, the average treatment
# effect, and Pretest x Treatment as well as the pretest effects; no delta
# line is given there.
.marginal_block <- function(x, digits = 3) {
  s <- x$pretest_sensitivity
  if (is.null(x$pretest_center) || is.null(s) || !isTRUE(s$differ)) return(character(0))
  n <- x$standardized
  used <- x$pretest_center
  wrap <- function(text) strwrap(text, width = 80, indent = 2, exdent = 2)
  lines <- c(
    if (n[["no_outcome"]] > 0L) {
      "Missing outcomes and the standardization of the pretested cells"
    } else {
      "Pretested participants not analyzed and the standardization of the pretested cells"
    },
    wrap(sprintf("Pretested cells standardized over %s (%s).", .marginal_set_phrase(n, used),
                 if (used == "analyzed") "default" else "pretest_center = \"observed\"")),
    wrap(sprintf("Over %s instead, the contrasts below move.",
                 .marginal_set_phrase(n, .other_center(used))))
  )
  d <- .marginal_moved(s$difference)
  for (sc in unique(d$scale)) {
    lines <- c(lines, paste0("  ", sc, ":"), .center_table(d[d$scale == sc, , drop = FALSE],
                                                           digits, indent = "    "))
  }
  dis <- .marginal_disagreement(s)
  c(lines, unlist(lapply(c(
    .attrition_line(x$attrition, "Outcomes missing"),
    "No delta line is given on the risk or rate scale.",
    .marginal_disagreement_text(dis, x$conf_level),
    s$note,
    "The data cannot tell which standardization is right. See ?marginal_solomon."
  ), wrap)))
}

# ---- pretested participants without a pretest (issue #142) -------------------------

# The printed block on the pretested participants who have a posttest but
# no pretest (issue #142; the decision record posted there, section 5), as
# lines. It is printed under either value of `missing_pretest`, whenever
# there is such a participant: how many there are and how they are
# handled, with what that assumes; the pretest effects with them kept and
# dropped, and the difference; their mean posttest beside that of the
# participants with a pretest; and the arm check, a description.
# character(0) when there is none, and for fits made before #142. `x` is a
# fit of fit_solomon_glm(), its summary, or, with `marginal`, a result of
# marginal_solomon(), whose table is given scale by scale.
.missing_pretest_block <- function(x, digits = 3, marginal = FALSE) {
  s <- x$missing_pretest_sensitivity
  n <- x$n_missing_pretest
  if (is.null(x$missing_pretest) || is.null(n) || !isTRUE(sum(n) > 0L)) return(character(0))
  k <- sum(n)
  kept <- identical(x$missing_pretest, "indicator")
  m <- max(digits - 1L, 1L)
  wrap <- function(text) strwrap(text, width = 80, indent = 2, exdent = 2)
  outcome <- if (marginal) "an outcome" else "a posttest"
  one <- k == 1L
  counts <- paste(sprintf("%d %s", n, names(n)), collapse = ", ")
  share <- .pct(x$missing_pretest_share)
  handling <- if (kept) {
    sprintf(paste(
      "%d pretested participant%s had %s but no pretest (%s; %s of the pretested participants",
      "analyzed). %s kept by the missing-indicator method (Groenwold et al., 2012;",
      "missing_pretest = \"indicator\"): their centered pretest is 0, and %s marks them. %s This",
      "assumes that losing the pretest did not depend on the treatment condition."
    ), k, if (one) "" else "s", outcome, counts, share, if (one) "That participant is" else "They are",
    .missing_pretest_column, if (marginal) {
      paste("The pretested cells average the predictions of the participants with and",
            "without a pretest, each at their own value of the indicator.")
    } else {
      paste("The coefficient of pretested compares the participants with a pretest; each",
            "pretest effect adds the indicator's mean times its coefficient, which standardizes",
            "it over all the pretested participants analyzed (solomonR's step).")
    })
  } else {
    sprintf(paste(
      "%d pretested participant%s had %s but no pretest (%s) and %s dropped (missing_pretest =",
      "\"drop\"). For the pretest effects, dropping them assumes that, in each pretested group,",
      "they would have had the same mean %s as those with a pretest. Kept by the",
      "missing-indicator method (Groenwold et al., 2012; missing_pretest = \"indicator\"), they",
      "would be %s of the pretested participants analyzed."
    ), k, if (one) "" else "s", outcome, counts, if (one) "was" else "were",
    if (marginal) "outcome" else "posttest", share)
  }
  lines <- c("Pretested participants without a pretest", wrap(handling))

  # The pretest effects kept and dropped.
  if (!is.null(s)) {
    e <- s$effects
    if (marginal) {
      for (sc in unique(e$scale)) {
        lines <- c(lines, paste0("  ", sc, ":"),
                   .handling_table(e[e$scale == sc, , drop = FALSE], digits, indent = "    "))
      }
    } else {
      link <- if (is.null(x$family)) "identity" else x$family$link
      if (!identical(link, "identity")) {
        lines <- c(lines, wrap(sprintf("Pretest effects on the %s scale:", .center_scale(x$family))))
      }
      lines <- c(lines, .handling_table(e, digits))
    }
  }
  notes <- character(0)

  # Their posttests beside those of the participants with a pretest.
  p <- s$posttests
  if (!is.null(p)) {
    p <- p[p$n_without > 0L, , drop = FALSE]
    notes <- c(notes, sprintf(
      "Mean %s without and with a pretest: %s.", if (marginal) "outcome" else "posttest",
      paste(sprintf("%s %s (%d) and %s (%d)", p$condition, .center_num(p$mean_without, m),
                    p$n_without, .center_num(p$mean_with, m), p$n_with), collapse = "; ")
    ))
  }

  # Separation: every kept participant has the same outcome, at a bound.
  if (!is.null(s) && !identical(s$separation, 0L)) {
    notes <- c(notes, .missing_pretest_separation_text(x, k, s$separation, marginal, kept))
  }

  # One such participant: the indicator fits them exactly.
  if (one && !marginal && (is.null(s) || identical(s$separation, 0L))) {
    notes <- c(notes, if (kept) {
      paste0(
        "With one such participant, the indicator fits that participant exactly: the ",
        "treatment contrasts equal those of missing_pretest = \"drop\", and each pretest ",
        "effect differs from that of \"drop\" by the indicator's mean times the participant's ",
        "residual.",
        if (identical(x$robust, "HC3")) {
          paste0(" Under HC3 the pretest effects then have no standard error (see \"HC3 at ",
                 "leverage 1\" in ?fit_solomon_glm); robust = \"none\" gives one.")
        } else if (identical(x$robust, "CR2")) {
          paste0(" Under CR2 the pretest effects keep a standard error, which leaves out the ",
                 "variance of that participant's posttest (see \"Pretested participants ",
                 "without a pretest\" in ?fit_solomon_glm); robust = \"none\" includes it.")
        }
      )
    } else {
      paste0(
        "With one such participant, missing_pretest = \"indicator\" gives the same treatment ",
        "contrasts, and pretest effects that differ by the difference shown; under HC3 they ",
        "would have no standard error."
      )
    })
  }

  # The arm check, a description.
  a <- x$missing_pretest_arm_check
  if (!is.null(a)) {
    check <- vapply(seq_len(nrow(a)), function(i) {
      r <- a[i, ]
      arms <- strsplit(r$comparison, " vs ", fixed = TRUE)[[1]]
      sprintf("%s %s, %s %s%s; mean observed pretest %s and %s%s",
              arms[2], .pct(r$share_control), arms[1], .pct(r$share_treatment),
              if (is.finite(r$fisher_p)) {
                paste0(" (Fisher's exact test, ", .apa_p(r$fisher_p, md = FALSE), ")")
              } else {
                ""
              },
              .center_num(r$mean_control, m), .center_num(r$mean_treatment, m),
              if (is.finite(r$p.value)) {
                sprintf(" (Hedges' g = %s, Welch %s)", .center_num(r$g, 2L),
                        .apa_p(r$p.value, md = FALSE))
              } else {
                ""
              })
    }, "")
    notes <- c(notes, sprintf(
      "Pretests missing: %s. A description; see ?check_solomon_missing.",
      paste(check, collapse = "; ")
    ))
  }

  if (identical(x$pretest_center, "observed")) {
    notes <- c(notes, paste(
      "The observed center assumes pretests were missing completely at random among the",
      "pretested participants."
    ))
  }

  # The bootstrap's refits (marginal_solomon()).
  refits <- x$indicator_refits
  if (marginal && !is.null(refits) && sum(refits) > 0L) {
    notes <- c(notes, sprintf(paste(
      "Bootstrap refits: %d resample%s held no participant kept without a pretest, and the",
      "indicator was left out; in %d, all those drawn had the same outcome, and their risks",
      "were set to it (see ?marginal_solomon)."
    ), refits[["none_kept"]], if (refits[["none_kept"]] == 1L) "" else "s",
    refits[["shared_outcome"]]))
  }
  c(lines, unlist(lapply(notes, wrap)))
}

# The sentence of the block on separation by the indicator: every kept
# participant has the same outcome, at a bound of the mean, so the
# indicator's coefficient has no finite estimate (Albert & Anderson, 1984).
# On the model's scale the pretest effects with them kept are infinite; on
# the risk or rate scale of marginal_solomon() their risks are set to that
# outcome, the limit of the fit.
.missing_pretest_separation_text <- function(x, k, direction, marginal, kept) {
  binary <- if (marginal) identical(x$outcome, "binary") else {
    x$family$family %in% c("binomial", "quasibinomial")
  }
  shared <- if (!binary) "a count of 0" else if (direction < 0) "the outcome 0" else "the outcome 1"
  who <- if (k == 1L) "That participant has" else sprintf("All %d have", k)
  if (marginal) {
    return(sprintf(paste(
      "%s %s, so the indicator's coefficient has no finite estimate (separation; Albert &",
      "Anderson, 1984). With them kept, their %s set to that outcome, the limit of the fit,",
      "and the other participants are fitted without them."
    ), who, shared, if (binary) "risks are" else "rates are"))
  }
  sprintf(paste(
    "%s %s, so the indicator's coefficient has no finite estimate on the %s scale",
    "(separation; Albert & Anderson, 1984). With them kept, the pretest effects are %sinfinity",
    "on that scale, with no standard error, test, or interval%s. marginal_solomon()",
    "estimates them from standardized %s, which are defined."
  ), who, shared, .center_scale(x$family), if (direction < 0) "minus " else "",
  if (kept) "; the treatment contrasts keep theirs" else "",
  if (binary) "risks" else "rates")
}

# What a print method adds for the pretested participants without a
# pretest: the block, preceded by an empty line, or nothing.
.missing_pretest_output <- function(x, digits = 3, marginal = FALSE) {
  block <- .missing_pretest_block(x, digits, marginal)
  if (length(block)) cat("\n", paste(block, collapse = "\n"), "\n", sep = "")
  invisible(NULL)
}
