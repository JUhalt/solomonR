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
# `pretest_sensitivity$difference`.
.center_differences <- function(d) {
  zero <- abs(d$difference) <= 1e-12 * pmax(abs(d$analyzed), abs(d$observed))
  ifelse(is.finite(d$difference) & zero, 0, d$difference)
}

# The attrition of a fit's Solomon groups: a data frame with `group`,
# `condition`, `pretested`, `n`, `posttest_missing`, and `share`, one row
# for each group. It is `pretest_sensitivity$attrition` when the fit has a
# pretest covariate, and is otherwise counted from the data the model was
# fitted to. NULL for fits made before #135 and for fits whose data do not
# hold the posttests (fit_solomon_mmrm(), which .mmrm_attrition() covers).
.fit_attrition <- function(x) {
  if (is.null(x$pretest_center)) return(NULL)
  # A summary carries the table it was given.
  if (!is.null(x$attrition)) return(x$attrition)
  if (!is.null(x$pretest_sensitivity$attrition)) return(x$pretest_sensitivity$attrition)
  d <- x$data
  if (!is.data.frame(d) || !all(c("y", "pretested") %in% names(d))) return(NULL)
  condition <- as.character(if ("condition" %in% names(d)) d$condition else d$treat)
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
      sprintf("Pretested participants analyzed: %s (%d participants).", analyzed,
              n[["analyzed"]]))
  } else {
    c(sprintf("Centered at %s, the mean pretest of the %d pretested participants analyzed (default).",
              analyzed, n[["analyzed"]]),
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
  if (is.null(fit$pretest_center) || is.null(s) || !isTRUE(s$differ)) {
    return(sprintf("Adjusted means: pretested groups at the mean pretest (%s)", num(pretest_mean)))
  }
  if (identical(fit$pretest_center, "observed")) {
    sprintf(paste0("Adjusted means: pretested groups at the mean of every observed pretest\n",
                   "(%s; pretested participants analyzed: %s)"),
            num(fit$pretest_mean_observed), num(fit$pretest_mean_analyzed))
  } else {
    sprintf(paste0("Adjusted means: pretested groups at the mean pretest of the %d pretested\n",
                   "participants analyzed (%s; every observed pretest: %s)"),
            fit$pretest_n[["analyzed"]], num(fit$pretest_mean_analyzed),
            num(fit$pretest_mean_observed))
  }
}

# The participants each standardization of marginal_solomon() runs over, as
# in "the 74 pretested participants analyzed" or "every pretested
# participant with a pretest (118, 44 without an outcome)", from
# `standardized` (`n`).
.marginal_set_phrase <- function(n, center) {
  if (center == "analyzed") {
    return(sprintf("the %d pretested participants analyzed", n[["analyzed"]]))
  }
  outside <- .center_outside_phrase(
    c(analyzed = n[["analyzed"]], observed = n[["observed"]], no_posttest = n[["no_outcome"]]),
    "without an outcome"
  )
  sprintf("every pretested participant with a pretest (%d, %s)", n[["observed"]], outside)
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
