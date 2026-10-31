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
  by_group <- .pct(s$range, 0L)
  spread <- if (by_group[1] == by_group[2]) {
    paste(by_group[1], "in every group")
  } else {
    paste(by_group[1], "to", by_group[2], "by group")
  }
  sprintf("%s: pretested %s, unpretested %s (%s).", what, .pct(s$pretested),
          .pct(s$unpretested), spread)
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
  "protected against them. See ?check_solomon_missing."
)

# The pretested participants outside the model, as in "44 without a
# posttest" or "5 not analyzed, 3 of them without a posttest", from
# `pretest_n` (`n`). `posttest` words the lack of a posttest.
.center_outside_phrase <- function(n, posttest = "without a posttest") {
  outside <- n[["observed"]] - n[["analyzed"]]
  none <- n[["no_posttest"]]
  if (none == outside) {
    sprintf("%d %s", outside, posttest)
  } else if (none == 0L) {
    sprintf("%d not analyzed", outside)
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
               .center_num(d$observed, digits), .center_num(d$difference, digits))
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
  if (!any(lost)) {
    text <- paste(
      "No unpretested group lost a posttest. If the data hold everyone randomized to",
      "those groups, this is the case for \"observed\": the default's assumption then",
      "holds only if the pretested groups' losses were unrelated to baseline."
    )
  } else if (!all(lost)) {
    none <- d$condition[!lost]
    text <- c(text, sprintf(paste(
      "The unpretested %s %s lost no posttests: if the data hold everyone randomized",
      "to %s, the default's assumption holds in %s only if the pretested %s",
      "losses were unrelated to baseline."
    ), .series_and(none), if (length(none) > 1L) "groups" else "group",
    if (length(none) > 1L) "them" else "it",
    if (length(none) > 1L) "those conditions" else "that condition",
    if (length(none) > 1L) "groups'" else "group's"))
  }
  text
}

# The pretest effects whose sign differs between the two centers, or whose
# interval excludes 0 at one center only: the record's rule for giving no
# verdict on a pretest effect (issue #135, section 6). `s` is
# `pretest_sensitivity`; returns the labels of each kind.
.center_disagreement <- function(s) {
  e <- s$effects
  a <- e[e$center == "analyzed", , drop = FALSE]
  o <- e[e$center == "observed", , drop = FALSE]
  labels <- .center_labels(a)
  sign_differs <- sign(a$estimate) != sign(o$estimate)
  excludes <- function(r) r$conf.low > 0 | r$conf.high < 0
  ea <- excludes(a)
  eo <- excludes(o)
  interval_differs <- !sign_differs & !is.na(ea) & !is.na(eo) & ea != eo
  list(sign = labels[sign_differs], interval = labels[interval_differs])
}

# The lines on the disagreement, if there is any.
.center_disagreement_text <- function(s, level) {
  dis <- .center_disagreement(s)
  parts <- c(
    if (length(dis$sign)) {
      sprintf("for %s, the sign differs between the centers", .series_and(dis$sign))
    },
    if (length(dis$interval)) {
      sprintf("for %s, only one center's %s%% interval excludes 0", .series_and(dis$interval),
              format(100 * level))
    }
  )
  if (!length(parts)) return(character(0))
  paste0(.capitalize(paste(parts, collapse = "; ")),
         ". Whether these effects differ from zero depends on the center.")
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
    "No delta line is derived for repeated posttests."
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

# The contrasts of a marginal_solomon() result whose sign differs between
# the two standardizations, or whose interval excludes 0 under one only, as
# "Risk difference: Pretest main effect". The bootstrap gives no interval
# for the other standardization, so only the sign is compared there.
.marginal_disagreement <- function(s) {
  e <- .marginal_moved(s$effects)
  a <- e[e$center == "analyzed", , drop = FALSE]
  o <- e[e$center == "observed", , drop = FALSE]
  o <- o[match(paste(a$scale, a$contrast), paste(o$scale, o$contrast)), , drop = FALSE]
  excludes <- function(r) r$conf.low > 0 | r$conf.high < 0
  ea <- excludes(a)
  eo <- excludes(o)
  differs <- (sign(a$estimate) != sign(o$estimate)) %in% TRUE |
    (!is.na(ea) & !is.na(eo) & ea != eo)
  if (!any(differs)) return(character(0))
  paste0(a$scale[differs], ": ", a$contrast[differs])
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
    "No delta line is given on the risk or rate scale.",
    if (length(dis)) {
      sprintf(paste("For %s, the sign differs between the standardizations or only one",
                    "interval excludes 0: whether %s from zero depends on the standardization."),
              .series_and(dis), if (length(dis) > 1L) "they differ" else "it differs")
    },
    s$note,
    "The data cannot tell which standardization is right. See ?marginal_solomon."
  ), wrap)))
}
