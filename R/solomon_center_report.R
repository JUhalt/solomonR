# The sentences report_solomon() adds for the two pretest centers (issue
# #135): the center used and its assumption, the pretest effects at the
# other center and the gap when the two centers differ, a sentence when
# whether a pretest effect differs from zero, or its direction, depends on
# the center, and an attrition sentence whenever a posttest is missing.
# Reports of complete data and of fits made before #135 are as before.

# The assumption under which each center gives the pretest effects (the
# decision record on issue #135, section 4), as report_solomon() states it.
.center_assumption <- c(
  analyzed = paste("assumes that, within each treatment condition, the pretested and the",
                   "unpretested groups lost the same kinds of participants"),
  observed = paste("assumes that the pretested participants' losses depended at most on",
                   "their pretest score and that the unpretested participants' losses were",
                   "unrelated to baseline and to the posttest")
)

# The center at which a fit's pretest effects were evaluated, for the
# report: "the mean pretest score of the 74 pretested participants analyzed
# (3.13)" or "the mean of every observed pretest score among the pretested
# participants (3.16; 118 participants, 44 of them without a posttest)".
.center_phrase <- function(x, center, digits, posttest = "without a posttest") {
  n <- x$pretest_n
  if (center == "analyzed") {
    return(sprintf("the mean pretest score of the %d pretested participants analyzed%s (%s)",
                   n[["analyzed"]], .analyzed_with_pretest(x), .apa_num(x$pretest_mean_analyzed, digits)))
  }
  sprintf(paste("the mean of every observed pretest score among the pretested participants",
                "(%s; %d participants, %s)"),
          .apa_num(x$pretest_mean_observed, digits), n[["observed"]],
          .center_outside_phrase(n, posttest, of_them = TRUE))
}

# One row of a pretest-effects table as a clause: the estimate, with its
# interval and test when it has a standard error. Without one (a CR2 fit may
# give none at the other center, and an HC3 standard error is undefined at
# leverage 1, issue #173), the estimate is followed by .no_inference_phrase,
# as in .effect_clause().
.center_clause <- function(r, level, digits, md, where = NULL, p = "p.value", p_label = NULL) {
  if (!is.finite(r$std.error)) {
    return(paste0(.apa_num(r$estimate, digits), if (!is.null(where)) paste0(" ", where),
                  if (is.finite(r$estimate)) .no_inference_phrase))
  }
  .effect_clause(r, level, digits, md, p = p, p_label = p_label, where = where)
}

# How much larger the pretest effects are at the analyzed center than at
# the observed center, as in "0.01 smaller", from the differences `d`
# (analyzed minus observed).
.center_larger <- function(d, digits) {
  paste(.apa_num(abs(d), digits), ifelse(d < 0, "smaller", "larger"))
}

# The sentence on the difference between the centers: one number when the
# differences are all equal (a four-group fit, or several treatments), as
# they are when the model has one pretest slope.
.center_difference_sentence <- function(d, digits) {
  if (length(unique(round(d, 10))) == 1L) {
    sprintf("Each pretest effect was %s at the analyzed center than at the observed center.",
            .center_larger(d[1], digits))
  } else {
    sprintf("The pretest effects were %s at the analyzed center than at the observed center.",
            .series_and(.center_larger(d, digits)))
  }
}

# The sentences report_solomon() adds for a fit whose two centers differ:
# the center used and its assumption; the other center, its assumption, and
# the pretest effects there (`other_text`, a clause such as "the pretest
# effect was ..."); the difference between the centers
# (`difference_sentence`); and, when only one center's interval excludes 0,
# that whether those effects differ from zero depends on the center, or,
# when only the sign differs, that their direction does (the record's rule
# for giving no verdict; .disagreement_report()).
# `names` gives the phrase of each pretest effect, in the order of the
# analyzed rows of `pretest_sensitivity$effects`. character(0) when the
# centers do not differ or the fit predates #135.
.center_sentences <- function(x, digits, other_text, difference_sentence, names,
                              posttest = "without a posttest") {
  s <- x$pretest_sensitivity
  if (is.null(x$pretest_center) || is.null(s) || !isTRUE(s$differ)) return(character(0))
  used <- x$pretest_center
  other <- .other_center(used)
  c(
    sprintf("The pretest effects were evaluated at %s, which %s.",
            .center_phrase(x, used, digits, posttest), .center_assumption[[used]]),
    sprintf("At %s, which %s, %s.",
            .center_phrase(x, other, digits, posttest), .center_assumption[[other]], other_text),
    difference_sentence,
    .disagreement_report(names, unname(.center_disagreement(s)), x$conf_level),
    s$note
  )
}

# The report's sentences on a disagreement between the two centers (or
# standardizations, `what`), one for each kind present
# (.disagreement_kind()): whether a contrast differs from zero depends on
# the center only when one center's interval excludes 0 and the other's
# does not; when only the sign differs, its direction does. `names` gives
# the phrase of each contrast and `kind` its kind.
.disagreement_report <- function(names, kind, level, what = "pretest center", null = 0) {
  pct <- format(100 * level)
  both <- if (what == "pretest center") "at both centers" else "under both standardizations"
  value <- if (null == 0) "zero" else "1"
  sentence <- function(k) {
    l <- names[kind %in% k]
    if (!length(l)) return(NULL)
    several <- length(l) > 1L
    direction <- sprintf("The %s of %s %s on the %s", if (several) "directions" else "direction",
                         .series_and(l), if (several) "depend" else "depends", what)
    interval <- sprintf("%s %s%% confidence %s", if (several) "their" else "its", pct,
                        if (several) "intervals" else "interval")
    switch(
      k,
      interval = sprintf("Whether %s %s from %s therefore depends on the %s.", .series_and(l),
                         if (several) "differ" else "differs", value, what),
      sign_none = sprintf("%s; %s %s %s %s.", direction, interval,
                          if (several) "include" else "includes", null, both),
      sign_both = sprintf("%s, although %s %s %s %s.", direction, interval,
                          if (several) "exclude" else "excludes", null, both),
      sign = paste0(direction, ".")
    )
  }
  unlist(lapply(c("interval", "sign_none", "sign_both", "sign"), sentence))
}

# The attrition sentence of a report, whenever a posttest is missing, as in
# "Posttests were missing for 37.3% of the pretested and 36.6% of the
# unpretested participants (29% to 46% by group).", with the sentence for
# losses in the unpretested groups only when the fit has a pretest
# covariate. `what` words repeated posttests.
.attrition_sentences <- function(a, pretest, what = "Posttests were missing") {
  s <- .attrition_shares(a)
  if (is.null(s) || !sum(s$missing)) return(character(0))
  c(
    sprintf("%s for %s of the pretested and %s of the unpretested participants (%s).",
            what, .pct(s$pretested), .pct(s$unpretested), .attrition_spread(s)),
    if (pretest && .unpretested_only(a)) {
      paste("Only the unpretested groups lost posttests, so the two pretest centers coincide;",
            "those losses cannot be checked against the pretest. The pretest effects assume,",
            "at either center, that the unpretested participants' losses were unrelated to",
            "baseline and to the posttest.")
    }
  )
}

# The center sentences of a four-group fit (fit_solomon_glm()).
.four_group_center_sentences <- function(x, digits, md) {
  s <- x$pretest_sensitivity
  if (is.null(s) || !isTRUE(s$differ)) return(character(0))
  level <- x$conf_level
  other <- s$effects[s$effects$center == .other_center(x$pretest_center), , drop = FALSE]
  rows <- other[match(.solomon_pretest_order, other$contrast), , drop = FALSE]
  other_text <- sprintf(
    "the pretest effect was %s, and %s; their average, the pretest main effect, was %s",
    .center_clause(rows[1, ], level, digits, md, "among control participants"),
    .center_clause(rows[2, ], level, digits, md, "among treated participants"),
    .center_clause(rows[3, ], level, digits, md)
  )
  .center_sentences(x, digits, other_text,
                    .center_difference_sentence(s$difference$difference, digits), c(
    "the pretest effect among control participants",
    "the pretest effect among treated participants",
    "the pretest main effect"
  ))
}

# The center sentences of a design with several treatments: the pretest
# effect in each condition and the pretest main effect, with the treatments'
# p-values adjusted as in the main sentence.
.ngroup_center_sentences <- function(x, digits, md) {
  s <- x$pretest_sensitivity
  if (is.null(s) || !isTRUE(s$differ)) return(character(0))
  level <- x$conf_level
  adjusted <- x$adjust != "none"
  label <- if (adjusted) switch(x$adjust, holm = "Holm-adjusted", bonferroni = "Bonferroni-adjusted")
  other <- s$effects[s$effects$center == .other_center(x$pretest_center), , drop = FALSE]
  by_condition <- other$contrast != .solomon_pretest_order[3]
  clauses <- vapply(which(by_condition), function(i) {
    r <- other[i, ]
    where <- sprintf("in the %s condition", r$comparison)
    if (adjusted && r$contrast == .solomon_pretest_order[2]) {
      .center_clause(r, level, digits, md, where, "p.adjusted", label)
    } else {
      .center_clause(r, level, digits, md, where)
    }
  }, "")
  main <- other[!by_condition, , drop = FALSE]
  n <- length(clauses)
  listed <- paste0(paste(clauses[-n], collapse = "; "), "; and ", clauses[n])
  other_text <- sprintf("the pretest effect was %s; averaged over the conditions, the pretest main effect was %s",
                        listed, .center_clause(main, level, digits, md))
  a <- s$effects[s$effects$center == "analyzed", , drop = FALSE]
  names <- ifelse(a$contrast == .solomon_pretest_order[3], "the pretest main effect",
                  sprintf("the pretest effect in the %s condition", a$comparison))
  .center_sentences(x, digits, other_text,
                    .center_difference_sentence(s$difference$difference, digits), names)
}

# The center sentences of a repeated-measures fit: the pretest effects at
# the other center at each occasion, and their difference at each occasion
# (the pretest slope differs between occasions).
.mmrm_center_sentences <- function(x, digits, md) {
  s <- x$pretest_sensitivity
  if (is.null(s) || !isTRUE(s$differ)) return(character(0))
  level <- x$conf_level
  other <- s$effects[s$effects$center == .other_center(x$pretest_center), , drop = FALSE]
  by_occasion <- vapply(unique(other$occasion), function(o) {
    rows <- other[other$occasion == o, , drop = FALSE]
    rows <- rows[match(.solomon_pretest_order, rows$contrast), , drop = FALSE]
    sprintf("at occasion %s, %s, and %s, with a pretest main effect of %s", o,
            .center_clause(rows[1, ], level, digits, md, "among control participants"),
            .center_clause(rows[2, ], level, digits, md, "among treated participants"),
            .center_clause(rows[3, ], level, digits, md))
  }, "")
  other_text <- paste("the pretest effect was", paste(by_occasion, collapse = "; "))
  # One pretest slope per occasion, so the difference is the same for the
  # three pretest effects of an occasion.
  d <- s$difference[s$difference$contrast == .solomon_pretest_order[3], , drop = FALSE]
  difference_sentence <- sprintf(
    "At the analyzed center, the pretest effects were %s than at the observed center.",
    .series_and(sprintf("%s at occasion %s", .center_larger(d$difference, digits), d$occasion))
  )
  a <- s$effects[s$effects$center == "analyzed", , drop = FALSE]
  phrase <- c(
    "Pretest effect | control" = "the pretest effect among control participants",
    "Pretest effect | treated" = "the pretest effect among treated participants",
    "Pretest main effect" = "the pretest main effect"
  )
  names <- sprintf("%s at occasion %s", phrase[a$contrast], a$occasion)
  .center_sentences(x, digits, other_text, difference_sentence, names, "without any posttest")
}

# The sentences of a marginal_solomon() report when the two standardizations
# of the pretested cells differ (issue #135): the participants standardized
# over and the assumption, and the moved contrasts under the other
# standardization, scale by scale. No delta line is given on these scales.
# The note on the pretested participants left out of the "observed"
# standardization (issue #171) is given whether or not the two differ; it
# says "pretest score", as the sets do here.
.marginal_center_sentences <- function(x, digits) {
  s <- x$pretest_sensitivity
  if (is.null(x$pretest_center) || is.null(s)) return(character(0))
  note <- sub("with a pretest leaves out", "with a pretest score leaves out", s$note, fixed = TRUE)
  if (!isTRUE(s$differ)) return(note)
  n <- x$standardized
  used <- x$pretest_center
  other <- .other_center(used)
  set <- function(center) {
    sub("with a pretest (", "with a pretest score (", .marginal_set_phrase(n, center), fixed = TRUE)
  }
  # Short names of the contrasts that move, in the order of the table.
  phrase <- c(
    "ATE (avg over pretest)" = "the average treatment effect",
    "Pretest x Treatment" = "the Pretest x Treatment interaction",
    "Treatment | pretested" = "the treatment effect among pretested participants",
    "Pretest effect | control" = "the pretest effect among control participants",
    "Pretest effect | treated" = "the pretest effect among treated participants",
    "Pretest main effect" = "the pretest main effect"
  )
  d <- .marginal_moved(s$difference)
  by_scale <- vapply(unique(d$scale), function(sc) {
    r <- d[d$scale == sc, , drop = FALSE]
    sprintf("On the %s scale, they were %s.", tolower(sc),
            .series_and(sprintf("%s for %s", .apa_num(r[[other]], digits), phrase[r$contrast])))
  }, "")
  dis <- .marginal_disagreement(s)
  c(
    sprintf("The pretested cells were standardized over %s, which %s.", set(used),
            .center_assumption[[used]]),
    sprintf(paste("Standardized over %s instead, which %s, the contrasts that depend on the",
                  "standardization were as follows."), set(other), .center_assumption[[other]]),
    by_scale,
    unlist(lapply(unique(dis$scale), function(sc) {
      r <- dis[dis$scale == sc, , drop = FALSE]
      .disagreement_report(sprintf("%s on the %s scale", phrase[r$contrast], tolower(sc)),
                           r$kind, x$conf_level, "standardization", r$null[1])
    })),
    note
  )
}

# The sentence of a fit_solomon_1949() report when the two inferred
# pretests differ (issue #135): the one used, from the participants with a
# posttest by default, and the other with the interaction and the
# unpretested groups' improvements at it.
.solomon_1949_center_sentences <- function(x, digits, md) {
  s <- x$pretest_sensitivity
  if (is.null(x$pretest_center) || is.null(s)) return(character(0))
  if (isTRUE(all.equal(s$inferred_pretest[1], s$inferred_pretest[2]))) return(character(0))
  used <- x$pretest_center
  other <- s[s$center == .other_center(used), , drop = FALSE]
  I <- if (md) "*I*" else "I"
  four <- identical(x$design, "four-group")
  source <- c(
    analyzed = "the pretests of the pretested participants with a posttest",
    observed = "every observed pretest of the pretested groups"
  )
  without <- sum(x$inferred_from$without_posttest)
  c(
    sprintf(paste("The inferred pretest mean was computed from %s; %d pretested participant%s",
                  "had a pretest and no posttest."), source[[used]], without,
            if (without == 1L) "" else "s"),
    sprintf("Computed from %s instead, it would be %s, with %s = %s%s.",
            source[[.other_center(used)]], .apa_num(other$inferred_pretest, digits), I,
            .apa_num(other$I, digits),
            if (four) {
              "; in the four-group design the inferred pretest cancels from the interaction"
            } else {
              ""
            })
  )
}

# ---- pretested participants without a pretest (issue #142) -------------------------

# The sentences report_solomon() adds on the pretested participants who
# have a posttest but no pretest (issue #142): how many there are and how
# they were handled, with what that assumes; the pretest effects under the
# other handling; and, where they apply, the sentences on one such
# participant and on separation. `other_text` gives the pretest effects
# under the other handling as a clause ("the pretest effect was ..."), and
# `outcome` words the posttest ("an outcome" for marginal_solomon()).
# character(0) when there is no such participant and for fits made before
# #142.
.missing_pretest_sentences <- function(x, digits, other_text, outcome = "a posttest",
                                       marginal = FALSE) {
  s <- x$missing_pretest_sensitivity
  n <- x$n_missing_pretest
  if (is.null(x$missing_pretest) || is.null(n) || !isTRUE(sum(n) > 0L)) return(character(0))
  k <- sum(n)
  kept <- identical(x$missing_pretest, "indicator")
  counts <- .series_and(sprintf("%d in the %s condition", n[n > 0L], names(n)[n > 0L]))
  share <- .pct(x$missing_pretest_share)
  handled <- if (kept) {
    sprintf(paste0(
      "Of the pretested participants analyzed, %d (%s) had %s but no pretest score (%s). ",
      "%s kept by the missing-indicator method (Groenwold et al., 2012): their ",
      "centered pretest score was set to 0, an indicator marked them, and %s. This assumes ",
      "that losing the pretest did not depend on the treatment condition."
    ), k, share, outcome, counts, if (k == 1L) "That participant was" else "They were",
    if (marginal) {
      "each was standardized over with their own indicator"
    } else {
      "the pretest effects were standardized over all the pretested participants analyzed"
    })
  } else {
    sprintf(paste0(
      "Pretested participants who had %s but no pretest score were excluded from the model ",
      "(%d: %s). For the pretest effects, that assumes that, within each pretested ",
      "condition, they would have had the same mean %s as the participants with a pretest ",
      "score."
    ), outcome, k, counts, if (marginal) "outcome" else "posttest")
  }
  other <- if (!is.null(other_text)) {
    if (kept) {
      sprintf("With %s excluded instead, %s.",
              if (k == 1L) "this participant" else "these participants", other_text)
    } else {
      sprintf(paste0("With %s kept by the missing-indicator method (Groenwold et al., 2012), ",
                     "%s."), if (k == 1L) "this participant" else "them", other_text)
    }
  }
  separation <- if (!is.null(s) && !identical(s$separation, 0L)) {
    binary <- if (marginal) identical(x$outcome, "binary") else {
      x$family$family %in% c("binomial", "quasibinomial")
    }
    shared <- if (!binary) "a count of 0" else if (s$separation < 0) "the outcome 0" else
      "the outcome 1"
    if (marginal) {
      sprintf(paste0(
        "Every such participant had %s, so the indicator's coefficient had no finite ",
        "estimate (separation; Albert & Anderson, 1984); with them kept, their %s set to ",
        "that outcome, the limit of the fit."
      ), shared, if (binary) "risks were" else "rates were")
    } else {
      sprintf(paste0(
        "Every such participant had %s, so the indicator's coefficient had no finite ",
        "estimate on the %s scale (separation; Albert & Anderson, 1984), and the pretest ",
        "effects with them kept are %sinfinite on that scale, without a standard error; ",
        "they can be estimated from standardized %s."
      ), shared, .center_scale(x$family), if (s$separation < 0) "negatively " else "",
      if (binary) "risks" else "rates")
    }
  }
  single <- if (k == 1L && kept && !marginal && is.null(separation)) {
    paste0(
      "With one such participant, the indicator fitted that participant exactly, so the ",
      "treatment contrasts equal those of the analysis that excludes the participant.",
      if (identical(x$robust, "HC3")) {
        " The pretest effects then have no HC3 standard error."
      } else if (identical(x$robust, "CR2")) {
        paste0(" The CR2 standard errors of the pretest effects leave out the variance of ",
               "that participant's posttest.")
      }
    )
  }
  c(handled, other, separation, single)
}

# The pretest effects under the other handling as a clause, from the rows
# of `missing_pretest_sensitivity$effects` (`e`) and the column of the
# other handling (`column`): for a four-group fit, among control and among
# treated participants and their average; for several treatments, in each
# condition and averaged over the conditions. A value that is not finite
# is worded.
.other_handling_clause <- function(e, column, digits) {
  value <- function(v) {
    if (is.na(v)) "not available" else if (is.finite(v)) .apa_num(v, digits) else
      if (v < 0) "minus infinity" else "infinity"
  }
  v <- vapply(e[[column]], value, "")
  if ("comparison" %in% names(e)) {
    main <- e$contrast == .solomon_pretest_order[3]
    return(sprintf(
      "the pretest effect was %s; averaged over the conditions, the pretest main effect was %s",
      .series_and(sprintf("%s in the %s condition", v[!main], e$comparison[!main])), v[main]
    ))
  }
  at <- match(.solomon_pretest_order, e$contrast)
  sprintf(paste0("the pretest effect was %s among control participants and %s among treated ",
                 "participants, and the pretest main effect was %s"), v[at[1]], v[at[2]], v[at[3]])
}

# The sentences of a four-group or several-treatment fit.
.fit_missing_pretest_sentences <- function(x, digits) {
  s <- x$missing_pretest_sensitivity
  if (is.null(s)) return(character(0))
  column <- if (identical(x$missing_pretest, "indicator")) "dropped" else "kept"
  .missing_pretest_sentences(x, digits, .other_handling_clause(s$effects, column, digits))
}

# The sentences of a marginal_solomon() result: the pretest effects under
# the other handling on each scale.
.marginal_missing_pretest_sentences <- function(x, digits) {
  s <- x$missing_pretest_sensitivity
  if (is.null(s)) return(character(0))
  column <- if (identical(x$missing_pretest, "indicator")) "dropped" else "kept"
  e <- s$effects
  by_scale <- vapply(unique(e$scale), function(sc) {
    r <- e[e$scale == sc, , drop = FALSE]
    r <- r[match(.solomon_pretest_order, r$contrast), , drop = FALSE]
    v <- ifelse(is.na(r[[column]]), "not available", .apa_num(r[[column]], digits))
    sprintf("on the %s scale, %s among control participants, %s among treated participants, and %s for the pretest main effect",
            tolower(sc), v[1], v[2], v[3])
  }, "")
  .missing_pretest_sentences(
    x, digits, paste("the pretest effects were as follows:", paste(by_scale, collapse = "; ")),
    outcome = "an outcome", marginal = TRUE
  )
}

# The references of those sentences: Groenwold et al. (2012) for the
# method, and Albert and Anderson (1984) under separation.
.missing_pretest_refs <- function(x) {
  s <- x$missing_pretest_sensitivity
  if (is.null(x$missing_pretest) || is.null(s)) return(character(0))
  c("groenwold2012", if (!identical(s$separation, 0L)) "albert1984")
}

# Whether the pretest effects of a fit are infinite because the indicator
# separates the outcome: the report then leaves out their sentence, which
# the sentences above replace.
.pretest_effects_separated <- function(x) {
  identical(x$missing_pretest, "indicator") &&
    !is.null(x$missing_pretest_sensitivity) &&
    !identical(x$missing_pretest_sensitivity$separation, 0L)
}
