# The sentences report_solomon() adds for the two pretest centers (issue
# #135): the center used and its assumption, the pretest effects at the
# other center and the gap when the two centers differ, a sentence when
# whether a pretest effect differs from zero depends on the center, and an
# attrition sentence whenever a posttest is missing. Reports of complete data
# and of fits made before #135 are as before.

# The assumption under which each center gives the pretest effects (the
# decision record on issue #135, section 4), as report_solomon() states it.
.center_assumption <- c(
  analyzed = paste("assumes that, within each treatment condition, the pretested and the",
                   "unpretested groups lost the same kinds of participants"),
  observed = paste("assumes that the pretested participants' losses depended at most on",
                   "their pretest score and that the unpretested participants' losses were",
                   "unrelated to baseline")
)

# The center at which a fit's pretest effects were evaluated, for the
# report: "the mean pretest score of the 74 pretested participants analyzed
# (3.13)" or "the mean of every observed pretest score among the pretested
# participants (3.16; 118 participants, 44 of them without a posttest)".
.center_phrase <- function(x, center, digits, posttest = "without a posttest") {
  n <- x$pretest_n
  if (center == "analyzed") {
    return(sprintf("the mean pretest score of the %d pretested participants analyzed (%s)",
                   n[["analyzed"]], .apa_num(x$pretest_mean_analyzed, digits)))
  }
  outside <- n[["observed"]] - n[["analyzed"]]
  none <- n[["no_posttest"]]
  who <- if (none == outside) {
    sprintf("%d of them %s", outside, posttest)
  } else if (none == 0L) {
    sprintf("%d of them not analyzed", outside)
  } else {
    sprintf("%d of them not analyzed, %d %s", outside, none, posttest)
  }
  sprintf(paste("the mean of every observed pretest score among the pretested participants",
                "(%s; %d participants, %s)"),
          .apa_num(x$pretest_mean_observed, digits), n[["observed"]], who)
}

# One row of a pretest-effects table as a clause: the estimate, with its
# interval and test when it has a standard error (a CR2 fit may give none at
# the other center).
.center_clause <- function(r, level, digits, md, where = NULL, p = "p.value", p_label = NULL) {
  if (!is.finite(r$std.error)) {
    return(paste0(.apa_num(r$estimate, digits), if (!is.null(where)) paste0(" ", where)))
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
# (`difference_sentence`); and, when the sign differs between the centers or
# only one center's interval excludes 0, that whether those effects differ
# from zero depends on the center (the record's rule for giving no verdict).
# `names` gives the phrase of each pretest effect, in the order of the
# analyzed rows of `pretest_sensitivity$effects`. character(0) when the
# centers do not differ or the fit predates #135.
.center_sentences <- function(x, digits, other_text, difference_sentence, names,
                              posttest = "without a posttest") {
  s <- x$pretest_sensitivity
  if (is.null(x$pretest_center) || is.null(s) || !isTRUE(s$differ)) return(character(0))
  used <- x$pretest_center
  other <- .other_center(used)
  dis <- .center_disagreement(s)
  a <- s$effects[s$effects$center == "analyzed", , drop = FALSE]
  flagged <- .center_labels(a) %in% c(dis$sign, dis$interval)
  c(
    sprintf("The pretest effects were evaluated at %s, which %s.",
            .center_phrase(x, used, digits, posttest), .center_assumption[[used]]),
    sprintf("At %s, which %s, %s.",
            .center_phrase(x, other, digits, posttest), .center_assumption[[other]], other_text),
    difference_sentence,
    if (any(flagged)) {
      sprintf("Whether %s %s from zero therefore depends on the pretest center.",
              .series_and(names[flagged]), if (sum(flagged) > 1L) "differ" else "differs")
    },
    s$note
  )
}

# The attrition sentence of a report, whenever a posttest is missing, as in
# "Posttests were missing for 37.3% of the pretested and 36.6% of the
# unpretested participants (29% to 46% by group).", with the sentence for
# losses in the unpretested groups only when the fit has a pretest
# covariate. `what` words repeated posttests.
.attrition_sentences <- function(a, pretest, what = "Posttests were missing") {
  s <- .attrition_shares(a)
  if (is.null(s) || !sum(s$missing)) return(character(0))
  by_group <- .pct(s$range, 0L)
  c(
    sprintf("%s for %s of the pretested and %s of the unpretested participants (%s).",
            what, .pct(s$pretested), .pct(s$unpretested),
            if (by_group[1] == by_group[2]) {
              paste(by_group[1], "in every group")
            } else {
              paste(by_group[1], "to", by_group[2], "by group")
            }),
    if (pretest && .unpretested_only(a)) {
      paste("Only the unpretested groups lost posttests, so the two pretest centers coincide;",
            "those losses cannot be checked against the pretest.")
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
.marginal_center_sentences <- function(x, digits) {
  s <- x$pretest_sensitivity
  if (is.null(x$pretest_center) || is.null(s) || !isTRUE(s$differ)) return(character(0))
  n <- x$standardized
  used <- x$pretest_center
  other <- .other_center(used)
  set <- function(center) {
    phrase <- .marginal_set_phrase(n, center)
    sub("^every pretested participant with a pretest", "every pretested participant with a pretest score",
        phrase)
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
    if (length(dis)) {
      scale <- sub(":.*$", "", dis)
      contrast <- sub("^[^:]+: ", "", dis)
      sprintf("Whether %s %s from zero therefore depends on the standardization.",
              .series_and(sprintf("%s on the %s scale", phrase[contrast], tolower(scale))),
              if (length(dis) > 1L) "differ" else "differs")
    },
    s$note
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
