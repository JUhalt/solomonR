# APA 7 results text with method-specific references (issue #52).

# ---- Reference registry ---------------------------------------------------------
#
# The references each exported analysis function rests on, as keys of
# .solomon_reference_text. report_solomon() cites these, plus the
# option-specific references of a particular fit. A unit test checks that
# every exported analysis function is registered here.
.solomon_function_refs <- list(
  fit_solomon_glm = character(0),
  fit_solomon_ml = "vanengelenburg1999",
  marginal_solomon = c("localio2007", "daniel2021"),
  perm_solomon = "phipson2010",
  equivalence_solomon = c("schuirmann1987", "lakens2017", "lakens2018"),
  fit_solomon_classic = "waltonbraver1988",
  fit_solomon_1949 = c("solomon1949", "campbell1966"),
  fisher_solomon = c("elkarkri2025b", "gelman2006"),
  stouffer_solomon = c("stouffer1949", "waltonbraver1988"),
  fit_solomon_sem = "rosseel2012",
  fit_solomon_sem_latent = c("rosseel2012", "meredith1993", "vandenberg2000"),
  invariance_solomon = c("rosseel2012", "meredith1993", "vandenberg2000", "satorra2001",
                         "chen2007"),
  solomon_from_summary = "waltonbraver1988",
  solomon_effect_sizes = c("morris2008", "hedges1981"),
  baseline_solomon = c("cumming2001", "kelley2007"),
  fit_solomon_mi = c("carpenter2023", "vanbuuren2018", "cro2019"),
  tipping_point_solomon = c("white2011", "little2012", "carpenter2023"),
  fit_solomon_mmrm = c("mallinckrodt2008", "laird1982", "sabanesbove2026"),
  fit_solomon_steyn = c("steyn2009", "waltonbraver1988", "scheffe1953")
)

# ---- What report_solomon() leaves out ----------------------------------------------
#
# Exported analysis functions whose results report_solomon() does not report,
# each with the reason (issue #111). ?report_solomon and the article "How to
# Cite solomonR and the Methods It Implements" list them, and a unit test
# checks that the result of every other analysis function has a handler in
# .report_handlers.
.report_exclusions <- c(
  stouffer_solomon = paste(
    "It combines the z statistics of tests computed elsewhere. Test I of the",
    "historical sequence, which is that combination, is reported with",
    "`fit_solomon_classic()`."
  ),
  solomon_effect_sizes = paste(
    "It computes effect sizes and their sampling variances to be combined in a",
    "meta-analysis of several studies, not the results of one study."
  )
)

# ---- Formatting helpers -----------------------------------------------------------

.apa_num <- function(x, digits) {
  out <- formatC(x, format = "f", digits = digits)
  sub("^-(0\\.?0*)$", "\\1", out)
}

# A number that cannot exceed 1 in absolute value, without the leading zero
# (APA 7): fit indexes and their changes. Three decimals by default, the
# precision of Chen's (2007) cutoffs for the change in fit.
.apa_unit <- function(x, digits = 3) {
  out <- formatC(x, format = "f", digits = digits)
  out <- sub("^(-?)0\\.", "\\1.", out)
  sub("^-(\\.0+)$", "\\1", out)
}

.apa_p <- function(p, md) {
  sym <- if (md) "*p*" else "p"
  if (is.na(p)) return(paste(sym, "= NA"))
  if (p < .001) return(paste(sym, "< .001"))
  paste(sym, "=", sub("^0", "", sprintf("%.3f", p)))
}

.apa_df <- function(df) {
  if (abs(df - round(df)) < 1e-8) format(round(df)) else sprintf("%.1f", df)
}

.apa_stat <- function(statistic, df, md, name = NULL, digits = 2) {
  if (is.null(name)) name <- if (is.finite(df)) "t" else "z"
  sym <- if (md) paste0("*", name, "*") else name
  if (name == "z" || !is.finite(df)) {
    sprintf("%s = %s", sym, .apa_num(statistic, digits))
  } else {
    sprintf("%s(%s) = %s", sym, .apa_df(df), .apa_num(statistic, digits))
  }
}

# An omnibus test: F(df1, df2) = x, or the chi-square statistic with df1
# degrees of freedom when `reference` is "chisq" (fixed dispersion).
.apa_omnibus <- function(statistic, df1, df2, reference, md, digits = 2) {
  if (identical(reference, "chisq")) {
    # Greek chi and a superscript two, as \u escapes (R code must be ASCII).
    sym <- if (md) "*\u03c7*\u00b2" else "\u03c7\u00b2"
    sprintf("%s(%s) = %s", sym, .apa_df(df1), .apa_num(statistic, digits))
  } else {
    sym <- if (md) "*F*" else "F"
    sprintf("%s(%s, %s) = %s", sym, .apa_df(df1), .apa_df(df2), .apa_num(statistic, digits))
  }
}

# Whole numbers below 10 as words, as APA style asks; larger ones as numerals.
.number_word <- function(n) {
  words <- c("one", "two", "three", "four", "five", "six", "seven", "eight", "nine")
  if (n >= 1 && n <= 9 && n == round(n)) words[n] else format(n)
}

.apa_ci <- function(lo, hi, level, digits) {
  sprintf("%s%% CI [%s, %s]", format(100 * level), .apa_num(lo, digits), .apa_num(hi, digits))
}

.capitalize <- function(x) paste0(toupper(substr(x, 1, 1)), substring(x, 2))

.contrast_phrase <- function(contrast) {
  phrase <- c(
    "ATE (avg over pretest)" = "the average treatment effect across pretest conditions",
    "ATE" = "the average treatment effect across pretest conditions",
    "Pretest x Treatment" = "the Pretest x Treatment interaction (pretest sensitization)",
    "Sens" = "the Pretest x Treatment interaction (pretest sensitization)",
    "Treatment | pretested" = "the treatment effect among pretested participants",
    "Pre_Eff" = "the treatment effect among pretested participants",
    "Treatment | unpretested" = "the treatment effect among unpretested participants",
    "Unpre_Eff" = "the treatment effect among unpretested participants"
  )
  out <- phrase[contrast]
  ifelse(is.na(out), contrast, out)
}

# Who a comparison of a design with several treatments compares. A
# comparison of one condition with another (weights +1 and -1) reads "RP
# relative to Control"; any other set of weights is named, as in "the
# comparison Lecture". Without weights, a name of the form "A vs B" is read
# as A relative to B.
.comparison_clause <- function(comparison, weights = NULL) {
  if (!is.null(weights)) {
    plus <- names(weights)[weights == 1]
    minus <- names(weights)[weights == -1]
    if (length(plus) == 1L && length(minus) == 1L && sum(weights != 0) == 2L) {
      return(paste(plus, "relative to", minus))
    }
    return(paste("the comparison", comparison))
  }
  parts <- regmatches(comparison, regexec("^(.+?) vs (.+)$", comparison))[[1]]
  if (length(parts) == 3L) paste(parts[2], "relative to", parts[3]) else paste("the comparison", comparison)
}

# The nonzero weights of a comparison, as "RP = 0.5, GS = 0.5, Control = -1":
# the treatments first, then the control. `conditions` lists the control
# first.
.weights_phrase <- function(weights, conditions) {
  conditions <- as.character(conditions)
  v <- weights[c(conditions[-1], conditions[1])]
  v <- v[v != 0]
  paste(sprintf("%s = %s", names(v), trimws(formatC(v, digits = 3, format = "fg"))),
        collapse = ", ")
}

# .contrast_phrase() for one comparison of a design with several
# treatments. The phrases keep the words "treatment effect", which
# .nonrandom_wording() rewrites for nonrandomized designs.
.ngroup_contrast_phrase <- function(contrast, comparison, weights = NULL) {
  who <- .comparison_clause(comparison, weights)
  of <- if (startsWith(who, "the comparison")) paste("for", who) else paste("of", who)
  phrase <- switch(
    contrast,
    "ATE (avg over pretest)" = sprintf("the average treatment effect %s across pretest conditions", of),
    "Pretest x Treatment" = sprintf("the Pretest x Treatment interaction (pretest sensitization) for %s", who),
    "Treatment | pretested" = sprintf("the treatment effect %s among pretested participants", of),
    "Treatment | unpretested" = sprintf("the treatment effect %s among unpretested participants", of),
    sprintf("%s for %s", contrast, who)
  )
  phrase
}

# One sentence per contrast: estimate, interval, test.
.contrast_sentences <- function(eff, level, digits, md, stat_name = NULL, scale = "") {
  vapply(seq_len(nrow(eff)), function(i) {
    e <- eff[i, ]
    test <- if (!is.null(e$statistic) && is.finite(e$statistic)) {
      df <- if (!is.null(e$df)) e$df else Inf
      paste0(.apa_stat(e$statistic, df, md, stat_name, digits), ", ")
    } else {
      ""
    }
    sprintf("%s was %s%s, %s, %s%s.",
            .capitalize(.contrast_phrase(e$contrast)), .apa_num(e$estimate, digits), scale,
            .apa_ci(e$conf.low, e$conf.high, level, digits), test, .apa_p(e$p.value, md))
  }, "")
}

# Cell sizes in the order pretested treated, pretested control, unpretested
# treated, unpretested control.
#
# With `conditions` (the control first, then the treatments), `treat` holds
# each participant's condition and the 2(k + 1) counts follow the order of
# .solomon_cells(conditions): the pretested treatments, the pretested
# control, the unpretested treatments, and the unpretested control.
.cell_counts <- function(treat, pretested, conditions = NULL) {
  keep <- !is.na(treat) & !is.na(pretested)
  treat <- treat[keep]
  pretested <- pretested[keep]
  if (!is.null(conditions)) {
    cells <- .solomon_cells(conditions)
    treat <- as.character(treat)
    return(vapply(seq_len(nrow(cells)), function(i) {
      sum(treat == cells$treat[i] & pretested == cells$pretested[i])
    }, integer(1)))
  }
  # A treatment with more than two conditions would give zero counts below.
  if (!all(treat %in% c(0, 1))) {
    stop(
      "The four-group cell counts need a 0/1 treatment indicator; `treat` has the values ",
      paste(utils::head(unique(as.character(treat)), 5L), collapse = ", "), ".",
      call. = FALSE
    )
  }
  c(sum(treat == 1 & pretested == 1), sum(treat == 0 & pretested == 1),
    sum(treat == 1 & pretested == 0), sum(treat == 0 & pretested == 0))
}

# ---- Per-class reporting -----------------------------------------------------------

.report_parts <- function(fit, digits, md) {
  cls <- intersect(class(fit), names(.report_handlers))
  if (!length(cls)) {
    stop("report_solomon() does not support objects of class ",
         paste(class(fit), collapse = "/"), ". See ?report_solomon for the analyses it ",
         "reports and the results it leaves out on purpose.", call. = FALSE)
  }
  .report_handlers[[cls[1]]](fit, digits, md)
}

# The phrases that describe a fit_solomon_glm() model (the model, the
# pretest adjustment, the offset, the covariance estimator, and the scale of
# the contrasts) and the references they cite. Shared by the four-group and
# the N-group reports.
.glm_method_phrases <- function(fit) {
  fam <- fit$family$family
  link <- fit$family$link
  pre <- "pre_obs" %in% names(fit$data)
  exposure <- "log_exposure" %in% names(fit$data)
  refs <- character(0)

  model <- if (startsWith(fam, "Negative Binomial")) {
    refs <- c(refs, "venables2002", "cameron2013")
    "a negative-binomial (NB2) regression fitted by maximum likelihood (Venables & Ripley, 2002; Cameron & Trivedi, 2013)"
  } else if (fam == "poisson") {
    refs <- c(refs, "cameron2013")
    "a Poisson regression (Cameron & Trivedi, 2013)"
  } else if (fam == "binomial") {
    "a logistic regression"
  } else {
    "a linear model"
  }
  scale <- if (fam == "binomial" && link == "logit") {
    " Contrasts are on the log-odds scale."
  } else if (link == "log") {
    " Contrasts are log rate ratios."
  } else {
    ""
  }

  covariance <- switch(
    fit$robust,
    HC3 = {
      refs <- c(refs, "mackinnon1985", "long2000")
      "HC3 heteroskedasticity-consistent standard errors (MacKinnon & White, 1985; Long & Ervin, 2000)"
    },
    CR2 = {
      refs <- c(refs, "bell2002", "pustejovsky2018")
      "CR2 cluster-robust standard errors with Satterthwaite degrees of freedom (Bell & McCaffrey, 2002; Pustejovsky & Tipton, 2018)"
    },
    "model-based standard errors"
  )
  adjust <- if (pre) {
    refs <- c(refs, "lin2013")
    ", adjusting for the pretest score among pretested participants (Lin, 2013)"
  } else {
    ""
  }
  if (pre && !link %in% c("identity", "log")) refs <- c(refs, "daniel2021")

  list(
    model = model,
    adjust = adjust,
    exposure = if (exposure) " and the log of exposure as an offset" else "",
    covariance = covariance,
    scale = scale,
    refs = refs
  )
}

# The rows of a fit_solomon_glm() data frame that entered the model.
.glm_rows_used <- function(fit) {
  used <- rep(TRUE, nrow(fit$data))
  if (!is.null(fit$model$na.action)) used[fit$model$na.action] <- FALSE
  used
}

.report_glm <- function(fit, digits, md) {
  p <- .glm_method_phrases(fit)

  method <- sprintf(
    "Posttest outcomes were analyzed with %s containing treatment, pretesting, and their interaction%s%s, with %s.%s",
    p$model, p$adjust, p$exposure, p$covariance, p$scale
  )

  used <- .glm_rows_used(fit)

  list(
    method = method,
    results = .contrast_sentences(fit$effects, fit$conf_level, digits, md),
    table = fit$effects,
    refs = p$refs,
    cells = .cell_counts(fit$data$treat[used], fit$data$pretested[used])
  )
}

# The design of a study with several treatments, for the design statement:
# the group labels in the order of .solomon_cells(conditions) and the
# sentence naming the design. `conditions` lists the control first.
.ngroup_design_parts <- function(conditions) {
  conditions <- as.character(conditions)
  control <- conditions[1]
  treatments <- conditions[-1]
  k <- length(treatments)
  cells <- .solomon_cells(conditions)
  list(
    groups = paste(ifelse(cells$pretested == 1L, "pretested", "unpretested"), cells$treat),
    design_text = sprintf(paste0(
      "The design was a Solomon N-group design, the extension of the four-group design ",
      "(Solomon, 1949) to several treatments (Steyn, 2009): %s treatment%s (%s) and a ",
      "control (%s), each with and without a pretest, giving %s groups"
    ), .number_word(k), if (k == 1L) "" else "s", .series_and(treatments), control,
    .number_word(2L * (k + 1L))),
    refs = "steyn2009"
  )
}

# The multiplicity statement of a design with several treatments: nothing
# for one comparison; otherwise whether, and how, the p-values of each
# contrast were adjusted across the `n_comp` comparisons. Returns the
# sentence, its references, and whether the p-values were adjusted.
.ngroup_multiplicity <- function(n_comp, adjust) {
  if (n_comp == 1L) {
    return(list(text = character(0), refs = character(0), adjusted = FALSE))
  }
  if (identical(adjust, "none")) {
    return(list(
      text = "The p-values and confidence intervals were not adjusted for multiple comparisons.",
      refs = character(0), adjusted = FALSE
    ))
  }
  list(
    text = sprintf(paste0(
      "Within each contrast, the p-values of the %s comparisons were adjusted with %s; ",
      "the confidence intervals were not adjusted."
    ), .number_word(n_comp), switch(
      adjust,
      holm = "Holm's (1979) procedure",
      bonferroni = "the Bonferroni procedure (see Holm, 1979)"
    )),
    refs = "holm1979",
    adjusted = TRUE
  )
}

# One results paragraph per comparison of a design with several treatments:
# each Solomon contrast with its estimate, confidence interval, test, and
# p-value, the adjusted p-value (labeled) when `adjusted`. `weights` holds
# one row of condition weights per comparison, or is NULL when every
# comparison is named "A vs B".
.ngroup_comparison_paragraphs <- function(effects, comparisons, weights, adjusted, adjust,
                                          conf_level, digits, md) {
  label <- if (adjusted) {
    switch(adjust, holm = "Holm-adjusted", bonferroni = "Bonferroni-adjusted")
  }
  vapply(comparisons, function(cmp) {
    rows <- effects[effects$comparison == cmp, , drop = FALSE]
    w <- if (is.null(weights)) NULL else weights[cmp, ]
    sentences <- vapply(seq_len(nrow(rows)), function(i) {
      r <- rows[i, ]
      p_text <- if (adjusted) {
        paste0(.apa_p(r$p.adjusted, md), ", ", label)
      } else {
        .apa_p(r$p.value, md)
      }
      sprintf("%s was %s, %s, %s, %s.",
              .capitalize(.ngroup_contrast_phrase(r$contrast, cmp, w)),
              .apa_num(r$estimate, digits),
              .apa_ci(r$conf.low, r$conf.high, conf_level, digits),
              .apa_stat(r$statistic, r$df, md, NULL, digits), p_text)
    }, "")
    paste(sentences, collapse = " ")
  }, "", USE.NAMES = FALSE)
}

# Report for a design with several treatments (class solomon_ngroup).
.report_ngroup <- function(fit, digits, md) {
  p <- .glm_method_phrases(fit)
  conditions <- fit$conditions$condition[order(fit$conditions$role != "control")]
  control <- conditions[1]
  treatments <- conditions[-1]
  design <- .ngroup_design_parts(conditions)
  refs <- c(p$refs, design$refs)
  w <- fit$weights
  comparisons <- rownames(w)
  n_comp <- length(comparisons)
  om <- fit$omnibus

  # The comparisons: each treatment against the control, every pair of
  # conditions, or weights the user gave.
  known <- function(type) {
    ref <- .ngroup_weights(type, control, treatments)
    identical(dim(ref), dim(w)) && identical(rownames(ref), comparisons) &&
      isTRUE(all.equal(unname(ref), unname(w[, colnames(ref), drop = FALSE])))
  }
  compared <- if (known("control")) {
    sprintf("Each treatment was compared with the control (%s)", .series_and(comparisons))
  } else if (known("pairwise")) {
    sprintf("Every pair of conditions was compared (%s)", .series_and(comparisons))
  } else {
    # Weights the user gave: state them, so that each comparison can be
    # reproduced (the treatments first, then the control; zeros left out).
    weight_text <- vapply(comparisons, function(cmp) .weights_phrase(w[cmp, ], conditions), "")
    sprintf("The comparison%s %s %s defined by weights over the conditions (%s)",
            if (n_comp > 1L) "s" else "", .series_and(comparisons),
            if (n_comp > 1L) "were" else "was",
            if (n_comp > 1L) {
              paste(sprintf("%s: %s", comparisons, weight_text), collapse = "; ")
            } else {
              weight_text
            })
  }

  omnibus <- if (identical(fit$robust, "CR2")) {
    "Omnibus Wald tests, using the small-sample test of Pustejovsky and Tipton (2018) for CR2 covariance,"
  } else if (all(om$reference == "chisq")) {
    "Omnibus Wald chi-square tests"
  } else {
    "Omnibus Wald F tests"
  }

  mult <- .ngroup_multiplicity(n_comp, fit$adjust)
  refs <- c(refs, mult$refs)

  method <- paste(c(
    sprintf(
      "Posttest outcomes of the %s groups were analyzed jointly with %s containing an indicator for each treatment, pretesting, and their interactions%s%s, with %s.%s",
      .number_word(2L * (length(treatments) + 1L)), p$model, p$adjust, p$exposure,
      p$covariance, p$scale
    ),
    sprintf(paste0(
      "%s examined whether the differences between the conditions depended on pretesting ",
      "(the Pretest x Condition interaction) and whether the conditions differed when ",
      "averaged over pretest conditions."
    ), omnibus),
    sprintf("%s, and the Solomon contrasts were estimated for %s.", compared,
            if (n_comp > 1L) "each comparison" else "this comparison"),
    mult$text
  ), collapse = " ")

  # Results: the omnibus tests of the interaction and of the conditions
  # averaged over pretest conditions, then each comparison's contrasts.
  om_sentence <- function(test) {
    r <- om[om$test == test, , drop = FALSE]
    paste0(.apa_omnibus(r$statistic, r$df1, r$df2, r$reference, md, digits), ", ",
           .apa_p(r$p.value, md))
  }
  results <- sprintf(
    "The omnibus test of the Pretest x Condition interaction (pretest sensitization) gave %s, and the omnibus test of the conditions, averaged over pretest conditions, gave %s.",
    om_sentence("Pretest x Condition"), om_sentence("Condition (avg over pretest)")
  )

  results <- c(results, .ngroup_comparison_paragraphs(
    fit$effects, comparisons, w, mult$adjusted, fit$adjust, fit$conf_level, digits, md
  ))

  used <- .glm_rows_used(fit)

  list(
    method = method,
    results = results,
    table = fit$effects,
    refs = refs,
    cells = .cell_counts(fit$data$condition[used], fit$data$pretested[used], conditions),
    groups = design$groups,
    design_text = design$design_text
  )
}

.report_ml <- function(fit, digits, md) {
  satt <- identical(fit$inference, "satterthwaite")
  refs <- .solomon_function_refs$fit_solomon_ml
  if (satt) refs <- c(refs, "satterthwaite1946", "welch1947")
  method <- paste0(
    "The Solomon contrasts were estimated by full-information maximum likelihood, ",
    "treating the absent pretests as structurally missing (van Engelenburg, 1999), with ",
    if (satt) {
      "t tests using Satterthwaite degrees of freedom (Satterthwaite, 1946; Welch, 1947)."
    } else {
      "Wald z tests."
    }
  )
  list(
    method = method,
    results = .contrast_sentences(fit$effects, fit$conf_level, digits, md),
    table = fit$effects,
    refs = refs,
    cells = .cell_counts(fit$data$treat, fit$data$pretested)
  )
}

.report_classic <- function(fit, digits, md) {
  flow <- if (is.null(fit$settings$flow)) "1988" else fit$settings$flow
  refs <- .solomon_function_refs$fit_solomon_classic
  method <- "The historical Tests A-I sequence was followed (Walton Braver & Braver, 1988)"
  if (flow == "1990") {
    refs <- c(refs, "braver1990")
    method <- paste0(method, ", as amended by Braver and Walton Braver (1990)")
  }
  if (flow == "1995") {
    refs <- c(refs, "sawilowsky1996")
    method <- paste0(method, ", as revised to omit Test D (Walton Braver & Braver, ",
                     "1995, as cited in Sawilowsky, 1996)")
  }
  allocation <- fit$settings$alpha_allocation
  if (!is.null(allocation) && allocation != "none") {
    refs <- c(refs, "sawilowsky1996")
    lv <- fit$settings$alpha_levels
    method <- paste0(method, sprintf(
      ", with the test-wise significance levels of Sawilowsky's (1996) Method %s under the %s robustness criterion (Test A at %s; Tests E, H, and I at %s)",
      substr(allocation, 7, 7), sub("^method[12]_", "", allocation),
      sub("^0", "", format(lv[["A"]])), sub("^0", "", format(lv[["S"]]))
    ))
  }
  method <- paste0(method, ".")

  label <- c(
    A = "Pretest x Treatment interaction (Test A)",
    B = "treatment effect among pretested groups (Test B)",
    C = "treatment effect among unpretested groups (Test C)",
    D = "treatment main effect (Test D)",
    E = "ANCOVA treatment effect in the pretested groups (Test E)",
    F = "gain-score treatment effect (Test F)",
    G = "Treatment x Time interaction (Test G)",
    H = "posttest-only treatment effect (Test H)"
  )
  results <- character(0)
  for (letter in intersect(fit$path, names(label))) {
    r <- fit$tests[[letter]]$result
    stat <- if (letter == "H") {
      .apa_stat(r$statistic, r$df, md, "t", digits)
    } else {
      sprintf("%s(1, %s) = %s", if (md) "*F*" else "F", .apa_df(r$df), .apa_num(r$F, digits))
    }
    results <- c(results, sprintf("The %s was %s, %s, %s.", label[[letter]],
                                  .apa_num(r$estimate, digits), stat, .apa_p(r$p.value, md)))
  }
  if ("E" %in% fit$path) refs <- c(refs, "huck1973")
  if ("I" %in% fit$path) {
    ii <- fit$tests$I$result
    refs <- c(refs, "stouffer1949", "sawilowsky1994")
    results <- c(results, sprintf(
      "Test I, the Stouffer combination (Stouffer et al., 1949; Walton Braver & Braver, 1988), gave %s, %s. The experiment-wise Type I error of this sequence exceeds its nominal level (Sawilowsky et al., 1994).",
      .apa_stat(ii$z, Inf, md, "z", digits), .apa_p(ii$p.value, md)
    ))
  }
  results <- c(results, sprintf(
    "The %s sequence ended at Test %s%s", flow, utils::tail(fit$path, 1),
    if (flow == "1990" && "I" %in% fit$path) {
      ", which the 1990 amendment regards as the most definitive test."
    } else {
      "."
    }
  ))

  if (!is.null(fit$history) && nrow(fit$history) > 0L) {
    refs <- c(refs, "campbell1966", "mai2020")
    h <- fit$history
    results <- c(results, sprintf(
      "In the history/maturation check (Campbell & Stanley, 1963/1966; Mai et al., 2020), the unpretested control posttest differed from the treated-group pretest by %s, %s, %s, and from the control-group pretest by %s, %s, %s.",
      .apa_num(h$estimate[1], digits), .apa_stat(h$statistic[1], h$df[1], md, "t", digits), .apa_p(h$p.value[1], md),
      .apa_num(h$estimate[2], digits), .apa_stat(h$statistic[2], h$df[2], md, "t", digits), .apa_p(h$p.value[2], md)
    ))
  }

  frame <- fit$tests$A$model$model
  table <- do.call(rbind, lapply(fit$path[fit$path %in% names(label)], function(l) fit$tests[[l]]$result[, c("test", "estimate", "p.value")]))
  list(method = method, results = results, table = table, refs = refs,
       cells = .cell_counts(frame$treat, frame$pretested))
}

.report_perm <- function(fit, digits, md) {
  cluster <- identical(fit$level, "cluster")
  refs <- character(0)
  statistic <- if (identical(fit$statistic, "difference")) {
    # The difference statistic tests only the sharp null hypothesis (#113).
    refs <- c(refs, "romano1990")
    sprintf(paste0(
      "the difference between groups as the statistic (a test of the sharp null ",
      "hypothesis of no treatment effect %s; Romano, 1990)"
    ), if (cluster) "in any cluster" else "for any participant")
  } else {
    if (cluster) {
      refs <- c(refs, "wu2021")
      "a studentized statistic (Wu & Ding, 2021)"
    } else {
      refs <- c(refs, "diciccio2017", "wu2021")
      "a studentized statistic (DiCiccio & Romano, 2017; Wu & Ding, 2021)"
    }
  }
  allocations <- if (isTRUE(fit$exact)) {
    sprintf("all %s possible allocations", format(fit$n_allocations, big.mark = ","))
  } else {
    refs <- c(refs, "phipson2010")
    sprintf("%s random allocations (Phipson & Smyth, 2010)", format(fit$reps, big.mark = ","))
  }
  method <- if (cluster) {
    refs <- c(refs, "gail1996", "hayes2017", "bennett2002")
    sprintf(
      "A cluster-level randomization test permuted the treatment labels of whole clusters (Gail et al., 1996; Hayes & Moulton, 2017), using covariate-adjusted cluster summaries (Bennett et al., 2002), %s, and %s.",
      statistic, allocations
    )
  } else {
    sprintf("A randomization test permuted treatment labels within pretest conditions, using %s and %s.",
            statistic, allocations)
  }
  results <- sprintf("For %s, the estimate was %s, randomization %s.",
                     .contrast_phrase(fit$contrast), .apa_num(fit$estimate, digits), .apa_p(fit$p_perm, md))
  list(method = method, results = results,
       table = data.frame(contrast = fit$contrast, estimate = fit$estimate, p.value = fit$p_perm),
       refs = refs, cells = NULL)
}

.report_marginal <- function(fit, digits, md) {
  if (identical(fit$method, "cluster_summary")) {
    method <- paste0(
      "Risk differences were estimated from unweighted cluster-level summaries (",
      fit$design, ") and compared with t intervals using separate variances and ",
      "Satterthwaite degrees of freedom (Hayes & Moulton, 2017)."
    )
    e <- fit$effects
    e$statistic <- NA_real_
    results <- paste0("On the risk difference scale: ",
                      paste(.contrast_sentences(e, fit$conf_level, digits, md), collapse = " "))
    return(list(method = method, results = results, table = fit$effects, refs = "hayes2017",
                cells = NULL))
  }
  count <- identical(fit$outcome, "count")
  refs <- .solomon_function_refs$marginal_solomon
  if (count) refs <- c(refs, "cameron2013")
  clustered <- !is.null(fit$effects$df) && any(is.finite(fit$effects$df))
  intervals <- if (identical(fit$method, "bootstrap")) {
    sprintf("percentile intervals from %d cell-stratified bootstrap resamples", fit$R)
  } else if (clustered) {
    refs <- c(refs, "bell2002", "pustejovsky2018")
    paste("delta-method intervals from CR2 cluster-robust standard errors, with t",
          "reference distributions using Satterthwaite degrees of freedom",
          "(Bell & McCaffrey, 2002; Pustejovsky & Tipton, 2018)")
  } else {
    "delta-method intervals"
  }
  method <- sprintf(
    "Marginal (standardized) %s were compared across the Solomon groups (Localio et al., 2007; Daniel et al., 2021)%s, with %s.",
    if (count) "rates" else "risks", if (count) " from the count model (Cameron & Trivedi, 2013)" else "",
    intervals
  )
  results <- unlist(lapply(unique(fit$effects$scale), function(sc) {
    e <- fit$effects[fit$effects$scale == sc, ]
    e$statistic <- NA_real_
    paste0("On the ", tolower(sc), " scale: ",
           paste(.contrast_sentences(e, fit$conf_level, digits, md), collapse = " "))
  }))
  list(method = method, results = results, table = fit$effects, refs = refs, cells = NULL)
}

.report_equivalence <- function(fit, digits, md) {
  # A test of one comparison of a design with several treatments carries
  # `comparison`, its `weights`, and the `conditions`.
  ngroup <- !is.null(fit$comparison)
  phrase <- if (ngroup) {
    .ngroup_contrast_phrase(fit$contrast, fit$comparison, fit$weights)
  } else {
    .contrast_phrase(fit$contrast)
  }
  method <- sprintf(
    "Equivalence of %s was tested with two one-sided tests (Schuirmann, 1987; Lakens, 2017; Lakens et al., 2018) against bounds of [%s, %s], using the standard error of the fitted model.",
    phrase, .apa_num(fit$bounds[1], digits), .apa_num(fit$bounds[2], digits)
  )
  results <- sprintf(
    "The estimate was %s, %s; the larger one-sided %s. %s",
    .apa_num(fit$estimate, digits), .apa_ci(fit$conf.low, fit$conf.high, fit$conf_level, digits),
    .apa_p(fit$p_equivalence, md), fit$interpretation
  )
  table <- data.frame(contrast = fit$contrast, estimate = fit$estimate,
                      p_equivalence = fit$p_equivalence, outcome = fit$outcome)
  out <- list(method = method, results = results, table = table,
              refs = .solomon_function_refs$equivalence_solomon, cells = NULL)
  if (ngroup) {
    conditions <- fit$conditions$condition[order(fit$conditions$role != "control")]
    design <- .ngroup_design_parts(conditions)
    # A comparison that is not one condition against another is named in
    # the text, so its weights are stated.
    defined <- if (startsWith(.comparison_clause(fit$comparison, fit$weights), "the comparison")) {
      sprintf("The comparison %s was defined by weights over the conditions (%s).",
              fit$comparison, .weights_phrase(fit$weights, conditions))
    }
    out$method <- paste(c(
      method, defined, "The test was not adjusted for the other comparisons of the design."
    ), collapse = " ")
    out$table <- cbind(data.frame(comparison = fit$comparison), table)
    out$refs <- c(out$refs, design$refs)
    out$groups <- design$groups
    out$design_text <- design$design_text
  }
  out
}

.report_fisher <- function(fit, digits, md) {
  t <- fit$tests
  method <- paste(
    "Treatment was compared within each pretest condition with Fisher's exact test,",
    "following the historical categorical analysis of El Karkri et al. (2025b)."
  )
  where <- ifelse(t$condition == "Combined", "Across both pretest conditions",
                  paste("Among", tolower(t$condition), "participants"))
  results <- sprintf(
    "%s, the risk was %s with treatment and %s without, Fisher's exact %s.",
    where, .apa_num(t$risk_treatment, digits), .apa_num(t$risk_control, digits),
    vapply(t$fisher_p, .apa_p, "", md = md)
  )
  results <- c(results, paste(
    "A difference in significance between the pretest conditions is not itself a test of",
    "sensitization (Gelman & Stern, 2006)."
  ))
  cells <- c(t$n_treatment[1], t$n_control[1], t$n_treatment[2], t$n_control[2])
  list(method = method, results = results, table = t, refs = .solomon_function_refs$fisher_solomon,
       cells = cells)
}

.report_summary_fit <- function(fit, digits, md) {
  a <- fit$anova
  f <- function(src) {
    r <- a[a$source == src, ]
    sprintf("%s(1, %s) = %s, %s", if (md) "*F*" else "F", .apa_df(fit$df_error),
            .apa_num(r$F, digits), .apa_p(r$p.value, md))
  }
  method <- paste(
    "The published cell statistics were reanalyzed with a two-way analysis of variance of the",
    "posttest (Type III sums of squares), the model behind Tests A-D (Walton Braver & Braver, 1988)."
  )
  results <- c(
    sprintf("The Pretest x Treatment interaction was %s.", f("Treatment x Pretest")),
    sprintf("The treatment main effect was %s, and the pretest main effect %s.", f("Treatment"), f("Pretest")),
    .contrast_sentences(fit$contrasts[fit$contrasts$test %in% c("B", "C"), ], fit$conf_level, digits, md)
  )
  list(method = method, results = results, table = fit$anova,
       refs = .solomon_function_refs$solomon_from_summary, cells = fit$cells$n)
}

# Report for the summary statistics of a design with several treatments
# (class solomon_summary_ngroup, issue #111): the omnibus F tests of the
# two-way analysis of variance, then the Holm-adjusted Solomon contrasts of
# each treatment against the control, worded as in .report_ngroup().
.report_summary_ngroup <- function(fit, digits, md) {
  conditions <- fit$conditions$condition[order(fit$conditions$role != "control")]
  design <- .ngroup_design_parts(conditions)
  comparisons <- unique(fit$contrasts$comparison)
  n_comp <- length(comparisons)
  mult <- .ngroup_multiplicity(n_comp, fit$adjust)
  a <- fit$anova
  omnibus <- function(src) {
    r <- a[a$source == src, , drop = FALSE]
    paste0(.apa_omnibus(r$F, r$df, fit$df_error, "F", md, digits), ", ", .apa_p(r$p.value, md))
  }
  method <- paste(c(
    sprintf(paste0(
      "The published cell statistics of the %s groups were reanalyzed with a two-way analysis ",
      "of variance of the posttest, Condition x Pretest, with Type III sums of squares and a ",
      "pooled error variance, which assumes equal variances in the groups."
    ), .number_word(nrow(fit$cells))),
    paste0(
      "Omnibus F tests examined whether the differences between the conditions depended on ",
      "pretesting (the Pretest x Condition interaction), whether the conditions differed when ",
      "averaged over pretest conditions, and whether pretesting had an effect when averaged ",
      "over the conditions."
    ),
    sprintf(paste0(
      "Each treatment was compared with the control (%s), and the Solomon contrasts were ",
      "estimated for %s with t tests on the pooled error variance."
    ), .series_and(comparisons), if (n_comp > 1L) "each comparison" else "this comparison"),
    mult$text
  ), collapse = " ")
  results <- c(
    sprintf(paste0(
      "The omnibus test of the Pretest x Condition interaction (pretest sensitization) gave %s; ",
      "the omnibus test of the conditions, averaged over pretest conditions, gave %s; and the ",
      "test of pretesting, averaged over the conditions, gave %s."
    ), omnibus("Pretest x Condition"), omnibus("Condition"), omnibus("Pretest")),
    .ngroup_comparison_paragraphs(fit$contrasts, comparisons, NULL, mult$adjusted, fit$adjust,
                                  fit$conf_level, digits, md)
  )
  list(
    method = method,
    results = results,
    table = fit$contrasts,
    refs = c(mult$refs, design$refs),
    cells = fit$cells$n,
    groups = design$groups,
    design_text = design$design_text
  )
}

.report_sem <- function(fit, digits, md) {
  refs <- .solomon_function_refs$fit_solomon_sem
  method <- if (identical(fit$mode, "ancova_pretested")) {
    refs <- c(refs, "huck1973")
    "A two-group structural equation model of the pretested groups, with the pretest as a covariate (Huck & Sandler, 1973; Rosseel, 2012), estimated the treatment effect among pretested participants."
  } else {
    "A four-group mean-structure structural equation model (Rosseel, 2012) estimated the posttest means of the Solomon groups; the model is saturated, so global fit is not reported."
  }
  eff <- as.data.frame(fit$effects)
  list(method = method,
       results = .contrast_sentences(eff, fit$conf_level, digits, md, "z"),
       table = eff, refs = refs, cells = NULL)
}

# ---- Latent SEM and measurement invariance (issue #112) ---------------------------
#
# What these reports give follows the reporting standards for structural
# equation models (Appelbaum et al., 2018, Table 7) and for measurement
# invariance (Putnick & Bornstein, 2016): the software and its version, the
# estimator, the handling of missing data, the group sizes, the
# identification of the latent scale, the global fit of each model, the
# difference tests between models with their criteria, and the parameters
# freed.

# The four Solomon groups as lavaan labels them in fit_solomon_sem_latent()
# and invariance_solomon().
.sem_group_names <- c(P1 = "pretested treated", P0 = "pretested control",
                      U1 = "unpretested treated", U0 = "unpretested control")

.need_lavaan <- function() {
  if (!requireNamespace("lavaan", quietly = TRUE)) {
    stop("Package 'lavaan' is needed to report this fit; please install.packages('lavaan').",
         call. = FALSE)
  }
}

# The numbers of participants in a multiple-group lavaan fit, in the order
# P1, P0, U1, U0 of .sem_group_names (the four-group order of the report).
.sem_group_sizes <- function(fit) {
  n <- lavaan::lavInspect(fit, "nobs")
  labels <- lavaan::lavInspect(fit, "group.label")
  as.integer(n[match(names(.sem_group_names), labels)])
}

# The software, the estimator, and its references for the method statement.
# `estimator` is the value given to fit_solomon_sem_latent() or
# invariance_solomon(); `scaled` says whether the fit has a scaled test
# statistic as well as the standard one.
.sem_estimator <- function(estimator, fit) {
  version <- tryCatch(as.character(fit@version), error = function(e) "")
  software <- if (length(version) == 1L && nzchar(version)) {
    sprintf("lavaan (Version %s; Rosseel, 2012)", version)
  } else {
    "lavaan (Rosseel, 2012)"
  }
  tests <- lavaan::lavInspect(fit, "options")$test
  scaled <- any(!tests %in% c("standard", "none"))
  est <- toupper(estimator)
  if (est == "MLR") {
    text <- paste(
      "maximum likelihood with Huber-White robust standard errors and a scaled test",
      "statistic asymptotically equal to the Yuan-Bentler statistic (MLR; Yuan & Bentler, 2000)"
    )
    refs <- "yuan2000"
  } else if (est == "ML" && !scaled) {
    text <- "maximum likelihood (ML)"
    refs <- character(0)
  } else {
    label <- if (scaled) {
      robust <- lavaan::lavInspect(fit, "test")
      robust <- robust[!names(robust) %in% c("standard", "none")]
      tolower(robust[[1]]$label)
    }
    text <- sprintf("the %s estimator%s", estimator,
                    if (scaled) sprintf(", with a scaled test statistic (%s)", label) else "")
    refs <- character(0)
  }
  list(software = software, text = text, refs = refs, scaled = scaled)
}

# One sentence on the global fit of a lavaan model: the chi-square test, CFI,
# RMSEA, and SRMR. With a scaled test the chi-square is the scaled one; the
# CFI and RMSEA, as fit_solomon_sem_latent() stores and prints them, are
# those of the unscaled maximum likelihood chi-square, which is also given.
.sem_fit_sentence <- function(subject, fit, scaled, digits, md) {
  df <- unname(lavaan::fitMeasures(fit, "df"))
  if (isTRUE(df == 0)) {
    return(sprintf("%s has no degrees of freedom, so its global fit cannot be tested.", subject))
  }
  fm <- lavaan::fitMeasures(fit, c("chisq", "df", "pvalue", "cfi", "rmsea", "srmr",
                                   if (scaled) c("chisq.scaled", "df.scaled", "pvalue.scaled")))
  chi <- function(stat, d) .apa_omnibus(stat, d, NA, "chisq", md, digits)
  indexes <- sprintf("CFI = %s, RMSEA = %s, and SRMR = %s", .apa_unit(fm[["cfi"]]),
                     .apa_unit(fm[["rmsea"]]), .apa_unit(fm[["srmr"]]))
  if (scaled) {
    sprintf(
      "%s gave a scaled %s, %s, %s; the CFI and RMSEA are computed from the unscaled maximum likelihood %s.",
      subject, chi(fm[["chisq.scaled"]], fm[["df.scaled"]]), .apa_p(fm[["pvalue.scaled"]], md),
      indexes, chi(fm[["chisq"]], fm[["df"]])
    )
  } else {
    sprintf("%s gave %s, %s, %s.", subject, chi(fm[["chisq"]], fm[["df"]]),
            .apa_p(fm[["pvalue"]], md), indexes)
  }
}

# A freed parameter in lavaan syntax, in words: "y4 ~ 1" is the intercept of
# y4 and "POST =~ y3" the loading of y3.
.freed_parameter <- function(par) {
  par <- trimws(par)
  if (grepl("^[^=~]+~\\s*1$", par)) {
    return(sprintf("the intercept of %s", trimws(sub("~\\s*1$", "", par))))
  }
  if (grepl("=~", par, fixed = TRUE)) {
    return(sprintf("the loading of %s", trimws(sub("^.*=~", "", par))))
  }
  if (grepl("~~", par, fixed = TRUE)) {
    v <- trimws(strsplit(par, "~~", fixed = TRUE)[[1]])
    return(if (length(v) == 2L && v[1] == v[2]) {
      sprintf("the residual variance of %s", v[1])
    } else {
      sprintf("the residual covariance of %s", paste(v, collapse = " and "))
    })
  }
  sprintf("the parameter %s", par)
}

# ", except the intercept of y4, which was estimated separately in each
# group (partial invariance; Byrne et al., 1989)", or "" with nothing freed.
.freed_clause <- function(freed) {
  if (!length(freed)) return("")
  sprintf(", except %s, which %s estimated separately in each group (partial invariance; Byrne et al., 1989)",
          .series_and(vapply(freed, .freed_parameter, "")),
          if (length(freed) > 1L) "were" else "was")
}

# Each invariance criterion's decision: "supported scalar invariance",
# "supported only metric invariance", or no decision when a statistic was
# unavailable. `target` is "scalar" or "partial scalar".
.invariance_verdicts <- function(supported, target) {
  vapply(supported, function(level) {
    if (level == "undetermined") {
      "reached no decision, because a statistic could not be computed"
    } else if (level == target) {
      sprintf("supported %s invariance", target)
    } else {
      sprintf("supported only %s invariance", level)
    }
  }, "")
}

.chisq_criterion <- function(scaled) {
  if (scaled) {
    "the scaled chi-square difference test (Satorra & Bentler, 2001)"
  } else {
    "the chi-square difference test"
  }
}

# The method sentence on the invariance check of fit_solomon_sem_latent(),
# with its references. `n_items` is the number of posttest indicators; the
# check needs at least three.
.latent_check_sentence <- function(fit, n_items) {
  inv <- fit$invariance
  if (is.null(inv)) {
    text <- if (n_items < 3L) {
      "Measurement invariance was not tested, because the test needs at least three indicators."
    } else {
      "Measurement invariance was not tested in this analysis."
    }
    return(list(text = text, refs = character(0)))
  }
  target <- if (is.null(inv$partial)) "scalar" else "partial scalar"
  verdict <- .invariance_verdicts(inv$supported, target)
  support <- inv$supported == target
  tested <- sprintf(
    "Measurement invariance was tested by comparing configural, metric, and scalar models in sequence%s (Vandenberg & Lance, 2000)",
    if (is.null(inv$partial)) "" else ", with the same parameters freed"
  )
  criteria <- if (identical(verdict[["chisq"]], verdict[["chen2007"]])) {
    sprintf("%s and the change in fit (Chen, 2007) both %s", .chisq_criterion(isTRUE(inv$scaled)),
            verdict[["chisq"]])
  } else {
    sprintf("%s %s, and the change in fit (Chen, 2007) %s", .chisq_criterion(isTRUE(inv$scaled)),
            verdict[["chisq"]], verdict[["chen2007"]])
  }
  text <- if (all(support)) {
    sprintf("%s; %s.", tested, criteria)
  } else {
    sprintf("%s: %s, so %s the %s invariance that the latent mean contrasts assume.",
            tested, criteria,
            if (any(support)) "only one criterion supported" else "neither criterion supported",
            target)
  }
  list(text = text, refs = c("chen2007", if (isTRUE(inv$scaled)) "satorra2001"))
}

# Report for fit_solomon_sem_latent() (issue #112). The invariance sentence
# follows the fit's own constraints and check rather than asserting scalar
# invariance.
.report_sem_latent <- function(fit, digits, md) {
  .need_lavaan()
  s <- fit$settings
  level <- if (is.null(s$conf_level)) 0.95 else s$conf_level
  fp <- fit$fit_post
  estimator <- if (is.null(s$estimator)) lavaan::lavInspect(fp, "options")$estimator else s$estimator
  est <- .sem_estimator(estimator, fp)
  refs <- c(.solomon_function_refs$fit_solomon_sem_latent, est$refs)
  pt <- lavaan::parTable(fp)
  items <- unique(pt$rhs[pt$op == "=~" & pt$lhs == "POST"])

  # The latent scale: a latent variance fixed at 1 (std_lv = TRUE fixes it
  # in the first group only once the loadings are equal) or a loading fixed
  # at 1.
  labels <- lavaan::lavInspect(fp, "group.label")
  var_fixed <- pt$op == "~~" & pt$lhs == "POST" & pt$rhs == "POST" & pt$free == 0L
  scale <- if (any(var_fixed)) {
    g <- unique(pt$group[var_fixed])
    if (length(g) == length(labels)) {
      "and the latent variance at 1 in every group"
    } else {
      where <- .sem_group_names[labels[g]]
      sprintf("and the latent variance at 1 in the %s group%s; the contrasts are therefore in standard deviations of the latent posttest in %s",
              .series_and(where), if (length(g) > 1L) "s" else "",
              if (length(g) > 1L) "those groups" else "that group")
    }
  } else {
    marker <- unique(pt$rhs[pt$op == "=~" & pt$lhs == "POST" & pt$free == 0L])[1]
    sprintf("and the loading of %s at 1, which puts the latent posttest on the scale of %s",
            marker, marker)
  }

  freed_post <- s$partial_post
  if (length(freed_post)) refs <- c(refs, "byrne1989")
  check <- .latent_check_sentence(fit, length(items))
  refs <- c(refs, check$refs)

  method <- c(
    sprintf(paste0(
      "Latent posttest means were compared across the four Solomon groups in a multiple-group ",
      "structural equation model in which %s posttest indicators (%s) measured one latent ",
      "posttest, fitted with %s by %s, using full-information maximum likelihood for missing ",
      "indicator values."
    ), .number_word(length(items)), .series_and(items), est$software, est$text),
    sprintf(paste0(
      "The loadings and intercepts of the indicators were constrained to be equal across the ",
      "four groups, the scalar measurement invariance that latent mean comparisons require ",
      "(Meredith, 1993; Vandenberg & Lance, 2000)%s."
    ), .freed_clause(freed_post)),
    sprintf(paste0(
      "For identification, the latent posttest mean of the unpretested control group was ",
      "fixed at 0, %s."
    ), scale),
    check$text,
    "The contrasts between latent means were tested with Wald z tests."
  )

  eff <- as.data.frame(fit$effects_post)
  results <- c(
    .sem_fit_sentence("The four-group model", fp, est$scaled, digits, md),
    paste(.contrast_sentences(eff, level, digits, md, "z"), collapse = " ")
  )

  if (!is.null(fit$fit_pre)) {
    refs <- c(refs, "huck1973")
    ptp <- lavaan::parTable(fit$fit_pre)
    pre_items <- unique(ptp$rhs[ptp$op == "=~" & ptp$lhs == "PRE"])
    freed_pre <- c(s$partial_post, s$partial_pre)
    if (length(freed_pre)) refs <- c(refs, "byrne1989")
    method <- c(method, sprintf(paste0(
      "A latent analysis of covariance in the two pretested groups (Huck & Sandler, 1973) ",
      "regressed the latent posttest on a latent pretest measured by %s indicators (%s), with ",
      "a common slope. The loadings and intercepts of the pretest and posttest indicators were ",
      "constrained to be equal in the two groups%s. The latent posttest intercept of the ",
      "pretested control group was fixed at 0, and the invariance of this model was not tested."
    ), .number_word(length(pre_items)), .series_and(pre_items), .freed_clause(freed_pre)))
    e <- fit$effects_pre[fit$effects_pre$contrast == "Pre_Eff", , drop = FALSE][1, ]
    results <- c(results, paste(
      .sem_fit_sentence("The latent analysis of covariance", fit$fit_pre, est$scaled, digits, md),
      sprintf(
        "In it, the treatment effect among pretested participants, adjusted for the latent pretest, was %s, %s, %s, %s.",
        .apa_num(e$estimate, digits), .apa_ci(e$conf.low, e$conf.high, level, digits),
        .apa_stat(e$statistic, Inf, md, "z", digits), .apa_p(e$p.value, md)
      )
    ))
  }

  list(method = paste(method, collapse = " "), results = results, table = eff, refs = refs,
       cells = .sem_group_sizes(fp))
}

# Report for invariance_solomon() (issue #112): the minimal information
# Putnick and Bornstein (2016) propose, namely the group sizes, the handling
# of missing data, the criteria, and each model's fit with the comparisons
# and decisions.
.report_invariance <- function(fit, digits, md) {
  .need_lavaan()
  configural <- fit$fits$configural
  est <- .sem_estimator(fit$estimator, configural)
  scaled_diff <- isTRUE(fit$scaled)
  partial <- fit$partial
  target <- if (is.null(partial)) "scalar" else "partial scalar"
  items <- fit$items
  cut <- fit$cutoffs
  refs <- c(.solomon_function_refs$invariance_solomon, est$refs,
            if (!is.null(partial)) "byrne1989")
  if (!scaled_diff) refs <- setdiff(refs, "satorra2001")

  n_total <- sum(fit$sizes)
  equal <- length(unique(fit$sizes)) == 1L
  setting <- if (n_total > 300 && equal) {
    "for a total N above 300 with equal group sizes"
  } else if (n_total <= 300 && !equal) {
    "for a total N of 300 or less with unequal group sizes"
  } else {
    paste("for a total N of 300 or less with unequal group sizes, the smaller of",
          "Chen's two sets, as Chen gives none for these sample sizes")
  }

  method <- paste(c(
    sprintf(paste0(
      "Measurement invariance (Meredith, 1993) of %s indicators (%s) across the four Solomon ",
      "groups was tested with multiple-group confirmatory factor analyses, fitted with %s by ",
      "%s, using full-information maximum likelihood for missing values."
    ), .number_word(length(items)), .series_and(items), est$software, est$text),
    paste0(
      "Configural, metric (equal loadings), and scalar (equal loadings and intercepts) models ",
      "were fitted in sequence, and each was compared with the one before it (Vandenberg & ",
      "Lance, 2000)."
    ),
    if (!is.null(partial)) {
      sprintf("In the metric and scalar models, %s %s estimated separately in each group (partial invariance; Byrne et al., 1989).",
              .series_and(vapply(partial, .freed_parameter, "")),
              if (length(partial) > 1L) "were" else "was")
    },
    sprintf(paste0(
      "Noninvariance was judged by %s at alpha = %s and by the change in fit: a drop in CFI of ",
      "at least %s together with a rise in RMSEA of at least %s or in SRMR of at least %s ",
      "(loadings) or %s (intercepts), the cutoffs of Chen (2007) %s."
    ), .chisq_criterion(scaled_diff), sub("^0", "", format(fit$alpha)), .apa_unit(cut$cfi),
    .apa_unit(cut$rmsea), .apa_unit(cut$srmr[["metric"]]), .apa_unit(cut$srmr[["scalar"]]),
    setting),
    if (est$scaled) {
      "The chi-square of each model, the CFI, the RMSEA, and the SRMR are those of the unscaled maximum likelihood fit."
    }
  ), collapse = " ")

  m <- fit$models
  model_phrase <- vapply(seq_len(nrow(m)), function(i) {
    if (isTRUE(m$df[i] == 0)) {
      sprintf("the %s model, with no degrees of freedom, fit exactly", m$model[i])
    } else {
      sprintf("the %s model gave %s, CFI = %s, RMSEA = %s, and SRMR = %s", m$model[i],
              .apa_omnibus(m$chisq[i], m$df[i], NA, "chisq", md, digits),
              .apa_unit(m$cfi[i]), .apa_unit(m$rmsea[i]), .apa_unit(m$srmr[i]))
    }
  }, "")
  models <- paste0(.capitalize(paste(model_phrase[-length(model_phrase)], collapse = "; ")),
                   "; and ", model_phrase[length(model_phrase)], ".")

  tt <- fit$tests
  step <- function(i, what) {
    chi <- if (is.na(tt$chisq_diff[i])) {
      sprintf("the %schi-square difference could not be computed", if (scaled_diff) "scaled " else "")
    } else {
      sprintf("%s\u0394%s, %s", if (scaled_diff) "a scaled " else "",
              .apa_omnibus(tt$chisq_diff[i], tt$df_diff[i], NA, "chisq", md, digits),
              .apa_p(tt$p.value[i], md))
    }
    sprintf("%s gave %s, \u0394CFI = %s, \u0394RMSEA = %s, and \u0394SRMR = %s", what, chi,
            .apa_unit(tt$delta_cfi[i]), .apa_unit(tt$delta_rmsea[i]), .apa_unit(tt$delta_srmr[i]))
  }
  steps <- paste0(
    step(1L, "Constraining the loadings to be equal (metric against configural)"), "; ",
    step(2L, "constraining the intercepts as well (scalar against metric)"), "."
  )

  verdict <- .invariance_verdicts(fit$supported, target)
  chisq_test <- if (scaled_diff) "The scaled chi-square difference test" else "The chi-square difference test"
  decisions <- if (identical(verdict[["chisq"]], verdict[["chen2007"]])) {
    sprintf("%s and the change in fit both %s.", chisq_test, verdict[["chisq"]])
  } else {
    sprintf("%s %s, and the change in fit %s.", chisq_test, verdict[["chisq"]],
            verdict[["chen2007"]])
  }

  na_row <- tt[1L, , drop = FALSE]
  na_row[] <- NA
  table <- cbind(m, compared_with = c(NA, "configural", "metric"),
                 rbind(na_row, tt)[, setdiff(names(tt), "comparison")])
  rownames(table) <- NULL

  list(method = method, results = c(models, steps, decisions), table = table, refs = refs,
       cells = .sem_group_sizes(configural))
}

# Baseline report for a design with several treatments: the pretests of
# each treatment's pretested group against those of the pretested control
# group (baseline_solomon() with `control`).
.report_baseline_ngroup <- function(fit, digits, md) {
  conditions <- fit$conditions$condition[order(fit$conditions$role != "control")]
  design <- .ngroup_design_parts(conditions)
  g <- fit$groups
  cmp <- fit[["comparisons"]]
  m_sym <- if (md) "*M*" else "M"
  sd_sym <- if (md) "*SD*" else "SD"
  method <- paste(
    "Baseline equivalence of the pretested arms was examined by comparing the pretest",
    "scores of each treatment's pretested group with those of the pretested control group,",
    "with Hedges's g and a noncentral-t interval (Cumming & Finch, 2001; Kelley, 2007).",
    "Each comparison used only the two groups it compared, and the p-values were not",
    "adjusted for the number of comparisons."
  )
  scores <- sprintf(
    "the %s group scored %s = %s (%s = %s)", sub("^Pretested, ", "", g$group),
    m_sym, .apa_num(g$mean, digits), sd_sym, .apa_num(g$sd, digits)
  )
  differences <- vapply(seq_len(nrow(cmp)), function(i) {
    r <- cmp[i, ]
    sprintf(
      "The difference of %s was %s, %s, %s, %s, g = %s, %s.",
      .comparison_clause(r$comparison), .apa_num(r$difference, digits),
      .apa_ci(r$conf.low, r$conf.high, fit$conf_level, digits),
      .apa_stat(r$statistic, r$df, md, "t", digits), .apa_p(r$p.value, md),
      .apa_num(r[["g"]], digits), .apa_ci(r$g.low, r$g.high, fit$conf_level, digits)
    )
  }, "")
  list(
    method = method,
    results = c(
      sprintf("At pretest, %s.", .series_and(scores)),
      differences,
      "The unpretested arms have no pretest, so their baseline could not be checked."
    ),
    table = cmp,
    refs = c(.solomon_function_refs$baseline_solomon, design$refs),
    cells = NULL,
    groups = design$groups,
    design_text = design$design_text
  )
}

.report_baseline <- function(fit, digits, md) {
  # A design with several treatments has one comparison for each treatment.
  # `[[` matches the name exactly; `$` would match it partially.
  if (!is.null(fit[["comparisons"]])) {
    return(.report_baseline_ngroup(fit, digits, md))
  }
  g <- fit$groups
  method <- paste(
    "Baseline equivalence of the pretested arms was examined by comparing their pretest",
    "scores, with Hedges's g and a noncentral-t interval (Cumming & Finch, 2001; Kelley, 2007)."
  )
  results <- c(
    sprintf(
      "At pretest, the treated group scored %s = %s (%s = %s) and the control group %s = %s (%s = %s), a difference of %s, %s, %s, %s, g = %s, %s.",
      if (md) "*M*" else "M", .apa_num(g$mean[1], digits), if (md) "*SD*" else "SD", .apa_num(g$sd[1], digits),
      if (md) "*M*" else "M", .apa_num(g$mean[2], digits), if (md) "*SD*" else "SD", .apa_num(g$sd[2], digits),
      .apa_num(fit$difference, digits), .apa_ci(fit$conf.low, fit$conf.high, fit$conf_level, digits),
      .apa_stat(fit$statistic, fit$df, md, "t", digits), .apa_p(fit$p.value, md),
      .apa_num(fit$g, digits), .apa_ci(fit$g.low, fit$g.high, fit$conf_level, digits)
    ),
    "The unpretested arms have no pretest, so their baseline could not be checked."
  )
  list(method = method, results = results, table = fit$groups,
       refs = .solomon_function_refs$baseline_solomon, cells = NULL)
}

# "a", "a and b", or "a, b, and c".
.series_and <- function(x) {
  if (length(x) <= 1L) return(paste(x))
  if (length(x) == 2L) return(paste(x, collapse = " and "))
  paste0(paste(x[-length(x)], collapse = ", "), ", and ", x[length(x)])
}

.mi_offset_phrase <- function(delta, digits) {
  groups <- .solomon_group_labels[delta != 0]
  values <- .apa_num(delta[delta != 0], digits)
  if (length(unique(values)) == 1L) {
    sprintf("%s in the %s group%s", values[1L], .series_and(groups),
            if (length(groups) > 1L) "s" else "")
  } else {
    paste(sprintf("%s in the %s group", values, groups), collapse = "; ")
  }
}

.mi_method <- function(x, pretest, robust, m) {
  missing <- sum(x$missing)
  total <- sum(x$n)
  sprintf(
    paste0(
      "Missing posttests (%d of %d, %.1f%%) were multiply imputed (m = %d) with a normal ",
      "linear regression fitted separately in each Solomon group%s (Carpenter et al., 2023). ",
      "Each completed data set was analyzed with a linear model containing treatment, pretesting, ",
      "and their interaction%s, with %s, and the estimates were combined with Rubin's rules, ",
      "using the small-sample degrees of freedom of Barnard and Rubin (1999, as cited in van Buuren, 2018)."
    ),
    missing, total, 100 * missing / total, m,
    if (pretest) ", on the pretest in the pretested groups" else "",
    if (pretest) ", adjusting for the pretest score among pretested participants" else "",
    if (robust == "HC3") "HC3 heteroskedasticity-consistent standard errors" else "model-based standard errors"
  )
}

.report_mi <- function(fit, digits, md) {
  refs <- .solomon_function_refs$fit_solomon_mi
  method <- .mi_method(fit$missing, fit$pretest, fit$robust, fit$m)
  method <- paste(method, if (all(fit$delta == 0)) {
    "The imputations assumed the posttests were missing at random."
  } else {
    refs <- c(refs, "white2011")
    sprintf(paste0(
      "As a sensitivity analysis, the imputed posttests were shifted by %s, a delta-adjusted ",
      "pattern-mixture analysis (Carpenter et al., 2023; White et al., 2011)."
    ), .mi_offset_phrase(fit$delta, digits))
  })
  if (fit$robust == "HC3") refs <- c(refs, "mackinnon1985", "long2000")
  if (fit$pretest) refs <- c(refs, "lin2013")
  list(
    method = method,
    results = .contrast_sentences(fit$effects, fit$conf_level, digits, md),
    table = fit$effects,
    refs = refs,
    cells = .cell_counts(fit$data$treat, fit$data$pretested)
  )
}

.report_tipping <- function(fit, digits, md) {
  refs <- c(.solomon_function_refs$tipping_point_solomon, "vanbuuren2018")
  if (fit$robust == "HC3") refs <- c(refs, "mackinnon1985", "long2000")
  r <- fit$results
  pretest <- isTRUE(fit$pretest)
  method <- paste(
    .mi_method(fit$missing, pretest, fit$robust, fit$m),
    sprintf(paste0(
      "In a tipping-point sensitivity analysis (White et al., 2011; Little et al., 2012), ",
      "the imputed posttests of the %s group%s were shifted by offsets from %s to %s ",
      "(%s to %s pooled within-group standard deviations of the observed posttests)."
    ), .series_and(fit$groups), if (length(fit$groups) > 1L) "s" else "",
    .apa_num(min(r$delta), digits), .apa_num(max(r$delta), digits),
    .apa_num(min(r$delta_sd), 2), .apa_num(max(r$delta_sd), 2))
  )
  base <- r[r$delta == 0, ]
  level <- 1 - fit$alpha
  results <- sprintf(
    "Assuming the posttests were missing at random, %s was %s, %s, %s.",
    .contrast_phrase(fit$contrast), .apa_num(base$estimate, digits),
    .apa_ci(base$conf.low, base$conf.high, level, digits), .apa_p(base$p.value, md)
  )
  sig <- if (md) "*p*" else "p"
  for (side in c("negative", "positive")) {
    tp <- fit$tipping[[side]]
    results <- c(results, if (is.na(tp)) {
      sprintf("Its statistical significance at alpha = %s did not change for any %s offset tried.",
              sub("^0", "", format(fit$alpha)), side)
    } else {
      row <- r[r$delta == tp, ]
      sprintf("Its statistical significance at alpha = %s changed at an offset of %s (%s SD), where the estimate was %s, %s.",
              sub("^0", "", format(fit$alpha)), .apa_num(tp, digits),
              .apa_num(fit$tipping_sd[[side]], 2), .apa_num(row$estimate, digits),
              .apa_p(row$p.value, md))
    })
  }
  if (pretest) refs <- c(refs, "lin2013")
  list(
    method = method,
    results = results,
    table = r,
    refs = refs,
    cells = .cell_counts(fit$data$treat, fit$data$pretested)
  )
}

.report_mmrm <- function(fit, digits, md) {
  refs <- .solomon_function_refs$fit_solomon_mmrm
  kr <- fit$df_method == "kenward-roger"
  refs <- c(refs, if (kr) "fitzmaurice2011" else "satterthwaite1946")
  if (fit$pretest) refs <- c(refs, "lin2013")
  method <- paste0(
    "Posttest outcomes at ", length(fit$occasions), " occasions were analyzed with a mixed model ",
    "for repeated measures (Mallinckrodt et al., 2008) containing occasion, treatment, pretesting, ",
    "and all their interactions",
    if (fit$pretest) ", adjusting for the pretest score among pretested participants separately at each occasion (Lin, 2013)" else "",
    ", with ", if (fit$covariance == "unstructured") "an unstructured" else paste("a", fit$covariance),
    " within-participant covariance",
    if (fit$grouped) " estimated separately for pretested and unpretested participants" else "",
    " by restricted maximum likelihood (Laird & Ware, 1982). Tests used ",
    if (kr) "Kenward-Roger degrees of freedom (Kenward & Roger, 1997, as cited in Fitzmaurice et al., 2011)"
    else "Satterthwaite (1946) degrees of freedom",
    ". The model was fitted with the mmrm package (Sabanes Bove et al., 2026) and assumes that ",
    "missing posttests are missing at random."
  )
  eff <- fit$effects
  results <- unlist(lapply(unique(eff$occasion), function(o) {
    e <- eff[eff$occasion == o, ]
    s <- .contrast_sentences(e, fit$conf_level, digits, md)
    change <- e$contrast == "Change in Pretest x Treatment"
    s[!change] <- paste0("At occasion ", o, ", ", sub("^(.)", "\\L\\1", s[!change], perl = TRUE))
    s[change] <- sub("^Change in Pretest x Treatment",
                     sprintf("The change in the Pretest x Treatment interaction from occasion %s to occasion %s",
                             fit$occasions[1], fit$occasions[length(fit$occasions)]), s[change])
    s
  }))
  list(method = method, results = results, table = eff, refs = refs,
       cells = .cell_counts(fit$data$treat, fit$data$pretested))
}

# Report for Solomon's (1949) improvement-score analysis (class
# solomon_1949, issue #111): his inferred pretest, the improvements, and the
# interaction I, with the fact that he gave no test for it. The three-group
# design has its own design statement and group labels.
.report_1949 <- function(fit, digits, md) {
  g <- fit$groups
  four <- identical(fit$design, "four-group")
  k <- nrow(g)
  I <- if (md) "*I*" else "I"
  d <- if (md) sprintf("*d*~%d~", seq_len(k)) else sprintf("d%d", seq_len(k))
  formula <- if (four) {
    sprintf("%s = %s - (%s + %s - %s)", I, d[1], d[2], d[3], d[4])
  } else {
    sprintf("%s = %s - (%s + %s)", I, d[1], d[2], d[3])
  }
  inferred <- if (identical(fit$inferred_method, "pooled")) {
    "the mean of the pooled pretested groups (Solomon, 1949, p. 141)"
  } else {
    "the average of the two pretested groups' pretest means, as in Solomon's tables"
  }
  method <- paste(
    sprintf("Solomon's (1949) improvement-score analysis of the %s design was reproduced.", fit$design),
    sprintf(paste0(
      "The pretest mean of the unpretested groups was inferred as %s, and each group's ",
      "improvement was its posttest mean minus its observed or inferred pretest mean."
    ), inferred),
    sprintf("The interaction is %s (Solomon, 1949, p. %d), where %s are the improvements of the experimental group and Control Groups %s.",
            formula, if (four) 147L else 143L, .series_and(d),
            if (four) "I, II, and III" else "I and II")
  )
  results <- c(
    sprintf("The inferred pretest mean was %s. The improvements were %s, and the interaction was %s = %s.",
            .apa_num(fit$inferred_pretest, digits),
            .series_and(sprintf("%s = %s", d, .apa_num(g$change, digits))), I,
            .apa_num(fit$I, digits)),
    if (four) {
      sprintf("The interaction equals the posttest interaction contrast, %s, less the pretest difference between the two pretested groups, %s.",
              .apa_num(fit$posttest_contrast, digits), .apa_num(fit$pretest_difference, digits))
    },
    sprintf(paste0(
      "Solomon (1949) gave no standard error or test for %s and judged the interaction from ",
      "the variability of the scores (p. 144); Campbell and Stanley (1963/1966, p. 25) later ",
      "judged his gain-score analysis unacceptable."
    ), I)
  )
  out <- list(
    method = method,
    results = results,
    table = g,
    refs = .solomon_function_refs$fit_solomon_1949,
    cells = if (anyNA(g$n)) NULL else as.integer(g$n)
  )
  if (!four) {
    out$groups <- c("experimental", "Control I", "Control II")
    out$design_text <- paste(
      "The design was Solomon's (1949) three-group design: an experimental group and",
      "Control Group I, both pretested, and Control Group II, not pretested; the",
      "experimental group and Control Group II were trained"
    )
  }
  out
}

# Reporting language for nonrandomized designs.
.nonrandom_wording <- function(x) {
  # A comparison of a design with several treatments names the conditions it
  # compares (see .ngroup_contrast_phrase()), which may both be treatments.
  x <- gsub("treatment effect (of|for) ", "difference \\1 ", x)
  x <- gsub("average treatment effect", "average treatment-control difference", x)
  x <- gsub("treatment effect", "treatment-control difference", x)
  x <- gsub("treatment main effect", "treatment-control main difference", x)
  x
}

# Reporting functions by result class. Each takes (fit, digits, md), where
# `md` is TRUE for markdown output, and returns a list with:
#   method       character: the sentences describing the analysis;
#   results      character: one element per paragraph of results;
#   table        data frame: the estimates reported;
#   refs         character: keys of .solomon_reference_text for every work
#                the method and results cite (report_solomon() adds
#                Solomon, 1949);
#   cells        integer counts analyzed per group, in the order of
#                `groups`, or NULL when the analysis has no group counts;
#   groups       optional character: the group labels for `cells` and for
#                `design$randomized` and `design$measurement`; NULL means
#                the four Solomon groups (pretested treatment, pretested
#                control, unpretested treatment, unpretested control);
#   design_text  optional character: the sentence naming the design,
#                without its final period; NULL means "The design was a
#                Solomon four-group design (Solomon, 1949)". For designs
#                with several treatments, .ngroup_design_parts(conditions)
#                gives `groups`, `design_text`, and its `refs`.
# The solomon_steyn handler is .report_steyn() in R/solomon_steyn.R; it is
# called through a wrapper because that file is collated after this one.
.report_handlers <- list(
  solomon_baseline = .report_baseline,
  solomon_glm = .report_glm,
  solomon_ngroup = .report_ngroup,
  solomon_steyn = function(fit, digits, md) .report_steyn(fit, digits, md),
  solomon_ml = .report_ml,
  solomon_classic = .report_classic,
  solomon_perm = .report_perm,
  solomon_marginal = .report_marginal,
  solomon_equivalence = .report_equivalence,
  solomon_fisher = .report_fisher,
  solomon_summary_fit = .report_summary_fit,
  solomon_summary_ngroup = .report_summary_ngroup,
  solomon_1949 = .report_1949,
  solomon_sem = .report_sem,
  solomon_sem_latent = .report_sem_latent,
  solomon_invariance = .report_invariance,
  solomon_mi = .report_mi,
  solomon_tipping = .report_tipping,
  solomon_mmrm = .report_mmrm
)

# ---- Design reporting ----------------------------------------------------------------

# `groups` and `design_text` come from the report handler (see
# .report_handlers); NULL gives the four-group design.
.report_design <- function(cells, design, groups = NULL, design_text = NULL) {
  four <- is.null(groups)
  if (four) {
    groups <- c("pretested treatment", "pretested control", "unpretested treatment",
                "unpretested control")
  }
  n_groups <- length(groups)
  series <- function(x) paste0(paste(x[-length(x)], collapse = ", "), ", and ", x[length(x)])
  out <- character(0)
  if (!is.null(design_text)) {
    out <- c(out, paste0(design_text, "."))
    if (!is.null(cells)) {
      out <- c(out, sprintf(
        "The numbers of participants analyzed in the %s groups were %s, respectively.",
        series(groups), series(cells)
      ))
    }
  } else if (!is.null(cells)) {
    out <- c(out, sprintf(
      "The design was a Solomon four-group design (Solomon, 1949), with %s participants analyzed in the %s groups, respectively.",
      series(cells), series(groups)
    ))
  } else {
    out <- c(out, "The design was a Solomon four-group design (Solomon, 1949).")
  }
  if (!is.null(design$randomized)) {
    randomized <- design$randomized
    if (!is.numeric(randomized) || length(randomized) != n_groups) {
      stop(sprintf(
        "`design$randomized` must give the number assigned to each of the %s groups%s.",
        .number_word(n_groups), if (four) "" else paste0(", in the order ", series(groups))
      ), call. = FALSE)
    }
    if (is.null(cells)) {
      out <- c(out, sprintf("The numbers assigned were %s.", paste(randomized, collapse = ", ")))
    } else {
      lost <- 100 * (randomized - cells) / randomized
      out <- c(out, sprintf(
        "Of %s participants assigned to these groups, %s were analyzed (attrition of %s, respectively).",
        series(randomized), series(cells), series(sprintf("%.1f%%", lost))
      ))
    }
  }
  if (identical(design$assignment, "random")) {
    out <- c(out, sprintf("Participants were randomly assigned to the %s groups.",
                          .number_word(n_groups)))
  } else if (identical(design$assignment, "nonrandom")) {
    out <- c(out, paste(
      "The groups were not formed by random assignment, so the contrasts below are",
      "differences between groups rather than treatment effects. Selection bias, the",
      "largest threat to internal validity in quasi-experimental research, and",
      "instrumentation are the threats most relevant to a nonrandomized Solomon design",
      "(Edmonds & Kennedy, 2017). Baseline differences can be examined only in the",
      "pretested arms; the unpretested arms, whose comparison isolates pretest",
      "sensitization, have no baseline and form a static-group comparison, whose",
      "groups cannot be shown to have been equivalent (Campbell & Stanley, 1963/1966)."
    ))
  }
  if (!is.null(design$plan)) {
    s <- design$plan$settings
    out <- c(out, sprintf(
      "The analysis followed an analysis plan dated %s, whose confirmatory contrast%s %s%s.",
      format(design$plan$date), if (length(s$confirmatory) > 1L) "s were" else " was",
      paste(.contrast_phrase(s$confirmatory), collapse = " and "),
      if (!is.null(s$primary_occasion)) sprintf(" at occasion %s", s$primary_occasion) else ""
    ))
  }
  if (!is.null(design$prespecified)) {
    out <- c(out, if (isTRUE(design$prespecified)) {
      "The analysis of pretest sensitization was pre-specified."
    } else {
      "The analysis of pretest sensitization was not pre-specified and is exploratory."
    })
  }
  if (!is.null(design$measurement)) {
    m <- as.character(design$measurement)
    out <- c(out, if (length(m) == 1L) {
      sprintf("Measurement: %s", m)
    } else if (length(m) == n_groups) {
      sprintf("Measurement by group: %s.", paste(sprintf("%s, %s", groups, m), collapse = "; "))
    } else {
      stop("`design$measurement` must be one description or one per group.", call. = FALSE)
    })
  }
  out
}

# APA order of the reference list: by the first author's surname and
# initials, then the other authors, then year (as in tools/check-references.R).
.apa_sort <- function(refs) {
  key <- vapply(refs, function(entry) {
    authors <- sub(" \\(\\d{4}.*$", "", entry)
    year <- sub("^[^()]+? \\((\\d{4}[a-z]?).*$", "\\1", entry)
    parts <- strsplit(authors, ", ", fixed = TRUE)[[1]]
    first <- paste(parts[seq_len(min(2L, length(parts)))], collapse = " ")
    rest <- paste(parts[-seq_len(min(2L, length(parts)))], collapse = " ")
    letters_only <- function(x) gsub("[^a-z]", "", tolower(gsub("&", "", x, fixed = TRUE)))
    paste(letters_only(first), letters_only(rest), year)
  }, "")
  refs[order(key)]
}

#' APA 7 results text for a Solomon analysis
#'
#' `r lifecycle::badge("stable")`
#' Turns a fitted Solomon analysis into APA 7 results sentences, a short
#' design statement, and the references for exactly the methods that analysis
#' used.
#'
#' @details
#' Supported objects come from [fit_solomon_glm()], [fit_solomon_ml()],
#' [fit_solomon_classic()], [fit_solomon_1949()], [perm_solomon()],
#' [marginal_solomon()], [equivalence_solomon()], [fisher_solomon()],
#' [solomon_from_summary()] (for four-group designs and for designs with
#' several treatments), [fit_solomon_sem()], [fit_solomon_sem_latent()],
#' [invariance_solomon()], [baseline_solomon()], [fit_solomon_mi()],
#' [tipping_point_solomon()], [fit_solomon_mmrm()], and [fit_solomon_steyn()].
#' The references depend
#' on the options the fit used: for example, a CR2 fit cites Bell and
#' McCaffrey (2002) and Pustejovsky and Tipton (2018), and the 1990 flow of
#' the classic analysis adds Braver and Walton Braver (1990). Every reference
#' matches the package's canonical APA 7 bibliography.
#'
#' **What the report leaves out.** Two analysis functions are not reported,
#' on purpose:
#'
#' - [stouffer_solomon()] combines the z statistics of tests computed
#'   elsewhere. Test I of the historical sequence, which is that
#'   combination, is reported with [fit_solomon_classic()].
#' - [solomon_effect_sizes()] computes effect sizes and their sampling
#'   variances to be combined in a meta-analysis of several studies, not the
#'   results of one study.
#'
#' Their help pages give the references for their methods. The design
#' checks ([validate_solomon()], [check_solomon_missing()], and
#' [check_solomon_assumptions()]), planning and power ([plan_solomon()],
#' [power_solomon()], and [analysis_plan_solomon()]), [simulate_solomon()],
#' and [compare_solomon_methods()] do not analyze a study's results and are
#' not reported either. Any other object is refused with an error.
#'
#' The design statement follows the MERIT recommendations on reporting
#' measurement in trials (French et al., 2021b): the numbers analyzed in each
#' group, attrition by group when `design$randomized` is given, whether the
#' sensitization analysis was pre-specified, and the measurement procedure
#' in each group (Recommendation 11 is to use identical measurement protocols
#' in all arms, p. 34).
#'
#' **Designs with several treatments.** For a [fit_solomon_glm()] fit with
#' `control` and three or more conditions (class `solomon_ngroup`), the
#' design statement names the treatments, the control, and the 2(k + 1)
#' groups of the design (Steyn, 2009), with the numbers analyzed in each.
#' The results give the omnibus tests of the Pretest x Condition interaction
#' and of the conditions averaged over pretest conditions, then the Solomon
#' contrasts of each comparison. Their p-values are adjusted within each
#' contrast across the comparisons, by Holm's (1979) procedure unless the
#' fit chose another adjustment; the confidence intervals are not adjusted.
#' Comparisons defined by weights are reported with their weights. An
#' [equivalence_solomon()] test of one comparison and the
#' [baseline_solomon()] comparisons of such a design are reported with the
#' same design statement. So is a [solomon_from_summary()] analysis of the
#' design's cell statistics, whose results give the omnibus F tests of the
#' two-way analysis of variance (including the test of pretesting) and the
#' Holm-adjusted contrasts of each treatment against the control. solomonR
#' follows a pre-publication draft of Steyn (2009), which the author
#' provided; see [fit_solomon_steyn()].
#'
#' **Latent structural equation models.** For [fit_solomon_sem_latent()],
#' the report gives what the reporting standards for structural equation
#' models ask for (Appelbaum et al., 2018, Table 7): the software and its
#' version, the estimator, the handling of missing data, the numbers in each
#' group, how the latent scale was identified (the latent posttest mean of
#' the unpretested control group is fixed at 0), the chi-square test, the
#' CFI, RMSEA, and SRMR, and any parameters freed for partial invariance.
#' The invariance statement follows the fit rather than asserting scalar
#' invariance: it states the constraints the model imposed, names the freed
#' parameters (`partial_post` and `partial_pre`), and gives the result of
#' the invariance check under each criterion, including when the check did
#' not support the invariance the contrasts assume or was not run. With the
#' default MLR estimator, the chi-square is the scaled one, and the CFI and
#' RMSEA are those of the unscaled maximum likelihood chi-square, as the fit
#' prints them. For [invariance_solomon()], the report gives the minimal
#' information Putnick and Bornstein (2016) propose for invariance tests:
#' the group sizes, the handling of missing data, the criteria, and the fit
#' of each model with the comparisons and the decision under each
#' criterion.
#'
#' **Nonrandomized designs.** With `design$assignment = "nonrandom"`, the
#' results describe differences between groups rather than treatment
#' effects, and the design statement names the threats that random
#' assignment would otherwise control: selection bias, the largest threat to
#' internal validity in quasi-experimental research (Edmonds & Kennedy, 2017,
#' p. 7), and instrumentation, which with selection bias Edmonds and Kennedy
#' name as the threats most common in quasi-experimental Solomon designs
#' (p. 94). It adds that baseline differences can be examined only in the
#' pretested arms ([baseline_solomon()]): without random assignment the
#' unpretested arms form a static-group comparison, whose groups cannot be
#' shown to have been equivalent (Campbell & Stanley, 1963/1966, pp. 12,
#' 25).
#'
#' **What the report does not decide.** It states results; it does not
#' interpret them. Whether the analysis was pre-specified must be supplied,
#' never inferred, and the choice of analysis, the reading of the results,
#' and the conclusions remain the researcher's responsibility.
#'
#' @param fit A fitted Solomon analysis (see Details).
#' @param design Optional list describing the design: `randomized`, the
#'   numbers assigned to the four groups (pretested treatment, pretested
#'   control, unpretested treatment, unpretested control); `prespecified`,
#'   `TRUE` or `FALSE` for whether the sensitization analysis was
#'   pre-specified; `plan`, the [analysis_plan_solomon()] result the study
#'   registered, which sets `prespecified` from the plan's confirmatory
#'   contrasts; `measurement`, one description of the measurement
#'   procedure or one per group; and `assignment`, `"random"` or
#'   `"nonrandom"`. For a design with k treatments, `randomized` (and
#'   `measurement`, when given per group) has one entry for each of the
#'   2(k + 1) groups, in this order: the pretested treatments (in the order
#'   of the levels of `treat`), the pretested control, the unpretested
#'   treatments, and the unpretested control. For `mai2020`, that is
#'   pretested RP, pretested GS, pretested Control, unpretested RP,
#'   unpretested GS, and unpretested Control.
#' @param digits Decimal places for estimates and statistics. Default 2.
#' @param format `"text"` (default) or `"markdown"`, which italicizes
#'   statistical symbols.
#'
#' @return An object of class `solomon_report` with `method`, `results`, and
#'   `design` (character vectors of sentences), `table` (the estimates), and
#'   `references` (APA 7 reference entries, in APA order).
#'
#' @references
#' Appelbaum, M., Cooper, H., Kline, R. B., Mayo-Wilson, E., Nezu, A. M., &
#' Rao, S. M. (2018). Journal article reporting standards for quantitative
#' research in psychology: The APA Publications and Communications Board task
#' force report. *American Psychologist, 73*(1), 3–25.
#' https://doi.org/10.1037/amp0000191
#'
#' Campbell, D. T., & Stanley, J. C. (1966). *Experimental and
#' quasi-experimental designs for research*. Rand McNally. (Original work
#' published 1963)
#'
#' Edmonds, W. A., & Kennedy, T. D. (2017). *An applied guide to research
#' designs: Quantitative, qualitative, and mixed methods* (2nd ed.). SAGE
#' Publications. https://doi.org/10.4135/9781071802779
#'
#' French, D. P., Miles, L. M., Elbourne, D., Farmer, A., Gulliford, M.,
#' Locock, L., Sutton, S., McCambridge, J., & MERIT Collaborative Group.
#' (2021b). Reducing bias in trials from reactions to measurement: The MERIT
#' study including developmental work and expert workshop. *Health Technology
#' Assessment, 25*(55), 1–72. https://doi.org/10.3310/hta25550
#'
#' Holm, S. (1979). A simple sequentially rejective multiple test procedure.
#' *Scandinavian Journal of Statistics, 6*(2), 65–70.
#' https://www.jstor.org/stable/4615733
#'
#' Putnick, D. L., & Bornstein, M. H. (2016). Measurement invariance
#' conventions and reporting: The state of the art and future directions for
#' psychological research. *Developmental Review, 41*, 71–90.
#' https://doi.org/10.1016/j.dr.2016.06.004
#'
#' Steyn, R. (2009). Re-designing the Solomon four-group: Can we improve on
#' this exemplary model? *Design Principles and Practices: An International
#' Journal—Annual Review, 3*(1), 383–394. https://doi.org/10.18848/1833-1874/CGP/v03i01/37588
#'
#' @examples
#' fit <- with(solomon_example, fit_solomon_glm(y_post, treat, pretested, y_pre))
#' report_solomon(fit, design = list(prespecified = TRUE))
#'
#' # A six-group design: two treatments and a control (Mai et al., 2020).
#' fit6 <- fit_solomon_glm(post_behavior, condition, pretested, pre_behavior,
#'                         control = "Control", data = mai2020)
#' report_solomon(fit6)
#'
#' @export
report_solomon <- function(fit, design = NULL, digits = 2, format = c("text", "markdown")) {
  format <- match.arg(format)
  if (!is.null(design) && !is.list(design)) {
    stop("`design` must be a list (see ?report_solomon).", call. = FALSE)
  }
  md <- format == "markdown"
  if (!is.null(design$plan)) {
    # Pre-specification comes from the plan itself, never from a guess.
    if (!inherits(design$plan, "solomon_analysis_plan")) {
      stop("`design$plan` must come from analysis_plan_solomon().", call. = FALSE)
    }
    planned <- "Pretest x Treatment" %in% design$plan$settings$confirmatory
    if (!is.null(design$prespecified) && !identical(isTRUE(design$prespecified), planned)) {
      stop("`design$prespecified` contradicts `design$plan`, in which the sensitization ",
           "analysis is ", if (planned) "confirmatory." else "exploratory.", call. = FALSE)
    }
    design$prespecified <- planned
  }
  parts <- .report_parts(fit, digits, md)
  nonrandom <- identical(design$assignment, "nonrandom")
  if (!is.null(design$assignment) && !design$assignment %in% c("random", "nonrandom")) {
    stop("`design$assignment` must be \"random\" or \"nonrandom\".", call. = FALSE)
  }
  if (nonrandom) {
    # Without random assignment, contrasts are differences between groups,
    # not treatment effects.
    parts$results <- .nonrandom_wording(parts$results)
    parts$method <- .nonrandom_wording(parts$method)
    parts$refs <- c(parts$refs, "campbell1966", "edmonds2017")
  }
  keys <- unique(c("solomon1949", parts$refs))
  refs <- .apa_sort(unname(.solomon_reference_text[keys]))
  refs <- if (md) gsub("(https?://[^ ]+)", "<\\1>", refs) else gsub("*", "", refs, fixed = TRUE)
  structure(
    list(
      method = parts$method,
      results = parts$results,
      design = .report_design(parts$cells, design, parts$groups, parts$design_text),
      table = parts$table,
      references = refs,
      format = format
    ),
    class = "solomon_report"
  )
}

#' @export
print.solomon_report <- function(x, ...) {
  wrap <- function(s) cat(strwrap(paste(s, collapse = " "), width = 78), sep = "\n")
  wrap(x$design)
  cat("\n")
  wrap(x$method)
  cat("\n")
  for (r in x$results) {
    wrap(r)
  }
  cat("\nReferences\n\n")
  for (r in x$references) {
    cat(strwrap(r, width = 78, exdent = 4), sep = "\n")
    cat("\n")
  }
  invisible(x)
}
