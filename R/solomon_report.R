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
  solomon_from_summary = "waltonbraver1988",
  solomon_effect_sizes = c("morris2008", "hedges1981"),
  baseline_solomon = c("cumming2001", "kelley2007"),
  fit_solomon_mi = c("carpenter2023", "vanbuuren2018", "cro2019"),
  tipping_point_solomon = c("white2011", "little2012", "carpenter2023"),
  fit_solomon_mmrm = c("mallinckrodt2008", "laird1982", "sabanesbove2026")
)

# ---- Formatting helpers -----------------------------------------------------------

.apa_num <- function(x, digits) {
  out <- formatC(x, format = "f", digits = digits)
  sub("^-(0\\.?0*)$", "\\1", out)
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
.cell_counts <- function(treat, pretested) {
  keep <- !is.na(treat) & !is.na(pretested)
  treat <- treat[keep]
  pretested <- pretested[keep]
  c(sum(treat == 1 & pretested == 1), sum(treat == 0 & pretested == 1),
    sum(treat == 1 & pretested == 0), sum(treat == 0 & pretested == 0))
}

# ---- Per-class reporting -----------------------------------------------------------

.report_parts <- function(fit, digits, md) {
  cls <- intersect(class(fit), names(.report_handlers))
  if (!length(cls)) {
    stop("report_solomon() does not yet support objects of class ",
         paste(class(fit), collapse = "/"), ".", call. = FALSE)
  }
  .report_handlers[[cls[1]]](fit, digits, md)
}

.report_glm <- function(fit, digits, md) {
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

  method <- sprintf(
    "Posttest outcomes were analyzed with %s containing treatment, pretesting, and their interaction%s%s, with %s.%s",
    model, adjust, if (exposure) " and the log of exposure as an offset" else "", covariance, scale
  )

  used <- rep(TRUE, nrow(fit$data))
  if (!is.null(fit$model$na.action)) used[fit$model$na.action] <- FALSE

  list(
    method = method,
    results = .contrast_sentences(fit$effects, fit$conf_level, digits, md),
    table = fit$effects,
    refs = refs,
    cells = .cell_counts(fit$data$treat[used], fit$data$pretested[used])
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
    "the difference between groups as the statistic"
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
  method <- sprintf(
    "Equivalence of %s was tested with two one-sided tests (Schuirmann, 1987; Lakens, 2017; Lakens et al., 2018) against bounds of [%s, %s], using the standard error of the fitted model.",
    .contrast_phrase(fit$contrast), .apa_num(fit$bounds[1], digits), .apa_num(fit$bounds[2], digits)
  )
  results <- sprintf(
    "The estimate was %s, %s; the larger one-sided %s. %s",
    .apa_num(fit$estimate, digits), .apa_ci(fit$conf.low, fit$conf.high, fit$conf_level, digits),
    .apa_p(fit$p_equivalence, md), fit$interpretation
  )
  list(method = method, results = results,
       table = data.frame(contrast = fit$contrast, estimate = fit$estimate,
                          p_equivalence = fit$p_equivalence, outcome = fit$outcome),
       refs = .solomon_function_refs$equivalence_solomon, cells = NULL)
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

.report_sem_latent <- function(fit, digits, md) {
  level <- if (is.null(fit$settings$conf_level)) 0.95 else fit$settings$conf_level
  method <- paste(
    "Latent posttest means were compared across the four Solomon groups in a multiple-group",
    "structural equation model (Rosseel, 2012) with scalar measurement invariance",
    "(Meredith, 1993; Vandenberg & Lance, 2000)."
  )
  eff <- as.data.frame(fit$effects_post)
  list(method = method,
       results = .contrast_sentences(eff, level, digits, md, "z"),
       table = eff, refs = .solomon_function_refs$fit_solomon_sem_latent, cells = NULL)
}

# Reporting functions by result class.
.report_baseline <- function(fit, digits, md) {
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

.mi_offset_phrase <- function(delta, digits) {
  groups <- .solomon_group_labels[delta != 0]
  values <- .apa_num(delta[delta != 0], digits)
  if (length(unique(values)) == 1L) {
    sprintf("%s in the %s group%s", values[1L], paste(groups, collapse = ", "),
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
    ), paste(fit$groups, collapse = ", "), if (length(fit$groups) > 1L) "s" else "",
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

# Reporting language for nonrandomized designs.
.nonrandom_wording <- function(x) {
  x <- gsub("average treatment effect", "average treatment-control difference", x)
  x <- gsub("treatment effect", "treatment-control difference", x)
  x <- gsub("treatment main effect", "treatment-control main difference", x)
  x
}

.report_handlers <- list(
  solomon_baseline = .report_baseline,
  solomon_glm = .report_glm,
  solomon_ml = .report_ml,
  solomon_classic = .report_classic,
  solomon_perm = .report_perm,
  solomon_marginal = .report_marginal,
  solomon_equivalence = .report_equivalence,
  solomon_fisher = .report_fisher,
  solomon_summary_fit = .report_summary_fit,
  solomon_sem = .report_sem,
  solomon_sem_latent = .report_sem_latent,
  solomon_mi = .report_mi,
  solomon_tipping = .report_tipping,
  solomon_mmrm = .report_mmrm
)

# ---- Design reporting ----------------------------------------------------------------

.report_design <- function(cells, design) {
  groups <- c("pretested treatment", "pretested control", "unpretested treatment",
              "unpretested control")
  series <- function(x) paste0(paste(x[-length(x)], collapse = ", "), ", and ", x[length(x)])
  out <- character(0)
  if (!is.null(cells)) {
    out <- c(out, sprintf(
      "The design was a Solomon four-group design (Solomon, 1949), with %s participants analyzed in the %s groups, respectively.",
      series(cells), series(groups)
    ))
  } else {
    out <- c(out, "The design was a Solomon four-group design (Solomon, 1949).")
  }
  if (!is.null(design$randomized)) {
    randomized <- design$randomized
    if (!is.numeric(randomized) || length(randomized) != 4L) {
      stop("`design$randomized` must give the number assigned to each of the four groups.",
           call. = FALSE)
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
    out <- c(out, "Participants were randomly assigned to the four groups.")
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
      "The analysis followed an analysis plan dated %s, whose confirmatory contrast%s %s.",
      format(design$plan$date), if (length(s$confirmatory) > 1L) "s were" else " was",
      paste(.contrast_phrase(s$confirmatory), collapse = " and ")
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
    } else if (length(m) == 4L) {
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
#' [fit_solomon_classic()], [perm_solomon()], [marginal_solomon()],
#' [equivalence_solomon()], [fisher_solomon()], [solomon_from_summary()],
#' [fit_solomon_sem()], [fit_solomon_sem_latent()], [baseline_solomon()],
#' [fit_solomon_mi()], [tipping_point_solomon()], and [fit_solomon_mmrm()].
#' The references depend
#' on the options the fit used: for example, a CR2 fit cites Bell and
#' McCaffrey (2002) and Pustejovsky and Tipton (2018), and the 1990 flow of
#' the classic analysis adds Braver and Walton Braver (1990). Every reference
#' matches the package's canonical APA 7 bibliography.
#'
#' The design statement follows the MERIT recommendations on reporting
#' measurement in trials (French et al., 2021b): the numbers analyzed in each
#' group, attrition by group when `design$randomized` is given, whether the
#' sensitization analysis was pre-specified, and the measurement procedure
#' in each group (Recommendation 11 is to use identical measurement protocols
#' in all arms, p. 34).
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
#'   `"nonrandom"`.
#' @param digits Decimal places for estimates and statistics. Default 2.
#' @param format `"text"` (default) or `"markdown"`, which italicizes
#'   statistical symbols.
#'
#' @return An object of class `solomon_report` with `method`, `results`, and
#'   `design` (character vectors of sentences), `table` (the estimates), and
#'   `references` (APA 7 reference entries, in APA order).
#'
#' @references
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
#' @examples
#' fit <- with(solomon_example, fit_solomon_glm(y_post, treat, pretested, y_pre))
#' report_solomon(fit, design = list(prespecified = TRUE))
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
      design = .report_design(parts$cells, design),
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
