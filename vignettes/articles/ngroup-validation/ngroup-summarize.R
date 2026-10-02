# Summaries for the issue #45 simulation: performance measures with Monte
# Carlo standard errors, the agreement check against the package, and the
# run information. Sourced by ngroup-simulation.R after the run; it can also
# be sourced on its own to rebuild the summaries from the checkpoints.

summarize_run <- function(results, commit, started, workers) {

  pkgload::load_all(".", quiet = TRUE)

  # ---- Performance ----------------------------------------------------------------

  by_scenario <- split(results, vapply(results, `[[`, numeric(1), "scenario"))
  perf <- do.call(rbind, lapply(by_scenario, function(parts) {
    sc <- parts[[1]]$scenario
    s <- scenarios[sc, ]
    k <- s$k
    R <- sum(vapply(parts, `[[`, numeric(1), "n"))
    total <- function(name) Reduce(`+`, lapply(parts, function(p) p$sums[[name]]))
    pairs <- pairs_of(k)
    np <- length(pairs)
    comp <- rep(pair_names(k), length(types))
    contrast <- rep(types, each = np)
    truth <- unlist(lapply(seq_along(types), function(ti) {
      vapply(pairs, function(pr) truth_pair(s, k, pr[1], pr[2])[ti], numeric(1))
    }))
    e <- effects[[s$effects]]
    delta <- e$delta(k); sens <- e$sens(k)
    omni_null <- c(all(delta + sens / 2 == 0), all(sens == 0), all(delta + sens == 0),
                   all(delta == 0))
    row <- function(method, measure, contrast, comparison, estimate, mcse) {
      data.frame(scenario = sc, k = k, n = s$n, rho = s$rho, spread = s$spread,
                 effects = s$effects, method = method, measure = measure,
                 contrast = contrast, comparison = comparison, estimate = estimate,
                 mcse = mcse, reps = R, stringsAsFactors = FALSE)
    }
    rate <- function(x) list(est = x / R, mcse = sqrt(pmax(x / R * (1 - x / R), 0) / R))

    out <- list()
    for (m in c("M1", "M2", "M3")) {
      bias <- total(paste(m, "bias_sum", sep = "|")) / R
      sq <- total(paste(m, "bias_sq", sep = "|")) / R
      emp_se <- sqrt(pmax(sq - bias^2, 0) * R / (R - 1))
      model_se <- total(paste(m, "se_sum", sep = "|")) / R
      cover <- rate(total(paste(m, "cover", sep = "|")))
      reject <- rate(total(paste(m, "reject", sep = "|")))
      out[[length(out) + 1]] <- row(m, "bias", contrast, comp, bias, emp_se / sqrt(R))
      out[[length(out) + 1]] <- row(m, "empirical_se", contrast, comp, emp_se, emp_se / sqrt(2 * (R - 1)))
      out[[length(out) + 1]] <- row(m, "model_se", contrast, comp, model_se, NA_real_)
      out[[length(out) + 1]] <- row(m, "coverage", contrast, comp, cover$est, cover$mcse)
      null_rows <- abs(truth) < 1e-12
      if (any(null_rows)) {
        out[[length(out) + 1]] <- row(m, "type1_unadjusted", contrast[null_rows],
                                      comp[null_rows], reject$est[null_rows], reject$mcse[null_rows])
      }
      omni <- rate(total(paste(m, "omni_reject", sep = "|")))
      omni_tests <- c("Condition (avg over pretest)", "Pretest x Condition",
                      "Condition | pretested", "Condition | unpretested")
      if (any(omni_null)) {
        out[[length(out) + 1]] <- row(m, "type1_omnibus", omni_tests[omni_null], NA,
                                      omni$est[omni_null], omni$mcse[omni_null])
      }
      for (fam in c("fwer_control", "fwer_pairwise", "fwer_unadjusted")) {
        x <- total(paste(m, fam, sep = "|"))
        ok <- !is.na(x)
        if (any(ok)) {
          fr <- rate(x[ok])
          label <- switch(fam, fwer_control = "fwer_holm_control",
                          fwer_pairwise = "fwer_holm_pairwise",
                          fwer_unadjusted = "fwer_unadjusted_control")
          out[[length(out) + 1]] <- row(m, label, types[ok], NA, fr$est, fr$mcse)
        }
      }
    }
    sens_null <- abs(truth[contrast == "Pretest x Treatment"]) < 1e-12
    if (any(sens_null[seq_len(k)])) {
      x <- rate(total("M4|any_control"))
      out[[length(out) + 1]] <- row("M4", "any_interaction_control", "Pretest x Treatment", NA, x$est, x$mcse)
    }
    if (any(sens_null)) {
      x <- rate(total("M4|any_pair"))
      out[[length(out) + 1]] <- row("M4", "any_interaction_pairs", "Pretest x Treatment", NA, x$est, x$mcse)
    }
    for (m5 in c("M5", "M5h")) {
      for (meas in c("e1", "null_effective", "flag", "differ")) {
        x <- rate(total(paste(m5, meas, sep = "|")))
        out[[length(out) + 1]] <- row(m5, meas, NA, NA, x$est, x$mcse)
      }
    }
    do.call(rbind, out)
  }))
  rownames(perf) <- NULL
  utils::write.csv(perf, file.path(out_dir, "performance.csv"), row.names = FALSE)

  # ---- Agreement with the package ------------------------------------------------

  agreement <- do.call(rbind, lapply(results, function(part) {
    if (!length(part$agree)) return(NULL)
    s <- scenarios[part$scenario, ]
    des <- design_of(s$k, s$n)
    lab <- ifelse(des$cond == 0, "Control", paste0("T", des$cond))
    do.call(rbind, lapply(part$agree, function(a) {
      df <- data.frame(y = a$d$y, treat = lab, pretested = des$pretested,
                       pre = ifelse(des$pretested == 1, a$d$pre, NA))
      if (a$method %in% c("M5", "M5h")) {
        return(agree_steyn_row(part$scenario, a, df, s$k))
      }
      y_pre <- if (a$method == "M1") df$pre else NULL
      robust <- if (a$method == "M3") "none" else "HC3"
      fp <- fit_solomon_glm(df$y, df$treat, df$pretested, y_pre, control = "Control",
                            contrasts = "pairwise", robust = robust)
      fc <- fit_solomon_glm(df$y, df$treat, df$pretested, y_pre, control = "Control",
                            robust = robust)
      vc <- rep(seq_len(length(pairs_of(s$k))) <= s$k, length(types))
      data.frame(
        scenario = part$scenario, rep = a$rep, method = a$method,
        max_diff_estimate = max(abs(fp$effects$estimate - a$x$est)),
        max_diff_se = max(abs(fp$effects$std.error - a$x$se)),
        max_diff_p = max(abs(fp$effects$p.value - a$x$p)),
        max_diff_p_holm_control = max(abs(fc$effects$p.adjusted - a$padj_vc[vc])),
        max_diff_p_holm_pairwise = max(abs(fp$effects$p.adjusted - a$padj_all)),
        max_diff_omnibus_p = max(abs(fc$omnibus$p.value - a$x$omni)),
        decisions_identical = NA,
        stringsAsFactors = FALSE
      )
    }))
  }))
  utils::write.csv(agreement, file.path(out_dir, "agreement.csv"), row.names = FALSE)

  # ---- Run information ------------------------------------------------------------

  info <- data.frame(
    item = c("commit", "started", "finished", "replications_per_scenario", "scenarios",
             "workers", "R_version", "platform", "seed", "rng"),
    value = c(commit, format(started, tz = "UTC", usetz = TRUE),
              format(Sys.time(), tz = "UTC", usetz = TRUE), reps, nrow(scenarios),
              workers, R.version.string, R.version$platform, "45045", "L'Ecuyer-CMRG"),
    stringsAsFactors = FALSE
  )
  utils::write.csv(info, file.path(out_dir, "run-information.csv"), row.names = FALSE)

  invisible(list(performance = perf, agreement = agreement, info = info))
}

# Compare the direct implementation of Steyn's sequence with
# fit_solomon_steyn() on one data set: the p-values of E1, E2, E4, and E5,
# and the decisions the performance measures are built from.
agree_steyn_row <- function(scenario, a, df, k) {
  posthoc <- if (a$method == "M5") "scheffe" else "holm"
  fit <- fit_solomon_steyn(df$y, df$treat, df$pretested, control = "Control",
                           posthoc = posthoc)
  tests <- fit$effects$tests
  st <- a$steyn
  short <- function(x) sub(" [(].*$", "", x)

  e1_p <- tests$p.value[tests$step == "E1"]
  e2 <- tests[tests$step == "E2", ]
  sides <- strsplit(e2$groups, " vs ", fixed = TRUE)
  key <- function(x, y) paste(pmin(x, y), pmax(x, y))
  pkg_key <- key(short(vapply(sides, `[`, "", 1)), short(vapply(sides, `[`, "", 2)))
  fast_key <- key(st$pairs$a, st$pairs$b)
  e2_diff <- max(abs(e2$p.adjusted[match(fast_key, pkg_key)] - st$pairs$p))

  e4_p <- tests$p.value[tests$step == "E4"]
  e5 <- tests[tests$step == "E5" & tests$test != "One-way ANOVA, groups combined", ]
  flag <- any(e4_p < alpha)
  differ <- !flag && any(e5$p.adjusted < alpha)
  groups <- fit$effects$groups
  effective <- vapply(seq_len(k), function(j) {
    e1_p < alpha && any(groups$differs_both[groups$condition == paste0("T", j)])
  }, logical(1))

  diffs <- c(abs(e1_p - st$e1_p), e2_diff, abs(e4_p - st$e4_p),
             if (!flag) abs(min(e5$p.adjusted) - st$e5_min_p))
  data.frame(
    scenario = scenario, rep = a$rep, method = a$method,
    max_diff_estimate = NA_real_, max_diff_se = NA_real_,
    max_diff_p = max(diffs),
    max_diff_p_holm_control = NA_real_, max_diff_p_holm_pairwise = NA_real_,
    max_diff_omnibus_p = NA_real_,
    decisions_identical = identical(flag, st$flag) && identical(differ, st$differ) &&
      identical(unname(effective), unname(st$effective)) &&
      identical(e1_p < alpha, st$e1_p < alpha),
    stringsAsFactors = FALSE
  )
}
