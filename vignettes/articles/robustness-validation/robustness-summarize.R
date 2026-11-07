# Summaries for the issue #116 simulation (robustness of the recommended
# inference): the agreement checks, and, in the study run only, the
# performance measures with their Monte Carlo standard errors (MCSEs), the
# decision rules, the failed fits and warnings, the scenario table, and the
# study's part of run-information.csv.
#
# Sourced by robustness-simulation.R, whose constants and scenario table it
# uses; it defines functions only. It is base R, because the harness does
# not load the package in its summary stage (tools/sim/contract.R).
#
# What may be reported depends on the run's purpose: the protocol's
# "Disclosed runs", and Amendment 2 of #116 (item 9). In code-testing runs
# (purpose "test", which includes the smoke mode of
# tools/check-validation-scripts.R) and in the timing gate,
# summarize_robustness() reports the agreement checks and the run
# information only: no performance measure, no decision, and no count of
# failed fits. To decide whether agreement checks 1, 3, 4, and 5 passed, it
# counts in memory the comparisons each check made, finds check 5's largest
# difference, and sums the warnings muffled in checks 1, 3, and 4, and it
# writes none of these. Before that it checks only that every status is
# one that replicate() writes and that the preflight's checks are present,
# which counts nothing. Outside the study run, agreement.csv is also kept
# from standing in for a count of failed fits: a comparison is skipped when
# a fit failed, so that run's rows are aggregated over the scenarios on a
# fixed plan of checks, methods, and quantities, without the number of
# comparisons and without the warnings counted in checks 1, 3, and 4 (check
# 1 says only whether #122's warning was muffled); the study run gives every
# scenario's row, with these counts. The contract stops if a summary
# outside the study run contains performance, decisions, or failures. The
# unit tests of these functions are in robustness-summarize-tests.R and use
# synthetic inputs only.
#
# Sources (APA 7):
# Bradley, J. V. (1978). Robustness? British Journal of Mathematical and
#   Statistical Psychology, 31(2), 144-152.
#   https://doi.org/10.1111/j.2044-8317.1978.tb00581.x
# Efron, B., & Stein, C. (1981). The jackknife estimate of variance. The
#   Annals of Statistics, 9(3), 586-596. https://doi.org/10.1214/aos/1176345462
# Morris, T. P., White, I. R., & Crowther, M. J. (2019). Using simulation
#   studies to evaluate statistical methods. Statistics in Medicine, 38(11),
#   2074-2102. https://doi.org/10.1002/sim.8086
# Romano, J. P. (1990). On the behavior of randomization tests without a
#   group invariance assumption. Journal of the American Statistical
#   Association, 85(411), 686-692. https://doi.org/10.1080/01621459.1990.10474928
# Schuirmann, D. J. (1987). A comparison of the two one-sided tests procedure
#   and the power approach for assessing the equivalence of average
#   bioavailability. Journal of Pharmacokinetics and Biopharmaceutics, 15(6),
#   657-680. https://doi.org/10.1007/BF01068419
# Siepe, B. S., Bartos, F., Morris, T. P., Boulesteix, A.-L., Heck, D. W., &
#   Pawel, S. (2024). Simulation studies for methodological research in
#   psychology: A standardized template for planning, preregistration, and
#   reporting. Psychological Methods. Advance online publication.
#   https://doi.org/10.1037/met0000695

# ---- Formulas shared with the replications ---------------------------------------

# The two one-sided tests of equivalence_solomon() (its internal .tost();
# Schuirmann, 1987), written out here so that agreement check 3 compares the
# package with an independent copy: H0 effect <= -delta and H0 effect >= delta.
tost_p <- function(estimate, std_error, df, delta) {
  p_lower <- stats::pt((estimate + delta) / std_error, df = df, lower.tail = FALSE)
  p_upper <- stats::pt((estimate - delta) / std_error, df = df, lower.tail = TRUE)
  list(p_lower = p_lower, p_upper = p_upper, equivalent = pmax(p_lower, p_upper) < alpha)
}

# P(X >= x) for X ~ Binomial(n, p0): the p-value of the one-sided exact
# binomial test of rate <= p0 (rules 1 and 2, and the material tier). NA
# when no fit succeeded.
binom_p_upper <- function(x, n, p0) {
  ifelse(n > 0, stats::pbinom(x - 1, n, p0, lower.tail = FALSE), NA_real_)
}

# The p-value of the two-sided z test of bias / MCSE (rule 3).
bias_p <- function(bias, mcse) {
  ifelse(is.finite(bias) & is.finite(mcse) & mcse > 0, 2 * stats::pnorm(-abs(bias / mcse)), NA_real_)
}

# The largest value, NA when there is none.
biggest <- function(v) if (!length(v) || all(is.na(v))) NA_real_ else max(v, na.rm = TRUE)

# ---- The judged cells and C --------------------------------------------------------

# One row per rate that a decision rule checks, from the scenario table
# alone (protocol, "Where each method is judged" and the table of rules).
# Rule 2 checks two rates per judged cell of rule 1 for M1 and M2.
judged_rates <- function(sc) {
  cells <- expand.grid(contrast = contrasts_all, method = c("M1", "M2", "M3", "M4"),
                       scenario = sc$scenario, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  x <- sc[match(cells$scenario, sc$scenario), , drop = FALSE]
  e0 <- x$pattern == "E0"
  main <- x$level == "main"
  skew <- x$level == "opposite skew"
  slopes <- x$level == "slopes"
  h1 <- main & !is.na(x$h) & x$h == 1
  hgt1 <- main & !is.na(x$h) & x$h > 1
  eq <- x$equal_arms
  treat <- cells$contrast %in% contrasts_treat
  m2_five <- cells$contrast %in% c(contrasts_treat, "Pretest main effect")
  pxt <- cells$contrast == "Pretest x Treatment"
  m <- cells$method
  rule1 <- (m == "M1" & e0) |
    (m == "M2" & e0 & (h1 | skew | (hgt1 & eq & m2_five))) |
    (m == "M3" & treat & (e0 | (x$pattern == "ET" & main & pxt))) |
    (m == "M4" & treat & e0 & (h1 | (hgt1 & eq) | skew | (slopes & eq)))
  rule2 <- rule1 & m %in% c("M1", "M2")
  rule3 <- m == "M1" & e0 & (main | skew)
  part <- function(keep, rule, measure) {
    if (!any(keep)) return(NULL)
    data.frame(rule = rule, scenario = cells$scenario[keep], method = m[keep],
               contrast = cells$contrast[keep], measure = measure, stringsAsFactors = FALSE)
  }
  out <- rbind(part(rule1, "1. error control", "rejection"),
               part(rule2, "2. one-sided error (TOST)", "upper_limit_below_truth"),
               part(rule2, "2. one-sided error (TOST)", "lower_limit_above_truth"),
               part(rule3, "3. bias", "bias"))
  rownames(out) <- NULL
  out
}

# The counts the protocol states, for the full table and for the fallback
# that drops the uniform distribution from the main grid. Stops if the
# table gives other counts.
check_family <- function(sc, fallback = FALSE) {
  j <- judged_rates(sc)
  n <- function(rule, method = NULL) {
    sum(startsWith(j$rule, rule) & (if (is.null(method)) TRUE else j$method == method))
  }
  got <- c(M1 = n("1", "M1"), M2 = n("1", "M2"), M3 = n("1", "M3"), M4 = n("1", "M4"),
           rule1 = n("1"), rule2 = n("2"), rule3 = n("3"))
  got <- c(got, C = sum(got[c("rule1", "rule2", "rule3")]))
  # Protocol: M1 2,142; M2 756 + 500; M3 1,224 + 270; M4 360 + 400 + 72 + 40;
  # rule 2 4,284 + 2,512; rule 3 2,016; C = 14,576. Fallback: C = 11,964
  # (4,738 + 5,588 + 1,638).
  want <- if (fallback) {
    c(M1 = 1764, M2 = 1030, M3 = 1224, M4 = 720, rule1 = 4738, rule2 = 5588, rule3 = 1638, C = 11964)
  } else {
    c(M1 = 2142, M2 = 1256, M3 = 1494, M4 = 872, rule1 = 5764, rule2 = 6796, rule3 = 2016, C = 14576)
  }
  if (!identical(as.numeric(got), as.numeric(want[names(got)]))) {
    stop("The judged cells do not give the protocol's counts: ",
         paste(names(got), got, sep = " = ", collapse = ", "), ".", call. = FALSE)
  }
  as.integer(got[["C"]])
}

# The family size of a summary's scenario set: the full table (C = 14,576)
# or the protocol's fallback without the uniform distribution in the main
# grid (C = 11,964). Any other set cannot be judged.
family_size <- function(sc) {
  if (nrow(sc) == nrow(scenarios) && identical(sc$scenario, scenarios$scenario)) {
    return(check_family(sc))
  }
  kept <- scenarios$scenario[study_fallback(scenarios)]
  if (identical(sort(sc$scenario), kept)) return(check_family(sc, fallback = TRUE))
  stop("A study run is judged on the full scenario table or on the protocol's fallback ",
       "(the main grid without the uniform distribution).", call. = FALSE)
}

# ---- Inputs ---------------------------------------------------------------------------

# The statuses replicate() writes (robustness-simulation.R): "ok", a failed
# fit of M1 or M2 ("error", "nonconverged"), or an undefined statistic
# ("undefined", "no_valid_permutations").
status_values <- c("ok", "error", "nonconverged", "undefined", "no_valid_permutations")

# Stops, naming the scenario and the replication, if a status is missing or
# is not one of these, so that an unknown spelling is never taken for a
# failed fit or a success. It counts nothing.
check_statuses <- function(results) {
  bad <- which(!results$status %in% status_values)
  if (length(bad)) {
    i <- bad[[1L]]
    stop("Status ", if (is.na(results$status[i])) "NA" else paste0("\"", results$status[i], "\""),
         " in scenario ", results$scenario[i], ", replication ", results$rep[i], " (", results$method[i],
         ", ", results$contrast[i], ") is not one of ", paste(status_values, collapse = ", "), ".",
         call. = FALSE)
  }
  invisible(TRUE)
}

# The preflight's checks (plan/preflight.csv, read as text): the protocol's
# two stops and agreement check 2 on its 20 fixed data sets, for M3 and M4
# and the four treatment contrasts. Every summary needs them: check 2 is
# reported in agreement.csv, and its failure is a defect.
check_preflight <- function(preflight) {
  if (is.null(preflight) || !all(c("dataset", "method", "contrast", "passed") %in% names(preflight))) {
    stop("The preflight's checks (plan/preflight.csv) are missing, so agreement check 2 cannot be ",
         "reported.", call. = FALSE)
  }
  ds <- as.character(preflight$dataset)
  p <- preflight[!is.na(ds) & nzchar(ds) & ds != "NA", , drop = FALSE]
  want <- expand.grid(contrast = contrasts_treat, method = c("M3", "M4"), dataset = seq_len(check2_sets),
                      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  key <- paste(as.character(p$dataset), p$method, p$contrast, sep = "\r")
  if (nrow(p) != nrow(want) || anyDuplicated(key) ||
      !setequal(key, paste(want$dataset, want$method, want$contrast, sep = "\r"))) {
    stop("The preflight does not hold agreement check 2's ", nrow(want), " rows (", check2_sets,
         " data sets, M3 and M4, four contrasts).", call. = FALSE)
  }
  p
}

# ---- Per-scenario blocks of the results ----------------------------------------------

# The rows of each scenario in `results`, which the contract binds in the
# order of scenario and replication.
scenario_ranges <- function(results) {
  r <- rle(results$scenario)
  if (anyDuplicated(r$values) || is.unsorted(r$values)) {
    stop("The results are not in the order of scenario and replication.", call. = FALSE)
  }
  ends <- cumsum(r$lengths)
  data.frame(scenario = r$values, start = ends - r$lengths + 1L, end = ends)
}

# The results of one scenario as matrices with one row per method and
# contrast (the layout replicate() writes) and one column per replication.
scenario_block <- function(results, from, to, layout) {
  i <- seq.int(from, to)
  k <- nrow(layout)
  if (length(i) %% k != 0L) stop("A scenario's results do not fill whole replications.", call. = FALSE)
  reps <- length(i) %/% k
  if (!identical(results$method[i], rep(layout$method, reps)) ||
      !identical(results$contrast[i], rep(layout$contrast, reps))) {
    stop("A scenario's rows are not in the layout replicate() writes.", call. = FALSE)
  }
  rep_no <- matrix(results$rep[i], k)
  if (any(rep_no != rep_no[rep(1L, k), , drop = FALSE])) {
    stop("A replication's rows are not together.", call. = FALSE)
  }
  m <- function(col) matrix(results[[col]][i], k)
  list(reps = reps, rep = rep_no[1L, ], status = m("status"), estimate = m("estimate"),
       std_error = m("std_error"), df = m("df"), statistic = m("statistic"), p_value = m("p_value"))
}

# Agreement check 5 (protocol): in every E0 replication, M2's seven
# estimates equal M1's within 1e-6, over the replications in which both
# fits of the contrast succeeded. One row per scenario and contrast, with
# the number of replications compared; with none, `passed` is NA (not
# compared).
check5_rows <- function(results, sc) {
  rg <- scenario_ranges(results)
  out <- list()
  for (j in seq_len(nrow(rg))) {
    s <- sc[match(rg$scenario[j], sc$scenario), , drop = FALSE]
    if (!identical(s$pattern, "E0")) next
    b <- scenario_block(results, rg$start[j], rg$end[j], layout_E0)
    i1 <- which(layout_E0$method == "M1")
    i2 <- which(layout_E0$method == "M2")
    both <- b$status[i1, , drop = FALSE] == "ok" & b$status[i2, , drop = FALSE] == "ok"
    d <- abs(b$estimate[i2, , drop = FALSE] - b$estimate[i1, , drop = FALSE])
    d[!both] <- NA_real_
    n <- as.integer(rowSums(both))
    mx <- apply(d, 1L, biggest)
    out[[length(out) + 1L]] <- data.frame(scenario = s$scenario, contrast = contrasts_all,
                                          max_abs_diff = unname(mx), comparisons = n,
                                          passed = ifelse(n > 0L, !is.na(mx) & mx <= tol_check5, NA),
                                          stringsAsFactors = FALSE)
  }
  if (!length(out)) return(NULL)
  do.call(rbind, out)
}

# ---- Agreement summary ---------------------------------------------------------------

# The comparisons that agreement checks 1, 3, and 4 plan, from the protocol
# alone: check 1 in every scenario for M3 and M4, checks 3 and 4 in every E0
# scenario for M1 and M2, each with these quantities and tolerances.
# agreement.csv has a row for each planned group whether or not a fit
# succeeded, so that the presence of a row never depends on a fit.
agreement_plan <- function() {
  data.frame(
    check = rep(c("1. permutation tests", "3. TOST", "4. equivariance"), c(8L, 6L, 6L)),
    method = c(rep(c("M3", "M4"), each = 4L), rep(c("M1", "M2"), each = 3L), rep(c("M1", "M2"), each = 3L)),
    quantity = c(rep(c("p_value", "valid_permutations", "observed_statistic", "permuted_statistics"), 2L),
                 rep(c("equivalent", "p_lower", "p_upper"), 2L), rep(c("estimate_shift", "std_error", "df"), 2L)),
    tolerance = c(rep(c(0, 0, tol_check1, tol_check1), 2L), rep(c(0, tol_check3, tol_check3), 2L),
                  rep(tol_check4[["M1"]], 3L), rep(tol_check4[["M2"]], 3L)),
    e0_only = rep(c(FALSE, TRUE, TRUE), c(8L, 6L, 6L)), stringsAsFactors = FALSE)
}

empty_agreement <- function() {
  data.frame(check = character(), scenario = integer(), method = character(), quantity = character(),
             abs_diff = numeric(), passed = logical(), warnings_muffled = integer(), warnings_122 = integer(),
             warning_classes = character(), stringsAsFactors = FALSE)
}

# Per group g (1 to n_groups) of the agreement rows: the number of
# comparisons, whether all passed (NA without one), the largest absolute
# difference, the warnings and their classes, and the scenarios in which a
# comparison failed.
group_summary <- function(rows, g, n_groups) {
  f <- factor(g, levels = seq_len(n_groups))
  zero <- function(v) { v <- as.integer(v); v[is.na(v)] <- 0L; v }
  text <- function(v) { v <- as.character(v); v[is.na(v)] <- ""; v }
  failed <- !rows$passed
  list(
    n = tabulate(g, nbins = n_groups),
    passed = unname(as.logical(tapply(rows$passed, f, all))),
    max_abs_diff = unname(as.numeric(tapply(rows$abs_diff, f, biggest))),
    warnings = zero(unname(tapply(rows$warnings_muffled, f, sum))),
    warnings_122 = zero(unname(tapply(rows$warnings_122, f, sum))),
    classes = text(unname(tapply(rows$warning_classes, f, function(v) {
      v <- unique(unlist(strsplit(v[nzchar(v)], " ", fixed = TRUE)))
      paste(v, collapse = " ")
    }))),
    failed_in = text(unname(tapply(rows$scenario[failed], f[failed], function(v) {
      paste(sort(unique(v)), collapse = " ")
    }))))
}

# agreement.csv. Check 2 comes from the preflight (one row per method and
# contrast, over its 20 fixed data sets). Checks 1, 3, and 4 come from the
# replications: in the study run one row per scenario, method, and quantity
# of the plan above, with the number of comparisons, the warnings muffled
# in the check, and their classes, and "not compared" where every fit
# failed; outside the study run one row per check, method, and quantity
# over all scenarios, with the scenarios in which a comparison failed, and,
# for check 1, whether #122's warning was muffled (yes or no, never a
# count). Check 5 comes from the stored estimates: per scenario and
# contrast in the study run, per contrast outside it, where its largest
# difference, one of the protocol's descriptive measures, is not shown.
# [Left to the script: the aggregation outside the study run, so that a
# code-testing or gate summary cannot stand in for a count of failed fits,
# which the protocol and Amendment 2 of #116 (item 9) keep out of those
# runs.]
agreement_summary <- function(agreement, results, sc, preflight, study_run) {
  out <- list()
  num <- function(v) suppressWarnings(as.integer(v))
  # Check 2 (ties), from the preflight, which runs before any task.
  p <- check_preflight(preflight)
  p$passed <- as.logical(p$passed)
  p$max_abs_diff <- suppressWarnings(as.numeric(p$max_abs_diff))
  pk <- paste(p$method, p$contrast, sep = "\r")
  for (key in unique(pk)) {
    q <- p[pk == key, , drop = FALSE]
    out[[length(out) + 1L]] <- data.frame(
      check = "2. ties (preflight)", scenario = NA_integer_, method = q$method[[1]],
      contrast = q$contrast[[1]], quantity = "p-value, valid and tie-labelled permutations, statistics",
      max_abs_diff = biggest(q$max_abs_diff), tolerance = tol_check1, passed = all(q$passed),
      comparisons = nrow(q), warnings_muffled = sum(num(q$warnings_muffled)),
      warnings_122 = if ("warnings_122" %in% names(q)) sum(num(q$warnings_122)) else NA_integer_,
      detail = paste0(nrow(q), " fixed data sets; tie-labelled permutations ", sum(num(q$tie_labelled)),
                      ", counted by the package ", sum(num(q$counted_package)),
                      ", by the direct computation ", sum(num(q$counted_direct)),
                      if ("n_tied_stored" %in% names(q)) paste0("; stored n_tied ", sum(num(q$n_tied_stored))) else ""),
      stringsAsFactors = FALSE)
  }
  # Checks 1, 3, and 4, from the replications, on the plan's grid. A row
  # outside the grid is a code fault.
  plan <- agreement_plan()
  rows <- if (is.null(agreement) || !nrow(agreement)) empty_agreement() else agreement
  e0 <- sc$scenario[sc$pattern == "E0"]
  grid <- do.call(rbind, lapply(seq_len(nrow(plan)), function(i) {
    scen <- if (plan$e0_only[i]) e0 else sc$scenario
    if (!length(scen)) return(NULL)
    data.frame(plan = i, check = plan$check[i], scenario = scen, method = plan$method[i],
               quantity = plan$quantity[i], tolerance = plan$tolerance[i], stringsAsFactors = FALSE)
  }))
  grid <- grid[order(match(grid$check, unique(plan$check)), grid$scenario, grid$plan), , drop = FALSE]
  rownames(grid) <- NULL
  g <- match(paste(rows$check, rows$scenario, rows$method, rows$quantity, sep = "\r"),
             paste(grid$check, grid$scenario, grid$method, grid$quantity, sep = "\r"))
  if (anyNA(g)) {
    i <- which(is.na(g))[[1]]
    stop("An agreement row is not one of the planned comparisons: ", rows$check[i], ", scenario ",
         rows$scenario[i], ", ", rows$method[i], ", ", rows$quantity[i], ".", call. = FALSE)
  }
  if (study_run) {
    s <- group_summary(rows, g, nrow(grid))
    out[[length(out) + 1L]] <- data.frame(
      check = grid$check, scenario = grid$scenario, method = grid$method, contrast = "all",
      quantity = grid$quantity, max_abs_diff = s$max_abs_diff, tolerance = grid$tolerance, passed = s$passed,
      comparisons = s$n, warnings_muffled = s$warnings, warnings_122 = s$warnings_122,
      detail = ifelse(s$n == 0L, "not compared",
                      ifelse(nzchar(s$classes), paste("warnings muffled:", s$classes), "")),
      stringsAsFactors = FALSE)
  } else {
    s <- group_summary(rows, grid$plan[g], nrow(plan))
    cm <- paste(plan$check, plan$method, sep = "\r")
    any122 <- tapply(s$warnings_122 > 0L, cm, any)[cm]
    detail <- mapply(function(n, failed_in, check, w) {
      parts <- c(if (n == 0L) "not compared in any scenario",
                 if (nzchar(failed_in)) paste("failed in scenarios", failed_in),
                 if (startsWith(check, "1")) paste("#122 warnings muffled:", if (isTRUE(w)) "yes" else "none"))
      paste(parts, collapse = "; ")
    }, s$n, s$failed_in, plan$check, any122, USE.NAMES = FALSE)
    out[[length(out) + 1L]] <- data.frame(
      check = plan$check, scenario = NA_integer_, method = plan$method, contrast = "all",
      quantity = plan$quantity, max_abs_diff = s$max_abs_diff, tolerance = plan$tolerance, passed = s$passed,
      comparisons = NA_integer_, warnings_muffled = NA_integer_, warnings_122 = NA_integer_,
      detail = detail, stringsAsFactors = FALSE)
  }
  # Check 5, from the stored estimates.
  c5 <- if (is.null(results)) NULL else check5_rows(results, sc)
  if (!is.null(c5)) {
    if (study_run) {
      out[[length(out) + 1L]] <- data.frame(
        check = "5. M2 equals M1", scenario = c5$scenario, method = "M2", contrast = c5$contrast,
        quantity = "estimate", max_abs_diff = c5$max_abs_diff, tolerance = tol_check5, passed = c5$passed,
        comparisons = c5$comparisons, warnings_muffled = 0L, warnings_122 = 0L,
        detail = ifelse(c5$comparisons == 0L, "not compared", ""), stringsAsFactors = FALSE)
    } else {
      f <- factor(c5$contrast, levels = contrasts_all)
      bad <- c5$passed %in% FALSE
      n <- as.integer(tapply(c5$comparisons, f, sum))
      passed <- as.logical(tapply(c5$passed, f, function(v) if (all(is.na(v))) NA else all(v, na.rm = TRUE)))
      failed_in <- as.character(tapply(c5$scenario[bad], f[bad], function(v) paste(sort(unique(v)), collapse = " ")))
      failed_in[is.na(failed_in)] <- ""
      detail <- vapply(seq_along(contrasts_all), function(k) {
        paste(c("largest difference shown in the study run only (a descriptive measure)",
                if (n[k] == 0L) "not compared in any scenario",
                if (nzchar(failed_in[k])) paste("failed in scenarios", failed_in[k])), collapse = "; ")
      }, "")
      out[[length(out) + 1L]] <- data.frame(
        check = "5. M2 equals M1", scenario = NA_integer_, method = "M2", contrast = contrasts_all,
        quantity = "estimate", max_abs_diff = NA_real_, tolerance = tol_check5, passed = passed,
        comparisons = NA_integer_, warnings_muffled = NA_integer_, warnings_122 = NA_integer_,
        detail = detail, stringsAsFactors = FALSE)
    }
  }
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  if (!study_run) {
    out$comparisons <- NULL
    out$warnings_122 <- NULL
  }
  out
}

# ---- Performance measures ------------------------------------------------------------

note_rate <- "proportion at 0 or 1: the binomial formula gives 0"
note_t3 <- "unavailable: the fourth moment of t (3 df) is infinite (protocol)"

# Proportions with their binomial MCSE, sqrt(r(1 - r)/n), n the number of
# successful fits (Morris et al., 2019, Table 6).
rate_rows <- function(hits, n) {
  x <- vapply(hits, function(h) sum(h), 0)
  r <- if (n > 0) x / n else rep(NA_real_, length(x))
  note <- if (n == 0) rep("no successful fit", length(x)) else ifelse(r %in% c(0, 1), note_rate, "")
  data.frame(measure = names(hits), estimate = unname(r), mcse = unname(sqrt(r * (1 - r) / n)),
             mcse_note = note, n_success = as.integer(n), events = as.integer(x), stringsAsFactors = FALSE,
             row.names = NULL)
}

# The leave-one-out jackknife estimate of the standard error of a statistic
# from its n leave-one-out values (Efron & Stein, 1981).
jackknife_se <- function(loo) {
  n <- length(loo)
  sqrt((n - 1) / n * sum((loo - mean(loo))^2))
}

# The leave-one-out empirical SEs, sqrt(SS_(j) / (n - 2)) with
# SS_(j) = SS - n/(n - 1) (x_j - mean)^2, computed in O(n).
loo_empirical_se <- function(est) {
  n <- length(est)
  dev2 <- (est - mean(est))^2
  sqrt(pmax(sum(dev2) - n / (n - 1) * dev2, 0) / (n - 2))
}

# Empirical SE, SD(estimate); model SE, sqrt(mean(SE^2)); and the relative
# error of the model SE, 100 (model SE / empirical SE - 1), with jackknife
# MCSEs over the replications (protocol). The leave-one-out values are
# computed in O(n): the empirical SE by loo_empirical_se(), the model SE as
# sqrt((sum SE^2 - SE_j^2)/(n - 1)).
# [Left to the script: the protocol names jackknife estimates without a
# formula; this is the standard leave-one-out jackknife.] Under t (3 df) the
# MCSEs are unavailable (protocol).
se_rows <- function(est, se, heavy_tailed) {
  n <- length(est)
  emp <- if (n >= 2) stats::sd(est) else NA_real_
  mod <- if (n >= 1) sqrt(mean(se^2)) else NA_real_
  rel <- 100 * (mod / emp - 1)
  mcse <- rep(NA_real_, 3)
  note <- rep("", 3)
  if (heavy_tailed) {
    note[] <- note_t3
  } else if (n < 3) {
    note[] <- "fewer than 3 successful fits"
  } else {
    emp_j <- loo_empirical_se(est)
    s2 <- se^2
    mod_j <- sqrt((sum(s2) - s2) / (n - 1))
    rel_j <- 100 * (mod_j / emp_j - 1)
    mcse <- c(jackknife_se(emp_j), jackknife_se(mod_j), jackknife_se(rel_j))
  }
  if (n < 2) note[c(1, 3)] <- ifelse(nzchar(note[c(1, 3)]), note[c(1, 3)], "fewer than 2 successful fits")
  data.frame(measure = c("empirical_se", "model_se", "relative_error_model_se"),
             estimate = c(emp, mod, rel), mcse = mcse, mcse_note = note, n_success = as.integer(n),
             events = NA_integer_, stringsAsFactors = FALSE)
}

# Bias of M1, mean(estimate) in E0, with MCSE SD(estimate)/sqrt(n)
# (protocol); and the standardized bias, bias / empirical SE.
# [Left to the script: the protocol says that each measure carries an MCSE
# but gives none for the standardized bias. It is the leave-one-out
# jackknife over the replications, as for the measures built on the
# empirical SE, and, like theirs, unavailable under t (3 df).]
bias_rows <- function(est, heavy_tailed) {
  n <- length(est)
  b <- if (n >= 1) mean(est) else NA_real_
  sdv <- if (n >= 2) stats::sd(est) else NA_real_
  mcse_std <- NA_real_
  note_std <- ""
  if (heavy_tailed) {
    note_std <- note_t3
  } else if (n < 3) {
    note_std <- "fewer than 3 successful fits"
  } else {
    mcse_std <- jackknife_se(((sum(est) - est) / (n - 1)) / loo_empirical_se(est))
  }
  data.frame(measure = c("bias", "standardized_bias"), estimate = c(b, b / sdv),
             mcse = c(sdv / sqrt(n), mcse_std),
             mcse_note = c(if (n < 2) "fewer than 2 successful fits" else "", note_std),
             n_success = as.integer(n), events = NA_integer_, stringsAsFactors = FALSE)
}

# Every measure of one fit (M1 or M2) and contrast in an E0 scenario, from
# the stored estimates, standard errors, and degrees of freedom of the
# successful fits. Every estimand is 0 in E0, so the two-sided test rejects
# when 2 pt(-|est/se|, df) < .05 (the package's p-value); the 90% interval's
# upper limit falls below the truth when est/se < -t(.95, df), and its lower
# limit above it when est/se > t(.95, df); the TOST size at the bound uses
# the estimate with Delta added (truth on the upper bound) or subtracted
# (lower bound); power at 0.5 uses the estimate with 0.5 added.
fit_measures <- function(method, est, se, df, heavy_tailed) {
  n <- length(est)
  t <- est / se
  reject <- 2 * stats::pt(-abs(t), df) < alpha
  crit90 <- stats::qt(1 - alpha, df)
  hits <- list(rejection = reject, rejection_below = reject & t < 0, rejection_above = reject & t > 0,
               coverage = !reject, upper_limit_below_truth = t < -crit90,
               lower_limit_above_truth = t > crit90)
  for (d in tost_deltas) {
    hits[[paste0("tost_size_upper_", d)]] <- tost_p(est + d, se, df, d)$equivalent
    hits[[paste0("tost_size_lower_", d)]] <- tost_p(est - d, se, df, d)$equivalent
  }
  hits$power <- 2 * stats::pt(-abs((est + power_effect) / se), df) < alpha
  out <- rate_rows(hits, n)
  if (method == "M1") out <- rbind(out, bias_rows(est, heavy_tailed))
  rbind(out, se_rows(est, se, heavy_tailed))
}

# Rejection rates of a permutation test (M3, M4) and contrast: p < .05, and
# the rejections with an observed statistic below and above 0. In ET the
# treatment contrasts have true value 0.5, so their rate is power.
perm_measures <- function(p, stat, null_true) {
  reject <- p < alpha
  hits <- if (null_true) {
    list(rejection = reject, rejection_below = reject & stat < 0, rejection_above = reject & stat > 0)
  } else {
    list(power = reject)
  }
  rate_rows(hits, length(p))
}

# The role of each performance row: "judged" (a rate or bias a decision rule
# checks), "described" (the protocol's "Described, not judged" cells),
# "descriptive" (the protocol's "Descriptive only" measures), or "reported".
performance_roles <- function(perf, sc) {
  j <- judged_rates(sc)
  key <- function(s, m, c, me) paste(s, m, c, me, sep = "\r")
  judged <- key(perf$scenario, perf$method, perf$contrast, perf$measure) %in%
    key(j$scenario, j$method, j$contrast, j$measure)
  x <- sc[match(perf$scenario, sc$scenario), , drop = FALSE]
  e0 <- x$pattern == "E0"
  hgt1 <- x$level == "main" & !is.na(x$h) & x$h > 1
  slopes <- x$level == "slopes"
  eq <- x$equal_arms
  rate <- perf$measure %in% c("rejection", "upper_limit_below_truth", "lower_limit_above_truth")
  described <- !judged & (
    (perf$method == "M2" & e0 & (hgt1 | slopes) & rate) |
      (perf$method == "M1" & slopes & perf$measure == "bias") |
      (perf$method == "M4" & e0 & (hgt1 | slopes) & !eq & perf$measure == "rejection") |
      (perf$method == "M4" & x$pattern == "ET" & perf$contrast == "Pretest x Treatment" &
         perf$measure == "rejection"))
  descriptive <- perf$measure %in% c("power", "max_abs_difference_from_M1")
  ifelse(judged, "judged", ifelse(described, "described", ifelse(descriptive, "descriptive", "reported")))
}

# Performance measures of every scenario, method, and contrast, plus, per E0
# scenario, the proportion of rejected contrasts of M4 in each replication
# (for the Romano check's MCSE).
performance_table <- function(results, sc, check5) {
  rg <- scenario_ranges(results)
  perf <- vector("list", nrow(rg))
  romano <- vector("list", nrow(rg))
  fails <- vector("list", nrow(rg))
  for (j in seq_len(nrow(rg))) {
    s <- sc[match(rg$scenario[j], sc$scenario), , drop = FALSE]
    e0 <- identical(s$pattern, "E0")
    layout <- if (e0) layout_E0 else layout_ET
    b <- scenario_block(results, rg$start[j], rg$end[j], layout)
    heavy <- identical(s$dist, "t3")
    rows <- vector("list", nrow(layout))
    for (k in seq_len(nrow(layout))) {
      m <- layout$method[k]
      ok <- b$status[k, ] == "ok"
      x <- if (m %in% c("M1", "M2")) {
        fit_measures(m, b$estimate[k, ok], b$std_error[k, ok], b$df[k, ok], heavy)
      } else {
        null_true <- e0 || layout$contrast[k] == "Pretest x Treatment"
        perm_measures(b$p_value[k, ok], b$statistic[k, ok], null_true)
      }
      rows[[k]] <- data.frame(scenario = s$scenario, method = m, contrast = layout$contrast[k], x,
                              stringsAsFactors = FALSE)
    }
    if (e0) {
      c5 <- check5[check5$scenario == s$scenario, , drop = FALSE]
      rows[[length(rows) + 1L]] <- data.frame(
        scenario = s$scenario, method = "M2", contrast = c5$contrast, measure = "max_abs_difference_from_M1",
        estimate = c5$max_abs_diff, mcse = NA_real_, mcse_note = "a maximum; no MCSE",
        n_success = as.integer(c5$comparisons), events = NA_integer_, stringsAsFactors = FALSE)
      i4 <- which(layout$method == "M4")
      ok4 <- b$status[i4, , drop = FALSE] == "ok"
      rej4 <- ok4 & b$p_value[i4, , drop = FALSE] < alpha
      defined <- colSums(ok4)
      q <- colSums(rej4)[defined > 0] / defined[defined > 0]
      romano[[j]] <- data.frame(scenario = s$scenario, mean_q = if (length(q)) mean(q) else NA_real_,
                                var_q = if (length(q) > 1L) stats::var(q) else NA_real_,
                                reps_q = length(q), failed = any(!ok4))
    }
    # Failed fits per method: the replications in which the method failed,
    # and, per status, the replications with a contrast in that status. The
    # protocol counts failures "per method and scenario"; a replication
    # counts as a failed fit of a method when any one of its contrasts
    # failed, as Amendment 2 of #116 (item 11) says.
    fr <- list()
    for (m in unique(layout$method)) {
      st <- b$status[layout$method == m, , drop = FALSE]
      bad <- colSums(st != "ok") > 0
      fr[[length(fr) + 1L]] <- data.frame(scenario = s$scenario, type = "failed_fit", method = m,
                                          status = "any", count = sum(bad), replications = b$reps,
                                          stringsAsFactors = FALSE)
      for (code in setdiff(sort(unique(as.vector(st))), "ok")) {
        fr[[length(fr) + 1L]] <- data.frame(scenario = s$scenario, type = "failed_fit", method = m,
                                            status = code, count = sum(colSums(st == code) > 0),
                                            replications = b$reps, stringsAsFactors = FALSE)
      }
    }
    fails[[j]] <- do.call(rbind, fr)
    perf[[j]] <- do.call(rbind, rows)
  }
  perf <- do.call(rbind, perf)
  rownames(perf) <- NULL
  perf$role <- performance_roles(perf, sc)
  romano <- do.call(rbind, romano)
  fails <- do.call(rbind, fails)
  rownames(fails) <- NULL
  list(performance = perf, romano = romano, failed_fits = fails)
}

# ---- Decision rules ------------------------------------------------------------------

# Which function a failure belongs to (protocol, "Which function"): a rule 1
# failure to the function that computes the test, a rule 2 failure to
# equivalence_solomon() for fits of that method, and a fit failure to the
# function that fits or tests. Rule 2's function is kept apart by the fit's
# method, and M3 and M4 are one function, perm_solomon(), as the protocol
# names it; Amendment 2 of #116 (item 11) fixes both.
method_function <- c(M1 = "fit_solomon_glm", M2 = "fit_solomon_ml", M3 = "perm_solomon", M4 = "perm_solomon")
rule_function <- function(rule, method) {
  ifelse(startsWith(rule, "2"), paste0("equivalence_solomon (", method_function[method], " fits)"),
         method_function[method])
}

# The columns of decisions.csv. `n_success` is the number of successful fits
# of a judged rate and `count` its events; elsewhere `count` and `total`
# are what `measure` names (failed replications of all replications, flagged
# rates of all rates, comparisons, material failures).
decision_frame <- function(n, ...) {
  x <- data.frame(rule = rep(NA_character_, max(n, 1L)), scenario = NA_integer_, level = NA_character_,
                  pattern = NA_character_, dist = NA_character_, alloc = NA_character_,
                  alloc_type = NA_character_, min_cell = NA_integer_, h = NA_real_, rho = NA_real_,
                  method = NA_character_, contrast = NA_character_, measure = NA_character_,
                  n_success = NA_integer_, count = NA_integer_, total = NA_integer_, estimate = NA_real_,
                  mcse = NA_real_, p_value = NA_real_, p_value_material = NA_real_, threshold = NA_real_,
                  verdict = "", class = "", `function` = NA_character_, stress_test = NA, feature = "",
                  stringsAsFactors = FALSE, check.names = FALSE)[seq_len(n), , drop = FALSE]
  values <- list(...)
  for (nm in names(values)) {
    if (!nm %in% names(x)) stop("decisions.csv has no column ", nm, ".", call. = FALSE)
    x[[nm]] <- rep_len(if (is.factor(values[[nm]])) as.character(values[[nm]]) else values[[nm]], n)
  }
  rownames(x) <- NULL
  x
}

# A row of decisions.csv that is not a judged rate.
decision_row <- function(rule, verdict, fun = NA_character_, ...) {
  decision_frame(1L, rule = rule, verdict = verdict, `function` = fun, ...)
}

# The scenario's factor columns of decision rows.
scenario_columns <- function(scenario, sc) {
  x <- sc[match(scenario, sc$scenario), , drop = FALSE]
  list(scenario = scenario, level = x$level, pattern = x$pattern, dist = x$dist, alloc = x$alloc,
       alloc_type = x$alloc_type, min_cell = x$min_cell, h = x$h, rho = x$rho,
       stress_test = x$dist %in% stress_dists)
}

# Rules 1 to 3 on every judged rate, at the Bonferroni level .05/C, with the
# material tier (a failing rate whose one-sided exact binomial test of rate
# <= .075 also has p < .05/C; Bradley's, 1978, liberal criterion) and the
# defects. A rate whose fits all failed counts as failed (protocol); so does
# a bias whose test cannot be computed (no successful fit, a single one, or
# an MCSE of 0; bias_p() is then NA). Such a rate or bias is classed as a
# documented limitation (a defect under rule 3, and under rule 1 of M3 and
# M4 at the sharp null), never as material, as Amendment 2 of #116
# (item 11) says.
judge_rates <- function(perf, sc, C) {
  level <- 0.05 / C
  j <- judged_rates(sc)
  key <- function(x) paste(x$scenario, x$method, x$contrast, x$measure, sep = "\r")
  at <- match(key(j), key(perf))
  if (anyNA(at)) stop("A judged rate has no performance row.", call. = FALSE)
  p <- perf[at, , drop = FALSE]
  x <- sc[match(j$scenario, sc$scenario), , drop = FALSE]
  bias <- startsWith(j$rule, "3")
  p_value <- ifelse(bias, bias_p(p$estimate, p$mcse), binom_p_upper(p$events, p$n_success, nominal))
  p_material <- ifelse(bias, NA_real_, binom_p_upper(p$events, p$n_success, material_rate))
  fail <- is.na(p_value) | p_value < level
  material <- !bias & fail & !is.na(p_material) & p_material < level
  sharp <- j$method %in% c("M3", "M4") & startsWith(j$rule, "1") & x$sharp_null
  class <- ifelse(!fail, "", ifelse(bias | sharp, "defect",
                                    ifelse(material, "material failure", "documented limitation")))
  do.call(decision_frame, c(list(nrow(j)), scenario_columns(j$scenario, sc),
                            list(rule = j$rule, method = j$method, contrast = j$contrast, measure = j$measure,
                                 n_success = p$n_success, count = p$events, estimate = p$estimate,
                                 mcse = p$mcse, p_value = p_value, p_value_material = p_material,
                                 threshold = level, verdict = ifelse(fail, "fail", "pass"), class = class,
                                 `function` = rule_function(j$rule, j$method))))
}

# The fit-failure rule (protocol): a method with more than 1% failed fits in
# a scenario has a documented limitation there, and more than 5% is a
# material failure. It is not part of C. The rows carry the scenario's
# factors, so that a material one enters (a) to (c) below like a material
# rate failure, and one under a stress-test distribution does not.
fit_failure_rows <- function(failed_fits, sc) {
  ff <- failed_fits[failed_fits$status == "any", , drop = FALSE]
  prop <- ff$count / ff$replications
  do.call(decision_frame, c(list(nrow(ff)), scenario_columns(ff$scenario, sc),
                            list(rule = "fit failures", method = ff$method,
                                 measure = "proportion of replications with a failed fit",
                                 count = as.integer(ff$count), total = as.integer(ff$replications),
                                 estimate = prop, threshold = 0.01,
                                 verdict = ifelse(prop > 0.01, "fail", "pass"),
                                 class = ifelse(prop > 0.05, "material failure",
                                                ifelse(prop > 0.01, "documented limitation", "")),
                                 `function` = unname(method_function[ff$method]))))
}

# Consequence (a) of the material failures (protocol): per function, n* is
# the smallest of 10, 20, and 60 such that the function has no material
# failure in a balanced allocation with n* or more per group (NA when one
# remains at 60 per group). Each material failure outside the stress tests
# counts toward the first feature that covers it: a small cell (n* of 20 or
# 60, and a smallest cell below n*), unequal arms (1:3 or 3:1), or unequal
# pretest conditions. [Left to the script: computed mechanically.] The
# protocol calls more than 5% failed fits "a material failure", and a
# material failure leads to whichever of (a), (b), and (c) applies; the
# script takes this literally, so the fit-failure rule's material failures
# count toward n*, the features, (b), and (c) like those of rules 1 and 2,
# as Amendment 2 of #116 (item 11) says.
assign_features <- function(dec) {
  dec$feature[dec$class == "material failure" & dec$stress_test %in% TRUE] <-
    "stress test: documented, not counted toward (a), (b), or (c)"
  nstar <- list()
  counts <- list()
  for (f in sort(unique(dec$`function`))) {
    mat <- dec$`function` == f & dec$class == "material failure" & dec$stress_test %in% FALSE
    bal <- dec$min_cell[mat & dec$alloc_type == "balanced"]
    n_star <- if (!length(bal)) 10L else if (max(bal) < 20L) 20L else if (max(bal) < 60L) 60L else NA_integer_
    nstar[[f]] <- n_star
    counts[[f]] <- sum(mat)
    small <- !is.na(n_star) & n_star > 10L & dec$min_cell < n_star
    feat <- ifelse(small, paste0("small cell (n* = ", n_star, ")"),
                   ifelse(dec$alloc_type %in% c("treated smaller", "treated larger"), "unequal arms",
                          ifelse(dec$alloc_type %in% c("unpretested smaller", "unpretested larger"),
                                 "unequal pretest conditions", "none: (b) is considered")))
    dec$feature[mat] <- feat[mat]
  }
  list(decisions = dec, n_star = unlist(nstar), material = unlist(counts))
}

# The Romano check (protocol): over the k M4 scenarios of each pairing
# (main grid, E0, h > 1; k = 40, or 32 in the fallback), the mean of the 4k
# M4 rejection rates (k scenarios, 4 contrasts; the protocol's "160 M4
# rates" of the full design) should exceed .05 by more than 2 MCSE where
# the treated arms are smaller and fall below .05 by more than 2 MCSE where
# they are larger. The MCSE comes from each replication's proportion q of
# rejected contrasts among those defined: (1/k) sqrt(sum_s var_r(q) / R_s),
# with R_s the replications of scenario s with a contrast defined, as
# Amendment 2 of #116 (item 11) fixes it; the protocol says only that the
# MCSE "is computed from each replication's proportion of rejected
# contrasts". Where no M4 contrast failed in the pairing, the mean of the
# rates equals (1/k) sum_s mean_r(q), and the summary stops if it does not.
# The row's label names no count (4k is 160, or 128 in the fallback); its
# `total` gives the number of rates averaged (4k, less any rate whose fits
# all failed).
romano_rows <- function(romano, perf, sc) {
  x <- sc[match(romano$scenario, sc$scenario), , drop = FALSE]
  pick <- x$level == "main" & !is.na(x$h) & x$h > 1
  out <- list()
  for (pairing in c("treated smaller", "treated larger")) {
    r <- romano[pick & x$alloc_type == pairing, , drop = FALSE]
    k <- nrow(r)
    rates <- perf$estimate[perf$method == "M4" & perf$measure == "rejection" & perf$scenario %in% r$scenario]
    mean_rate <- if (any(!is.na(rates))) mean(rates, na.rm = TRUE) else NA_real_
    if (k && !any(r$failed) && !isTRUE(abs(mean_rate - mean(r$mean_q)) <= 1e-12)) {
      stop("The Romano check's two means of the M4 rates disagree (", pairing, ").", call. = FALSE)
    }
    mcse <- sqrt(sum(r$var_q / r$reps_q)) / k
    as_predicted <- if (pairing == "treated smaller") mean_rate - nominal > 2 * mcse else nominal - mean_rate > 2 * mcse
    as_predicted <- isTRUE(as_predicted)
    out[[pairing]] <- decision_row(
      rule = paste0("Romano check (", pairing, " arms, h > 1)"),
      verdict = if (as_predicted) "as predicted" else "contradicted",
      class = if (as_predicted) "" else "defect", method = "M4",
      measure = "mean of the pairing's M4 rejection rates", estimate = mean_rate, mcse = mcse,
      threshold = nominal, total = as.integer(sum(!is.na(rates))), fun = "perm_solomon",
      contrast = "the four treatment contrasts")
  }
  do.call(rbind, out)
}

# Patterns (protocol): for each rule and method, the number of rates flagged
# at an unadjusted .05 against the 5% expected by chance, overall and broken
# down by distribution, allocation, spread, correlation, and contrast, one
# factor at a time, as in #45 (the pattern tables of the article on designs
# with several treatments). Each row gives the flagged rates (`count`), the
# rates (`total`), their share (`estimate`), and .05 (`threshold`); the
# level is in the factor's own column. In the slopes level h and rho are NA.
pattern_factors <- c(distribution = "dist", allocation = "alloc", spread = "h", correlation = "rho",
                     contrast = "contrast")
pattern_rows <- function(dec) {
  judged <- dec[dec$rule %in% c("1. error control", "2. one-sided error (TOST)", "3. bias"), , drop = FALSE]
  out <- list()
  one <- function(rule, m, label, d, extra = list()) {
    flagged <- !is.na(d$p_value) & d$p_value < 0.05
    do.call(decision_row, c(list(rule = paste("patterns,", rule), verdict = "described", method = m,
                                 measure = paste("share of rates with p < .05 (unadjusted),", label),
                                 estimate = mean(flagged), count = as.integer(sum(flagged)),
                                 total = nrow(d), threshold = 0.05), extra))
  }
  for (rule in unique(judged$rule)) {
    for (m in unique(judged$method[judged$rule == rule])) {
      d <- judged[judged$rule == rule & judged$method == m, , drop = FALSE]
      out[[length(out) + 1L]] <- one(rule, m, "all rates", d)
      for (f in names(pattern_factors)) {
        col <- pattern_factors[[f]]
        lev <- d[[col]]
        key <- ifelse(is.na(lev), "NA", as.character(lev))
        for (k in unique(key)) {
          i <- key == k
          extra <- stats::setNames(list(lev[i][[1]]), col)
          out[[length(out) + 1L]] <- one(rule, m, paste("by", f), d[i, , drop = FALSE], extra)
        }
      }
    }
  }
  do.call(rbind, out)
}

# The agreement checks (protocol): a failure of check 1, 2, 3, or 5 is a
# defect; if check 4 fails, the E0 shortcut is abandoned. One row per check,
# method, and quantity, with that quantity's own tolerance and the number of
# comparisons.
agreement_decision_rows <- function(agreement) {
  key <- paste(agreement$check, agreement$method, agreement$quantity, sep = "\r")
  out <- list()
  for (k in unique(key)) {
    a <- agreement[key == k, , drop = FALSE]
    compared <- !is.na(a$passed)
    verdict <- if (!any(compared)) "not compared" else if (all(a$passed[compared])) "pass" else "fail"
    four <- startsWith(a$check[[1]], "4")
    out[[length(out) + 1L]] <- decision_row(
      rule = paste("agreement check", a$check[[1]]), verdict = verdict, method = a$method[[1]],
      class = if (verdict != "fail") "" else if (four) "E0 shortcut abandoned (an amendment adds a pattern)" else "defect",
      measure = paste("largest absolute difference:", a$quantity[[1]]), estimate = biggest(a$max_abs_diff),
      threshold = a$tolerance[[1]], count = as.integer(sum(a$comparisons, na.rm = TRUE)))
  }
  do.call(rbind, out)
}

decisions_table <- function(perf, romano, failed_fits, agreement, sc, C) {
  dec <- rbind(judge_rates(perf, sc, C), fit_failure_rows(failed_fits, sc))
  feat <- assign_features(dec)
  dec <- feat$decisions
  level <- 0.05 / C
  outside <- dec$class == "material failure" & dec$stress_test %in% FALSE
  more <- list()
  # n* per function (the level is that of the material tier), and
  # consequences (b) and (c).
  for (f in names(feat$n_star)) {
    n_star <- feat$n_star[[f]]
    more[[length(more) + 1L]] <- decision_row(
      rule = "(a) small-cell warning: n*", fun = f, measure = "n_star", estimate = n_star,
      count = as.integer(feat$material[[f]]), threshold = level,
      verdict = if (is.na(n_star)) "a material failure remains at 60 per group: no small-cell warning; (b) is considered"
                else if (n_star == 10L) "n* = 10: no small-cell warning"
                else paste0("n* = ", n_star, ": warning when a cell has fewer than ", n_star))
  }
  # (b) counts the material failures of M1's own function,
  # fit_solomon_glm() (rule 1 and fit failures); rule 2 failures of M1's
  # intervals belong to equivalence_solomon() ("each function from its own
  # failures") and are reported on the next row, not counted, as
  # Amendment 2 of #116 (item 11) says.
  m1 <- sum(dec$`function` == "fit_solomon_glm" & outside)
  more[[length(more) + 1L]] <- decision_row(
    rule = "(b) changed default", fun = "fit_solomon_glm", method = "M1",
    measure = "material failures outside the stress tests", estimate = m1, count = as.integer(m1),
    threshold = level,
    verdict = if (m1 > 0) "judgment required (see the M1 material failures)" else "not applicable")
  m1_rule2 <- sum(dec$`function` == "equivalence_solomon (fit_solomon_glm fits)" & outside)
  more[[length(more) + 1L]] <- decision_row(
    rule = "(b) changed default: rule 2 of fit_solomon_glm fits", fun = "equivalence_solomon (fit_solomon_glm fits)",
    method = "M1", measure = "material failures outside the stress tests (reported, not counted toward (b))",
    estimate = m1_rule2, count = as.integer(m1_rule2), threshold = level, verdict = "reported")
  m4 <- sum(dec$method == "M4" & outside & sc$equal_arms[match(dec$scenario, sc$scenario)], na.rm = TRUE)
  more[[length(more) + 1L]] <- decision_row(
    rule = "(c) difference statistic", fun = "perm_solomon", method = "M4",
    measure = "material failures in equal-arm allocations outside the stress tests", estimate = m4,
    count = as.integer(m4), threshold = level,
    verdict = if (m4 > 0) "warning extended to every call with statistic = \"difference\"" else "not triggered")
  # The Romano check.
  more[[length(more) + 1L]] <- romano_rows(romano, perf, sc)
  # The agreement checks.
  more[[length(more) + 1L]] <- agreement_decision_rows(agreement)
  # Patterns.
  more[[length(more) + 1L]] <- pattern_rows(dec)
  out <- rbind(dec, do.call(rbind, more))
  rownames(out) <- NULL
  out
}

# ---- Failures and warnings -----------------------------------------------------------

# failures.csv: the failed fits per method and scenario (every scenario and
# method, zeros included), and the warnings the contract muffled and
# counted in the replications, by scenario and class (protocol: "counted by
# class and reported"), with every distinct message of the class
# (`messages`, separated by " | ", in the order they first occur).
failures_table <- function(failed_fits, failures, ranges, sc) {
  ff <- failed_fits
  ff$class <- NA_character_
  ff$messages <- NA_character_
  ff$n_messages <- NA_integer_
  ff$proportion <- ff$count / ff$replications
  cols <- c("scenario", "type", "method", "status", "class", "messages", "n_messages", "count", "replications",
            "proportion")
  out <- ff[, cols]
  if (!is.null(failures) && nrow(failures)) {
    reps <- ranges$reps[match(failures$scenario, ranges$scenario)]
    key <- paste(failures$scenario, failures$class, sep = "\r")
    grp <- match(key, unique(key))
    first <- !duplicated(grp)
    w <- data.frame(scenario = failures$scenario[first], type = "warning", method = NA_character_,
                    status = NA_character_, class = failures$class[first],
                    messages = unname(as.character(tapply(failures$message, grp, function(v) {
                      paste(unique(v), collapse = " | ")
                    }))),
                    n_messages = unname(as.integer(tapply(failures$message, grp, function(v) length(unique(v))))),
                    count = unname(as.numeric(tapply(failures$count, grp, sum))),
                    replications = reps[first], stringsAsFactors = FALSE)
    w$proportion <- unname(as.numeric(tapply(failures$rep, grp, function(r) length(unique(r))))) /
      w$replications
    out <- rbind(out, w[, cols])
  }
  out <- out[order(out$scenario, out$type, out$method, out$status), , drop = FALSE]
  rownames(out) <- NULL
  out
}

# ---- Run information -----------------------------------------------------------------

# The study's part of run-information.csv: what the protocol lists that the
# harness does not record itself. The harness (tools/sim/finalize-results.R)
# adds, after these, the commit, the package and R versions, the platform,
# BLAS and LAPACK, the RNG kinds, the base seed, R, the workers, and the
# times; sessionInfo() is saved for each job in timing/sessionInfo-<job>.txt.
# Package versions come from the Version field only; broom, which left the
# Imports in #153, is "not installed" on the runners, as Amendment 2 of
# #116 (item 5) says; elsewhere its installed version, if any, is recorded.
version_of <- function(pkg) {
  v <- suppressWarnings(utils::packageDescription(pkg, fields = "Version"))
  if (is.null(v) || length(v) != 1L || is.na(v)) "not installed" else as.character(v)
}

run_information <- function(C) {
  data.frame(item = c(
    "permutations", "permutations_check1", "permutations_check2", "check1_seed",
    "check2_seeds", "tie_rule", "tost_bounds", "tost_alpha", "test_alpha", "family_size_C",
    "family_level", "material_rate", "sandwich_version", "broom_version", "withr_version",
    "summarize_sha256"),
    value = c(
      n_perm, n_perm_check1, n_perm_check2,
      "1000 * scenario + replication, through perm_solomon(seed =) under the task's RNG kind",
      paste0("data ", check2_data_seed + 1L, " to ", check2_data_seed + check2_sets, ", permutations ",
             check2_perm_seed + 1L, " to ", check2_perm_seed + check2_sets,
             ", under L'Ecuyer-CMRG, Inversion, Rejection"),
      "|s_b| >= |s| - 1e-10 * max(|s|, u); u = 1 (studentized) or the largest finite |s_b| (difference)",
      paste(tost_deltas, collapse = " "), format(alpha), format(alpha), format(C),
      format(0.05 / C, digits = 6), format(material_rate),
      version_of("sandwich"), version_of("broom"), version_of("withr"),
      sim_sha256_text(file.path(robustness_dir, "robustness-summarize.R"))),
    stringsAsFactors = FALSE)
}

# ---- The summary -----------------------------------------------------------------------

summarize_robustness <- function(results, agreement, failures, scenarios, info, preflight = NULL) {
  purpose <- info[["purpose"]]
  study_run <- identical(purpose, "study")
  sc <- scenarios
  if (!is.null(results)) check_statuses(results)
  out <- list()
  out$agreement <- agreement_summary(agreement, results, sc, preflight, study_run)
  # Outside the study run, nothing else is computed: no performance measure,
  # no decision, and no count of failed fits, as the protocol and
  # Amendment 2 of #116 (item 9) say. The family size reported is the
  # protocol's.
  if (!study_run) {
    out$`run-information` <- run_information(family_full)
    bad <- out$agreement[out$agreement$passed %in% FALSE, , drop = FALSE]
    if (nrow(bad) && identical(purpose, "test")) {
      stop("Agreement check failed: ",
           paste(unique(paste0(bad$check, " (", bad$method, ", ", bad$quantity, ")",
                               ifelse(nzchar(bad$detail), paste0(": ", bad$detail), ""))), collapse = "; "),
           ".", call. = FALSE)
    }
    return(out)
  }
  C <- family_size(sc)
  out$`run-information` <- run_information(C)
  ranges <- scenario_ranges(results)
  ranges$reps <- vapply(seq_len(nrow(ranges)), function(j) {
    k <- if (identical(sc$pattern[match(ranges$scenario[j], sc$scenario)], "E0")) nrow(layout_E0) else nrow(layout_ET)
    as.integer((ranges$end[j] - ranges$start[j] + 1L) %/% k)
  }, 0L)
  c5 <- check5_rows(results, sc)
  pt <- performance_table(results, sc, c5)
  out$performance <- pt$performance
  out$decisions <- decisions_table(pt$performance, pt$romano, pt$failed_fits, out$agreement, sc, C)
  out$failures <- failures_table(pt$failed_fits, failures, ranges, sc)
  # The scenario table, so that the rows of performance.csv and decisions.csv
  # can be read by their factors (the protocol's tables by factor).
  out$scenarios <- sc
  out
}
