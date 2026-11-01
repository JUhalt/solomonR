# Summaries for the issue #116 simulation (robustness of the recommended
# inference): the agreement checks, and, in the study run only, the
# performance measures with their Monte Carlo standard errors (MCSEs), the
# decision rules, the failed fits and warnings, and the study's part of
# run-information.csv.
#
# Sourced by robustness-simulation.R, whose constants and scenario table it
# uses; it defines functions only. It is base R, because the harness does
# not load the package in its summary stage (tools/sim/contract.R).
#
# What may be computed depends on the run's purpose (protocol, "Disclosed
# runs"; draft Amendment 2, item 9). In code-testing runs (purpose "test",
# which includes the smoke mode of tools/check-validation-scripts.R) and in
# the timing gate, summarize_robustness() computes the agreement checks and
# the run information only: no performance measure, no decision, and no
# count of failed fits. It checks the purpose before it computes anything
# else, and the contract stops if such a summary contains performance,
# decisions, or failures. The unit tests of these functions are in
# robustness-summarize-tests.R and use synthetic inputs only.
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
# fits of the contrast succeeded. One row per scenario and contrast.
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
    mx <- apply(d, 1L, function(v) if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE))
    out[[length(out) + 1L]] <- data.frame(scenario = s$scenario, contrast = contrasts_all,
                                          max_abs_diff = unname(mx), passed = is.na(mx) | mx <= tol_check5,
                                          stringsAsFactors = FALSE)
  }
  if (!length(out)) return(NULL)
  do.call(rbind, out)
}

# ---- Agreement summary ---------------------------------------------------------------

# One row per check, scenario, method, and quantity: whether every
# comparison passed and the largest absolute difference. In code-testing
# runs and the gate no row says how many comparisons were made, because a
# comparison is skipped when a fit failed, so a count would reveal failed
# fits; the study run adds the count (column `comparisons`). Check 5's
# largest difference is also one of the protocol's descriptive measures, so
# outside the study run only its verdict is shown.
agreement_summary <- function(agreement, results, sc, preflight, study_run) {
  rows <- list()
  biggest <- function(v) if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE)
  # Check 2 (ties), from the preflight, which runs before any task.
  if (!is.null(preflight) && "dataset" %in% names(preflight)) {
    p <- preflight[nzchar(preflight$dataset) & preflight$dataset != "NA", , drop = FALSE]
    p$passed <- as.logical(p$passed)
    p$max_abs_diff <- as.numeric(p$max_abs_diff)
    for (key in unique(paste(p$method, p$contrast, sep = "\r"))) {
      q <- p[paste(p$method, p$contrast, sep = "\r") == key, , drop = FALSE]
      rows[[length(rows) + 1L]] <- data.frame(
        check = "2. ties (preflight)", scenario = NA_integer_, method = q$method[[1]],
        contrast = q$contrast[[1]], quantity = "p-value, valid and tie-labelled permutations, statistics",
        max_abs_diff = biggest(q$max_abs_diff), tolerance = tol_check1,
        passed = all(q$passed), warnings_muffled = sum(as.integer(q$warnings_muffled)),
        detail = paste0(nrow(q), " fixed data sets; tie-labelled permutations ", sum(as.integer(q$tie_labelled)),
                        ", counted by the package ", sum(as.integer(q$counted_package)),
                        ", by the direct computation ", sum(as.integer(q$counted_direct))),
        comparisons = nrow(q), stringsAsFactors = FALSE)
    }
  }
  # Checks 1, 3, and 4, from the replications.
  if (!is.null(agreement) && nrow(agreement)) {
    key <- paste(agreement$check, agreement$scenario, agreement$method, agreement$quantity, sep = "\r")
    grp <- match(key, unique(key))
    first <- !duplicated(grp)
    classes <- tapply(agreement$warning_classes, grp, function(v) {
      v <- unique(unlist(strsplit(v[nzchar(v)], " ", fixed = TRUE)))
      paste(v, collapse = " ")
    })
    rows[[length(rows) + 1L]] <- data.frame(
      check = agreement$check[first], scenario = agreement$scenario[first],
      method = agreement$method[first], contrast = "all", quantity = agreement$quantity[first],
      max_abs_diff = unname(tapply(agreement$abs_diff, grp, biggest)),
      tolerance = agreement$tolerance[first],
      passed = unname(tapply(agreement$passed, grp, all)),
      warnings_muffled = unname(as.integer(tapply(agreement$warnings_muffled, grp, sum))),
      detail = ifelse(nzchar(classes), paste("warnings muffled:", classes), ""),
      comparisons = unname(as.integer(tabulate(grp))), stringsAsFactors = FALSE)
  }
  # Check 5, from the stored estimates.
  c5 <- if (is.null(results)) NULL else check5_rows(results, sc)
  if (!is.null(c5)) {
    rows[[length(rows) + 1L]] <- data.frame(
      check = "5. M2 equals M1", scenario = c5$scenario, method = "M2", contrast = c5$contrast,
      quantity = "estimate", max_abs_diff = if (study_run) c5$max_abs_diff else NA_real_,
      tolerance = tol_check5, passed = c5$passed, warnings_muffled = 0L,
      detail = if (study_run) "" else "largest difference shown in the study run only (a descriptive measure)",
      comparisons = NA_integer_, stringsAsFactors = FALSE)
  }
  if (!length(rows)) stop("No agreement check was run.", call. = FALSE)
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  if (!study_run) out$comparisons <- NULL
  out
}

# ---- Performance measures ------------------------------------------------------------

note_rate <- "proportion at 0 or 1: the binomial formula gives 0"

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

# Empirical SE, SD(estimate); model SE, sqrt(mean(SE^2)); and the relative
# error of the model SE, 100 (model SE / empirical SE - 1), with jackknife
# MCSEs over the replications (protocol). The leave-one-out values are
# computed in O(n) from centered sums:
#   SS_(j) = SS - n/(n - 1) (x_j - mean)^2, empirical SE_(j) = sqrt(SS_(j)/(n - 2));
#   model SE_(j) = sqrt((sum SE^2 - SE_j^2)/(n - 1)).
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
    note[] <- "unavailable: the fourth moment of t (3 df) is infinite (protocol)"
  } else if (n < 3) {
    note[] <- "fewer than 3 successful fits"
  } else {
    dev2 <- (est - mean(est))^2
    ss <- sum(dev2)
    emp_j <- sqrt(pmax(ss - n / (n - 1) * dev2, 0) / (n - 2))
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

# Bias of M1, mean(estimate) in E0, with MCSE SD(estimate)/sqrt(n); and the
# standardized bias, bias / empirical SE. [Left to the script: the protocol
# gives no MCSE for the standardized bias.]
bias_rows <- function(est) {
  n <- length(est)
  b <- if (n >= 1) mean(est) else NA_real_
  sdv <- if (n >= 2) stats::sd(est) else NA_real_
  data.frame(measure = c("bias", "standardized_bias"), estimate = c(b, b / sdv),
             mcse = c(sdv / sqrt(n), NA_real_),
             mcse_note = c(if (n < 2) "fewer than 2 successful fits" else "",
                           "not specified by the protocol; the bias carries its MCSE"),
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
  if (method == "M1") out <- rbind(out, bias_rows(est))
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
# (for the Romano check).
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
        n_success = NA_integer_, events = NA_integer_, stringsAsFactors = FALSE)
      i4 <- which(layout$method == "M4")
      ok4 <- b$status[i4, , drop = FALSE] == "ok"
      rej4 <- ok4 & b$p_value[i4, , drop = FALSE] < alpha
      defined <- colSums(ok4)
      q <- colSums(rej4)[defined > 0] / defined[defined > 0]
      romano[[j]] <- data.frame(scenario = s$scenario, mean_q = if (length(q)) mean(q) else NA_real_,
                                var_q = if (length(q) > 1L) stats::var(q) else NA_real_,
                                reps_q = length(q))
    }
    # Failed fits per method: the replications in which the method failed for
    # at least one contrast, and, per status, the replications with a
    # contrast in that status. [Left to the script: the protocol counts
    # failures "per method and scenario".]
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
# equivalence_solomon() for fits of that method. [Left to the script: rule
# 2's function is kept apart by the fit's method, and M3 and M4 are one
# function, perm_solomon(), as the protocol names it.]
rule_function <- function(rule, method) {
  fit <- c(M1 = "fit_solomon_glm", M2 = "fit_solomon_ml", M3 = "perm_solomon", M4 = "perm_solomon")
  ifelse(startsWith(rule, "2"), paste0("equivalence_solomon (", fit[method], " fits)"), fit[method])
}

# Rules 1 to 3 on every judged rate, at the Bonferroni level .05/C, with the
# material tier (a failing rate whose one-sided exact binomial test of rate
# <= .075 also has p < .05/C; Bradley's, 1978, liberal criterion) and the
# defects. A rate whose fits all failed counts as failed (protocol); so does
# a bias that cannot be computed [left to the script].
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
  stress <- x$dist %in% stress_dists
  data.frame(rule = j$rule, scenario = j$scenario, level = x$level, pattern = x$pattern, dist = x$dist,
             alloc = x$alloc, alloc_type = x$alloc_type, min_cell = x$min_cell, h = x$h, rho = x$rho,
             method = j$method, contrast = j$contrast, measure = j$measure, n_success = p$n_success,
             estimate = p$estimate, mcse = p$mcse, p_value = p_value, p_value_material = p_material,
             threshold = level, verdict = ifelse(fail, "fail", "pass"), class = class,
             `function` = rule_function(j$rule, j$method), stress_test = stress, feature = "",
             stringsAsFactors = FALSE, check.names = FALSE)
}

# Consequence (a) of the material failures (protocol): per function, n* is
# the smallest of 10, 20, and 60 such that the function has no material
# failure in a balanced allocation with n* or more per group (NA when one
# remains at 60 per group). Each material failure outside the stress tests
# counts toward the first feature that covers it: a small cell (n* of 20 or
# 60, and a smallest cell below n*), unequal arms (1:3 or 3:1), or unequal
# pretest conditions. [Left to the script: computed mechanically from the
# rate failures of rules 1 and 2; the fit-failure rule's material failures
# are not counted toward n* (a question for the maintainer).]
assign_features <- function(dec) {
  dec$feature[dec$class == "material failure" & dec$stress_test] <-
    "stress test: documented, not counted toward (a), (b), or (c)"
  nstar <- list()
  for (f in sort(unique(dec$`function`))) {
    mat <- dec$`function` == f & dec$class == "material failure" & !dec$stress_test
    bal <- dec$min_cell[mat & dec$alloc_type == "balanced"]
    n_star <- if (!length(bal)) 10L else if (max(bal) < 20L) 20L else if (max(bal) < 60L) 60L else NA_integer_
    nstar[[f]] <- n_star
    small <- !is.na(n_star) & n_star > 10L & dec$min_cell < n_star
    feat <- ifelse(small, paste0("small cell (n* = ", n_star, ")"),
                   ifelse(dec$alloc_type %in% c("treated smaller", "treated larger"), "unequal arms",
                          ifelse(dec$alloc_type %in% c("unpretested smaller", "unpretested larger"),
                                 "unequal pretest conditions", "none: (b) is considered")))
    dec$feature[mat] <- feat[mat]
  }
  list(decisions = dec, n_star = unlist(nstar))
}

# A row of decisions.csv that is not a judged rate.
decision_row <- function(rule, verdict, class = "", method = NA_character_, scenario = NA_integer_,
                         measure = NA_character_, estimate = NA_real_, mcse = NA_real_,
                         threshold = NA_real_, n_success = NA_integer_, fun = NA_character_,
                         contrast = NA_character_, feature = "") {
  data.frame(rule = rule, scenario = scenario, level = NA_character_, pattern = NA_character_,
             dist = NA_character_, alloc = NA_character_, alloc_type = NA_character_,
             min_cell = NA_integer_, h = NA_real_, rho = NA_real_, method = method, contrast = contrast,
             measure = measure, n_success = n_success, estimate = estimate, mcse = mcse,
             p_value = NA_real_, p_value_material = NA_real_, threshold = threshold, verdict = verdict,
             class = class, `function` = fun, stress_test = NA, feature = feature,
             stringsAsFactors = FALSE, check.names = FALSE)
}

# The Romano check (protocol): over the 40 M4 scenarios of each pairing
# (main grid, E0, h > 1), the mean rate should exceed .05 by more than 2
# MCSE where the treated arms are smaller and fall below .05 by more than 2
# MCSE where they are larger. The mean and its MCSE come from each
# replication's proportion q of rejected contrasts among those defined:
# mean = (1/40) sum_s mean_r(q), MCSE = (1/40) sqrt(sum_s var_r(q) / R_s).
# [Left to the script: this formula; the protocol says only that the MCSE
# "is computed from each replication's proportion of rejected contrasts".]
romano_rows <- function(romano, sc) {
  x <- sc[match(romano$scenario, sc$scenario), , drop = FALSE]
  pick <- x$level == "main" & !is.na(x$h) & x$h > 1
  out <- list()
  for (pairing in c("treated smaller", "treated larger")) {
    r <- romano[pick & x$alloc_type == pairing, , drop = FALSE]
    k <- nrow(r)
    mean_rate <- mean(r$mean_q)
    mcse <- sqrt(sum(r$var_q / r$reps_q)) / k
    as_predicted <- if (pairing == "treated smaller") mean_rate - nominal > 2 * mcse else nominal - mean_rate > 2 * mcse
    as_predicted <- isTRUE(as_predicted)
    out[[pairing]] <- decision_row(
      rule = paste0("Romano check (", pairing, " arms, h > 1)"),
      verdict = if (as_predicted) "as predicted" else "contradicted",
      class = if (as_predicted) "" else "defect", method = "M4", measure = "mean rejection rate",
      estimate = mean_rate, mcse = mcse, threshold = nominal, n_success = as.integer(k),
      fun = "perm_solomon", contrast = "the four treatment contrasts")
  }
  do.call(rbind, out)
}

decisions_table <- function(perf, romano, failed_fits, agreement, sc, C) {
  dec <- judge_rates(perf, sc, C)
  feat <- assign_features(dec)
  dec <- feat$decisions
  level <- 0.05 / C
  more <- list()
  # n* per function, and consequences (b) and (c).
  for (f in names(feat$n_star)) {
    n_star <- feat$n_star[[f]]
    more[[length(more) + 1L]] <- decision_row(
      rule = "(a) small-cell warning: n*", fun = f, measure = "n_star", estimate = n_star,
      threshold = level,
      verdict = if (is.na(n_star)) "a material failure remains at 60 per group: no small-cell warning; (b) is considered"
                else if (n_star == 10L) "n* = 10: no small-cell warning"
                else paste0("n* = ", n_star, ": warning when a cell has fewer than ", n_star))
  }
  m1 <- sum(dec$`function` == "fit_solomon_glm" & dec$class == "material failure" & !dec$stress_test)
  more[[length(more) + 1L]] <- decision_row(
    rule = "(b) changed default", fun = "fit_solomon_glm", method = "M1",
    measure = "material failures outside the stress tests", estimate = m1, threshold = level,
    verdict = if (m1 > 0) "judgment required (see the M1 material failures)" else "not applicable")
  m4 <- sum(dec$method == "M4" & dec$class == "material failure" & !dec$stress_test &
              sc$equal_arms[match(dec$scenario, sc$scenario)])
  more[[length(more) + 1L]] <- decision_row(
    rule = "(c) difference statistic", fun = "perm_solomon", method = "M4",
    measure = "material failures in equal-arm allocations outside the stress tests", estimate = m4,
    threshold = level,
    verdict = if (m4 > 0) "warning extended to every call with statistic = \"difference\"" else "not triggered")
  # The Romano check.
  more[[length(more) + 1L]] <- romano_rows(romano, sc)
  # The fit-failure rule: more than 1% is a documented limitation, more than
  # 5% a material failure (protocol). It is not part of C.
  ff <- failed_fits[failed_fits$status == "any", , drop = FALSE]
  prop <- ff$count / ff$replications
  ff_class <- ifelse(prop > 0.05, "material failure", ifelse(prop > 0.01, "documented limitation", ""))
  fun <- c(M1 = "fit_solomon_glm", M2 = "fit_solomon_ml", M3 = "perm_solomon", M4 = "perm_solomon")
  more[[length(more) + 1L]] <- do.call(rbind, lapply(seq_len(nrow(ff)), function(i) {
    decision_row(rule = "fit failures", scenario = ff$scenario[i], method = ff$method[i],
                 measure = "proportion of replications with a failed fit", estimate = prop[i],
                 threshold = 0.01, n_success = as.integer(ff$replications[i]),
                 verdict = if (prop[i] > 0.01) "fail" else "pass", class = ff_class[i],
                 fun = unname(fun[ff$method[i]]))
  }))
  # The agreement checks: a failure of check 1, 2, 3, or 5 is a defect; if
  # check 4 fails, the E0 shortcut is abandoned (protocol).
  for (chk in unique(agreement$check)) {
    a <- agreement[agreement$check == chk, , drop = FALSE]
    ok <- all(a$passed)
    four <- startsWith(chk, "4")
    biggest <- if (all(is.na(a$max_abs_diff))) NA_real_ else max(a$max_abs_diff, na.rm = TRUE)
    more[[length(more) + 1L]] <- decision_row(
      rule = paste("agreement check", chk), verdict = if (ok) "pass" else "fail",
      class = if (ok) "" else if (four) "E0 shortcut abandoned (an amendment adds a pattern)" else "defect",
      measure = "largest absolute difference", estimate = biggest,
      threshold = max(a$tolerance, na.rm = TRUE))
  }
  # Patterns (protocol): for each rule and method, the share of rates
  # flagged at an unadjusted .05, against the 5% expected by chance. The
  # breakdown by factor is made from the rows above, which carry every
  # rate's p-value and factor columns.
  for (rule in unique(dec$rule)) {
    for (m in unique(dec$method[dec$rule == rule])) {
      d <- dec[dec$rule == rule & dec$method == m, , drop = FALSE]
      more[[length(more) + 1L]] <- decision_row(
        rule = paste("patterns,", rule), method = m, measure = "share of rates with p < .05 (unadjusted)",
        estimate = mean(!is.na(d$p_value) & d$p_value < 0.05), n_success = nrow(d), threshold = 0.05,
        verdict = "described")
    }
  }
  out <- rbind(dec, do.call(rbind, more))
  rownames(out) <- NULL
  out
}

# ---- Failures and warnings -----------------------------------------------------------

# failures.csv: the failed fits per method and scenario (every scenario and
# method, zeros included), and the warnings the contract muffled and
# counted in the replications, by scenario and class.
failures_table <- function(failed_fits, failures, ranges, sc) {
  ff <- failed_fits
  ff$class <- NA_character_
  ff$message <- NA_character_
  ff$proportion <- ff$count / ff$replications
  out <- ff[, c("scenario", "type", "method", "status", "class", "message", "count", "replications",
                "proportion")]
  if (!is.null(failures) && nrow(failures)) {
    reps <- ranges$reps[match(failures$scenario, ranges$scenario)]
    key <- paste(failures$scenario, failures$class, sep = "\r")
    grp <- match(key, unique(key))
    first <- !duplicated(grp)
    w <- data.frame(scenario = failures$scenario[first], type = "warning", method = NA_character_,
                    status = NA_character_, class = failures$class[first],
                    message = failures$message[first],
                    count = unname(as.numeric(tapply(failures$count, grp, sum))),
                    replications = reps[first], stringsAsFactors = FALSE)
    w$proportion <- unname(as.numeric(tapply(failures$rep, grp, function(r) length(unique(r))))) /
      w$replications
    out <- rbind(out, w)
  }
  out <- out[order(out$scenario, out$type, out$method, out$status), , drop = FALSE]
  rownames(out) <- NULL
  out
}

# ---- Run information -----------------------------------------------------------------

# The study's part of run-information.csv: what the protocol lists that the
# harness does not record itself (the harness adds the commit, the package
# and R versions, BLAS and LAPACK, the RNG kinds, the seed, R, the workers,
# and the times). Package versions come from the Version field only; broom
# is "not installed" once it has left the Imports (draft Amendment 2, item 5).
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
  out <- list()
  out$agreement <- agreement_summary(agreement, results, sc, preflight, study_run)
  # Outside the study run, nothing else is computed: no performance measure,
  # no decision, and no count of failed fits (protocol; draft Amendment 2,
  # item 9). The family size reported is the protocol's.
  if (!study_run) {
    out$`run-information` <- run_information(family_full)
    bad <- out$agreement[!out$agreement$passed, , drop = FALSE]
    if (nrow(bad) && identical(purpose, "test")) {
      stop("Agreement check failed: ",
           paste(unique(paste0(bad$check, " (", bad$method, ", scenario ",
                               ifelse(is.na(bad$scenario), "-", bad$scenario), ")")), collapse = "; "),
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
  out
}
