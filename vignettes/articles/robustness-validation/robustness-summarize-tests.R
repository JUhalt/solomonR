# Unit tests of the summary of the issue #116 simulation
# (robustness-summarize.R) and of the study's gate projection, on synthetic
# inputs with known answers. Run from the package root:
#
#   Rscript vignettes/articles/robustness-validation/robustness-summarize-tests.R
#
# No input comes from the study's data-generating mechanisms: the results,
# agreement rows, and timings below are constructed by hand (fixed numbers,
# or sin() and cos() of an index), so every expected value is known before
# the code runs. Nothing is drawn at random and the package is not loaded.
# The thresholds come from the protocol's table (exact binomial counts at R
# = 5,000 to 3,000) and its counts of judged rates.
#
# tools/check-validation-scripts.R runs every .R file of this folder; this
# one has nothing to scale down, so SOLOMONR_SMOKE_DIR and
# SOLOMONR_SMOKE_REPS do not change it, and it writes no file. It stops
# with status 1 if any test fails.

# The study script's definitions, without its final sim_main() call.
local({
  exprs <- parse(file.path("vignettes", "articles", "robustness-validation", "robustness-simulation.R"))
  for (e in exprs) {
    if (is.call(e) && identical(e[[1]], as.name("sim_main"))) next
    eval(e, globalenv())
  }
})

test_failures <- character()
n_tests <- 0L
expect <- function(name, ok) {
  n_tests <<- n_tests + 1L
  if (!isTRUE(ok)) test_failures <<- c(test_failures, name)
  invisible(ok)
}
near <- function(a, b, tol = 1e-12) isTRUE(all(abs(a - b) <= tol))

# ---- The family of rates and the Bonferroni level -------------------------------------------

j <- judged_rates(scenarios)
sc_of <- scenarios[match(j$scenario, scenarios$scenario), ]
r1 <- startsWith(j$rule, "1")
h_main <- ifelse(sc_of$level == "main", sc_of$h, NA)
expect("C is 14,576", identical(check_family(scenarios), 14576L))
expect("the fallback's C is 11,964",
       identical(check_family(scenarios[study_fallback(scenarios), ], fallback = TRUE), 11964L))
expect("M2: 756 rates with h = 1 and opposite skew",
       sum(r1 & j$method == "M2" & (sc_of$level == "opposite skew" | (!is.na(h_main) & h_main == 1))) == 756L)
expect("M2: 500 rates with h > 1 in equal-arm allocations",
       sum(r1 & j$method == "M2" & !is.na(h_main) & h_main > 1) == 500L)
expect("M3: 1,224 E0 and 270 ET rates",
       sum(r1 & j$method == "M3" & sc_of$pattern == "E0") == 1224L &&
         sum(r1 & j$method == "M3" & sc_of$pattern == "ET") == 270L)
expect("M4: 360, 400, 72, and 40 rates",
       identical(as.integer(c(sum(r1 & j$method == "M4" & !is.na(h_main) & h_main == 1),
                              sum(r1 & j$method == "M4" & !is.na(h_main) & h_main > 1),
                              sum(r1 & j$method == "M4" & sc_of$level == "opposite skew"),
                              sum(r1 & j$method == "M4" & sc_of$level == "slopes"))),
                 c(360L, 400L, 72L, 40L)))
expect("M4 in the slopes level only with equal arms", all(sc_of$equal_arms[r1 & j$method == "M4" & sc_of$level == "slopes"]))
expect("rule 3: M1 in the main grid and opposite skew only",
       all(j$method[startsWith(j$rule, "3")] == "M1") && !any(sc_of$level[startsWith(j$rule, "3")] == "slopes"))
expect(".05/C = 3.4303e-6", near(0.05 / family_full, 3.4303e-6, 5e-11))
expect("family_size() refuses a partial table",
       inherits(tryCatch(family_size(scenarios[1:10, ]), error = function(e) e), "error"))

# ---- Rules 1 and 2: the exact binomial test at and around its thresholds ----------------------

expect("binom_p_upper() is P(X >= x)",
       near(binom_p_upper(300, 5000, 0.05), sum(stats::dbinom(300:5000, 5000, 0.05)), 1e-12))
expect("binom_p_upper() is NA without a successful fit", is.na(binom_p_upper(0, 0, 0.05)))
level <- 0.05 / 14576
# The protocol's table: a rate fails if at least .0646 (323 of 5,000), ...;
# it is material if at least .0924 (462 of 5,000), ...
table_counts <- data.frame(R = c(5000, 4500, 4000, 3500, 3000), fail = c(323, 295, 266, 237, 207),
                           material = c(462, 421, 379, 336, 294))
for (i in seq_len(nrow(table_counts))) {
  R <- table_counts$R[i]
  f <- table_counts$fail[i]
  m <- table_counts$material[i]
  expect(paste("rule 1 fails at", f, "of", R), binom_p_upper(f, R, 0.05) < level)
  expect(paste("rule 1 passes at", f - 1, "of", R), binom_p_upper(f - 1, R, 0.05) >= level)
  expect(paste("material at", m, "of", R), binom_p_upper(m, R, 0.075) < level)
  expect(paste("not material at", m - 1, "of", R), binom_p_upper(m - 1, R, 0.075) >= level)
}

# ---- Rule 3: the z test of bias / MCSE -------------------------------------------------------

z_crit <- stats::qnorm(1 - level / 2)
expect("bias_p() is 2 Phi(-|z|)", near(bias_p(0.02, 0.01), 2 * stats::pnorm(-2)))
expect("rule 3 fails just above the critical z", bias_p(z_crit * 1.0001, 1) < level)
expect("rule 3 passes just below the critical z", bias_p(z_crit * 0.9999, 1) >= level)
expect("bias_p() is NA without an MCSE", is.na(bias_p(0.1, NA)))

# ---- MCSE formulas ----------------------------------------------------------------------------

x <- rate_rows(list(a = rep(c(TRUE, FALSE), c(250, 4750)), b = rep(FALSE, 10), c = rep(TRUE, 10)[1:5]),
               5000)
expect("rate and its binomial MCSE", near(x$estimate[1], 0.05) && near(x$mcse[1], sqrt(0.05 * 0.95 / 5000)))
expect("a rate of 0 carries a note", x$estimate[2] == 0 && x$mcse[2] == 0 && nzchar(x$mcse_note[2]))
y <- rate_rows(list(a = logical(0)), 0)
expect("a rate without fits is NA with a note", is.na(y$estimate) && nzchar(y$mcse_note))

est <- sin(seq_len(30) * 1.3) + 0.2 * cos(seq_len(30) * 3.1)
se <- 1 + cos(seq_len(30)) / 3
loo <- function(f) vapply(seq_along(est), function(i) f(est[-i], se[-i]), 0)
brute <- function(v) { n <- length(v); sqrt((n - 1) / n * sum((v - mean(v))^2)) }
emp_j <- loo(function(e, s) stats::sd(e))
mod_j <- loo(function(e, s) sqrt(mean(s^2)))
rel_j <- loo(function(e, s) 100 * (sqrt(mean(s^2)) / stats::sd(e) - 1))
z <- se_rows(est, se, heavy_tailed = FALSE)
expect("empirical and model SEs", near(z$estimate[1], stats::sd(est)) && near(z$estimate[2], sqrt(mean(se^2))) &&
         near(z$estimate[3], 100 * (sqrt(mean(se^2)) / stats::sd(est) - 1)))
expect("jackknife MCSE of the empirical SE equals brute-force leave-one-out", near(z$mcse[1], brute(emp_j), 1e-12))
expect("jackknife MCSE of the model SE equals brute-force leave-one-out", near(z$mcse[2], brute(mod_j), 1e-12))
expect("jackknife MCSE of the relative error equals brute-force leave-one-out", near(z$mcse[3], brute(rel_j), 1e-10))
expect("jackknife_se() of the mean is SD/sqrt(n)",
       near(jackknife_se(vapply(seq_along(est), function(i) mean(est[-i]), 0)), stats::sd(est) / sqrt(30)))
z3 <- se_rows(est, se, heavy_tailed = TRUE)
expect("under t (3 df) the SE MCSEs are unavailable, with a note", all(is.na(z3$mcse)) && all(nzchar(z3$mcse_note)))
b <- bias_rows(est)
expect("bias and its MCSE SD/sqrt(n)", near(b$estimate[1], mean(est)) && near(b$mcse[1], stats::sd(est) / sqrt(30)))
expect("standardized bias, without an MCSE, with a note",
       near(b$estimate[2], mean(est) / stats::sd(est)) && is.na(b$mcse[2]) && nzchar(b$mcse_note[2]))

# ---- The TOST and the measures of a fit ----------------------------------------------------------

tp <- tost_p(0.3, 0.1, 20, 0.5)
expect("tost_p() p_lower and p_upper", near(tp$p_lower, stats::pt(8, 20, lower.tail = FALSE)) &&
         near(tp$p_upper, stats::pt(-2, 20)) && isTRUE(tp$equivalent))
expect("tost_p() is not equivalent when one test fails", !tost_p(0.45, 0.1, 20, 0.5)$equivalent)

# Five fits with t = -3, -1.8, 0, 1.8, 3 and df = 1,000: two reject at .05
# (one each way); one upper limit of the 90% interval is below 0 (t = -3 and
# t = -1.8 are below -1.646: two), and two lower limits are above it.
t5 <- c(-3, -1.8, 0, 1.8, 3)
fm <- fit_measures("M1", est = t5 * 0.1, se = rep(0.1, 5), df = rep(1000, 5), heavy_tailed = FALSE)
val <- function(name) fm$estimate[fm$measure == name]
expect("two-sided rejection", near(val("rejection"), 2 / 5) && near(val("rejection_below"), 1 / 5) &&
         near(val("rejection_above"), 1 / 5) && near(val("coverage"), 3 / 5))
expect("one-sided error rates", near(val("upper_limit_below_truth"), 2 / 5) && near(val("lower_limit_above_truth"), 2 / 5))
# TOST size at the upper bound with Delta = 0.8: estimate + 0.8, bounds
# +-0.8: equivalent when t < -1.646 and t + 16 > 1.646, i.e. t = -3 or -1.8.
expect("TOST size at the upper bound", near(val("tost_size_upper_0.8"), 2 / 5) && near(val("tost_size_lower_0.8"), 2 / 5))
# With Delta = 0.2: equivalent when t < -1.646 and t + 4 > 1.646, i.e. t = -1.8
# only (t = -3 gives t + 4 = 1).
expect("TOST size with a narrow bound", near(val("tost_size_upper_0.2"), 1 / 5) &&
         near(val("tost_size_lower_0.2"), 1 / 5))
# Power at 0.5: (estimate + 0.5) / 0.1 = t + 5 = 2, 3.2, 5, 6.8, 8: all reject.
expect("power at 0.5", near(val("power"), 1))
expect("bias of M1 only", "bias" %in% fm$measure &&
         !"bias" %in% fit_measures("M2", t5 * 0.1, rep(0.1, 5), rep(1000, 5), FALSE)$measure)
pm <- perm_measures(p = c(0.01, 0.03, 0.2, 0.5), stat = c(-2, 3, 1, -1), null_true = TRUE)
expect("permutation rejection and its directions",
       identical(pm$events, c(2L, 1L, 1L)) && identical(pm$measure, c("rejection", "rejection_below", "rejection_above")))
expect("permutation power in ET", identical(perm_measures(c(0.01, 0.2), c(1, 1), FALSE)$measure, "power"))

# ---- Decisions at and around the thresholds -----------------------------------------------------

# A synthetic performance table with every judged rate at 250 of 5,000
# (passes) and every bias at 0 (MCSE 1), then a few cells moved to the
# thresholds.
perf <- data.frame(j[, c("scenario", "method", "contrast", "measure")], estimate = 250 / 5000,
                   mcse = sqrt(0.05 * 0.95 / 5000), mcse_note = "", n_success = 5000L, events = 250L,
                   stringsAsFactors = FALSE)
perf$estimate[perf$measure == "bias"] <- 0
perf$mcse[perf$measure == "bias"] <- 1
perf$events[perf$measure == "bias"] <- NA_integer_
pick <- function(method, measure, cond) which(perf$method == method & perf$measure == measure & cond)
sc_p <- scenarios[match(perf$scenario, scenarios$scenario), ]
i_lim <- pick("M1", "rejection", sc_p$scenario == 1 & perf$contrast == contrasts_all[1])
i_pass <- pick("M1", "rejection", sc_p$scenario == 1 & perf$contrast == contrasts_all[2])
i_mat <- pick("M2", "upper_limit_below_truth", sc_p$scenario == 1 & perf$contrast == contrasts_all[1])
i_nomat <- pick("M2", "upper_limit_below_truth", sc_p$scenario == 1 & perf$contrast == contrasts_all[2])
i_sharp <- pick("M3", "rejection", sc_p$scenario == 1 & perf$contrast == contrasts_all[1])
i_none <- pick("M1", "rejection", sc_p$scenario == 2 & perf$contrast == contrasts_all[1])
i_bias <- pick("M1", "bias", sc_p$scenario == 1 & perf$contrast == contrasts_all[1])
i_bias_ok <- pick("M1", "bias", sc_p$scenario == 1 & perf$contrast == contrasts_all[2])
perf$events[i_lim] <- 323L
perf$events[i_pass] <- 322L
perf$events[i_mat] <- 462L
perf$events[i_nomat] <- 461L
perf$events[i_sharp] <- 323L
perf$events[i_none] <- 0L
perf$n_success[i_none] <- 0L
perf$estimate[i_bias] <- z_crit * 1.0001
perf$estimate[i_bias_ok] <- z_crit * 0.9999
dec <- judge_rates(perf, scenarios, 14576L)
rowkey <- function(i) which(dec$scenario == perf$scenario[i] & dec$method == perf$method[i] &
                              dec$contrast == perf$contrast[i] & dec$measure == perf$measure[i])
cls <- function(i) dec$class[rowkey(i)]
expect("decisions has one row per judged rate", nrow(dec) == 14576L)
expect("323 of 5,000 fails, not material: a documented limitation", cls(i_lim) == "documented limitation")
expect("322 of 5,000 passes", dec$verdict[rowkey(i_pass)] == "pass" && cls(i_pass) == "")
expect("462 of 5,000 is material", cls(i_mat) == "material failure")
expect("461 of 5,000 fails but is not material", cls(i_nomat) == "documented limitation")
expect("M3 failing under the sharp null is a defect", cls(i_sharp) == "defect")
expect("a rate whose fits all failed counts as failed", dec$verdict[rowkey(i_none)] == "fail")
expect("rule 3 just above the critical z is a defect", cls(i_bias) == "defect")
expect("rule 3 just below the critical z passes", dec$verdict[rowkey(i_bias_ok)] == "pass")
expect("the threshold is .05/C in every rate row", all(dec$threshold == level))
expect("rule 2 belongs to equivalence_solomon()", startsWith(dec$`function`[rowkey(i_mat)], "equivalence_solomon"))

# ---- The material tier's consequences: n* and the warning features ---------------------------------

alloc_rows <- function(alloc) allocations[match(alloc, allocations$alloc), ]
feat_case <- function(fun, allocs, stress = rep(FALSE, length(allocs))) {
  a <- alloc_rows(allocs)
  data.frame(`function` = fun, class = "material failure", stress_test = stress,
             min_cell = pmin(a$n1, a$n2, a$n3, a$n4), alloc_type = a$alloc_type, feature = "",
             stringsAsFactors = FALSE, check.names = FALSE)
}
fd <- rbind(
  feat_case("f20", c("10-10-10-10", "10-30-10-30", "60-20-60-20", "30-30-10-10", "20-60-20-60")),
  feat_case("f20", "10-30-10-30", stress = TRUE),
  feat_case("f10", c("10-30-10-30", "10-10-30-30")),
  feat_case("fNA", c("60-60-60-60", "30-30-10-10")))
af <- assign_features(fd)
ff <- af$decisions$feature
expect("n* = 20 after a balanced failure at 10", identical(af$n_star[["f20"]], 20L))
expect("n* = 10 without a balanced failure", identical(af$n_star[["f10"]], 10L))
expect("no n* after a failure at 60 per group", is.na(af$n_star[["fNA"]]))
expect("features with n* = 20",
       identical(ff[1:5], c("small cell (n* = 20)", "small cell (n* = 20)", "unequal arms", "small cell (n* = 20)",
                            "unequal arms")))
expect("a stress test does not count", startsWith(ff[6], "stress test"))
expect("features with n* = 10", identical(ff[7:8], c("unequal arms", "unequal pretest conditions")))
expect("features without n*", identical(ff[9:10], c("none: (b) is considered", "unequal pretest conditions")))

# ---- The Romano check ---------------------------------------------------------------------------

rom_sc <- scenarios$scenario[scenarios$level == "main" & scenarios$pattern == "E0" & scenarios$h > 1 &
                               scenarios$alloc_type %in% c("treated smaller", "treated larger")]
expect("40 scenarios in each pairing", length(rom_sc) == 80L)
romano <- data.frame(scenario = rom_sc, mean_q = ifelse(scenarios$alloc_type[rom_sc] == "treated smaller", 0.06, 0.04),
                     var_q = 0.01, reps_q = 5000L)
rr <- romano_rows(romano, scenarios)
expect("Romano MCSE = (1/40) sqrt(sum var/R)", near(rr$mcse, rep(sqrt(40 * 0.01 / 5000) / 40, 2)))
expect("Romano means", near(rr$estimate, c(0.06, 0.04)))
expect("Romano check as predicted", identical(rr$verdict, c("as predicted", "as predicted")))
romano$mean_q[scenarios$alloc_type[rom_sc] == "treated smaller"] <- 0.05 + 1.9 * sqrt(40 * 0.01 / 5000) / 40
expect("Romano check contradicted within 2 MCSE", romano_rows(romano, scenarios)$verdict[1] == "contradicted")

# ---- Synthetic results: layout, check 5, and the purposes ---------------------------------------------

# Results of every scenario with `reps` replications in replicate()'s
# layout; estimates are fixed functions of the index, every status "ok".
synthetic_results <- function(sc, reps, m2_offset = 0) {
  parts <- lapply(seq_len(nrow(sc)), function(i) {
    lay <- if (sc$pattern[i] == "E0") layout_E0 else layout_ET
    k <- nrow(lay)
    n <- k * reps
    idx <- seq_len(n) + 1000 * sc$scenario[i]
    est <- sin(idx * 0.37) / 3
    out <- data.frame(scenario = sc$scenario[i], rep = rep(seq_len(reps), each = k),
                      method = rep(lay$method, reps), contrast = rep(lay$contrast, reps), status = "ok",
                      estimate = est, std_error = 0.2 + abs(cos(idx)) / 10, df = 30,
                      statistic = NA_real_, p_value = NA_real_, n_valid = NA_integer_, n_tied = NA_integer_,
                      convergence = NA_integer_, stringsAsFactors = FALSE)
    m1 <- out$method == "M1"
    m2 <- out$method == "M2"
    out$estimate[m2] <- out$estimate[m1] + m2_offset
    pm <- out$method %in% c("M3", "M4")
    out$estimate[pm] <- NA_real_
    out$std_error[pm] <- NA_real_
    out$df[pm] <- NA_real_
    out$statistic[pm] <- cos(idx[pm])
    out$p_value[pm] <- (seq_len(sum(pm)) %% 20) / 20 + 0.01
    out$n_valid[pm] <- 999L
    out$n_tied[pm] <- 0L
    out
  })
  do.call(rbind, parts)
}
small <- scenarios[scenarios$scenario %in% c(1, 271, 541, 559), ]
c5_ok <- check5_rows(synthetic_results(small, 3, 1e-7), small)
c5_bad <- check5_rows(synthetic_results(small, 3, 2e-6), small)
expect("check 5 passes within 1e-6", all(c5_ok$passed) && nrow(c5_ok) == 3L * 7L)
expect("check 5 fails beyond 1e-6", !any(c5_bad$passed))

agr <- data.frame(scenario = 1L, rep = 1L, check = "1. permutation tests", method = "M3",
                  contrast = contrasts_treat[1], statistic = "studentized", delta = NA_real_, bound = NA_character_,
                  quantity = "p_value", package = 0.5, direct = 0.5, abs_diff = 0, tolerance = 0, passed = TRUE,
                  warnings_muffled = 0L, warning_classes = "", stringsAsFactors = FALSE)
info <- function(purpose) c(purpose = purpose, protocol = "TRUE", design = if (purpose == "gate") "gate" else "run")
for (purpose in c("test", "gate")) {
  s <- summarize_robustness(synthetic_results(small, 3), agr, NULL, small, info(purpose))
  expect(paste("a", purpose, "summary has agreement and run information only"),
         setequal(names(s), c("agreement", "run-information")))
  expect(paste("a", purpose, "summary shows no check-5 difference"),
         all(is.na(s$agreement$max_abs_diff[startsWith(s$agreement$check, "5")])))
  expect(paste("a", purpose, "summary has no count of comparisons"), !"comparisons" %in% names(s$agreement))
}
bad <- agr
bad$passed <- FALSE
expect("a test summary stops on a failed agreement check",
       inherits(tryCatch(summarize_robustness(synthetic_results(small, 3), bad, NULL, small, info("test")),
                         error = function(e) e), "error"))
gate_bad <- summarize_robustness(synthetic_results(small, 3), bad, NULL, small, info("gate"))
expect("a gate summary reports a failed agreement check", any(!gate_bad$agreement$passed))
expect("a study summary refuses a partial table",
       inherits(tryCatch(summarize_robustness(synthetic_results(small, 3), agr, NULL, small, info("study")),
                         error = function(e) e), "error"))
ri <- run_information(14576L)
expect("run information names no harness item",
       !any(ri$item %in% c("study", "issue", "purpose", "package", "package_version", "r_version", "base_seed",
                           "reps", "rng_kinds", "workers_per_chunk", "started_utc", "finished_utc")))
expect("broom is \"not installed\" when absent", identical(version_of("a.package.that.is.not.installed"), "not installed"))
expect("the permutations are recorded", ri$value[ri$item == "permutations"] == "999")

# ---- A whole study summary on synthetic results ---------------------------------------------------------

full <- summarize_robustness(synthetic_results(scenarios, 4), agr,
                             data.frame(scenario = 1L, rep = 2L, type = "warning", class = "simpleWarning",
                                        message = "x", count = 3L, stringsAsFactors = FALSE),
                             scenarios, info("study"))
pf <- full$performance
dc <- full$decisions
expect("a study summary has every file", setequal(names(full), c("agreement", "run-information", "performance",
                                                                 "decisions", "failures")))
expect("every performance row has an MCSE above 0 or a note",
       all((is.finite(pf$mcse) & pf$mcse > 0) | nzchar(pf$mcse_note)))
expect("performance rows of every scenario", setequal(unique(pf$scenario), scenarios$scenario))
num <- cbind(is.finite(dc$mcse), is.finite(dc$threshold), is.finite(dc$p_value))
expect("every decision row has an MCSE, threshold, or p-value", all(rowSums(num) > 0))
expect("decisions hold the 14,576 judged rates",
       sum(dc$rule %in% c("1. error control", "2. one-sided error (TOST)", "3. bias")) == 14576L)
expect("the run information gives C", full$`run-information`$value[full$`run-information`$item == "family_size_C"] == "14576")
expect("the warnings are in failures.csv", any(full$failures$type == "warning" & full$failures$count == 3))
expect("every scenario and method has a failed-fit row",
       sum(full$failures$type == "failed_fit" & full$failures$status == "any") == 306L * 4L + 270L * 2L)
expect("M2's largest difference from M1 is reported", any(pf$measure == "max_abs_difference_from_M1"))
expect("judged and described rates are labelled", all(c("judged", "described", "descriptive") %in% pf$role))
expect("M2 in the slopes level is described",
       all(pf$role[pf$method == "M2" & pf$measure == "rejection" &
                     pf$scenario %in% scenarios$scenario[scenarios$level == "slopes"]] == "described"))

# ---- The gate's projection ---------------------------------------------------------------------------

# One task of 25 replications per scenario: 0.1 s per replication and 10 s
# of agreement checks each; half the tasks start 10 minutes into the job.
gate_times <- function(seconds_per_rep) {
  start <- as.POSIXct("2026-11-01 00:00:00", tz = "UTC")
  n <- nrow(scenarios)
  offset <- ifelse(seq_len(n) %% 2 == 0, 600, 0)
  data.frame(task_id = sprintf("s%04d-t0001", seq_len(n)), scenario = scenarios$scenario, reps = 25L,
             started = sim_utc(start + offset), finished = sim_utc(start + offset + 1),
             job_started = sim_utc(start), elapsed = 25 * seconds_per_rep + 10, cpu = 1,
             replication_seconds = 25 * seconds_per_rep, agreement_seconds = 10,
             cost = vapply(seq_len(n), function(i) cost_robustness(scenarios[i, ], 25L), 0),
             stringsAsFactors = FALSE)
}
gi <- c(purpose = "test")
g1 <- gate_project(gate_times(0.1), scenarios, gi)
p1 <- g1$`gate-projection`
expect("three readings of the projection", setequal(unique(p1$reading),
       c("runner jobs (decides R)", "runner worker-hours", "posted rule (10 workers, 24 hours)")))
wh <- p1[p1$reading == "runner worker-hours" & p1$R == 5000, ]
expect("worker-hours = (R x sum of seconds per replication + agreement seconds) / 3600",
       near(wh$projection_hours, (5000 * 0.1 * 576 + 10 * 576) / 3600, 1e-9))
jb <- p1[p1$reading == "runner jobs (decides R)" & p1$R == 5000, ]
expect("job hours are worker-hours / 1,024", near(jb$projection_hours, wh$projection_hours / 1024, 1e-12))
expect("R = 5,000 is chosen when it fits", p1$chosen[p1$reading == "runner jobs (decides R)" & p1$R == 5000])
summ <- g1$`gate-summary`
expect("the gate summary gives the decision and the posted rule's reading",
       summ$value[summ$item == "chosen_R"] == "5000" && "posted_rule_chosen_R" %in% summ$item)
expect("the gate summary gives seconds per unit of cost", "seconds_per_cost" %in% summ$item)
expect("the warm-up of 5 minutes is applied", all(g1$`gate-times`$warmup_applied))
# 2 s per replication: 1,601.6 worker-hours, 1.56 hours per job (2.03 with
# the allowance) fits at R = 5,000; the posted rule (160 hours) does not, not
# even with the fallback.
g2 <- gate_project(gate_times(2), scenarios, gi)
s2 <- g2$`gate-summary`
expect("the runner rule keeps R = 5,000 where the posted rule stops",
       s2$value[s2$item == "chosen_R"] == "5000" && s2$value[s2$item == "posted_rule_decision"] == "stop")
# 30 s per replication: 5,000 needs 30.5 job-hours with the allowance; even
# 3,000 does not fit, nor the fallback.
g3 <- gate_project(gate_times(30), scenarios, gi)
s3 <- g3$`gate-summary`
expect("with no rung within the cap the fallback is tried, then stop", s3$value[s3$item == "decision"] == "stop" &&
         any(g3$`gate-projection`$set == "main grid without the uniform distribution"))

# ---- Result ----------------------------------------------------------------------------------------

if (length(test_failures)) {
  cat("FAILED ", length(test_failures), " of ", n_tests, " tests:\n", paste0("  ", test_failures, "\n"), sep = "")
  quit(save = "no", status = 1L)
}
cat("All", n_tests, "tests of robustness-summarize.R and the gate projection passed.\n")
