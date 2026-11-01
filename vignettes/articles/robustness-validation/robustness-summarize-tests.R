# Unit tests of the summary of the issue #116 simulation
# (robustness-summarize.R) and of the study's gate projection, on synthetic
# inputs with known answers. Run from the package root:
#
#   Rscript vignettes/articles/robustness-validation/robustness-summarize-tests.R
#
# No input comes from the study's data-generating mechanisms: the results,
# agreement rows, preflight rows, and timings below are constructed by hand
# (fixed numbers, or sin() and cos() of an index), so every expected value
# is known before the code runs. Nothing is drawn at random and the package
# is not loaded. The thresholds come from the protocol's table (exact
# binomial counts at R = 5,000 to 3,000) and its counts of judged rates.
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
fails_with <- function(expr, pattern) {
  e <- tryCatch({ force(expr); NULL }, error = function(e) e)
  inherits(e, "error") && grepl(pattern, conditionMessage(e), fixed = TRUE)
}
# The bytes sim_write_csv() would write, without writing a file.
csv_bytes <- function(x) {
  con <- rawConnection(raw(0), "wb")
  on.exit(close(con))
  utils::write.table(x, con, sep = ",", qmethod = "double", row.names = FALSE, col.names = TRUE,
                     eol = "\n", na = "NA")
  rawConnectionValue(con)
}
same_bytes <- function(a, b) {
  identical(names(a), names(b)) && all(vapply(names(a), function(n) identical(csv_bytes(a[[n]]), csv_bytes(b[[n]])), TRUE))
}

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
std_j <- loo(function(e, s) mean(e) / stats::sd(e))
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
b <- bias_rows(est, heavy_tailed = FALSE)
expect("bias and its MCSE SD/sqrt(n)", near(b$estimate[1], mean(est)) && near(b$mcse[1], stats::sd(est) / sqrt(30)))
expect("standardized bias and its jackknife MCSE equal brute-force leave-one-out",
       near(b$estimate[2], mean(est) / stats::sd(est)) && near(b$mcse[2], brute(std_j), 1e-12) &&
         !nzchar(b$mcse_note[2]))
b3 <- bias_rows(est, heavy_tailed = TRUE)
expect("under t (3 df) the standardized bias has no MCSE, with a note, and the bias keeps its MCSE",
       is.na(b3$mcse[2]) && nzchar(b3$mcse_note[2]) && near(b3$mcse[1], stats::sd(est) / sqrt(30)))

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
expect("a judged rate's events are in `count`", dec$count[rowkey(i_lim)] == 323L && dec$n_success[rowkey(i_lim)] == 5000L)

# ---- The fit-failure rule ---------------------------------------------------------------------

sid <- function(level_, pattern_, dist_, alloc_, h_ = NA, rho_ = NA) {
  scenarios$scenario[scenarios$level == level_ & scenarios$pattern == pattern_ & scenarios$dist == dist_ &
                       scenarios$alloc == alloc_ & (is.na(h_) | scenarios$h %in% h_) &
                       (is.na(rho_) | scenarios$rho %in% rho_)]
}
s_bal20 <- sid("main", "E0", "normal", "20-20-20-20", 1, 0.5)
s_t3 <- sid("main", "E0", "t3", "10-10-10-10", 1, 0.5)
ffit <- data.frame(scenario = c(s_bal20, s_bal20, s_t3, 1L), type = "failed_fit", method = c("M1", "M2", "M1", "M1"),
                   status = "any", count = c(6L, 2L, 30L, 1L), replications = 100L, stringsAsFactors = FALSE)
fr <- fit_failure_rows(ffit, scenarios)
expect("more than 5% failed fits is a material failure, more than 1% a documented limitation, 1% passes",
       identical(fr$class, c("material failure", "documented limitation", "material failure", "")) &&
         identical(fr$verdict, c("fail", "fail", "fail", "pass")))
expect("fit-failure rows carry the scenario's factors and the counts",
       identical(fr$alloc_type[1], "balanced") && fr$min_cell[1] == 20L && identical(fr$dist[3], "t3") &&
         identical(fr$stress_test, c(FALSE, FALSE, TRUE, FALSE)) && identical(fr$count, c(6L, 2L, 30L, 1L)) &&
         all(fr$total == 100L) && all(is.na(fr$n_success)))
expect("fit failures belong to the fitting function", identical(fr$`function`[1:2], c("fit_solomon_glm", "fit_solomon_ml")))
af_fit <- assign_features(rbind(dec, fr))
expect("a material fit failure at 20 per group, balanced, sets n* = 60 (literal reading)",
       identical(af_fit$n_star[["fit_solomon_glm"]], 60L) &&
         identical(af_fit$decisions$feature[nrow(dec) + 1L], "small cell (n* = 60)"))
expect("a material fit failure under t (3 df) is a stress test",
       startsWith(af_fit$decisions$feature[nrow(dec) + 3L], "stress test"))

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
expect("material failures counted per function outside the stress tests",
       identical(as.integer(af$material[c("f20", "f10", "fNA")]), c(5L, 2L, 2L)))

# ---- The Romano check ---------------------------------------------------------------------------

rom_sc <- scenarios$scenario[scenarios$level == "main" & scenarios$pattern == "E0" & scenarios$h > 1 &
                               scenarios$alloc_type %in% c("treated smaller", "treated larger")]
expect("40 scenarios in each pairing", length(rom_sc) == 80L)
rom_perf <- function(romano) {
  data.frame(scenario = rep(romano$scenario, each = 4), method = "M4", contrast = rep(contrasts_treat, nrow(romano)),
             measure = "rejection", estimate = rep(romano$mean_q, each = 4), stringsAsFactors = FALSE)
}
romano <- data.frame(scenario = rom_sc, mean_q = ifelse(scenarios$alloc_type[rom_sc] == "treated smaller", 0.06, 0.04),
                     var_q = 0.01, reps_q = 5000L, failed = FALSE)
rr <- romano_rows(romano, rom_perf(romano), scenarios)
expect("Romano MCSE = (1/40) sqrt(sum var/R)", near(rr$mcse, rep(sqrt(40 * 0.01 / 5000) / 40, 2)))
expect("Romano means", near(rr$estimate, c(0.06, 0.04)) && identical(rr$total, c(160L, 160L)))
expect("Romano check as predicted", identical(rr$verdict, c("as predicted", "as predicted")))
romano$mean_q[scenarios$alloc_type[rom_sc] == "treated smaller"] <- 0.05 + 1.9 * sqrt(40 * 0.01 / 5000) / 40
expect("Romano check contradicted within 2 MCSE", romano_rows(romano, rom_perf(romano), scenarios)$verdict[1] == "contradicted")
# With a failed M4 contrast the mean of the 160 rates is used, not the mean
# of the replications' proportions; without one, the two must agree.
rp <- rom_perf(romano)
rp$estimate[1] <- 0.5
romano$failed[1] <- TRUE
pairing1 <- scenarios$alloc_type[romano$scenario[1]]
same_pairing <- romano$scenario[scenarios$alloc_type[romano$scenario] == pairing1]
expect("the Romano mean is the mean of the 160 rates",
       near(romano_rows(romano, rp, scenarios)$estimate[match(pairing1, c("treated smaller", "treated larger"))],
            mean(rp$estimate[rp$scenario %in% same_pairing])))
romano$failed[1] <- FALSE
expect("the Romano check stops if the two means disagree without a failed fit",
       fails_with(romano_rows(romano, rp, scenarios), "disagree"))

# ---- Synthetic inputs: results, agreement rows, and the preflight ---------------------------------

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
# Agreement rows as replicate() writes them (one per quantity and call;
# here one contrast per call), every comparison passing with a difference
# of 1e-15.
synthetic_agreement <- function(sc, reps) {
  plan <- agreement_plan()
  rows <- list()
  for (i in seq_len(nrow(sc))) {
    e0 <- sc$pattern[i] == "E0"
    for (r in seq_len(reps)) {
      use <- (startsWith(plan$check, "1") & r <= check1_reps) | (!startsWith(plan$check, "1") & e0 & r <= check34_reps)
      if (!any(use)) next
      p <- plan[use, , drop = FALSE]
      rows[[length(rows) + 1L]] <- data.frame(
        scenario = sc$scenario[i], rep = r, check = p$check, method = p$method, contrast = contrasts_all[1],
        statistic = NA_character_, delta = NA_real_, bound = NA_character_, quantity = p$quantity,
        package = 1, direct = 1, abs_diff = 1e-15, tolerance = p$tolerance, passed = TRUE,
        warnings_muffled = 0L, warnings_122 = 0L, warning_classes = "", stringsAsFactors = FALSE)
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
# plan/preflight.csv as the contract reads it (every column text, NA as
# "NA"): the two stops and check 2's 160 rows, all passed.
synthetic_preflight <- function(passed = TRUE) {
  w <- expand.grid(contrast = contrasts_treat, method = c("M3", "M4"), dataset = seq_len(check2_sets),
                   KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  ties <- data.frame(check = paste("agreement check 2 (ties): data set", w$dataset, w$method, w$contrast),
                     passed = passed, detail = "", dataset = w$dataset, method = w$method, contrast = w$contrast,
                     statistic = ifelse(w$method == "M3", "studentized", "difference"), max_abs_diff = 1e-15,
                     tie_labelled = 2L, counted_package = 2L, counted_direct = 2L, n_tied_stored = 2L,
                     warnings_muffled = 0L, warnings_122 = 0L, stringsAsFactors = FALSE)
  stops <- ties[1:2, ]
  stops$check <- c("the default inference", "the seven contrasts")
  stops$passed <- TRUE
  stops[, c("dataset", "method", "contrast", "statistic", "max_abs_diff", "tie_labelled", "counted_package",
            "counted_direct", "n_tied_stored", "warnings_muffled", "warnings_122")] <- NA
  x <- rbind(stops, ties)
  x[] <- lapply(x, function(v) ifelse(is.na(v), "NA", as.character(v)))
  rownames(x) <- NULL
  x
}
pre <- synthetic_preflight()
info <- function(purpose) c(purpose = purpose, protocol = "TRUE", design = if (purpose == "gate") "gate" else "run")

small <- scenarios[scenarios$scenario %in% c(1, 2, 271, 541, 559), ]
c5_ok <- check5_rows(synthetic_results(small, 3, 1e-7), small)
c5_bad <- check5_rows(synthetic_results(small, 3, 2e-6), small)
expect("check 5 passes within 1e-6", all(c5_ok$passed) && nrow(c5_ok) == 4L * 7L && all(c5_ok$comparisons == 3L))
expect("check 5 fails beyond 1e-6", !any(c5_bad$passed))

# ---- The purposes: what a test or gate summary may show ----------------------------------------------

res_a <- synthetic_results(small, 4)
agr_a <- synthetic_agreement(small, 4)
# The largest differences sit in scenario 541, replication 1, which both
# inputs below keep; #122's warning is muffled in scenario 541's check 1.
top <- agr_a$scenario == 541 & agr_a$rep == 1
agr_a$abs_diff[top] <- 5e-14
agr_a$warnings_122[top & agr_a$check == "1. permutation tests" & agr_a$method == "M4" & agr_a$quantity == "p_value"] <- 1L
agr_a$warning_classes[top & agr_a$check == "1. permutation tests" & agr_a$method == "M4" & agr_a$quantity == "p_value"] <- warning_class_122
# Input B differs only in fit statuses (M2 not converged in every
# replication of scenario 1; M1 failed in replication 1 of scenario 2), in
# the agreement rows that replicate() then skips, and in the warnings
# muffled in checks 1, 3, and 4.
res_b <- res_a
res_b$status[res_b$scenario == 1 & res_b$method == "M2"] <- "nonconverged"
gone <- res_b$scenario == 2 & res_b$rep == 1 & res_b$method == "M1"
res_b$status[gone] <- "error"
res_b$estimate[gone] <- res_b$std_error[gone] <- res_b$df[gone] <- NA_real_
agr_b <- agr_a[!(agr_a$scenario == 1 & agr_a$method == "M2") &
                 !(agr_a$scenario == 2 & agr_a$rep == 1 & (agr_a$check == "1. permutation tests" | agr_a$method == "M1")), ]
w4 <- agr_b$check == "4. equivariance" & agr_b$scenario == 541 & agr_b$quantity == "estimate_shift"
agr_b$warnings_muffled[w4] <- 3L
agr_b$warning_classes[w4] <- "simpleWarning"
w1 <- agr_b$check == "1. permutation tests" & agr_b$scenario == 559 & agr_b$quantity == "p_value"
agr_b$warnings_muffled[w1] <- 2L
agr_b$warning_classes[w1] <- "simpleWarning"
for (purpose in c("test", "gate")) {
  sa <- summarize_robustness(res_a, agr_a, NULL, small, info(purpose), pre)
  sb <- summarize_robustness(res_b, agr_b, NULL, small, info(purpose), pre)
  expect(paste("a", purpose, "summary has agreement and run information only"),
         setequal(names(sa), c("agreement", "run-information")))
  expect(paste("a", purpose, "summary is byte-identical whatever the fit statuses and the checks' warnings"),
         same_bytes(sa, sb))
  expect(paste("a", purpose, "summary shows no check-5 difference"),
         all(is.na(sa$agreement$max_abs_diff[startsWith(sa$agreement$check, "5")])))
  expect(paste("a", purpose, "summary has no count of comparisons and no #122 count"),
         !any(c("comparisons", "warnings_122") %in% names(sa$agreement)))
  expect(paste("a", purpose, "summary counts warnings in check 2 only"),
         all(is.na(sa$agreement$warnings_muffled[!startsWith(sa$agreement$check, "2")])) &&
           all(!is.na(sa$agreement$warnings_muffled[startsWith(sa$agreement$check, "2")])))
  expect(paste("a", purpose, "summary has one row per check, method, and quantity, plus checks 2 and 5"),
         nrow(sa$agreement) == 8L + nrow(agreement_plan()) + 7L)
  expect(paste("a", purpose, "summary says whether #122's warning was muffled"),
         all(grepl("#122 warnings muffled: yes", sa$agreement$detail[sa$agreement$check == "1. permutation tests" &
                                                                      sa$agreement$method == "M4"])) &&
           all(grepl("#122 warnings muffled: none", sa$agreement$detail[sa$agreement$check == "1. permutation tests" &
                                                                          sa$agreement$method == "M3"])))
  expect(paste("a", purpose, "summary gives the largest difference of the checks"),
         near(sa$agreement$max_abs_diff[sa$agreement$check == "4. equivariance"], rep(5e-14, 6), 1e-20))
}
bad <- agr_a
bad$passed[bad$scenario == 2 & bad$check == "3. TOST" & bad$quantity == "p_upper"] <- FALSE
expect("a test summary stops on a failed agreement check and names the scenario",
       fails_with(summarize_robustness(res_a, bad, NULL, small, info("test"), pre), "failed in scenarios 2"))
gate_bad <- summarize_robustness(res_a, bad, NULL, small, info("gate"), pre)
expect("a gate summary reports a failed agreement check", any(gate_bad$agreement$passed %in% FALSE))
expect("a summary without the preflight stops, whatever the purpose",
       all(vapply(c("test", "gate", "study"), function(p) {
         fails_with(summarize_robustness(res_a, agr_a, NULL, small, info(p), NULL), "preflight")
       }, TRUE)))
expect("a summary with an incomplete preflight stops",
       fails_with(summarize_robustness(res_a, agr_a, NULL, small, info("gate"), pre[-5, ]), "160 rows"))
na_status <- res_a
na_status$status[10] <- NA
expect("an NA status stops the summary and names the scenario and replication",
       fails_with(summarize_robustness(na_status, agr_a, NULL, small, info("test"), pre), "in scenario 1, replication 1"))
odd_status <- res_a
odd_status$status[10] <- "OK"
expect("an unknown status stops the summary", fails_with(summarize_robustness(odd_status, agr_a, NULL, small, info("gate"), pre), "\"OK\""))
stray <- agr_a[1, ]
stray$check <- "3. TOST"
stray$scenario <- 271L
expect("an agreement row outside the plan (check 3 in ET) stops the summary",
       fails_with(summarize_robustness(res_a, rbind(agr_a, stray), NULL, small, info("test"), pre), "not one of the planned"))
expect("a study summary refuses a partial table",
       inherits(tryCatch(summarize_robustness(res_a, agr_a, NULL, small, info("study"), pre),
                         error = function(e) e), "error"))
ri <- run_information(14576L)
expect("run information names no harness item",
       !any(ri$item %in% c("study", "issue", "purpose", "package", "package_version", "r_version", "base_seed",
                           "reps", "rng_kinds", "workers_per_chunk", "started_utc", "finished_utc")))
expect("broom is \"not installed\" when absent", identical(version_of("a.package.that.is.not.installed"), "not installed"))
expect("the permutations are recorded", ri$value[ri$item == "permutations"] == "999")

# The study run's agreement rows, on the same inputs (agreement_summary()
# directly, because a study summary needs the full table).
ag_a <- agreement_summary(agr_a, res_a, small, pre, study_run = TRUE)
ag_b <- agreement_summary(agr_b, res_b, small, pre, study_run = TRUE)
grid_rows <- function(a) a[a$check %in% c("1. permutation tests", "3. TOST", "4. equivariance"), ]
expect("the study run's agreement rows are the plan's grid, whatever the fits",
       nrow(grid_rows(ag_a)) == nrow(small) * 8L + sum(small$pattern == "E0") * 12L &&
         identical(grid_rows(ag_a)[, c("check", "scenario", "method", "quantity")],
                   grid_rows(ag_b)[, c("check", "scenario", "method", "quantity")]))
nc <- ag_b[ag_b$check == "4. equivariance" & ag_b$scenario == 1 & ag_b$method == "M2", ]
expect("a group whose fits all failed is \"not compared\", with no verdict",
       all(nc$comparisons == 0L) && all(is.na(nc$passed)) && all(nc$detail == "not compared"))
expect("the study run counts comparisons and the warnings of every check",
       all(ag_a$comparisons[ag_a$check == "4. equivariance" & ag_a$scenario == 1] == 4L) &&
         sum(ag_b$warnings_muffled[ag_b$check == "4. equivariance"]) == 24L &&
         sum(ag_a$warnings_122[ag_a$check == "1. permutation tests"]) == 1L)
c5b <- ag_b[ag_b$check == "5. M2 equals M1", ]
expect("check 5 is \"not compared\" where M2 failed in every replication, and counts the others",
       all(is.na(c5b$passed[c5b$scenario == 1])) && all(c5b$comparisons[c5b$scenario == 2] == 3L) &&
         all(c5b$passed[c5b$scenario == 2]))

# ---- A whole study summary on synthetic results ---------------------------------------------------------

full_res <- synthetic_results(scenarios, 4)
# Injected failed fits: M1 in replication 1 of a balanced normal scenario
# with 20 per group (25% of 4: material), and of a t (3 df) scenario (a
# stress test); M4 in replication 2 of the first Romano scenario.
inject <- function(res, s, r, m, status) {
  i <- res$scenario == s & res$rep == r & res$method == m
  res$status[i] <- status
  if (m %in% c("M1", "M2")) res$estimate[i] <- res$std_error[i] <- res$df[i] <- NA_real_
  else res$p_value[i] <- res$statistic[i] <- NA_real_
  res
}
full_res <- inject(full_res, s_bal20, 1, "M1", "error")
full_res <- inject(full_res, s_t3, 1, "M1", "error")
rom1 <- min(rom_sc[scenarios$alloc_type[rom_sc] == "treated smaller"])
full_res <- inject(full_res, rom1, 2, "M4", "undefined")
full_warn <- data.frame(scenario = c(1L, 1L), rep = c(2L, 3L), type = "warning", class = "simpleWarning",
                        message = c("first message", "second message"), count = c(3L, 1L), stringsAsFactors = FALSE)
full <- summarize_robustness(full_res, synthetic_agreement(scenarios[scenarios$scenario %in% c(1, 271), ], 3),
                             full_warn, scenarios, info("study"), pre)
pf <- full$performance
dc <- full$decisions
expect("a study summary has every file", setequal(names(full), c("agreement", "run-information", "performance",
                                                                 "decisions", "failures", "scenarios")))
expect("the scenario table is written with the results", identical(full$scenarios, scenarios))
expect("every performance row has an MCSE above 0 or a note",
       all((is.finite(pf$mcse) & pf$mcse > 0) | nzchar(pf$mcse_note)))
expect("performance rows of every scenario", setequal(unique(pf$scenario), scenarios$scenario))
num <- cbind(is.finite(dc$mcse), is.finite(dc$threshold), is.finite(dc$p_value))
expect("every decision row has an MCSE, threshold, or p-value", all(rowSums(num) > 0))
expect("decisions hold the 14,576 judged rates",
       sum(dc$rule %in% c("1. error control", "2. one-sided error (TOST)", "3. bias")) == 14576L)
expect("the run information gives C", full$`run-information`$value[full$`run-information`$item == "family_size_C"] == "14576")
fw <- full$failures[full$failures$type == "warning", ]
expect("warnings by class, with every distinct message",
       nrow(fw) == 1L && fw$count == 4 && fw$n_messages == 2L && fw$messages == "first message | second message" &&
         near(fw$proportion, 2 / 4))
expect("every scenario and method has a failed-fit row",
       sum(full$failures$type == "failed_fit" & full$failures$status == "any") == 306L * 4L + 270L * 2L)
expect("M2's largest difference from M1 is reported, with the replications compared",
       any(pf$measure == "max_abs_difference_from_M1") &&
         all(pf$n_success[pf$measure == "max_abs_difference_from_M1" & pf$scenario == s_bal20] == 3L))
expect("judged and described rates are labelled", all(c("judged", "described", "descriptive") %in% pf$role))
expect("M2 in the slopes level is described",
       all(pf$role[pf$method == "M2" & pf$measure == "rejection" &
                     pf$scenario %in% scenarios$scenario[scenarios$level == "slopes"]] == "described"))
expect("the standardized bias carries an MCSE outside t (3 df)",
       all(is.finite(pf$mcse[pf$measure == "standardized_bias" & !pf$scenario %in% scenarios$scenario[scenarios$dist == "t3"]])))
ffr <- dc[dc$rule == "fit failures", ]
expect("fit-failure rows: count and total, factors, and class",
       nrow(ffr) == 306L * 4L + 270L * 2L &&
         ffr$count[ffr$scenario == s_bal20 & ffr$method == "M1"] == 1L &&
         ffr$total[ffr$scenario == s_bal20 & ffr$method == "M1"] == 4L &&
         ffr$class[ffr$scenario == s_bal20 & ffr$method == "M1"] == "material failure" &&
         ffr$alloc_type[ffr$scenario == s_bal20 & ffr$method == "M1"] == "balanced" &&
         isTRUE(ffr$stress_test[ffr$scenario == s_t3 & ffr$method == "M1"]))
nstar <- dc[dc$rule == "(a) small-cell warning: n*", ]
expect("the material fit failure sets fit_solomon_glm's n* to 60; the stress test does not count",
       nstar$estimate[nstar$`function` == "fit_solomon_glm"] == 60 && nstar$count[nstar$`function` == "fit_solomon_glm"] == 1L)
expect("(b) counts M1's material fit failure; rule 2 of M1 fits is reported on its own row",
       dc$estimate[dc$rule == "(b) changed default"] == 1 &&
         dc$verdict[dc$rule == "(b) changed default: rule 2 of fit_solomon_glm fits"] == "reported")
expect("(c) counts M4's material fit failure only in equal-arm allocations",
       dc$estimate[dc$rule == "(c) difference statistic"] == 0)
rom_row <- dc[dc$rule == "Romano check (treated smaller arms, h > 1)", ]
rates <- pf$estimate[pf$method == "M4" & pf$measure == "rejection" &
                       pf$scenario %in% rom_sc[scenarios$alloc_type[rom_sc] == "treated smaller"]]
expect("the Romano mean is the mean of the 160 M4 rates", length(rates) == 160L && near(rom_row$estimate, mean(rates)))
ag4 <- dc[dc$rule == "agreement check 4. equivariance", ]
expect("agreement decisions: one row per check, method, and quantity, with its own tolerance",
       nrow(ag4) == 6L && all(ag4$threshold[ag4$method == "M1"] == 1e-10) && all(ag4$threshold[ag4$method == "M2"] == 1e-6))
expect("agreement decisions count the comparisons",
       dc$count[dc$rule == "agreement check 1. permutation tests" & dc$method == "M3" &
                  dc$measure == "largest absolute difference: p_value"] == 6L)
pat <- dc[dc$rule == "patterns, 1. error control" & dc$method == "M1", ]
by_f <- function(f) pat[pat$measure == paste("share of rates with p < .05 (unadjusted), by", f), ]
expect("patterns: all rates and each factor's levels, whose rates add up",
       pat$total[pat$measure == "share of rates with p < .05 (unadjusted), all rates"] == 2142L &&
         nrow(by_f("distribution")) == 5L && nrow(by_f("allocation")) == 9L && nrow(by_f("contrast")) == 7L &&
         nrow(by_f("spread")) == 4L && nrow(by_f("correlation")) == 3L &&
         all(vapply(names(pattern_factors), function(f) sum(by_f(f)$total) == 2142L, TRUE)) &&
         setequal(by_f("distribution")$dist, dists) && all(by_f("allocation")$threshold == 0.05))

# ---- The gate's projection ---------------------------------------------------------------------------

# One task of 25 replications per scenario (task 1, which carries the
# agreement checks and their cost): 0.1 s per replication and 10 s of
# agreement checks each; half the tasks start 10 minutes into the job.
gate_times <- function(seconds_per_rep) {
  start <- as.POSIXct("2026-11-01 00:00:00", tz = "UTC")
  n <- nrow(scenarios)
  offset <- ifelse(seq_len(n) %% 2 == 0, 600, 0)
  data.frame(task_id = sprintf("s%04d-t0001", seq_len(n)), scenario = scenarios$scenario, reps = 25L,
             started = sim_utc(start + offset), finished = sim_utc(start + offset + 1),
             job_started = sim_utc(start), elapsed = 25 * seconds_per_rep + 10, cpu = 1,
             replication_seconds = 25 * seconds_per_rep, agreement_seconds = 10,
             cost = vapply(seq_len(n), function(i) {
               cost_robustness(scenarios[i, ], 25L) + agreement_cost_robustness(scenarios[i, ])
             }, 0),
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
# Seconds per unit of cost from replication time and the replications' cost
# alone (the tasks after the warm-up: the even scenarios), and the agreement
# seconds per unit of agreement cost (all tasks).
even <- scenarios$scenario %% 2 == 0
rep_cost <- vapply(seq_len(nrow(scenarios)), function(i) cost_robustness(scenarios[i, ], 25L), 0)
agree_cost <- vapply(seq_len(nrow(scenarios)), function(i) agreement_cost_robustness(scenarios[i, ]), 0)
expect("seconds per unit of cost leave out the agreement checks and their cost",
       near(as.numeric(summ$value[summ$item == "seconds_per_cost"]), sum(even) * 2.5 / sum(rep_cost[even]), 1e-12))
expect("agreement seconds per unit of agreement cost",
       near(as.numeric(summ$value[summ$item == "agreement_seconds_per_cost"]), 576 * 10 / sum(agree_cost), 1e-12))
expect("task 1 of an E0 scenario carries check 1 on 3 replications and checks 3 and 4 on 20",
       agreement_cost_robustness(scenarios[1, ]) == 3 * (9000 + 9 * 40) + 20 * (2000 + 40) &&
         agreement_cost_robustness(scenarios[271, ]) == 3 * (9000 + 9 * 40))
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
