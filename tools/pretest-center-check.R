# Checks of the decision record on issue #135, "The pretest center when
# posttests are missing", run with the package's two centers. Run from the
# package root:
#
#   Rscript tools/pretest-center-check.R [out.csv] [parts]
#
# `out.csv` defaults to tools/pretest-center-check.csv; `parts` is a
# comma-separated list of the parts below, by default "1,2,3,4,5,6,8" (part
# 7 is the validation named under NOT HERE). The script exits with status 1
# if any check fails.
#
# A check of derivations, not a simulation study. The record (section 3)
# compares the two centers on data with no effects of any kind, and its
# conclusions rest on the direction and size of the biases, which are large
# against their Monte Carlo error. Nothing here was pre-registered.
#
# THE CENTERS
#   A, "analyzed": the default, the mean pretest of the pretested
#      participants in the model, read from `fit$effects`.
#   B, "observed": the mean of every observed pretest among the pretested
#      participants, read from `fit$pretest_sensitivity` of the same fit.
#      Every data set is also fitted with pretest_center = "observed", and
#      its pretest rows must equal those of `pretest_sensitivity` to 1e-10:
#      a check of the computation without a refit on every data set.
#
# THE RECORD'S SCRIPTS. The scripts that produced the record's table
# (JUDGE/judge_check.R and JUDGE2/judge2_check.R) were in a scratch folder
# that no longer exists. This script regenerates the scenarios from the
# record's description: 200 participants per group, a pretest-posttest
# correlation of .5, no effects, HC3 standard errors, and loss with
# probability plogis(-1 - 1.5 z), about 33%, with set.seed(20261015) before
# the draw of the loss rate p0 and set.seed(1350 + k) or set.seed(2026 + k)
# before scenario k of each block, as the record gives them. The order of
# the random draws within a data set is the description's (baseline, then
# posttest, then the uniforms of the loss), but the record's data sets are
# not reproduced exactly, so each value is compared with the record's within
# Monte Carlo error (see PASS RULES).
#
# PART 1: the record's table (section 3). R data sets per scenario; the
# pretest main effect (true value 0) at each center, with the share of data
# sets in which its 95% interval excludes 0.
#   Scenario  R    Loss case                                         Center whose
#                                                                    case it is
#   S1       500  by baseline, plogis(-1 - 1.5 x), all four groups   A
#   L1       300  by the attribute T that the pretest measures with  A
#                 reliability .5, plogis(-1 - 1.5 T), all groups
#   S2       500  by baseline, pretested groups only                 B
#   S2e      500  pretested by baseline; unpretested at random at    B
#                 the same rate p0
#   L2       300  by the attribute T, pretested groups only          neither
#   R3       300  by baseline in all groups, the pretested at a rate  neither
#                 of .41 (intercept solved numerically), the
#                 unpretested at .33
#   R4       300  by baseline in all groups, plus 15% at random in    A
#                 the pretested groups
#   S3       500  by baseline in the control arms, at random (rate p0) the main
#                 in the treated arms, alike in both pretest           effect: A
#                 conditions
#   P1       300  no posttest lost; pretests of pretested             neither
#                 participants missing by baseline
# Data, for x the baseline (the pretest of the pretested groups):
#   x ~ N(0, 1), y = 0.5 x + e, e ~ N(0, 0.75); in L1 and L2, T ~ N(0, 1),
#   x = sqrt(.5) T + sqrt(.5) e1, y = sqrt(.5) T + sqrt(.5) e2, so that the
#   pretest-posttest correlation is .5 and the pretest's reliability .5.
#   Groups 1 to 4 of 200 each: pretested treated, pretested control,
#   unpretested treated, unpretested control.
# Also reported: the fit without y_pre (L, R, and P scenarios), and in S3 the
# two conditional pretest effects, Pretest x Treatment, and the raw
# comparison of the control arms' posttest means (Welch t test).
#
# PART 2: the record's derivations, on one large sample per case
# (2,000,000 per group), as algebra on the completers' means and the
# pooled within-group slope of the pretested completers, without the
# package's fitting functions:
#   D1 (critic's derive2.R): loss by the posttest itself, plogis(-1 - 1.5 y),
#      in all groups, with a true pretest effect of 0, 0.3, and 0.5. Record:
#      A -0.001 / 0.234 / 0.389; B -0.073 / 0.174 / 0.337.
#   D2 (derive.R, R3 analogue): A +0.040, B -0.180.
#   D3 (derive.R, R4 analogue): A 0.000, B -0.178.
# The record's S3 case with unequal arm rates (main effect -0.0183) is not
# rerun: its arm rates are not given in the record.
#
# PART 3: the mai2020 block of the record's section 4, value by value, and
# the gender line from check_solomon_missing(covariates = "gender").
#
# PART 4: the full-information maximum likelihood comparison of section 2
# (lavaan; skipped when it is not installed). The record's data set is not
# reproducible; this part fits one data set of its own (2,000 per group,
# set.seed(135135), true pretest effects of 0.2 among controls and 0.4 among
# treated participants, pretested losses by the pretest, unpretested losses
# at random) by normal-theory FIML with the pretest, the slope equal in the
# two pretested groups and their pretest means equal, and compares its
# pretest effects with those at center B. Loss case: B's.
#
# PART 5: the repeated-measures caution of ?fit_solomon_mmrm (the record's
# section 6, the A advocate's check, which the record did not rerun). Its
# script is lost too; these scenarios are this script's. Two posttest
# occasions; x ~ N(0, 1); y_t = 0.5 x + u_t, (u_1, u_2) normal with
# variances .75 and correlation .5; no effects; 1,500 per group, 40 data
# sets per scenario, Satterthwaite degrees of freedom, set.seed(1354) once
# before the scenarios, in this order:
#   M1  dropout before occasion 2 by the occasion-1 posttest,
#       plogis(-1 - 1.5 y_1), all groups: missing at random; A = B
#       (everyone has a posttest).
#   M4  no posttest at either occasion, pretested groups by baseline,
#       unpretested at random (rate p0): B's case.
#   M2  dropout before occasion 2 by baseline, plogis(-1 - 1.5 x), all
#       groups: dependence on baseline beyond the earlier posttest, which
#       the unpretested groups' model cannot see; A = B.
#   M3  no posttest at either occasion, by baseline, all groups: A's case.
#
# PART 6: the article example ("Missing and Repeated Posttests"), whose
# data the article simulates with fixed seeds: the pretest main effect at
# each center (the record: 2.96 at A and 1.969 at B, for a true 2).
#
# PART 8: the rule by which check_solomon_missing() says that the pretested
# conditions' selection shifts differ (the record's section 5, item 3, whose
# rule the record leaves open): a Wald test of equal shifts at .05, with
# standard errors from the shifts' influence functions. Not a check of the
# record but of the package's rule, added after review of the pull request
# found that the earlier rule (a difference at two printed decimals) fired
# in nearly every data set. Four-group data as in part 1 (x ~ N(0, 1),
# y = 0.5 x + e), 1,000 data sets per setting, set.seed(1358) once:
#   S1    loss by baseline, plogis(-1 - 1.5 x), in all four groups: equal
#         selection, so the rule should fire in about 5% of data sets;
#   MCAR  33% lost completely at random in all groups: likewise;
#   S3    loss by baseline in the control arms and at random (rate .329)
#         in the treated arms: the arms select differently, so this is the
#         rule's power; no pass rule.
# Each at 200 and 50 per group. Pass rule: under S1 and MCAR, |rate - .05|
# <= z sqrt(.05 x .95 / 1000), z Bonferroni over the four settings; and in
# every data set the test statistic is the same, to 1e-8, with the pretest
# in other units (50 + 10 x).
#
# NOT HERE: the branch's validation of the "observed" center under its own
# assumption (branch fix/135-pretest-center, tools/pretest-center-check.R at
# 553ca52: designs A to F, about 1 hour 45 minutes on 8 workers). Its result
# is in that file's header; it is not rerun here, because it needs an idle
# machine.
#
# PASS RULES. A regenerated value is compared with the record's printed
# value. Means: |difference| <= z sqrt(2) MCSE, with MCSE the Monte Carlo
# standard error of this run's mean, since both runs carry Monte Carlo
# error. Rejection rates: |difference| <= z sqrt(2 p (1 - p) / R), with p
# the record's rate. z = qnorm(1 - .05 / (2 K)) for the K comparisons of a
# part (Bonferroni at a familywise .05). Each data set's equality of the
# observed rows with the refit at "observed" must hold to 1e-10. The
# derivations of part 2 must match to 0.006: with about 1,340,000
# completers per group, a difference of completers' means has a standard
# error of about 0.0012 on each sample, the record's and this one, and
# 0.006 is 3.5 sqrt(2) times that. Part 3 must match to the record's printed
# precision; part 4, FIML and B, to 0.002; part 6 exactly to the record's
# precision. Part 5 has no value from the record except the M2 occasion-2
# bias, -0.076 at both centers.
#
# RESULT: see the summary this script prints and the csv it writes. The
# result of the run made for the pull request is added below this line
# after the run.
#
# RESULT (added after the run; the header above is as it was before it)
# October 31, 2026, R 4.6.1, one process, 8.9 minutes; every value is in
# tools/pretest-center-check.csv. All checks passed.
#   Part 1: 54 comparisons with the record (z = 3.312); every one within
#   its tolerance. The largest departures, in Monte Carlo standard errors
#   of a difference: S3, the rejection rate at center A, .062 against .038
#   (2.0); L1, center A, 0.014 against 0.003 (1.8). In every one of
#   the 3,500 data sets, the rows at "observed" from pretest_sensitivity
#   equalled the refit at pretest_center = "observed" to within 1e-13.
#   S1: A 0.007 (rejection .058), B -0.172 (.508); S2: A 0.181 (.638), B
#   0.003 (.038); S2e: A 0.180 (.554), B 0.003 (.036); L1: A 0.014 (.047),
#   B -0.096 (.217); L2: A 0.257 (.947), B 0.147 (.497); R3: A 0.038
#   (.063), B -0.178 (.507); R4: A -0.005 (.057), B -0.184 (.533); S3: A
#   -0.007 (.062), B -0.095 (.182), by condition A -0.098 / 0.084, B
#   -0.185 / -0.004, Pretest x Treatment 0.181 (.176), raw control
#   contrast -0.010 (.042); P1: A = B 0.182 (.627), without y_pre 0.006
#   (.047).
#   Part 2: D1 A 0.000 / 0.232 / 0.392, B -0.071 / 0.173 / 0.340; D2 A
#   0.042, B -0.176; D3 A 0.001, B -0.177.
#   Part 3: every value of the record's block and gender line, to its
#   printed precision. The gender coefficient (-0.033, p = .761) is that of
#   the analysis without the pretest, with model-based standard errors;
#   with the pretest and HC3 it is -0.033, p = .810.
#   Part 4: A 0.4043 and 0.6300; B 0.2287 and 0.4544; FIML 0.2286 and
#   0.4543 (control and treated).
#   Part 5 (pretest main effect, mean over 40 data sets, Monte Carlo
#   standard error about 0.005): M1 A = B -0.005 and -0.002 at occasions 1
#   and 2; M4 A 0.184 and 0.181, B 0.004 and -0.000; M2 A = B -0.000 and
#   -0.067; M3 A 0.002 and 0.002, B -0.177 and -0.179.
#   Part 6: A 2.9583, B 1.9694.

args <- commandArgs(trailingOnly = TRUE)
out_csv <- if (length(args) >= 1L) args[1] else file.path("tools", "pretest-center-check.csv")
parts <- if (length(args) >= 2L) as.integer(strsplit(args[2], ",")[[1]]) else c(1:6, 8L)

suppressMessages(pkgload::load_all(".", quiet = TRUE))

rows <- list()
failures <- character()
add_row <- function(...) rows[[length(rows) + 1L]] <<- data.frame(..., stringsAsFactors = FALSE)
fail <- function(what) {
  failures <<- c(failures, what)
  message("FAIL: ", what)
}

main_label <- "Pretest main effect"
N <- 200L
lgt <- stats::plogis

# ---- part 1 ---------------------------------------------------------------------

part1 <- function() {
  cat("\n== Part 1: the record's table ==\n")
  set.seed(20261015)
  p0 <- mean(lgt(-1 - 1.5 * stats::rnorm(2e6)))
  cat(sprintf("p0 = %.3f\n", p0))
  # The intercept at which the expected loss rate is .41 (R3).
  a41 <- stats::uniroot(function(a) {
    stats::integrate(function(z) lgt(a - 1.5 * z) * stats::dnorm(z), -Inf, Inf)$value - 0.41
  }, c(-3, 3), tol = 1e-10)$root

  g <- rep(1:4, each = N)
  treat <- as.integer(g %in% c(1L, 3L))
  pretested <- as.integer(g <= 2L)
  pre_rows <- pretested == 1L

  gen <- function(s) {
    n <- 4L * N
    if (s %in% c("L1", "L2")) {
      t_attr <- stats::rnorm(n)
      x <- sqrt(0.5) * t_attr + sqrt(0.5) * stats::rnorm(n)
      y <- sqrt(0.5) * t_attr + sqrt(0.5) * stats::rnorm(n)
      driver <- t_attr
    } else {
      x <- stats::rnorm(n)
      y <- 0.5 * x + stats::rnorm(n, 0, sqrt(0.75))
      driver <- x
    }
    u <- stats::runif(n)
    by_base <- u < lgt(-1 - 1.5 * driver)
    lost <- switch(
      s,
      S1 = by_base, L1 = by_base,
      S2 = by_base & pre_rows, L2 = by_base & pre_rows,
      S2e = ifelse(pre_rows, by_base, u < p0),
      R3 = ifelse(pre_rows, u < lgt(a41 - 1.5 * x), by_base),
      R4 = by_base | (pre_rows & stats::runif(n) < 0.15),
      S3 = ifelse(treat == 0L, by_base, u < p0),
      P1 = rep(FALSE, n)
    )
    y_pre <- ifelse(pre_rows, x, NA_real_)
    if (s == "P1") y_pre[pre_rows & by_base] <- NA_real_
    data.frame(y_post = ifelse(lost, NA_real_, y), treat = treat, pretested = pretested,
               y_pre = y_pre)
  }

  record <- list(
    S1 = c(A = 0.001, A_rej = 0.052, B = -0.178, B_rej = 0.538),
    L1 = c(A = 0.003, A_rej = 0.040, B = -0.106, B_rej = 0.247, raw = 0.003, raw_rej = 0.040),
    S2 = c(A = 0.179, A_rej = 0.626, B = 0.002, B_rej = 0.050),
    S2e = c(A = 0.181, A_rej = 0.570, B = 0.004, B_rej = 0.050),
    L2 = c(A = 0.252, A_rej = 0.913, B = 0.142, B_rej = 0.470, raw = 0.252, raw_rej = 0.913),
    R3 = c(A = 0.042, A_rej = 0.093, B = -0.178, B_rej = 0.503, raw = 0.042, raw_rej = 0.093),
    R4 = c(A = -0.000, A_rej = 0.063, B = -0.181, B_rej = 0.500, raw = -0.000, raw_rej = 0.063),
    S3 = c(A = -0.002, A_rej = 0.038, B = -0.090, B_rej = 0.164,
           A_control = -0.094, A_treated = 0.090, B_control = -0.182, B_treated = 0.002,
           pxt = 0.184, pxt_rej = 0.184, raw_control = -0.003, raw_control_rej = 0.050),
    P1 = c(A = 0.184, A_rej = 0.640, B = 0.184, B_rej = 0.640, raw = 0.004, raw_rej = 0.037)
  )
  seeds <- c(S1 = 1351, S2 = 1352, S2e = 1353, S3 = 1354,
             L1 = 2027, L2 = 2028, R3 = 2029, R4 = 2030, P1 = 2031)
  reps <- c(S1 = 500, S2 = 500, S2e = 500, S3 = 500, L1 = 300, L2 = 300, R3 = 300, R4 = 300,
            P1 = 300)
  order <- c("S1", "S2", "S2e", "S3", "L1", "L2", "R3", "R4", "P1")

  results <- list()
  for (s in order) {
    set.seed(seeds[[s]])
    R <- reps[[s]]
    keep <- matrix(NA_real_, R, 14L, dimnames = list(NULL, c(
      "A", "A_rej", "B", "B_rej", "raw", "raw_rej", "A_control", "A_treated", "B_control",
      "B_treated", "pxt", "pxt_rej", "raw_control", "raw_control_rej")))
    mismatch <- 0
    lost_pre <- lost_un <- numeric(R)
    for (r in seq_len(R)) {
      d <- gen(s)
      lost_pre[r] <- mean(is.na(d$y_post[d$pretested == 1L]))
      lost_un[r] <- mean(is.na(d$y_post[d$pretested == 0L]))
      fit <- suppressWarnings(fit_solomon_glm(y_post, treat, pretested, y_pre, data = d))
      obs <- suppressWarnings(fit_solomon_glm(y_post, treat, pretested, y_pre, data = d,
                                              pretest_center = "observed"))
      e <- fit$effects
      sens <- fit$pretest_sensitivity$effects
      at_b <- sens[sens$center == "observed", ]
      refit <- obs$effects[obs$effects$contrast %in% .solomon_pretest_order, ]
      cols <- c("estimate", "std.error", "df", "p.value", "conf.low", "conf.high")
      mismatch <- max(mismatch, max(abs(as.matrix(at_b[, cols]) - as.matrix(refit[, cols]))))
      pick <- function(tab, label, what = "estimate") tab[[what]][tab$contrast == label]
      rej <- function(tab, label) as.numeric(pick(tab, label, "p.value") < 0.05)
      keep[r, c("A", "A_rej", "B", "B_rej")] <- c(pick(e, main_label), rej(e, main_label),
                                                   pick(at_b, main_label), rej(at_b, main_label))
      if (s %in% c("L1", "L2", "R3", "R4", "P1")) {
        raw <- fit_solomon_glm(y_post, treat, pretested, data = d)$effects
        keep[r, c("raw", "raw_rej")] <- c(pick(raw, main_label), rej(raw, main_label))
      }
      if (s == "S3") {
        keep[r, c("A_control", "A_treated")] <- c(pick(e, "Pretest effect | control"),
                                                  pick(e, "Pretest effect | treated"))
        keep[r, c("B_control", "B_treated")] <- c(pick(at_b, "Pretest effect | control"),
                                                  pick(at_b, "Pretest effect | treated"))
        keep[r, c("pxt", "pxt_rej")] <- c(pick(e, "Pretest x Treatment"),
                                          rej(e, "Pretest x Treatment"))
        a <- d$y_post[d$treat == 0L & d$pretested == 1L]
        b <- d$y_post[d$treat == 0L & d$pretested == 0L]
        tt <- stats::t.test(a, b)
        keep[r, c("raw_control", "raw_control_rej")] <- c(unname(diff(rev(tt$estimate))),
                                                          as.numeric(tt$p.value < 0.05))
      }
    }
    if (mismatch > 1e-10) {
      fail(sprintf("%s: the observed rows differ from the refit by %.3g", s, mismatch))
    }
    cat(sprintf("%-4s R = %d; posttests lost: pretested %.1f%%, unpretested %.1f%%; refit check %.1e\n",
                s, R, 100 * mean(lost_pre), 100 * mean(lost_un), mismatch))
    results[[s]] <- list(keep = keep, R = R, record = record[[s]])
  }

  # Compare every value the record prints.
  checks <- list()
  for (s in order) {
    rec <- results[[s]]$record
    keep <- results[[s]]$keep
    R <- results[[s]]$R
    for (nm in names(rec)) {
      if (grepl("_rej$", nm)) next
      v <- keep[, nm]
      checks[[length(checks) + 1L]] <- list(s = s, what = nm, kind = "mean", value = mean(v),
                                            mcse = stats::sd(v) / sqrt(R), record = rec[[nm]],
                                            R = R, sd = stats::sd(v))
      rej_name <- paste0(nm, "_rej")
      if (rej_name %in% names(rec)) {
        p <- rec[[rej_name]]
        checks[[length(checks) + 1L]] <- list(s = s, what = rej_name, kind = "rejection",
                                              value = mean(keep[, rej_name]),
                                              mcse = sqrt(p * (1 - p) / R), record = p, R = R,
                                              sd = NA_real_)
      }
    }
  }
  K <- length(checks)
  z <- stats::qnorm(1 - 0.05 / (2 * K))
  cat(sprintf("\n%d comparisons with the record; z = %.3f\n", K, z))
  cat(sprintf("%-4s %-16s %9s %8s %9s %9s  %s\n", "", "value", "this run", "MCSE", "record",
              "tolerance", "pass"))
  for (ch in checks) {
    tol <- z * sqrt(2) * ch$mcse
    ok <- abs(ch$value - ch$record) <= tol
    cat(sprintf("%-4s %-16s %9.3f %8.4f %9.3f %9.3f  %s\n", ch$s, ch$what, ch$value, ch$mcse,
                ch$record, tol, if (ok) "yes" else "NO"))
    if (!ok) fail(sprintf("part 1, %s %s: %.3f against the record's %.3f", ch$s, ch$what,
                          ch$value, ch$record))
    center <- if (startsWith(ch$what, "A")) "analyzed" else if (startsWith(ch$what, "B")) "observed" else
      "none"
    add_row(part = 1L, scenario = ch$s, center = center, value = ch$what, kind = ch$kind,
            reps = ch$R, estimate = ch$value, mcse = ch$mcse, sd = ch$sd, record = ch$record,
            tolerance = tol, pass = ok)
  }
}

# ---- part 2 ---------------------------------------------------------------------

# The pretest main effect at each center on one large sample of each group,
# without the package: the pretested completers' pooled regression of y on
# x, and the unpretested completers' mean, in the treated and control arms
# alike (no treatment effect).
large_sample_centers <- function(x, y, lost, pretested) {
  pc <- pretested == 1L & !lost
  fit <- stats::lm.fit(cbind(1, x[pc]), y[pc])
  b <- unname(fit$coefficients[2])
  y_pre_comp <- mean(y[pc])
  x_comp <- mean(x[pc])
  x_all <- mean(x[pretested == 1L])
  y_un <- mean(y[pretested == 0L & !lost])
  c(A = y_pre_comp - y_un, B = y_pre_comp + b * (x_all - x_comp) - y_un)
}

part2 <- function() {
  cat("\n== Part 2: the record's derivations ==\n")
  tol2 <- 0.006
  n <- 2e6
  pretested <- rep(c(1L, 0L), each = n)
  set.seed(1352)
  record <- list(`0` = c(A = -0.001, B = -0.073), `0.3` = c(A = 0.234, B = 0.174),
                 `0.5` = c(A = 0.389, B = 0.337))
  x <- stats::rnorm(2 * n)
  e <- stats::rnorm(2 * n, 0, sqrt(0.75))
  u <- stats::runif(2 * n)
  for (delta in c(0, 0.3, 0.5)) {
    y <- delta * pretested + 0.5 * x + e
    lost <- u < lgt(-1 - 1.5 * y)
    v <- large_sample_centers(x, y, lost, pretested)
    rec <- record[[format(delta)]]
    for (center in c("A", "B")) {
      ok <- abs(v[[center]] - rec[[center]]) <= tol2
      cat(sprintf("D1 loss by the posttest, true effect %.1f, %s: %.3f (record %.3f) %s\n", delta,
                  center, v[[center]], rec[[center]], if (ok) "yes" else "NO"))
      if (!ok) fail(sprintf("part 2, D1 delta %.1f %s", delta, center))
      add_row(part = 2L, scenario = sprintf("D1, true %.1f", delta),
              center = if (center == "A") "analyzed" else "observed", value = main_label,
              kind = "mean", reps = 1L, estimate = v[[center]], mcse = NA_real_, sd = NA_real_,
              record = rec[[center]], tolerance = tol2, pass = ok)
    }
  }
  a41 <- stats::uniroot(function(a) {
    stats::integrate(function(z) lgt(a - 1.5 * z) * stats::dnorm(z), -Inf, Inf)$value - 0.41
  }, c(-3, 3), tol = 1e-10)$root
  y <- 0.5 * x + e
  cases <- list(
    D2 = list(lost = ifelse(pretested == 1L, u < lgt(a41 - 1.5 * x), u < lgt(-1 - 1.5 * x)),
              record = c(A = 0.040, B = -0.180), label = "R3 analogue"),
    D3 = list(lost = (u < lgt(-1 - 1.5 * x)) | (pretested == 1L & stats::runif(2 * n) < 0.15),
              record = c(A = 0.000, B = -0.178), label = "R4 analogue")
  )
  for (nm in names(cases)) {
    v <- large_sample_centers(x, y, cases[[nm]]$lost, pretested)
    for (center in c("A", "B")) {
      rec <- cases[[nm]]$record[[center]]
      ok <- abs(v[[center]] - rec) <= tol2
      cat(sprintf("%s %s, %s: %.3f (record %.3f) %s\n", nm, cases[[nm]]$label, center,
                  v[[center]], rec, if (ok) "yes" else "NO"))
      if (!ok) fail(sprintf("part 2, %s %s", nm, center))
      add_row(part = 2L, scenario = paste(nm, cases[[nm]]$label),
              center = if (center == "A") "analyzed" else "observed", value = main_label,
              kind = "mean", reps = 1L, estimate = v[[center]], mcse = NA_real_, sd = NA_real_,
              record = rec, tolerance = tol2, pass = ok)
    }
  }
}

# ---- part 3 ---------------------------------------------------------------------

part3 <- function() {
  cat("\n== Part 3: the mai2020 block ==\n")
  fit <- fit_solomon_glm(post_behavior, condition, pretested, pre_behavior, control = "Control",
                         data = mai2020)
  print(fit)
  s <- fit$pretest_sensitivity
  a <- fit$pretest_sensitivity$attrition
  pre <- a$pretested == 1L
  m <- check_solomon_missing(post_behavior, condition, pretested, pre_behavior,
                             control = "Control", data = mai2020, covariates = "gender")
  g <- m$covariates$gender
  share <- function(p, level) g$share[g$pretested == p & g$level == level]
  # The record's gender coefficient (-0.03, p = .76) is that of the
  # posttest analysis without the pretest, with model-based standard errors.
  fit_g <- fit_solomon_glm(post_behavior, condition, pretested, covariates = "gender",
                           robust = "none", control = "Control", data = mai2020)
  gender <- fit_g$coefficients[fit_g$coefficients$term == "genderMale", ]
  values <- list(
    list("center analyzed", fit$pretest_mean_analyzed, 3.13, 2),
    list("participants analyzed", fit$pretest_n[["analyzed"]], 74, 0),
    list("center observed", fit$pretest_mean_observed, 3.16, 2),
    list("participants observed", fit$pretest_n[["observed"]], 118, 0),
    list("without a posttest", fit$pretest_n[["no_posttest"]], 44, 0),
    list("Control, analyzed", s$difference$analyzed[1], 0.097, 3),
    list("Control, observed", s$difference$observed[1], 0.111, 3),
    list("RP, analyzed", s$difference$analyzed[2], -0.193, 3),
    list("RP, observed", s$difference$observed[2], -0.179, 3),
    list("GS, analyzed", s$difference$analyzed[3], 0.005, 3),
    list("GS, observed", s$difference$observed[3], 0.020, 3),
    list("main, analyzed", s$difference$analyzed[4], -0.030, 3),
    list("main, observed", s$difference$observed[4], -0.016, 3),
    list("difference", s$difference$difference[1], -0.014, 3),
    list("posttests missing, pretested %", 100 * sum(a$posttest_missing[pre]) / sum(a$n[pre]), 37.3, 1),
    list("posttests missing, unpretested %", 100 * sum(a$posttest_missing[!pre]) / sum(a$n[!pre]), 36.6, 1),
    list("lowest group %", 100 * min(a$share), 29, 0),
    list("highest group %", 100 * max(a$share), 46, 0),
    list("mean pretest without a posttest", s$pretest_means[["not_analyzed"]], 3.21, 2),
    list("mean pretest with one", s$pretest_means[["analyzed"]], 3.13, 2),
    list("delta, Control", s$delta$delta[1], 0.04, 2),
    list("delta, RP", s$delta$delta[2], 0.05, 2),
    list("delta, GS", s$delta$delta[3], 0.03, 2),
    list("predicted", s$predicted, 0.04, 2),
    list("men lost, pretested %", 100 * share(1, "Male"), 57, 0),
    list("women lost, pretested %", 100 * share(1, "Female"), 33, 0),
    list("men lost, unpretested %", 100 * share(0, "Male"), 64, 0),
    list("women lost, unpretested %", 100 * share(0, "Female"), 32, 0),
    list("Fisher p, pretested", unique(g$p.value[g$pretested == 1L]), 0.048, 3),
    list("Fisher p, unpretested", unique(g$p.value[g$pretested == 0L]), 0.033, 3),
    list("gender coefficient", gender$estimate, -0.03, 2),
    list("gender p", gender$p.value, 0.76, 2)
  )
  for (v in values) {
    ok <- isTRUE(round(v[[2]], v[[4]]) == round(v[[3]], v[[4]]))
    cat(sprintf("%-34s %10.4f  record %s  %s\n", v[[1]], v[[2]], format(v[[3]]),
                if (ok) "yes" else "NO"))
    if (!ok) fail(sprintf("part 3, %s: %.4f against %s", v[[1]], v[[2]], format(v[[3]])))
    add_row(part = 3L, scenario = "mai2020", center = NA_character_, value = v[[1]],
            kind = "value", reps = NA_integer_, estimate = v[[2]], mcse = NA_real_,
            sd = NA_real_, record = v[[3]], tolerance = 0.5 * 10^(-v[[4]]), pass = ok)
  }
}

# ---- part 4 ---------------------------------------------------------------------

part4 <- function() {
  cat("\n== Part 4: FIML against the observed center ==\n")
  if (!requireNamespace("lavaan", quietly = TRUE)) {
    cat("lavaan is not installed; part 4 skipped.\n")
    return(invisible(NULL))
  }
  set.seed(135135)
  n <- 2000L
  g <- rep(1:4, each = n)
  treat <- as.integer(g %in% c(1L, 3L))
  pretested <- as.integer(g <= 2L)
  x <- stats::rnorm(4L * n)
  y <- 0.3 * treat + pretested * (0.2 + 0.2 * treat) + 0.5 * x + stats::rnorm(4L * n, 0, sqrt(0.75))
  lost <- ifelse(pretested == 1L, stats::runif(4L * n) < lgt(-1 - 1.5 * x),
                 stats::runif(4L * n) < 0.329)
  d <- data.frame(y_post = ifelse(lost, NA_real_, y), treat = treat, pretested = pretested,
                  y_pre = ifelse(pretested == 1L, x, NA_real_), group = g)
  fit <- fit_solomon_glm(y_post, treat, pretested, y_pre, data = d)
  s <- fit$pretest_sensitivity$effects
  b_rows <- s[s$center == "observed", ]
  a_rows <- s[s$center == "analyzed", ]

  # Normal-theory FIML: the pretested groups with the pretest, one slope and
  # one mean pretest for the two; the unpretested groups with the posttest
  # alone.
  pre_d <- d[d$pretested == 1L, ]
  pre_d$grp <- ifelse(pre_d$treat == 1L, "treated", "control")
  m_pre <- lavaan::sem(
    "y_post ~ c(i_t, i_c) * 1 + c(b, b) * y_pre
     y_pre ~ c(mx, mx) * 1",
    data = pre_d, group = "grp", group.label = c("treated", "control"), missing = "ml",
    fixed.x = FALSE
  )
  est <- lavaan::parameterEstimates(m_pre)
  value <- function(label) est$est[est$label == label][1]
  mean_at <- function(intercept) value(intercept) + value("b") * value("mx")
  un <- d[d$pretested == 0L, ]
  mu <- tapply(un$y_post, un$treat, mean, na.rm = TRUE)
  fiml <- c(control = mean_at("i_c") - mu[["0"]], treated = mean_at("i_t") - mu[["1"]])
  observed <- c(control = b_rows$estimate[b_rows$contrast == "Pretest effect | control"],
                treated = b_rows$estimate[b_rows$contrast == "Pretest effect | treated"])
  analyzed <- c(control = a_rows$estimate[a_rows$contrast == "Pretest effect | control"],
                treated = a_rows$estimate[a_rows$contrast == "Pretest effect | treated"])
  for (k in names(fiml)) {
    ok <- abs(fiml[[k]] - observed[[k]]) <= 0.002
    cat(sprintf("Pretest effect | %-8s A %.4f  B %.4f  FIML %.4f  %s\n", k, analyzed[[k]],
                observed[[k]], fiml[[k]], if (ok) "yes" else "NO"))
    if (!ok) fail(sprintf("part 4, %s: FIML %.4f against B %.4f", k, fiml[[k]], observed[[k]]))
    add_row(part = 4L, scenario = "FIML", center = "observed",
            value = paste("Pretest effect |", k), kind = "value", reps = 1L,
            estimate = observed[[k]], mcse = NA_real_, sd = NA_real_, record = fiml[[k]],
            tolerance = 0.002, pass = ok)
  }
}

# ---- part 5 ---------------------------------------------------------------------

part5 <- function() {
  cat("\n== Part 5: repeated posttests ==\n")
  if (!requireNamespace("mmrm", quietly = TRUE)) {
    cat("mmrm is not installed; part 5 skipped.\n")
    return(invisible(NULL))
  }
  n <- 1500L
  R <- 40L
  g <- rep(1:4, each = n)
  treat <- as.integer(g %in% c(1L, 3L))
  pretested <- as.integer(g <= 2L)
  set.seed(1354)
  for (s in c("M1", "M4", "M2", "M3")) {
    keep <- matrix(NA_real_, R, 4L, dimnames = list(NULL, c("A1", "B1", "A2", "B2")))
    lost_any <- lost_two <- numeric(R)
    for (r in seq_len(R)) {
      x <- stats::rnorm(4L * n)
      u1 <- stats::rnorm(4L * n, 0, sqrt(0.75))
      u2 <- 0.5 * u1 + stats::rnorm(4L * n, 0, sqrt(0.75 * 0.75))
      y1 <- 0.5 * x + u1
      y2 <- 0.5 * x + u2
      u <- stats::runif(4L * n)
      none <- switch(
        s,
        M1 = rep(FALSE, 4L * n), M2 = rep(FALSE, 4L * n),
        M4 = ifelse(pretested == 1L, u < lgt(-1 - 1.5 * x), u < 0.329),
        M3 = u < lgt(-1 - 1.5 * x)
      )
      drop2 <- switch(
        s,
        M1 = u < lgt(-1 - 1.5 * y1), M2 = u < lgt(-1 - 1.5 * x),
        M4 = none, M3 = none
      )
      y1[none] <- NA_real_
      y2[drop2 | none] <- NA_real_
      lost_any[r] <- mean(none)
      lost_two[r] <- mean(is.na(y2))
      long <- data.frame(
        id = rep(seq_len(4L * n), 2L), occasion = rep(1:2, each = 4L * n),
        y_post = c(y1, y2), treat = rep(treat, 2L), pretested = rep(pretested, 2L),
        y_pre = rep(ifelse(pretested == 1L, x, NA_real_), 2L)
      )
      fit <- suppressWarnings(fit_solomon_mmrm(y_post, treat, pretested, id, occasion, y_pre,
                                               df = "satterthwaite", data = long))
      e <- fit$pretest_sensitivity$effects
      pick <- function(center, o) {
        e$estimate[e$center == center & e$occasion == o & e$contrast == main_label]
      }
      keep[r, ] <- c(pick("analyzed", "1"), pick("observed", "1"), pick("analyzed", "2"),
                     pick("observed", "2"))
    }
    mcse <- apply(keep, 2, stats::sd) / sqrt(R)
    cat(sprintf(paste("%s  no posttest %.1f%%, missing at occasion 2 %.1f%%; occasion 1: A %+.3f",
                      "(%.3f), B %+.3f (%.3f); occasion 2: A %+.3f (%.3f), B %+.3f (%.3f)\n"),
                s, 100 * mean(lost_any), 100 * mean(lost_two), mean(keep[, 1]), mcse[1],
                mean(keep[, 2]), mcse[2], mean(keep[, 3]), mcse[3], mean(keep[, 4]), mcse[4]))
    for (j in seq_len(4L)) {
      add_row(part = 5L, scenario = s, center = c("analyzed", "observed")[(j + 1L) %% 2L + 1L],
              value = sprintf("%s, occasion %d", main_label, if (j <= 2L) 1L else 2L),
              kind = "mean", reps = R, estimate = mean(keep[, j]), mcse = mcse[j],
              sd = stats::sd(keep[, j]), record = if (s == "M2" && j >= 3L) -0.076 else NA_real_,
              tolerance = NA_real_, pass = NA)
    }
  }
}

# ---- part 8 ---------------------------------------------------------------------

part8 <- function() {
  cat("\n== Part 8: the selection-shift rule of check_solomon_missing() ==\n")
  # The rule as check_solomon_missing() applies it, from its internal
  # function, for a four-group design.
  shift_rule <- function(y, treat, pretested, x) {
    cond <- ifelse(treat == 1L, "treatment", "control")
    .missing_pretest_scores(y, ifelse(pretested == 1L, x, NA_real_), cond,
                            c("control", "treatment"), pretested)
  }
  draw <- function(n, case) {
    g <- rep(1:4, each = n)
    treat <- as.integer(g %in% c(1L, 3L))
    pretested <- as.integer(g <= 2L)
    x <- stats::rnorm(4L * n)
    y <- 0.5 * x + stats::rnorm(4L * n, 0, sqrt(0.75))
    u <- stats::runif(4L * n)
    lost <- switch(case,
                   S1 = u < lgt(-1 - 1.5 * x),
                   MCAR = u < 0.33,
                   S3 = ifelse(treat == 0L, u < lgt(-1 - 1.5 * x), u < 0.329))
    y[lost] <- NA
    list(y = y, treat = treat, pretested = pretested, x = x)
  }
  settings <- list(list("S1", 200L, 1000L), list("S1", 50L, 1000L), list("MCAR", 200L, 1000L),
                   list("MCAR", 50L, 1000L), list("S3", 200L, 1000L), list("S3", 50L, 1000L))
  K <- sum(vapply(settings, function(s) s[[1]] != "S3", TRUE))
  z <- stats::qnorm(1 - 0.05 / (2 * K))
  set.seed(1358)
  for (s in settings) {
    case <- s[[1]]
    n <- s[[2]]
    R <- s[[3]]
    fires <- logical(R)
    unit_free <- TRUE
    for (r in seq_len(R)) {
      d <- draw(n, case)
      a <- shift_rule(d$y, d$treat, d$pretested, d$x)
      fires[r] <- a$differ
      # The same data with the pretest in other units (50 + 10 x): the same
      # statistic.
      b <- shift_rule(d$y, d$treat, d$pretested, 50 + 10 * d$x)
      unit_free <- unit_free &&
        isTRUE(all.equal(a$test[["statistic"]], b$test[["statistic"]], tolerance = 1e-8))
    }
    rate <- mean(fires)
    null <- case != "S3"
    tolerance <- if (null) z * sqrt(0.05 * 0.95 / R) else NA_real_
    ok <- if (null) abs(rate - 0.05) <= tolerance && unit_free else unit_free
    cat(sprintf("%-4s %4d per group, %d data sets: the rule fires in %.3f%s; free of the units: %s %s\n",
                case, n, R, rate, if (null) sprintf(" (nominal .05, tolerance %.3f)", tolerance) else
                  " (power: the arms lose differently)", unit_free, if (ok) "yes" else "NO"))
    if (!ok) fail(sprintf("part 8, %s at %d per group: rate %.3f", case, n, rate))
    add_row(part = 8L, scenario = sprintf("%s, %d per group", case, n), center = NA_character_,
            value = "selection shifts differ (Wald test, p < .05)", kind = "rejection",
            reps = R, estimate = rate, mcse = sqrt(rate * (1 - rate) / R), sd = NA_real_,
            record = if (null) 0.05 else NA_real_, tolerance = tolerance, pass = ok)
  }
}

# ---- part 6 ---------------------------------------------------------------------

part6 <- function() {
  cat("\n== Part 6: the article example ==\n")
  d <- simulate_solomon(n = 50, delta = 5, sens = 0, pretest_effect = 2, rho = 0.6,
                        sigma = 10, mean = 50, digits = 0, seed = 2082)
  set.seed(2083)
  p_missing <- ifelse(d$pretested == 1, stats::plogis(-1.4 - 0.08 * (d$y_pre - 50)), 0.2)
  d$y_post[stats::runif(nrow(d)) < p_missing] <- NA
  fit <- fit_solomon_glm(y_post, treat, pretested, y_pre, data = d)
  s <- fit$pretest_sensitivity$difference
  main <- s[s$contrast == main_label, ]
  for (v in list(list("analyzed", main$analyzed, 2.96, 2), list("observed", main$observed, 1.969, 3))) {
    ok <- round(v[[2]], v[[4]]) == round(v[[3]], v[[4]])
    cat(sprintf("Pretest main effect at %s: %.4f (record %s) %s\n", v[[1]], v[[2]], format(v[[3]]),
                if (ok) "yes" else "NO"))
    if (!ok) fail(sprintf("part 6, %s: %.4f against %s", v[[1]], v[[2]], format(v[[3]])))
    add_row(part = 6L, scenario = "article", center = v[[1]], value = main_label, kind = "value",
            reps = 1L, estimate = v[[2]], mcse = NA_real_, sd = NA_real_, record = v[[3]],
            tolerance = 0.5 * 10^(-v[[4]]), pass = ok)
  }
}

start <- Sys.time()
cat("solomonR", as.character(utils::packageVersion("solomonR")), "from", getwd(), "\n")
cat("parts:", paste(parts, collapse = ", "), "\n")
for (p in parts) {
  switch(as.character(p), "1" = part1(), "2" = part2(), "3" = part3(), "4" = part4(),
         "5" = part5(), "6" = part6(), "8" = part8(), stop("Unknown part: ", p))
}
if (length(rows)) {
  out <- do.call(rbind, rows)
  utils::write.csv(out, out_csv, row.names = FALSE)
  cat("\nWrote", nrow(out), "rows to", out_csv, "\n")
}
cat(sprintf("Run time: %.1f minutes\n", as.numeric(difftime(Sys.time(), start, units = "mins"))))
if (length(failures)) {
  cat("\nFAILED CHECKS:\n", paste0("  ", failures, "\n"), sep = "")
  quit(status = 1L)
}
cat("\nALL CHECKS PASSED\n")
