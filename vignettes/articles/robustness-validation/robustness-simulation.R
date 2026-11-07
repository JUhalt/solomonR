# Simulation study for issue #116: robustness of the recommended inference.
#
# Protocol (ADEMP; Morris et al., 2019), posted on issue #116 before any run:
# https://github.com/JUhalt/solomonR/issues/116#issuecomment-6072371455
# Scope addition (the pretest effects of #104):
# https://github.com/JUhalt/solomonR/issues/116#issuecomment-5977354149
# Amendment 1 (the tie rule as merged in #132):
# https://github.com/JUhalt/solomonR/issues/116#issuecomment-6077966894
# Note on the pinned version:
# https://github.com/JUhalt/solomonR/issues/116#issuecomment-6083095621
# Amendment 2 of #116 (the timing gate and the run on GitHub-hosted
# runners, R = 5,000 unless the workflow's job limits force a lower value;
# cited below by item number, and its text wins where the script differs):
# TODO(#116): link Amendment 2 once posted
#
# Run it through the simulation harness (tools/sim; issue #117), from the
# root of a checkout of the pinned commit, for example a code-testing run:
#
#   Rscript tools/sim/run-local.R --study robustness --purpose test --reps 2 --out <empty folder>
#
# or on GitHub Actions with the workflow simulation-run (study "robustness");
# see tools/sim/README.md. The harness writes the plan (one task per
# scenario and block of replications, with its L'Ecuyer-CMRG substream),
# runs the protocol's stops before any task (study$preflight), runs the
# tasks, checks and merges them, and only then calls the summary in
# robustness-summarize.R. A study run writes performance.csv, decisions.csv,
# agreement.csv, failures.csv, scenarios.csv, and the study's part of
# run-information.csv; a code-testing run or the timing gate writes
# agreement.csv (aggregated over the scenarios) and run-information.csv
# only, and reports no performance measure, no decision, and no count of
# failed fits, as Amendment 2 of #116 (item 9) says. To decide whether the
# agreement checks passed, its summary counts in memory the comparisons
# each check made, finds check 5's largest difference, and sums the
# warnings muffled in checks 1, 3, and 4, and writes none of these
# (robustness-summarize.R).
#
# Smoke mode (issue #140): with SOLOMONR_SMOKE_DIR set, the contract plans,
# runs, and summarizes the scenarios in `smoke_scenarios` in that folder, in
# this process, with SOLOMONR_SMOKE_REPS replications (default 2), the
# code-testing base seed, the preflight, and every agreement check; it stops
# if an agreement check fails. tools/check-validation-scripts.R runs it.
#
# ADEMP in brief (the protocol has the details):
# - Aims: do the HC3 GLM (M1), fit_solomon_ml() with Satterthwaite degrees
#   of freedom (M2), and participant-level perm_solomon() with the
#   studentized (M3) and the difference statistic (M4) keep their nominal
#   Type I error and coverage with non-normal errors, unequal group sizes,
#   and treatment that changes the spread; and what is the size of the TOST
#   of equivalence_solomon() at an equivalence bound?
# - Data: 306 conditions and 576 scenarios (main grid of 5 distributions x 9
#   allocations x 3 SD ratios x 2 correlations, in E0 and ET; opposite skew
#   and heterogeneous slopes in E0).
# - Estimands: the seven contrasts of fit$effects (0 in E0).
# - Methods: M1 and M2 in E0; M3 and M4, computed directly in closed form,
#   in E0 and ET; the TOST from M1's and M2's stored values.
# - Performance: rejection rates, coverage, one-sided error rates, TOST
#   size, bias, empirical and model SEs, with MCSEs; rules 1 to 3 at
#   .05/C with C = 14,576.
#
# Choices the protocol leaves to the script are marked "[Left to the
# script]" below and listed in the results comment.
#
# References (APA 7):
# DiCiccio, C. J., & Romano, J. P. (2017). Robust permutation tests for
#   correlation and regression coefficients. Journal of the American
#   Statistical Association, 112(519), 1211-1220.
#   https://doi.org/10.1080/01621459.2016.1202117
# L'Ecuyer, P., Simard, R., Chen, E. J., & Kelton, W. D. (2002). An
#   object-oriented random-number package with many long streams and
#   substreams. Operations Research, 50(6), 1073-1075.
#   https://doi.org/10.1287/opre.50.6.1073.358
# MacKinnon, J. G., & White, H. (1985). Some heteroskedasticity-consistent
#   covariance matrix estimators with improved finite sample properties.
#   Journal of Econometrics, 29(3), 305-325.
#   https://doi.org/10.1016/0304-4076(85)90158-7
# Micceri, T. (1989). The unicorn, the normal curve, and other improbable
#   creatures. Psychological Bulletin, 105(1), 156-166.
#   https://doi.org/10.1037/0033-2909.105.1.156
# Morris, T. P., White, I. R., & Crowther, M. J. (2019). Using simulation
#   studies to evaluate statistical methods. Statistics in Medicine, 38(11),
#   2074-2102. https://doi.org/10.1002/sim.8086
# Phipson, B., & Smyth, G. K. (2010). Permutation p-values should never be
#   zero: Calculating exact p-values when permutations are randomly drawn.
#   Statistical Applications in Genetics and Molecular Biology, 9(1),
#   Article 39. https://doi.org/10.2202/1544-6115.1585
# Sawilowsky, S. S. (1996, June 23). Controlling experiment-wise Type I
#   error of meta-analysis in the Solomon four-group design [Paper
#   presentation]. First International Conference on Multiple Comparisons,
#   Tel Aviv, Israel. https://digitalcommons.wayne.edu/coe_tbf/29/
# Wu, J., & Ding, P. (2021). Randomization tests for weak null hypotheses in
#   randomized experiments. Journal of the American Statistical Association,
#   116(536), 1898-1913. https://doi.org/10.1080/01621459.2020.1750415

source(file.path("tools", "sim", "contract.R"))

robustness_dir <- file.path("vignettes", "articles", "robustness-validation")

# ---- Constants of the protocol ---------------------------------------------------------

contrasts_all <- c("ATE (avg over pretest)", "Pretest x Treatment", "Treatment | pretested",
                   "Treatment | unpretested", "Pretest effect | control", "Pretest effect | treated",
                   "Pretest main effect")
# perm_solomon() tests the four treatment contrasts; their weights on the
# pretested and the unpretested stratum's treatment effect.
contrasts_treat <- contrasts_all[1:4]
perm_weights <- list(c(0.5, 0.5), c(1, -1), c(1, 0), c(0, 1))
alpha <- 0.05                  # every test, and the TOST
nominal <- 0.05                # rules 1 and 2: rate <= .05
material_rate <- 0.075         # the material tier (Bradley, 1978)
tost_deltas <- c(0.2, 0.5, 0.8)
power_effect <- 0.5
stress_dists <- c("chisq1", "t3")
n_perm <- 999L                 # M3 and M4
n_perm_check1 <- 99L
n_perm_check2 <- 199L
tie_tol <- 1e-10               # Amendment 1, for the Gaussian identity model
check1_reps <- 3L              # agreement check 1: the first 3 replications
check34_reps <- 20L            # agreement checks 3 and 4: the first 20 (E0)
tol_check1 <- 1e-8
tol_check3 <- 1e-12
tol_check4 <- c(M1 = 1e-10, M2 = 1e-6)
tol_check5 <- 1e-6
# Agreement check 4 adds the group shifts of (delta, s, pi) = (0.5, 0.5,
# 0.5); the estimates must shift by the contrasts' true values.
check4_shift <- c(0.75, 0.5, 1.0, 0.5, 0.5, 1.0, 0.75)
# Agreement check 2's 20 fixed data sets (groups of 4, 4, 3, and 3): under
# L'Ecuyer-CMRG, data set d uses set.seed(9118600 + d) and its permutations
# seed 9118700 + d (d = 1 to 20), as Amendment 2 of #116 (item 8) fixes
# them. [Left to the script] The data sets come from the normal generator
# with rho = .5 and h = 1. These 40 seeds belong to check 2 in
# every run, the gate and the study included. They are not the study's base
# seed (1160116), the gate's (9118000), or the code-testing base seed
# (9118500); code-testing runs use other seeds of 9118501 to 9118999 for
# their other checks.
check2_sets <- 20L
check2_sizes <- c(4L, 4L, 3L, 3L)
check2_data_seed <- 9118600L
check2_perm_seed <- 9118700L
# The class of the warning that #122 added (perm_solomon() with the
# difference statistic when the arms differ in number); agreement checks 1
# and 2 count it ("Warnings from #122 are muffled and counted").
warning_class_122 <- "solomonR_unbalanced_arms_warning"

# ---- Scenarios -------------------------------------------------------------------------

dists <- c("normal", "uniform", "exponential", "chisq1", "t3")
allocations <- data.frame(
  alloc = c("10-10-10-10", "20-20-20-20", "60-60-60-60", "10-30-10-30", "20-60-20-60",
            "30-10-30-10", "60-20-60-20", "30-30-10-10", "10-10-30-30"),
  alloc_type = c(rep("balanced", 3), rep("treated smaller", 2), rep("treated larger", 2),
                 "unpretested smaller", "unpretested larger"),
  n1 = c(10L, 20L, 60L, 10L, 20L, 30L, 60L, 30L, 10L),
  n2 = c(10L, 20L, 60L, 30L, 60L, 10L, 20L, 30L, 10L),
  n3 = c(10L, 20L, 60L, 10L, 20L, 30L, 60L, 10L, 30L),
  n4 = c(10L, 20L, 60L, 30L, 60L, 10L, 20L, 10L, 30L),
  stringsAsFactors = FALSE)

# The order of the scenario table, which Amendment 2 of #116 (item 8) makes
# the order of plan/scenarios.csv, is slowest factor first:
# level (main grid, opposite skew, heterogeneous slopes); in the main grid,
# the pattern (E0, then ET); then the distribution, the allocation, h, and
# rho, each in the protocol's order. Scenarios 1-270 are the main grid in
# E0, 271-540 the main grid in ET, 541-558 opposite skew, and 559-576 the
# heterogeneous slopes (distribution, then allocation). `condition` is
# shared by an E0 scenario and its ET twin. In the slopes level rho and h
# play no part (the slopes are 0.3 and 0.8 and the residual SD 0.6) and are
# NA.
build_scenarios <- function() {
  main <- expand.grid(rho = c(0.5, 0.8), h = c(1, 2, 4), a = seq_len(nrow(allocations)), dist = dists,
                      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  skew <- expand.grid(rho = c(0.5, 0.8), a = seq_len(nrow(allocations)), KEEP.OUT.ATTRS = FALSE)
  slope <- expand.grid(a = seq_len(nrow(allocations)), dist = c("normal", "exponential"),
                       KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  part <- function(level, pattern, dist, a, h, rho, delta, condition) {
    al <- allocations[a, , drop = FALSE]
    data.frame(condition = as.integer(condition), level = level, pattern = pattern, dist = dist,
               alloc = al$alloc, alloc_type = al$alloc_type, n1 = al$n1, n2 = al$n2, n3 = al$n3,
               n4 = al$n4, h = as.numeric(h), rho = as.numeric(rho), delta = delta, sens = 0, pre_eff = 0,
               stringsAsFactors = FALSE)
  }
  k <- nrow(main)
  sc <- rbind(
    part("main", "E0", main$dist, main$a, main$h, main$rho, 0, seq_len(k)),
    part("main", "ET", main$dist, main$a, main$h, main$rho, 0.5, seq_len(k)),
    part("opposite skew", "E0", "exponential", skew$a, 1, skew$rho, 0, k + seq_len(nrow(skew))),
    part("slopes", "E0", slope$dist, slope$a, NA_real_, NA_real_, 0, k + nrow(skew) + seq_len(nrow(slope))))
  sc$equal_arms <- sc$n1 == sc$n2 & sc$n3 == sc$n4
  sc$min_cell <- pmin(sc$n1, sc$n2, sc$n3, sc$n4)
  sc$sharp_null <- sc$level == "main" & sc$pattern == "E0" & !is.na(sc$h) & sc$h == 1
  rownames(sc) <- NULL
  data.frame(scenario = seq_len(nrow(sc)), sc, stringsAsFactors = FALSE)
}
scenarios <- build_scenarios()
if (nrow(scenarios) != 576L || max(scenarios$condition) != 306L || sum(scenarios$pattern == "E0") != 306L ||
    sum(scenarios$pattern == "ET") != 270L) {
  stop("The scenario table does not have the protocol's 306 conditions and 576 scenarios.", call. = FALSE)
}

# The protocol's fallback if even R = 3,000 does not fit: the uniform
# distribution is dropped from the main grid.
study_fallback <- function(sc) !(sc$level == "main" & sc$dist == "uniform")

# The rows a replication stores, one per method and contrast (protocol,
# "Stored per replication"): E0 M1, M2 (seven contrasts each), M3, and M4
# (the four treatment contrasts each); ET M3 and M4.
layout_E0 <- data.frame(method = rep(c("M1", "M2", "M3", "M4"), c(7, 7, 4, 4)),
                        contrast = c(contrasts_all, contrasts_all, contrasts_treat, contrasts_treat),
                        stringsAsFactors = FALSE)
layout_ET <- layout_E0[layout_E0$method %in% c("M3", "M4"), , drop = FALSE]
rownames(layout_ET) <- NULL

# [Left to the script] The smoke scenarios: one of each level and pattern,
# covering unequal arms, unequal pretest conditions, and h > 1.
pick_scenario <- function(level, pattern, dist, alloc, h, rho) {
  i <- which(scenarios$level == level & scenarios$pattern == pattern & scenarios$dist == dist &
               scenarios$alloc == alloc & (is.na(h) & is.na(scenarios$h) | scenarios$h %in% h) &
               (is.na(rho) & is.na(scenarios$rho) | scenarios$rho %in% rho))
  if (length(i) != 1L) stop("No unique smoke scenario.", call. = FALSE)
  scenarios$scenario[[i]]
}
smoke_scenarios <- c(
  pick_scenario("main", "E0", "exponential", "10-30-10-30", 4, 0.5),
  pick_scenario("main", "ET", "t3", "30-30-10-10", 2, 0.8),
  pick_scenario("opposite skew", "E0", "exponential", "30-10-30-10", 1, 0.8),
  pick_scenario("slopes", "E0", "normal", "10-10-30-30", NA, NA))

source(file.path(robustness_dir, "robustness-summarize.R"))

# C = 14,576 rates for the full table (and 11,964 for the fallback); this
# stops if the judged cells give other counts.
family_full <- check_family(scenarios)
invisible(check_family(scenarios[study_fallback(scenarios), , drop = FALSE], fallback = TRUE))

# ---- Data ------------------------------------------------------------------------------

# Standardized generators: the shapes of classic-simulation.R (#51), which
# follow the notes to Sawilowsky's (1996) Table 1 (the protocol's code).
std_draw <- function(m, dist) switch(dist,
  normal      = stats::rnorm(m),
  uniform     = (stats::runif(m) - 0.5) * sqrt(12),
  exponential = stats::rexp(m) - 1,
  chisq1      = (stats::rchisq(m, df = 1) - 1) / sqrt(2),
  t3          = stats::rt(m, df = 3) / sqrt(3),
  stop("Unknown distribution: ", dist, call. = FALSE))

# One data set of scenario `s`, as in the protocol's code: the latent
# baseline A first, then the errors e; participants ordered G1 (pretested
# treated), G2 (pretested control), G3 (unpretested treated), G4
# (unpretested control).
make_data <- function(s) {
  n <- c(s$n1, s$n2, s$n3, s$n4)
  g <- rep(1:4, times = n)
  treat <- c(1, 0, 1, 0)[g]
  pretested <- c(1, 1, 0, 0)[g]
  A <- std_draw(sum(n), s$dist)
  e <- std_draw(sum(n), s$dist)
  if (s$level == "opposite skew") e[treat == 1] <- -e[treat == 1]
  slope <- if (s$level == "slopes") ifelse(treat == 1, 0.3, 0.8) else s$rho
  sd_e <- if (s$level == "slopes") 0.6 else sqrt(1 - s$rho^2) * ifelse(treat == 1, s$h, 1)
  y_post <- s$delta * treat + s$sens * treat * pretested + s$pre_eff * pretested + slope * A + sd_e * e
  y_pre <- ifelse(pretested == 1, A, NA)
  data.frame(y_post = y_post, treat = treat, pretested = pretested, y_pre = y_pre)
}

# ---- Permutations ----------------------------------------------------------------------

# [Left to the script] B permutations of the treatment labels within the
# pretest strata, drawn as perm_solomon()'s .perm_treat() draws them:
# stats::ave(treat, pretested, FUN = function(x) sample(x, length(x))),
# which samples the unpretested stratum (level "0") first and then the
# pretested one, each by sample.int() of the stratum's size. Returns 0/1
# matrices of treated labels, one column per permutation, for the pretested
# and the unpretested stratum. Agreement checks 1 and 2 give the package and
# this function the same seed, so a difference in the draws would show
# there.
draw_permutations <- function(treat, pretested, B) {
  tp <- treat[pretested == 1]
  tu <- treat[pretested == 0]
  np <- length(tp)
  nu <- length(tu)
  mp <- matrix(0, np, B)
  mu <- matrix(0, nu, B)
  for (b in seq_len(B)) {
    mu[, b] <- tu[sample.int(nu)]
    mp[, b] <- tp[sample.int(np)]
  }
  list(pre = mp, unpre = mu)
}

# ---- M3 and M4, computed directly --------------------------------------------------------

# One pretest stratum's treatment effect and its HC3 variance for every
# labelling at once (protocol, "Direct computation of M3 and M4"): `m` is a
# 0/1 matrix of treated labels, one column per labelling. The model is that
# of fit_solomon_glm(), y ~ treat * pretested + pre_obs, which within the
# pretested stratum has a mean per cell and a common slope on the pretest,
# and within the unpretested stratum a mean per cell. With x the pretest:
#   tau = (ybar_T - ybar_C) - b (xbar_T - xbar_C), b = Sxy / Sxx pooled
#     within the two cells (b and the x terms are absent without a pretest);
#   residual e_i = y_i - ybar_c - b (x_i - xbar_c);
#   leverage h_i = 1/n_c + (x_i - xbar_c)^2 / Sxx;
#   weight a_i = +-1/n_c - (xbar_T - xbar_C)(x_i - xbar_c) / Sxx (+ treated);
#   HC3 variance = sum a_i^2 e_i^2 / (1 - h_i)^2 (MacKinnon & White, 1985).
# Every sum is over centered values and done by colSums(), column by
# column, so that two equal labellings give identical results.
stratum_effects <- function(y, x, m) {
  n <- length(y)
  y <- y - mean(y)
  cm <- 1 - m
  n_t <- colSums(m)
  n_c <- colSums(cm)
  cell_n <- m * rep(n_t, each = n) + cm * rep(n_c, each = n)
  y_t <- colSums(m * y) / n_t
  y_c <- colSums(cm * y) / n_c
  dy <- y - (m * rep(y_t, each = n) + cm * rep(y_c, each = n))
  if (is.null(x)) {
    tau <- y_t - y_c
    e <- dy
    h <- 1 / cell_n
    a <- (2 * m - 1) / cell_n
  } else {
    x <- x - mean(x)
    x_t <- colSums(m * x) / n_t
    x_c <- colSums(cm * x) / n_c
    dx <- x - (m * rep(x_t, each = n) + cm * rep(x_c, each = n))
    sxx <- colSums(dx^2)
    b <- colSums(dx * dy) / sxx
    gap <- x_t - x_c
    tau <- (y_t - y_c) - b * gap
    e <- dy - dx * rep(b, each = n)
    h <- 1 / cell_n + dx^2 / rep(sxx, each = n)
    a <- (2 * m - 1) / cell_n - dx * rep(gap / sxx, each = n)
  }
  list(tau = tau, var = colSums(a^2 * e^2 / (1 - h)^2))
}

# The two statistics of the four treatment contrasts for the observed
# labels and every permutation. [Left to the script] The observed statistic
# comes from the same routine as the permuted ones (column 1 holds the
# observed labels), as perm_solomon() computes both with one refit(). The
# studentized statistic is NA where the variance is not finite or not
# positive, as in perm_solomon().
perm_statistics <- function(d, perms) {
  p <- d$pretested == 1
  pre <- stratum_effects(d$y_post[p], d$y_pre[p], cbind(d$treat[p], perms$pre))
  unpre <- stratum_effects(d$y_post[!p], NULL, cbind(d$treat[!p], perms$unpre))
  out <- list(studentized = vector("list", 4), difference = vector("list", 4))
  for (k in 1:4) {
    w <- perm_weights[[k]]
    est <- if (w[2] == 0) w[1] * pre$tau else if (w[1] == 0) w[2] * unpre$tau else w[1] * pre$tau + w[2] * unpre$tau
    v <- if (w[2] == 0) w[1]^2 * pre$var else if (w[1] == 0) w[2]^2 * unpre$var else w[1]^2 * pre$var + w[2]^2 * unpre$var
    z <- rep(NA_real_, length(est))
    defined <- is.finite(v) & v > 0
    z[defined] <- est[defined] / sqrt(v[defined])
    out$difference[[k]] <- est
    out$studentized[[k]] <- z
  }
  out
}

# The tie rule of Amendment 1, written out here rather than taken from the
# package (.perm_at_least() with .perm_unit() and .perm_tie_tolerance()), so
# that a change in the package shows in agreement checks 1 and 2: s_b counts
# as at least as extreme as s when |s_b| >= |s| - 1e-10 max(|s|, u), with
# u = 1 for the studentized statistic and the largest finite |s_b| for the
# difference statistic.
tie_unit <- function(sb, statistic) if (identical(statistic, "difference")) max(abs(sb)) else 1
at_least <- function(sb, s, u) abs(sb) >= abs(s) - tie_tol * max(abs(s), u)

# The test of one contrast from its statistics (observed first): non-finite
# permuted statistics are dropped, and p = (1 + the number at least as
# extreme) / (1 + the number of finite ones) (Phipson & Smyth, 2010; as in
# perm_solomon()). A test rejects when p < .05. [Left to the script] n_tied
# counts the finite s_b with ||s_b| - |s|| <= 1e-10 max(|s|, u), the
# permutations the tolerance treats as ties. An undefined observed
# statistic, or no finite permuted one (where perm_solomon() stops), is a
# failed fit, as Amendment 2 of #116 (item 11) says.
perm_test <- function(stat, statistic) {
  s <- stat[[1L]]
  sb <- stat[-1L]
  fin <- is.finite(sb)
  if (!is.finite(s)) {
    return(list(status = "undefined", statistic = NA_real_, p_value = NA_real_,
                n_valid = sum(fin), n_tied = NA_integer_))
  }
  if (!any(fin)) {
    return(list(status = "no_valid_permutations", statistic = s, p_value = NA_real_,
                n_valid = 0L, n_tied = NA_integer_))
  }
  sb <- sb[fin]
  u <- tie_unit(sb, statistic)
  band <- tie_tol * max(abs(s), u)
  list(status = "ok", statistic = s, p_value = (sum(at_least(sb, s, u)) + 1) / (length(sb) + 1),
       n_valid = length(sb), n_tied = sum(abs(abs(sb) - abs(s)) <= band))
}

# ---- M1 and M2 ---------------------------------------------------------------------------

# The protocol's calls, with every other argument at its default.
fit_m1 <- function(d) fit_solomon_glm(y_post, treat, pretested, y_pre, data = d)
fit_m2 <- function(d) fit_solomon_ml(y_post, treat, pretested, y_pre, data = d)
# Only the package's fitting calls are wrapped: an error there is a failed
# fit; any other error stops the task as a code fault (the contract).
try_fit <- function(f, d) tryCatch(f(d), error = function(e) e)

# A fit's seven contrasts, selected by label. A failed fit is an error (M1,
# M2) or convergence != 0 (M2) (protocol). A contrast whose estimate,
# standard error, or degrees of freedom is not finite, or whose standard
# error is not positive, is an undefined statistic, as
# Amendment 2 of #116 (item 11) says: this extends the protocol's
# "undefined observed statistic" to M1 and M2.
fit_values <- function(method, fit) {
  if (inherits(fit, "error")) {
    return(list(status = rep("error", 7), estimate = rep(NA_real_, 7), std_error = rep(NA_real_, 7),
                df = rep(NA_real_, 7), convergence = NA_integer_))
  }
  i <- match(contrasts_all, fit$effects$contrast)
  if (anyNA(i)) stop("fit$effects lacks a contrast of the protocol.", call. = FALSE)
  est <- as.numeric(fit$effects$estimate[i])
  se <- as.numeric(fit$effects$std.error[i])
  df <- as.numeric(fit$effects$df[i])
  status <- ifelse(is.finite(est) & is.finite(se) & se > 0 & is.finite(df), "ok", "undefined")
  conv <- if (method == "M2") as.integer(fit$convergence) else NA_integer_
  if (method == "M2" && !identical(conv, 0L)) status[] <- "nonconverged"
  list(status = status, estimate = est, std_error = se, df = df, convergence = conv)
}

# ---- Agreement checks ----------------------------------------------------------------------

agreement_rows <- function(check, method, contrast, quantity, package, direct, tolerance,
                           passed = NULL, statistic = NA_character_, delta = NA_real_,
                           bound = NA_character_, warnings = 0L, classes = "") {
  abs_diff <- abs(as.numeric(package) - as.numeric(direct))
  if (is.null(passed)) passed <- !is.na(abs_diff) & abs_diff <= tolerance
  data.frame(check = check, method = method, contrast = contrast, statistic = statistic,
             delta = as.numeric(delta), bound = bound, quantity = quantity,
             package = as.numeric(package), direct = as.numeric(direct), abs_diff = abs_diff,
             tolerance = as.numeric(tolerance), passed = as.logical(passed),
             warnings_muffled = as.integer(warnings), warnings_122 = 0L, warning_classes = classes,
             stringsAsFactors = FALSE)
}

# Evaluates `expr` with its warnings muffled and counted (agreement checks
# 1 and 2: "Warnings from #122 are muffled and counted"; the same is done in
# checks 3 and 4, so that the checks' warnings stay apart from the methods'
# warnings that the contract counts). `n122` counts the warnings of #122's
# class among them. The summary shows these counts in the study run only,
# because outside it a count could stand in for failed fits (see
# agreement_summary() in robustness-summarize.R).
with_counted_warnings <- function(expr) {
  n <- 0L
  n122 <- 0L
  classes <- character()
  value <- withCallingHandlers(expr, warning = function(w) {
    n <<- n + 1L
    if (inherits(w, warning_class_122)) n122 <<- n122 + 1L
    classes <<- union(classes, class(w)[[1L]])
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = n, n122 = n122, classes = paste(classes, collapse = " "))
}

# The union of two space-separated lists of warning classes.
union_classes <- function(a, b) {
  v <- unique(c(strsplit(a, " ", fixed = TRUE)[[1]], strsplit(b, " ", fixed = TRUE)[[1]]))
  paste(v[nzchar(v)], collapse = " ")
}

# perm_solomon(fit, contrast, reps = B, seed, statistic, return_dist = TRUE)
# against the direct computation given the same permutations (`direct`,
# from perm_statistics() on draw_permutations() under the same seed). The
# p-values and the numbers of valid permutations must be identical, and
# the observed and the permuted statistics (finite ones, by position) agree
# within 1e-8 (absolute). [Left to the script] If perm_solomon() stops and
# the direct observed statistic is undefined (or no permuted one is
# finite), the two agree; if only one fails, they disagree.
compare_perm <- function(check, fit, direct, k, statistic, B, seed) {
  method <- if (statistic == "studentized") "M3" else "M4"
  w <- with_counted_warnings(tryCatch(
    perm_solomon(fit, contrast = contrasts_treat[[k]], reps = B, seed = seed, statistic = statistic,
                 return_dist = TRUE),
    error = function(e) e))
  pk <- w$value
  own <- perm_test(direct[[statistic]][[k]], statistic)
  pk_ok <- !inherits(pk, "error")
  own_ok <- identical(own$status, "ok")
  row <- function(quantity, package, direct_value, tolerance, passed = NULL) {
    agreement_rows(check, method, contrasts_treat[[k]], quantity, package, direct_value, tolerance,
                   passed, statistic = statistic, warnings = 0L)
  }
  if (!pk_ok || !own_ok) {
    agree <- !pk_ok && !own_ok
    out <- row("p_value", NA_real_, NA_real_, 0, passed = agree)
  } else {
    sb <- direct[[statistic]][[k]][-1L]
    sb <- sb[is.finite(sb)]
    same_n <- length(pk$z_perm) == length(sb)
    perm_diff <- if (same_n) max(abs(pk$z_perm - sb)) else NA_real_
    out <- rbind(
      row("p_value", pk$p_perm, own$p_value, 0, passed = identical(pk$p_perm, own$p_value)),
      row("valid_permutations", pk$valid_reps, own$n_valid, 0,
          passed = identical(as.integer(pk$valid_reps), as.integer(own$n_valid))),
      row("observed_statistic", pk$z_obs, own$statistic, tol_check1),
      agreement_rows(check, method, contrasts_treat[[k]], "permuted_statistics", NA_real_, NA_real_,
                     tol_check1, passed = same_n && perm_diff <= tol_check1, statistic = statistic))
    out$abs_diff[out$quantity == "permuted_statistics"] <- perm_diff
  }
  out$warnings_muffled[1L] <- w$warnings
  out$warnings_122[1L] <- w$n122
  out$warning_classes[1L] <- w$classes
  list(rows = out, package = pk, own = own)
}

# Agreement check 1: on the first 3 replications of every scenario, the
# four treatment contrasts with both statistics, 99 permutations, seed s.
# s = 1000 * scenario + replication, as Amendment 2 of #116 (item 8) fixes
# it. perm_solomon() sets it with withr::local_seed() under the task's RNG
# kind and restores the task's stream on exit; the direct computation draws
# the same permutations inside withr::with_seed(s), so neither call moves
# the task's stream.
check1 <- function(s, rep, d, fit) {
  seed <- 1000L * as.integer(s$scenario) + as.integer(rep)
  perms <- withr::with_seed(seed, draw_permutations(d$treat, d$pretested, n_perm_check1))
  direct <- perm_statistics(d, perms)
  rows <- list()
  for (statistic in c("studentized", "difference")) {
    for (k in 1:4) {
      rows[[length(rows) + 1L]] <- compare_perm("1. permutation tests", fit, direct, k, statistic,
                                                n_perm_check1, seed)$rows
    }
  }
  do.call(rbind, rows)
}

# Agreement check 3: for each contrast of M1 and M2, each Delta, and both
# bounds, equivalence_solomon() on a copy of the fit whose estimate of that
# contrast is moved by +Delta or -Delta, against tost_p() on the stored
# values: the same `equivalent`, and p_lower and p_upper within 1e-12
# (absolute). Failed fits and undefined contrasts are not compared.
check3 <- function(fits, values) {
  rows <- list()
  for (m in names(fits)) {
    v <- values[[m]]
    if (inherits(fits[[m]], "error") || all(v$status != "ok")) next
    for (k in which(v$status == "ok")) {
      row <- match(contrasts_all[[k]], fits[[m]]$effects$contrast)
      for (d in tost_deltas) {
        for (bound in c("upper", "lower")) {
          shifted <- fits[[m]]
          moved <- v$estimate[[k]] + if (bound == "upper") d else -d
          shifted$effects$estimate[row] <- moved
          w <- with_counted_warnings(equivalence_solomon(shifted, bounds = d, contrast = contrasts_all[[k]]))
          pk <- w$value
          own <- tost_p(moved, v$std_error[[k]], v$df[[k]], d)
          r <- rbind(
            agreement_rows("3. TOST", m, contrasts_all[[k]], "equivalent", pk$equivalent, own$equivalent, 0,
                           passed = identical(isTRUE(pk$equivalent), isTRUE(own$equivalent)),
                           delta = d, bound = bound),
            agreement_rows("3. TOST", m, contrasts_all[[k]], "p_lower", pk$p_lower, own$p_lower, tol_check3,
                           delta = d, bound = bound),
            agreement_rows("3. TOST", m, contrasts_all[[k]], "p_upper", pk$p_upper, own$p_upper, tol_check3,
                           delta = d, bound = bound))
          r$warnings_muffled[1L] <- w$warnings
          r$warnings_122[1L] <- w$n122
          r$warning_classes[1L] <- w$classes
          rows[[length(rows) + 1L]] <- r
        }
      }
    }
  }
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}

# Agreement check 4: M1 and M2 refitted to the posttests with the group
# shifts of (delta, s, pi) = (0.5, 0.5, 0.5) added (G1 1.5, G2 0.5, G3 0.5,
# G4 0). The estimates must shift by the true values, and the standard
# errors and degrees of freedom stay unchanged, within 1e-10 for M1 and
# 1e-6 for M2 (absolute). [Left to the script] A fit or refit that fails is
# not compared. The refits' warnings are muffled and counted here, apart
# from the methods' warnings; because fit_solomon_ml() warns when optim()
# does not converge, the summary shows them in the study run only.
check4 <- function(d, values) {
  shifted <- d
  shifted$y_post <- d$y_post + 0.5 * d$treat + 0.5 * d$treat * d$pretested + 0.5 * d$pretested
  rows <- list()
  for (m in names(values)) {
    v <- values[[m]]
    if (all(v$status != "ok")) next
    w <- with_counted_warnings(try_fit(if (m == "M1") fit_m1 else fit_m2, shifted))
    v2 <- fit_values(m, w$value)
    k <- which(v$status == "ok" & v2$status == "ok")
    if (!length(k)) next
    tol <- tol_check4[[m]]
    r <- rbind(
      agreement_rows("4. equivariance", m, contrasts_all[k], "estimate_shift",
                     v2$estimate[k] - v$estimate[k], check4_shift[k], tol),
      agreement_rows("4. equivariance", m, contrasts_all[k], "std_error", v2$std_error[k], v$std_error[k], tol),
      agreement_rows("4. equivariance", m, contrasts_all[k], "df", v2$df[k], v$df[k], tol))
    r$warnings_muffled[1L] <- w$warnings
    r$warnings_122[1L] <- w$n122
    r$warning_classes[1L] <- w$classes
    rows[[length(rows) + 1L]] <- r
  }
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}

# ---- One replication ----------------------------------------------------------------------

result_frame <- function(method, contrast, status, estimate = NA_real_, std_error = NA_real_, df = NA_real_,
                         statistic = NA_real_, p_value = NA_real_, n_valid = NA_integer_,
                         n_tied = NA_integer_, convergence = NA_integer_) {
  data.frame(method = method, contrast = contrast, status = status, estimate = as.numeric(estimate),
             std_error = as.numeric(std_error), df = as.numeric(df), statistic = as.numeric(statistic),
             p_value = as.numeric(p_value), n_valid = as.integer(n_valid), n_tied = as.integer(n_tied),
             convergence = as.integer(convergence), stringsAsFactors = FALSE)
}

# [Left to the script] Within a task, each replication draws its data (A,
# then e) and then its 999 permutations from the task's substream, and
# nothing else draws from it: the agreement checks set their own seeds in
# withr::with_seed() and perm_solomon(seed =), and the replication stops as
# a code fault if the fits, the permutation statistics, or the checks moved
# the stream. The same 999 permutations serve the four contrasts and both
# statistics (protocol).
# [Left to the script] M3 and M4 are computed directly whether or not the
# M1 fit of the replication succeeded: the direct computation does not use
# that fit, while perm_solomon(fit_M1) could not run without it.
replicate_robustness <- function(s, rep) {
  d <- make_data(s)
  perms <- draw_permutations(d$treat, d$pretested, n_perm)
  # The state after the replication's own draws; nothing below may move it.
  stream <- get(".Random.seed", envir = globalenv())
  e0 <- identical(s$pattern, "E0")
  rows <- list()
  fits <- list()
  values <- list()
  if (e0) {
    fits <- list(M1 = try_fit(fit_m1, d), M2 = try_fit(fit_m2, d))
    values <- lapply(names(fits), function(m) fit_values(m, fits[[m]]))
    names(values) <- names(fits)
    for (m in names(values)) {
      v <- values[[m]]
      rows[[m]] <- result_frame(m, contrasts_all, v$status, v$estimate, v$std_error, v$df,
                                convergence = rep(v$convergence, 7))
    }
  }
  stats <- perm_statistics(d, perms)
  for (statistic in c("studentized", "difference")) {
    m <- if (statistic == "studentized") "M3" else "M4"
    tests <- lapply(stats[[statistic]], perm_test, statistic = statistic)
    rows[[m]] <- result_frame(m, contrasts_treat, vapply(tests, `[[`, "", "status"),
                              statistic = vapply(tests, `[[`, 0, "statistic"),
                              p_value = vapply(tests, `[[`, 0, "p_value"),
                              n_valid = vapply(tests, function(x) as.integer(x$n_valid), 0L),
                              n_tied = vapply(tests, function(x) as.integer(x$n_tied), 0L))
  }
  results <- do.call(rbind, unname(rows))
  rownames(results) <- NULL
  if (!identical(stream, get(".Random.seed", envir = globalenv()))) {
    stop("A fit or the permutation statistics moved the task's random-number stream.", call. = FALSE)
  }
  agreement <- NULL
  if (rep <= check1_reps || (e0 && rep <= check34_reps)) {
    agreement <- sim_time_agreement({
      parts <- list()
      if (rep <= check1_reps) {
        # In ET, M1 is fitted for this check only (protocol). Its warnings
        # are muffled and counted with the check's, not as warnings of the
        # scenario's methods; if that fit fails, check 1 cannot run in this
        # replication, and its warnings have no row to be counted in.
        fit1 <- if (e0) list(value = fits$M1, warnings = 0L, n122 = 0L, classes = "") else
          with_counted_warnings(try_fit(fit_m1, d))
        if (!inherits(fit1$value, "error")) {
          c1 <- check1(s, rep, d, fit1$value)
          c1$warnings_muffled[1L] <- c1$warnings_muffled[1L] + fit1$warnings
          c1$warnings_122[1L] <- c1$warnings_122[1L] + fit1$n122
          c1$warning_classes[1L] <- union_classes(c1$warning_classes[1L], fit1$classes)
          parts$check1 <- c1
        }
      }
      if (e0 && rep <= check34_reps) {
        parts$check3 <- check3(fits, values)
        parts$check4 <- check4(d, values)
      }
      if (length(parts)) do.call(rbind, unname(parts)) else NULL
    })
    if (!identical(stream, get(".Random.seed", envir = globalenv()))) {
      stop("An agreement check moved the task's random-number stream.", call. = FALSE)
    }
    if (!is.null(agreement)) rownames(agreement) <- NULL
  }
  list(results = results, agreement = agreement)
}

# ---- The protocol's stops before any task: Amendment 2 of #116 (item 10) -----------------------

# A fixed data set without random numbers, 10 per group, for the first two
# stops. [Left to the script] The default inference is checked on a fit made
# without `inference`, so that the default actually used is tested, and the
# seven labels in the effects of both an M1 and an M2 fit.
fixed_data <- function() {
  g <- rep(1:4, times = c(10L, 10L, 10L, 10L))
  treat <- c(1, 0, 1, 0)[g]
  pretested <- c(1, 1, 0, 0)[g]
  i <- seq_along(g)
  a <- sin(1.7 * i)
  data.frame(y_post = 0.6 * a + cos(2.3 * i) + 0.3 * treat, treat = treat, pretested = pretested,
             y_pre = ifelse(pretested == 1, a, NA))
}

# Agreement check 2 on fixed data set `ds`: the comparison of check 1 with
# 199 permutations, and, for every permutation whose labels, throughout the
# strata the contrast uses, all equal the observed labels or all differ
# from them (a full swap, possible because the arms are equal), both the
# package and the direct computation must count it as at least as extreme.
# The package's count applies the package's own rule to its own statistics
# (its z_perm against its z_obs, with .perm_at_least(), .perm_unit(), and
# .perm_tie_tolerance() of the fit's family, as perm_solomon() calls them;
# pkgload::load_all() makes these internal functions visible); the direct
# count applies the script's copy of the rule (at_least()) to the direct
# statistics. The stored n_tied of the direct computation (perm_test()) is
# reported beside the tie-labelled count, as information: the protocol's
# pass rule does not include it.
check2_set <- function(ds) {
  s <- list(level = "main", dist = "normal", n1 = check2_sizes[[1]], n2 = check2_sizes[[2]],
            n3 = check2_sizes[[3]], n4 = check2_sizes[[4]], h = 1, rho = 0.5, delta = 0, sens = 0,
            pre_eff = 0)
  set.seed(check2_data_seed + ds)
  d <- make_data(s)
  fit <- fit_m1(d)
  seed <- check2_perm_seed + ds
  perms <- withr::with_seed(seed, draw_permutations(d$treat, d$pretested, n_perm_check2))
  direct <- perm_statistics(d, perms)
  p <- d$pretested == 1
  same_pre <- colSums(perms$pre == d$treat[p]) == sum(p)
  swap_pre <- colSums(perms$pre != d$treat[p]) == sum(p)
  same_un <- colSums(perms$unpre == d$treat[!p]) == sum(!p)
  swap_un <- colSums(perms$unpre != d$treat[!p]) == sum(!p)
  out <- list()
  for (statistic in c("studentized", "difference")) {
    for (k in 1:4) {
      w <- perm_weights[[k]]
      tied <- if (w[2] == 0) same_pre | swap_pre else if (w[1] == 0) same_un | swap_un else
        (same_pre & same_un) | (swap_pre & swap_un)
      cmp <- compare_perm("2. ties", fit, direct, k, statistic, n_perm_check2, seed)
      r <- cmp$rows
      stat <- direct[[statistic]][[k]]
      sb <- stat[-1L]
      fin <- is.finite(sb)
      counted_direct <- counted_package <- 0L
      if (identical(cmp$own$status, "ok")) {
        counted_direct <- sum(tied & fin & at_least(sb, stat[[1L]], tie_unit(sb[fin], statistic)))
      }
      pk <- cmp$package
      if (!inherits(pk, "error") && length(pk$z_perm) == sum(fin)) {
        z <- rep(NA_real_, length(sb))
        z[fin] <- pk$z_perm
        counted_package <- sum(tied & fin & .perm_at_least(z, pk$z_obs, .perm_unit(pk$z_perm, statistic),
                                                             .perm_tie_tolerance(stats::family(fit$model))),
                               na.rm = TRUE)
      }
      n_tied <- sum(tied & fin)
      stored <- if (identical(cmp$own$status, "ok")) as.integer(cmp$own$n_tied) else NA_integer_
      ok <- all(r$passed) && counted_direct == n_tied && counted_package == n_tied
      diff <- suppressWarnings(max(r$abs_diff[r$quantity %in% c("observed_statistic", "permuted_statistics")],
                                   na.rm = TRUE))
      out[[length(out) + 1L]] <- data.frame(
        check = paste0("agreement check 2 (ties): data set ", ds, ", ", r$method[[1]], ", ",
                       contrasts_treat[[k]]),
        passed = ok,
        detail = paste0("p-values identical: ", r$passed[r$quantity == "p_value"][1],
                        "; tie-labelled permutations ", n_tied, ", counted by the package ", counted_package,
                        ", by the direct computation ", counted_direct, "; stored n_tied ", stored),
        dataset = ds, method = r$method[[1]], contrast = contrasts_treat[[k]], statistic = statistic,
        max_abs_diff = if (is.finite(diff)) diff else NA_real_, tie_labelled = n_tied,
        counted_package = counted_package, counted_direct = counted_direct, n_tied_stored = stored,
        warnings_muffled = sum(r$warnings_muffled), warnings_122 = sum(r$warnings_122),
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, out)
}

preflight_robustness <- function() {
  # The kind is set here and restored at the end: the plan job's preflight
  # is a new R process with the default Mersenne-Twister, while in smoke
  # mode it follows the plan in the same process (contract), and check 2's
  # data sets must be the same in both.
  kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = globalenv())
  on.exit({
    do.call(RNGkind, as.list(kind))
    if (had_seed) assign(".Random.seed", old_seed, envir = globalenv())
    else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
  })
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  d <- fixed_data()
  f1 <- fit_m1(d)
  f2 <- fit_m2(d)
  stops <- data.frame(
    check = c("the default inference of fit_solomon_ml() is \"satterthwaite\"",
              "fit$effects has the seven contrasts (M1 and M2)"),
    passed = c(identical(f2$inference, "satterthwaite"),
               all(contrasts_all %in% f1$effects$contrast) && all(contrasts_all %in% f2$effects$contrast)),
    detail = c(paste0("fit$inference of a fit without `inference`: ", f2$inference),
               paste0("M1: ", sum(contrasts_all %in% f1$effects$contrast), " of 7; M2: ",
                      sum(contrasts_all %in% f2$effects$contrast), " of 7")),
    dataset = NA_integer_, method = NA_character_, contrast = NA_character_, statistic = NA_character_,
    max_abs_diff = NA_real_, tie_labelled = NA_integer_, counted_package = NA_integer_,
    counted_direct = NA_integer_, n_tied_stored = NA_integer_, warnings_muffled = NA_integer_,
    warnings_122 = NA_integer_, stringsAsFactors = FALSE)
  ties <- do.call(rbind, lapply(seq_len(check2_sets), check2_set))
  out <- rbind(stops, ties)
  rownames(out) <- NULL
  out
}

# ---- The timing gate: the protocol and Amendment 2 of #116 (item 4) -------------------------------

# The gate's constants: 100 replications of every scenario in tasks of 25,
# the first 5 minutes of each job left out, times averaged by allocation and
# pattern, a 30% allowance, and the fallback without the uniform
# distribution. `divisor` and `cap_hours` are those of the reading that
# decides R (see gate_project()).
gate_constants <- list(
  reps = 100L, task_size = 25L, ladder = c(5000L, 4500L, 4000L, 3500L, 3000L),
  divisor = 1024, cap_hours = 4, allowance = 1.3, warmup_minutes = 5,
  group = function(sc) paste(sc$alloc, sc$pattern),
  fallback = study_fallback,
  fallback_label = "main grid without the uniform distribution")

# The projection of Amendment 2 of #116 (item 4), computed by
# sim_gate_default() of the contract three times, once for each of the
# item's three readings:
#   1. "runner jobs" decides R: R = 5,000 unless the run would not fit in at
#      most 256 jobs of 4 workers, each expected to take at most 4 hours
#      with the 30% allowance; a perfectly balanced run then takes
#      worker-hours / 1,024 per job (divisor 1,024, cap 4 hours).
#   2. "runner worker-hours" reports the same in runner worker-hours
#      (divisor 1, cap 256 x 4 x 4 = 4,096 worker-hours: the same decision).
#   3. "posted rule" is the protocol's reading as posted (divided by the 10
#      workers, against 24 hours), reported beside it.
# The uniform fallback and the stop rule apply only if R = 3,000 does not
# fit reading 1. gate-summary.csv gives reading 1's decision and its
# seconds per unit of cost, which size the run's jobs.
gate_project <- function(times, scenarios, info) {
  reading <- function(label, divisor, cap) {
    g <- gate_constants
    g$divisor <- divisor
    g$cap_hours <- cap
    x <- sim_gate_default(list(gate = g, agreement_cost = agreement_cost_robustness), times, scenarios, info)
    x$`gate-projection` <- data.frame(reading = label, x$`gate-projection`, stringsAsFactors = FALSE)
    x
  }
  r1 <- reading("runner jobs (decides R)", 1024, 4)
  r2 <- reading("runner worker-hours", 1, 4096)
  r3 <- reading("posted rule (10 workers, 24 hours)", 10, 24)
  posted <- r3$`gate-summary`
  extra <- data.frame(
    item = c("reading", "posted_rule_decision", "posted_rule_chosen_R"),
    value = c("runner jobs: 256 jobs of 4 workers, at most 4 hours each with the 30% allowance",
              posted$value[posted$item == "decision"], posted$value[posted$item == "chosen_R"]),
    stringsAsFactors = FALSE)
  list(`gate-times` = r1$`gate-times`,
       `gate-projection` = rbind(r1$`gate-projection`, r2$`gate-projection`, r3$`gate-projection`),
       `gate-summary` = rbind(r1$`gate-summary`, extra))
}

# ---- The study ----------------------------------------------------------------------------

# [Left to the script] The relative cost of a task, which only balances the
# jobs and sizes them from the gate's seconds per unit of cost (it changes
# no result): replications x (participants + 60), plus 200 per E0
# replication for the fits of M1 and M2. In code-testing timings on the
# laptop (one R process, seed 9118503, normal errors), a replication took
# about 0.00035 s per unit of this cost in every allocation and both
# patterns (E0 0.107 to 0.174 s, ET 0.034 to 0.106 s for 40 to 240
# participants); t (3 df) was slower in E0. In a run, a job starts its
# tasks in order of expected cost, largest first, ties by task id, as
# Amendment 2 of #116 (item 2) says. The agreement checks of task 1 are
# priced by agreement_cost_robustness() below, so the costliest task is a
# task 1 of an E0 scenario with the most participants, which the replay
# runs again as its task (a), as Amendment 2 of #116 (item 3) says.
cost_robustness <- function(s, reps) {
  reps * (s$n1 + s$n2 + s$n3 + s$n4 + 60 + if (identical(s$pattern, "E0")) 200 else 0)
}

# [Left to the script] The cost, in the units of cost_robustness(), of the
# agreement checks that task 1 of a scenario carries: check 1 on
# replications 1 to 3, and checks 3 and 4 on replications 1 to 20 in E0. The
# contract adds it to the cost of task 1, so that the gate's seconds per
# unit of cost come from replication time alone and the run's jobs are
# sized with each task 1's checks (it changes no result). From the
# 1-replication test gate on the laptop (seed 9118500, run rl-gate): check 1
# took about 2.5 + 0.0024 n seconds per replication (with n participants;
# ET, which also fits M1 for it), and checks 3 and 4 together about
# 0.55 + 0.0003 n more in E0, at about 0.00028 s per unit of cost: about
# 9,000 + 9 n and 2,000 + n units per replication. This assumes that task 1
# holds replications 1 to 20, as the study's tasks of 500 and the gate's of
# 25 do. A test with a task size below 20 (the workflow's task_size) also
# runs checks 3 and 4 in later tasks of an E0 scenario, and below 3 check 1
# too; task 1 is priced for all of them, the later tasks for none. Only
# that test's job sizes change; no result does.
agreement_cost_robustness <- function(s) {
  n <- s$n1 + s$n2 + s$n3 + s$n4
  check1_reps * (9000 + 9 * n) + if (identical(s$pattern, "E0")) check34_reps * (2000 + n) else 0
}

study <- list(
  key = "robustness", issue = 116L, protocol = TRUE,
  seeds = c(study = 1160116L, gate = 9118000L, test = 9118500L),
  scenarios = scenarios,
  task_size = 500L,
  reps_allowed = c(3000L, 3500L, 4000L, 4500L, 5000L),
  gate = c(gate_constants, list(project = gate_project)),
  rows_per_rep = function(s) if (identical(s$pattern, "E0")) nrow(layout_E0) else nrow(layout_ET),
  cost = cost_robustness,
  agreement_cost = agreement_cost_robustness,
  replicate = replicate_robustness,
  preflight = preflight_robustness,
  summarize = summarize_robustness,
  smoke_scenarios = smoke_scenarios,
  # robustness-summarize.R is sourced above, and agreement check 3 uses its
  # tost_p(); its SHA-256 is part of the plan and of every task's envelope.
  files = file.path(robustness_dir, "robustness-summarize.R"))

sim_main(study)
