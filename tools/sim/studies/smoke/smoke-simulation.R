# Smoke study of the simulation harness (tools/sim; issue #117).
#
# Not a validation protocol and not evidence about the package: a small study
# that exercises every stage of the contract in tools/sim/contract.R, seeded
# as a protocol is seeded (one L'Ecuyer-CMRG stream per scenario, in the
# order of the scenario table; one substream per task of 5 replications). Its
# base seeds, 1170901 (study), 1170902 (gate), and 1170903 (test), appear in
# no protocol and nowhere else in the repository.
#
# Run it through the harness, from the package root:
#
#   Rscript tools/sim/run-local.R --study smoke --chunks 3 --out <empty folder>
#
# or on GitHub Actions with the workflow simulation-run (study "smoke"); see
# tools/sim/README.md. With SOLOMONR_SMOKE_DIR set it runs its two smoke
# scenarios in that folder, as the validation scripts do (issue #140).
#
# Each replication draws a data set with simulate_solomon() and stores three
# fits: M1, fit_solomon_glm() with HC3 standard errors and the pretest as a
# covariate; M2, fit_solomon_glm() with conventional standard errors on the
# posttest only (both with their seven contrasts); and T0, Welch's two-sample
# t test of the treated against the control posttest means from
# stats::t.test(), which is not part of the package or of any protocol.
#
# Measures. To keep the harness's tests from producing evidence about the
# package's inference, no measure is computed for M1 or M2: they are stored
# so that the task files have the shape of a protocol's, and M2 feeds the
# agreement check. The measures that exercise the harness's Monte Carlo
# standard error check (bias, empirical SE, coverage of the 95% interval, and
# rejection rate at .05; Morris et al., 2019) are computed for T0 only.
#
# Agreement check: in the first 2 replications of each scenario, M2's seven
# contrasts and their standard errors are recomputed from the cell means and
# the pooled variance; every difference must be below 1e-10. Its time is
# measured apart from the replications' (sim_time_agreement()).
#
# The agreement check and the summary are in smoke-summarize.R, which this
# script sources and lists in study$files, as #116's script does with its
# summary: that file's SHA-256 is part of the plan and of every task's
# envelope, so the harness's own runs exercise that check.
#
# Reference:
# Morris, T. P., White, I. R., & Crowther, M. J. (2019). Using simulation
# studies to evaluate statistical methods. Statistics in Medicine, 38(11),
# 2074-2102. https://doi.org/10.1002/sim.8086

source(file.path("tools", "sim", "contract.R"))

smoke_dir <- file.path("tools", "sim", "studies", "smoke")

contrasts <- c("ATE (avg over pretest)", "Pretest x Treatment", "Treatment | pretested",
               "Treatment | unpretested", "Pretest effect | control",
               "Pretest effect | treated", "Pretest main effect")
methods <- c(M1 = "HC3, pretest as covariate", M2 = "conventional, posttest only")
measured <- "T0"
alpha <- 0.05
agree_reps <- 2L
agree_tol <- 1e-10

scenarios <- expand.grid(n = c(10L, 20L), delta = c(0, 0.5), KEEP.OUT.ATTRS = FALSE)
scenarios <- data.frame(scenario = seq_len(nrow(scenarios)), scenarios, sens = 0, rho = 0.5)

# The agreement check (direct_check()) and the summary (summarize_smoke()).
source(file.path(smoke_dir, "smoke-summarize.R"))

# ---- One replication -------------------------------------------------------------

fit_method <- function(m, d) {
  if (m == "M1") {
    fit_solomon_glm(d$y_post, d$treat, d$pretested, y_pre = d$y_pre, robust = "HC3")
  } else {
    fit_solomon_glm(d$y_post, d$treat, d$pretested, robust = "none")
  }
}

result_rows <- function(m, contrast, truth, estimate, std_error, conf_low, conf_high, p_value, status) {
  data.frame(method = m, contrast = contrast, truth = truth, estimate = estimate,
             std_error = std_error, conf_low = conf_low, conf_high = conf_high,
             p_value = p_value, status = status, stringsAsFactors = FALSE)
}

replicate_smoke <- function(s, rep) {
  d <- simulate_solomon(n = s$n, delta = s$delta, sens = s$sens, rho = s$rho)
  truth <- attr(d, "truth")
  truth <- truth$true_value[match(contrasts, truth$estimand)]
  rows <- lapply(names(methods), function(m) {
    fit <- tryCatch(fit_method(m, d), error = function(e) e)
    if (inherits(fit, "error")) {
      return(result_rows(m, contrasts, truth, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_,
                         paste("error:", substr(conditionMessage(fit), 1, 200))))
    }
    e <- fit$effects[match(contrasts, fit$effects$contrast), ]
    result_rows(m, contrasts, truth, e$estimate, e$std.error, e$conf.low, e$conf.high, e$p.value, "ok")
  })
  welch <- tryCatch(stats::t.test(d$y_post[d$treat == 1], d$y_post[d$treat == 0]), error = function(e) e)
  rows[[length(rows) + 1L]] <- if (inherits(welch, "error")) {
    result_rows(measured, contrasts[[1]], truth[[1]], NA_real_, NA_real_, NA_real_, NA_real_, NA_real_,
                paste("error:", substr(conditionMessage(welch), 1, 200)))
  } else {
    result_rows(measured, contrasts[[1]], truth[[1]], unname(welch$estimate[[1]] - welch$estimate[[2]]),
                unname(welch$stderr), welch$conf.int[[1]], welch$conf.int[[2]], welch$p.value, "ok")
  }
  results <- do.call(rbind, rows)
  agreement <- if (rep <= agree_reps) {
    sim_time_agreement(direct_check(d, results[results$method == "M2", ]))
  } else NULL
  list(results = results, agreement = agreement)
}

# ---- The protocol's stops before any task ------------------------------------------

preflight_smoke <- function() {
  # A fixed data set (no random numbers): 5 participants per group.
  d <- data.frame(treat = rep(c(1, 0, 1, 0), each = 5), pretested = rep(c(1, 1, 0, 0), each = 5))
  d$y_post <- 50 + 10 * sin(seq_len(20)) + 3 * d$treat
  fit <- fit_solomon_glm(d$y_post, d$treat, d$pretested, robust = "none")
  check <- direct_check(d, data.frame(contrast = fit$effects$contrast, estimate = fit$effects$estimate,
                                      std_error = fit$effects$std.error))
  data.frame(
    check = c("the default standard errors of fit_solomon_glm() are HC3",
              "fit$effects has the seven contrasts",
              "the agreement check passes on a fixed data set"),
    passed = c(identical(eval(formals(fit_solomon_glm)$robust)[[1]], "HC3"),
               all(contrasts %in% fit$effects$contrast),
               isTRUE(all(check$abs_diff < agree_tol))),
    detail = c("", "", format(max(check$abs_diff), digits = 3)),
    stringsAsFactors = FALSE)
}

study <- list(
  key = "smoke", issue = 117L, protocol = FALSE,
  seeds = c(study = 1170901L, gate = 1170902L, test = 1170903L),
  scenarios = scenarios,
  task_size = 5L,
  reps_allowed = c(100L, 500L),
  # The gate's rule is the one #116 posted, with small constants: groups of
  # scenarios by n (as #116 groups by allocation and pattern), tasks that
  # start in the first 0.5 seconds of their job left out, and a fallback that
  # keeps n = 10 only.
  gate = list(reps = 10L, task_size = 5L, ladder = c(1000L, 500L, 100L), divisor = 4,
              cap_hours = 0.05, allowance = 1.3, warmup_minutes = 0.5 / 60,
              group = function(sc) paste0("n", sc$n),
              fallback = function(sc) sc$n == 10L, fallback_label = "n = 10 only"),
  rows_per_rep = function(s) length(methods) * length(contrasts) + 1L,
  cost = function(s, reps) reps * s$n,
  # A nominal price for the agreement check that task 1 carries (M2 against
  # the cell means in its first 2 replications). It changes no result; it
  # is there so that the harness's own runs exercise the contract's
  # agreement cost, which #116's study gives.
  agreement_cost = function(s) agree_reps * s$n / 4,
  replicate = replicate_smoke,
  preflight = preflight_smoke,
  summarize = summarize_smoke,
  smoke_scenarios = c(1L, 4L),
  # smoke-summarize.R is sourced above, and the replications' agreement
  # check uses its direct_check(); its SHA-256 is part of the plan and of
  # every task's envelope.
  files = file.path(smoke_dir, "smoke-summarize.R"))

sim_main(study)
