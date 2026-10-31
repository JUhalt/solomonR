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
# Reference:
# Morris, T. P., White, I. R., & Crowther, M. J. (2019). Using simulation
# studies to evaluate statistical methods. Statistics in Medicine, 38(11),
# 2074-2102. https://doi.org/10.1002/sim.8086

source(file.path("tools", "sim", "contract.R"))

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

# ---- One replication -------------------------------------------------------------

fit_method <- function(m, d) {
  if (m == "M1") {
    fit_solomon_glm(d$y_post, d$treat, d$pretested, y_pre = d$y_pre, robust = "HC3")
  } else {
    fit_solomon_glm(d$y_post, d$treat, d$pretested, robust = "none")
  }
}

# Contrasts of the four cell means (TP, CP, TU, CU: treated or control,
# pretested or not) in the order of `contrasts`.
cell_contrasts <- rbind(
  c(0.5, -0.5, 0.5, -0.5), c(1, -1, -1, 1), c(1, -1, 0, 0), c(0, 0, 1, -1),
  c(0, 1, 0, -1), c(1, 0, -1, 0), c(0.5, 0.5, -0.5, -0.5))

direct_check <- function(d, m2) {
  cells <- list(d$treat == 1 & d$pretested == 1, d$treat == 0 & d$pretested == 1,
                d$treat == 1 & d$pretested == 0, d$treat == 0 & d$pretested == 0)
  mu <- vapply(cells, function(i) mean(d$y_post[i]), 0)
  n <- vapply(cells, sum, 0)
  sse <- sum(vapply(cells, function(i) sum((d$y_post[i] - mean(d$y_post[i]))^2), 0))
  mse <- sse / (length(d$y_post) - 4)
  est <- drop(cell_contrasts %*% mu)
  se <- sqrt(mse * drop(cell_contrasts^2 %*% (1 / n)))
  pkg <- m2[match(contrasts, m2$contrast), ]
  out <- data.frame(method = "M2", contrast = rep(contrasts, 2),
                    quantity = rep(c("estimate", "std_error"), each = length(contrasts)),
                    package = c(pkg$estimate, pkg$std_error), direct = c(est, se),
                    stringsAsFactors = FALSE)
  out$abs_diff <- abs(out$package - out$direct)
  out
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

# ---- Summaries ---------------------------------------------------------------------

mcse_row <- function(measure, estimate, mcse, n, note = "") {
  data.frame(measure = measure, estimate = estimate, mcse = mcse, mcse_note = note,
             n_success = n, stringsAsFactors = FALSE)
}

measures <- function(x) {
  n <- nrow(x)
  if (n < 2L) {
    note <- "fewer than 2 successful fits"
    return(rbind(mcse_row("bias", NA_real_, NA_real_, n, note),
                 mcse_row("empirical_se", NA_real_, NA_real_, n, note),
                 mcse_row("coverage", NA_real_, NA_real_, n, note),
                 mcse_row("rejection", NA_real_, NA_real_, n, note)))
  }
  truth <- x$truth[[1]]
  emp <- stats::sd(x$estimate)
  prop <- function(name, hit) {
    p <- mean(hit)
    mcse_row(name, p, sqrt(p * (1 - p) / n), n,
             if (p %in% c(0, 1)) "proportion at 0 or 1: the formula gives 0" else "")
  }
  rbind(mcse_row("bias", mean(x$estimate - truth), emp / sqrt(n), n),
        mcse_row("empirical_se", emp, emp / sqrt(2 * (n - 1)), n),
        prop("coverage", x$conf_low <= truth & truth <= x$conf_high),
        prop("rejection", x$p_value < alpha))
}

summarize_smoke <- function(results, agreement, failures, scenarios, info) {
  out <- list()
  if (is.null(agreement) || !nrow(agreement)) stop("No agreement check was run.", call. = FALSE)
  agr <- do.call(rbind, lapply(split(agreement, agreement$scenario), function(a) {
    data.frame(scenario = a$scenario[[1]], method = "M2", replications = length(unique(a$rep)),
               values = nrow(a), max_abs_diff = max(a$abs_diff), tolerance = agree_tol,
               passed = all(a$abs_diff < agree_tol), stringsAsFactors = FALSE)
  }))
  rownames(agr) <- NULL
  out$agreement <- agr
  # The study's part of run-information.csv.
  out$`run-information` <- data.frame(
    item = c("methods_stored", "methods_measured", "agreement_replications"),
    value = c("M1 M2 T0", measured, agree_reps), stringsAsFactors = FALSE)
  # A protocol computes performance measures only in the study run; this
  # smoke study is not a protocol, so its tests compute them too, for T0
  # only, which lets the harness check the Monte Carlo standard errors.
  if (identical(info[["purpose"]], "study") || identical(info[["protocol"]], "FALSE")) {
    x <- results[results$method == measured, , drop = FALSE]
    perf <- do.call(rbind, lapply(sort(unique(x$scenario)), function(s) {
      y <- x[x$scenario == s & x$status == "ok", , drop = FALSE]
      data.frame(scenario = s, method = measured, contrast = contrasts[[1]], measures(y),
                 stringsAsFactors = FALSE)
    }))
    rownames(perf) <- NULL
    out$performance <- perf
    # Failed T0 tests, and warnings caught in the replications, by scenario.
    fails <- data.frame(scenario = sort(unique(results$scenario)))
    fails$errors_T0 <- vapply(fails$scenario, function(s) sum(x$scenario == s & x$status != "ok"), 0L)
    fails$warnings <- vapply(fails$scenario, function(s) {
      if (is.null(failures)) 0 else sum(failures$count[failures$scenario == s])
    }, 0)
    out$failures <- fails
  }
  if (!all(agr$passed)) stop("Agreement check failed; see agreement.csv.", call. = FALSE)
  out
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
  replicate = replicate_smoke,
  preflight = preflight_smoke,
  summarize = summarize_smoke,
  smoke_scenarios = c(1L, 4L))

sim_main(study)
