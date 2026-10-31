# Smoke study of the simulation harness (tools/sim; issue #117).
#
# Not a validation protocol and not evidence about the package: a small study
# that exercises every stage of the contract in tools/sim/contract.R with the
# package's own simulate_solomon() and fit_solomon_glm(), seeded as a
# protocol is seeded (one L'Ecuyer-CMRG stream per scenario, in the order of
# the scenario table; one substream per task of 5 replications). Its base
# seeds, 1170901 (study), 1170902 (gate), and 1170903 (test), appear in no
# protocol and nowhere else in the repository.
#
# Run it through the harness, from the package root:
#
#   Rscript tools/sim/run-local.R --study smoke --chunks 3
#
# or on GitHub Actions with the workflow simulation-run (study "smoke"); see
# tools/sim/README.md. With SOLOMONR_SMOKE_DIR set it runs its two smoke
# scenarios in that folder, as the validation scripts do (issue #140).
#
# Methods: M1, fit_solomon_glm() with HC3 standard errors and the pretest as
# a covariate; M2, fit_solomon_glm() with conventional standard errors on the
# posttest only. Measures (Morris et al., 2019): bias, empirical SE, coverage
# of the 95% interval, and rejection rate at .05, each with its Monte Carlo
# standard error. Agreement check: in the first 2 replications of each
# scenario, M2's seven contrasts and their standard errors are recomputed
# from the cell means and the pooled variance; every difference must be below
# 1e-10.
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

replicate_smoke <- function(s, rep) {
  d <- simulate_solomon(n = s$n, delta = s$delta, sens = s$sens, rho = s$rho)
  truth <- attr(d, "truth")
  truth <- truth$true_value[match(contrasts, truth$estimand)]
  rows <- lapply(names(methods), function(m) {
    fit <- tryCatch(fit_method(m, d), error = function(e) e)
    if (inherits(fit, "error")) {
      return(data.frame(method = m, contrast = contrasts, truth = truth,
                        estimate = NA_real_, std_error = NA_real_, conf_low = NA_real_,
                        conf_high = NA_real_, p_value = NA_real_,
                        status = paste("error:", substr(conditionMessage(fit), 1, 200)),
                        stringsAsFactors = FALSE))
    }
    e <- fit$effects[match(contrasts, fit$effects$contrast), ]
    data.frame(method = m, contrast = contrasts, truth = truth, estimate = e$estimate,
               std_error = e$std.error, conf_low = e$conf.low, conf_high = e$conf.high,
               p_value = e$p.value, status = "ok", stringsAsFactors = FALSE)
  })
  results <- do.call(rbind, rows)
  agreement <- if (rep <= agree_reps) direct_check(d, results[results$method == "M2", ]) else NULL
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
  # A protocol computes performance measures only in the study run; this
  # smoke study is not a protocol, so its tests compute them too, which lets
  # the harness check the Monte Carlo standard errors.
  if (identical(info[["purpose"]], "study") || identical(info[["protocol"]], "FALSE")) {
    keys <- unique(results[c("scenario", "method", "contrast")])
    perf <- do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
      k <- keys[i, ]
      x <- results[results$scenario == k$scenario & results$method == k$method &
                     results$contrast == k$contrast & results$status == "ok", ]
      data.frame(scenario = k$scenario, method = k$method, contrast = k$contrast,
                 measures(x), stringsAsFactors = FALSE)
    }))
    rownames(perf) <- NULL
    out$performance <- perf
    # Fits that failed, by method, and warnings caught in the replications
    # (both methods together), by scenario.
    fails <- data.frame(scenario = sort(unique(results$scenario)))
    for (m in names(methods)) {
      fails[[paste0("fit_errors_", m)]] <- vapply(fails$scenario, function(s) {
        sum(results$scenario == s & results$method == m & results$status != "ok") / length(contrasts)
      }, 0)
    }
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
  gate = list(reps = 10L, task_size = 5L, ladder = c(1000L, 500L, 100L), divisor = 4,
              cap_hours = 0.05, allowance = 1.3),
  rows_per_rep = function(s) length(methods) * length(contrasts),
  cost = function(s, reps) reps * s$n,
  replicate = replicate_smoke,
  preflight = preflight_smoke,
  summarize = summarize_smoke,
  smoke_scenarios = c(1L, 4L))

sim_main(study)
