# Agreement check and summary of the harness's smoke study (issue #117),
# sourced by tools/sim/studies/smoke/smoke-simulation.R, whose constants
# (contrasts, measured, alpha, agree_reps, agree_tol) they use. This file
# defines the contrast matrix and functions only.
#
# The script lists this file in study$files, so its SHA-256 is part of the
# plan and of every task's envelope, and the merge compares it: the
# replications' agreement check calls direct_check() from here. #116's
# script treats its robustness-summarize.R the same way, so the harness's
# own runs test that path of the contract.
#
# Measures (for T0 only; see the script): bias, empirical SE, coverage of
# the 95% interval, and rejection rate at .05, each with its Monte Carlo
# standard error (Morris et al., 2019).
#
# Reference:
# Morris, T. P., White, I. R., & Crowther, M. J. (2019). Using simulation
# studies to evaluate statistical methods. Statistics in Medicine, 38(11),
# 2074-2102. https://doi.org/10.1002/sim.8086

# ---- Agreement check -------------------------------------------------------------

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

summarize_smoke <- function(results, agreement, failures, scenarios, info, preflight = NULL) {
  out <- list()
  if (is.null(agreement) || !nrow(agreement)) stop("No agreement check was run.", call. = FALSE)
  agr <- do.call(rbind, lapply(split(agreement, agreement$scenario), function(a) {
    data.frame(check = "M2 against the cell means", scenario = a$scenario[[1]], method = "M2",
               replications = length(unique(a$rep)), values = nrow(a), max_abs_diff = max(a$abs_diff),
               tolerance = agree_tol, passed = all(a$abs_diff < agree_tol), stringsAsFactors = FALSE)
  }))
  # The preflight's agreement check, on the fixed data set, as a protocol
  # reports a check that runs before any task.
  if (!is.null(preflight)) {
    p <- preflight[preflight$check == "the agreement check passes on a fixed data set", , drop = FALSE]
    agr <- rbind(data.frame(check = "M2 against the cell means, fixed data set (preflight)", scenario = NA_integer_,
                            method = "M2", replications = 0L, values = 2L * length(contrasts),
                            max_abs_diff = as.numeric(p$detail), tolerance = agree_tol,
                            passed = identical(p$passed, "TRUE"), stringsAsFactors = FALSE), agr)
  }
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
