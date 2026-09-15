# Simulation validation of fit_solomon_ml() inference options (issues #10, #22)
#
# Structured by aims, data-generating mechanisms, estimands, methods, and
# performance measures (ADEMP; Morris, White & Crowther, 2019). The protocol
# and tolerances were posted in issue #10 before the first run. Amendment for
# #22, posted before this run: 30 and 40 participants per cell are added to
# place the small-sample warning threshold, and the exploratory separate-OLS
# method is replaced by the package's implemented `inference = "satterthwaite"`
# option.
#
# Usage:
#   Rscript ml-simulation.R <repo> <output.rds> [nsim] [workers]

args <- commandArgs(trailingOnly = TRUE)
repo <- if (length(args) >= 1) args[[1]] else "."
output <- if (length(args) >= 2) args[[2]] else "ml-simulation-results.rds"
nsim <- if (length(args) >= 3) as.integer(args[[3]]) else 2000L
workers <- if (length(args) >= 4) as.integer(args[[4]]) else max(1L, parallel::detectCores() - 1L)

base_seed <- 20260914L
alpha <- 0.05
tau_unpretested <- 0.5

# ---- Data-generating mechanisms ------------------------------------------

scenarios <- expand.grid(
  n_cell = c(6L, 10L, 20L, 30L, 40L, 50L, 100L),
  sensitization = c(0, 0.3),
  rho = c(0, 0.5, 0.8),
  sd_ratio = c(1, 1.5),
  KEEP.OUT.ATTRS = FALSE
)
scenarios$scenario <- seq_len(nrow(scenarios))

# Every participant has a latent baseline X ~ N(0, 1); pretested participants
# observe it. The posttest is
#   Y = tau_U * T + delta * T * P + rho * X + sqrt(1 - rho^2) * s_T * e,
# with s_T = sd_ratio for treated and 1 for control participants, so the
# population treatment effects are tau_U (unpretested), tau_U + delta
# (pretested), delta (sensitization), and tau_U + delta / 2 (ATE).
simulate_solomon <- function(n_cell, sensitization, rho, sd_ratio) {
  treat <- rep(c(1L, 0L, 1L, 0L), each = n_cell)
  pretested <- rep(c(1L, 1L, 0L, 0L), each = n_cell)
  x <- stats::rnorm(4L * n_cell)
  e <- stats::rnorm(4L * n_cell)
  s <- ifelse(treat == 1L, sd_ratio, 1)
  y <- tau_unpretested * treat + sensitization * treat * pretested +
    rho * x + sqrt(1 - rho^2) * s * e
  data.frame(
    y_post = y,
    treat = treat,
    pretested = pretested,
    y_pre = ifelse(pretested == 1L, x, NA_real_)
  )
}

true_values <- function(sensitization) {
  c(
    "ATE (avg over pretest)" = tau_unpretested + sensitization / 2,
    "Pretest x Treatment" = sensitization,
    "Treatment | pretested" = tau_unpretested + sensitization,
    "Treatment | unpretested" = tau_unpretested
  )
}

contrasts <- names(true_values(0))

# ---- Methods -----------------------------------------------------------------

ml_effects <- function(d, inference) {
  fit <- fit_solomon_ml(d$y_post, d$treat, d$pretested, d$y_pre, inference = inference)
  if (fit$convergence != 0L) stop("optim convergence code ", fit$convergence)
  e <- fit$effects[match(contrasts, fit$effects$contrast), ]
  data.frame(estimate = e$estimate, se = e$std.error, df = e$df, p = e$p.value)
}

method_ml_wald <- function(d) ml_effects(d, "wald")

method_ml_satterthwaite <- function(d) ml_effects(d, "satterthwaite")

method_glm_hc3 <- function(d) {
  fit <- fit_solomon_glm(d$y_post, d$treat, d$pretested, d$y_pre, robust = "HC3")
  e <- fit$effects[match(contrasts, fit$effects$contrast), ]
  data.frame(estimate = e$estimate, se = e$std.error, df = e$df, p = e$p.value)
}

methods <- list(
  "ML, Wald (default)" = method_ml_wald,
  "ML, Satterthwaite (small-sample option)" = method_ml_satterthwaite,
  "GLM HC3 (t)" = method_glm_hc3
)

# ---- One scenario ------------------------------------------------------------

run_scenario <- function(i) {
  sc <- scenarios[i, ]
  set.seed(base_seed + i, kind = "L'Ecuyer-CMRG")
  truth <- true_values(sc$sensitization)

  rows <- vector("list", nsim * length(methods))
  failures <- list()
  k <- 0L

  for (r in seq_len(nsim)) {
    d <- simulate_solomon(sc$n_cell, sc$sensitization, sc$rho, sc$sd_ratio)
    for (m in names(methods)) {
      k <- k + 1L
      res <- tryCatch(
        withCallingHandlers(
          methods[[m]](d),
          warning = function(w) invokeRestart("muffleWarning")
        ),
        error = function(e) e
      )
      if (inherits(res, "error")) {
        failures[[length(failures) + 1L]] <- data.frame(
          scenario = i, replicate = r, method = m, message = conditionMessage(res)
        )
        next
      }
      rows[[k]] <- data.frame(
        scenario = i, replicate = r, method = m, contrast = contrasts,
        truth = unname(truth), res
      )
    }
  }

  reps <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  list(
    replicates = reps,
    failures = if (length(failures)) do.call(rbind, failures) else NULL
  )
}

# ---- Performance measures (Morris, White & Crowther, 2019) -------------------

summarise_performance <- function(reps, nsim_planned) {
  groups <- split(reps, list(reps$scenario, reps$method, reps$contrast), drop = TRUE)
  out <- lapply(groups, function(g) {
    n_ok <- nrow(g)
    crit <- stats::qt(1 - alpha / 2, g$df)
    covered <- abs(g$estimate - g$truth) <= crit * g$se
    emp_se <- stats::sd(g$estimate)
    mod_se <- sqrt(mean(g$se^2))
    coverage <- mean(covered)
    rejection <- mean(g$p < alpha)
    data.frame(
      scenario = g$scenario[1],
      method = g$method[1],
      contrast = g$contrast[1],
      truth = g$truth[1],
      n_successful = n_ok,
      n_failed = nsim_planned - n_ok,
      bias = mean(g$estimate) - g$truth[1],
      mcse_bias = emp_se / sqrt(n_ok),
      empirical_se = emp_se,
      mcse_empirical_se = emp_se / sqrt(2 * (n_ok - 1)),
      model_se = mod_se,
      relative_se_error = 100 * (mod_se / emp_se - 1),
      mcse_relative_se_error = 100 * (mod_se / emp_se) *
        sqrt(stats::var(g$se^2) / (4 * n_ok * mod_se^4) + 1 / (2 * (n_ok - 1))),
      coverage = coverage,
      mcse_coverage = sqrt(coverage * (1 - coverage) / n_ok),
      rejection = rejection,
      mcse_rejection = sqrt(rejection * (1 - rejection) / n_ok)
    )
  })
  res <- do.call(rbind, out)
  rownames(res) <- NULL
  merge(scenarios, res, by = "scenario")
}

# ---- Run -----------------------------------------------------------------------

started <- Sys.time()

cl <- parallel::makeCluster(workers)
on.exit(parallel::stopCluster(cl), add = TRUE)
parallel::clusterExport(
  cl,
  c("repo", "nsim", "base_seed", "alpha", "tau_unpretested", "scenarios",
    "simulate_solomon", "true_values", "contrasts", "methods", "ml_effects",
    "method_ml_wald", "method_ml_satterthwaite", "method_glm_hc3", "run_scenario")
)
invisible(parallel::clusterEvalQ(cl, {
  suppressMessages(pkgload::load_all(repo, quiet = TRUE))
  NULL
}))

results <- parallel::parLapplyLB(cl, seq_len(nrow(scenarios)), run_scenario)

replicates <- do.call(rbind, lapply(results, `[[`, "replicates"))
failures <- do.call(rbind, lapply(results, `[[`, "failures"))
performance <- summarise_performance(replicates, nsim)

saveRDS(
  list(
    performance = performance,
    failures = failures,
    scenarios = scenarios,
    settings = list(nsim = nsim, base_seed = base_seed, alpha = alpha,
                    tau_unpretested = tau_unpretested, workers = workers),
    session = utils::sessionInfo(),
    git_commit = Sys.getenv("SOLOMONR_COMMIT", NA_character_),
    started = started,
    finished = Sys.time()
  ),
  output
)

cat("Scenarios:", nrow(scenarios), " replications:", nsim,
    " failures:", if (is.null(failures)) 0L else nrow(failures),
    " minutes:", round(as.numeric(difftime(Sys.time(), started, units = "mins")), 1), "\n")
