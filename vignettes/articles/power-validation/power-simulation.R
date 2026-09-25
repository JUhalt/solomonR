# Simulation validation of power_solomon() (issue #18)
#
# Structured by aims, data-generating mechanisms, estimands, methods, and
# performance measures (ADEMP; Morris, White & Crowther, 2019). The protocol
# and its amendment were posted to issue #18 before any results were examined.
#
# The study drives power_solomon() itself, so what is validated is the
# function users call, and compares each rejection rate with a normal-theory
# analytic benchmark where one exists.
#
# Usage:
#   Rscript power-simulation.R <repo> <output.rds> [workers] [smoke]
#
# The optional fourth argument "smoke" runs three scenarios with 20
# replications each, to check the machinery before the full study.

args <- commandArgs(trailingOnly = TRUE)
repo <- if (length(args) >= 1) args[[1]] else "."
output <- if (length(args) >= 2) args[[2]] else "power-simulation-results.rds"
workers <- if (length(args) >= 3) as.integer(args[[3]]) else max(1L, parallel::detectCores() - 1L)
smoke <- length(args) >= 4 && identical(args[[4]], "smoke")

base_seed <- 20260916L
alpha <- 0.05
sigma <- 1
sims_null <- 5000L
sims_alternative <- 2000L

# ---- Data-generating mechanisms ----------------------------------------------

allocations <- list(
  "10 per cell" = c(n1 = 10L, n2 = 10L, n3 = 10L, n4 = 10L),
  "20 per cell" = c(n1 = 20L, n2 = 20L, n3 = 20L, n4 = 20L),
  "30 per cell" = c(n1 = 30L, n2 = 30L, n3 = 30L, n4 = 30L),
  "50 per cell" = c(n1 = 50L, n2 = 50L, n3 = 50L, n4 = 50L),
  "100 per cell" = c(n1 = 100L, n2 = 100L, n3 = 100L, n4 = 100L),
  "60/60/20/20" = c(n1 = 60L, n2 = 60L, n3 = 20L, n4 = 20L),
  "20/20/60/60" = c(n1 = 20L, n2 = 20L, n3 = 60L, n4 = 60L)
)

scenarios <- expand.grid(
  allocation = names(allocations),
  delta = c(0, 0.3, 0.6),
  sens = c(0, 0.3),
  rho = c(0, 0.5, 0.8),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
scenarios$scenario <- seq_len(nrow(scenarios))
scenarios$sims <- ifelse(scenarios$delta == 0 & scenarios$sens == 0,
                         sims_null, sims_alternative)

invariance_sims <- 4000L
if (smoke) {
  scenarios <- scenarios[c(1L, 2L, nrow(scenarios)), ]
  scenarios$sims <- 20L
  invariance_sims <- 20L
}

# ---- Analytic benchmarks -----------------------------------------------------

# Two-sided rejection probability of a t test, counting both tails.
power_t <- function(effect, se, df, alpha) {
  crit <- stats::qt(1 - alpha / 2, df)
  ncp <- effect / se
  stats::pt(-crit, df, ncp) + stats::pt(crit, df, ncp, lower.tail = FALSE)
}

# The package's own normal-theory benchmark for the four Solomon contrasts,
# plus the 2x2 ANOVA interaction, which pools all four cells without adjusting
# for the pretest.
benchmark <- function(cells, delta, rho, sens) {
  contrasts <- solomonR:::.solomon_power_analytic(
    n = cells, delta = delta, rho = rho, sens = sens, sigma = sigma, alpha = alpha
  )
  se_anova <- sigma * sqrt(sum(1 / cells))
  df_anova <- sum(cells) - 4
  rbind(
    data.frame(
      estimand = contrasts$estimand,
      test = "GLM (HC3, t)",
      analytic_power = contrasts$power,
      analytic_se = contrasts$std.error,
      analytic_df = contrasts$df
    ),
    data.frame(
      estimand = "Pretest x Treatment",
      test = "2x2 ANOVA interaction",
      analytic_power = power_t(sens, se_anova, df_anova, alpha),
      analytic_se = se_anova,
      analytic_df = df_anova
    ),
    data.frame(
      estimand = "Treatment (one-sided)",
      test = "Test I (Braver & Braver, 1988)",
      analytic_power = NA_real_,
      analytic_se = NA_real_,
      analytic_df = NA_real_
    )
  )
}

# ---- One scenario ------------------------------------------------------------

run_scenario <- function(i) {
  sc <- scenarios[i, ]
  cells <- allocations[[sc$allocation]]

  simulated <- suppressWarnings(
    power_solomon(
      n = cells,
      delta = sc$delta,
      rho = sc$rho,
      sens = sc$sens,
      sigma = sigma,
      sims = sc$sims,
      stouffer = TRUE,
      alpha = alpha,
      seed = base_seed + i
    )
  )

  bench <- benchmark(cells, sc$delta, sc$rho, sc$sens)
  out <- merge(simulated, bench, by = c("estimand", "test"), sort = FALSE)
  out$scenario <- i
  out$allocation <- sc$allocation
  out$n1 <- cells[["n1"]]; out$n2 <- cells[["n2"]]
  out$n3 <- cells[["n3"]]; out$n4 <- cells[["n4"]]
  out$delta <- sc$delta
  out$sens <- sc$sens
  out$rho <- sc$rho
  out$sigma <- sigma
  out$difference <- out$power - out$analytic_power
  out
}

# ---- Scale invariance check --------------------------------------------------

scale_invariance <- function() {
  unit <- suppressWarnings(
    power_solomon(n = 30, delta = 0.5, rho = 0.5, sens = 0.3, sigma = 1,
                  sims = invariance_sims, alpha = alpha, seed = base_seed + 9001L)
  )
  doubled <- suppressWarnings(
    power_solomon(n = 30, delta = 1.0, rho = 0.5, sens = 0.6, sigma = 2,
                  sims = invariance_sims, alpha = alpha, seed = base_seed + 9001L)
  )
  data.frame(
    estimand = unit$estimand,
    test = unit$test,
    power_sigma_1 = unit$power,
    power_sigma_2 = doubled$power,
    difference = doubled$power - unit$power,
    mcse = sqrt(unit$mcse^2 + doubled$mcse^2)
  )
}

# ---- Run ---------------------------------------------------------------------

started <- Sys.time()

cl <- parallel::makeCluster(workers)
on.exit(parallel::stopCluster(cl), add = TRUE)
parallel::clusterExport(
  cl,
  c("repo", "alpha", "sigma", "base_seed", "scenarios", "allocations",
    "power_t", "benchmark", "run_scenario")
)
invisible(parallel::clusterEvalQ(cl, {
  suppressMessages(pkgload::load_all(repo, quiet = TRUE))
  NULL
}))

results <- parallel::parLapplyLB(cl, seq_len(nrow(scenarios)), run_scenario)
performance <- do.call(rbind, results)
rownames(performance) <- NULL

suppressMessages(pkgload::load_all(repo, quiet = TRUE))
invariance <- scale_invariance()

saveRDS(
  list(
    performance = performance,
    invariance = invariance,
    scenarios = scenarios,
    allocations = allocations,
    settings = list(alpha = alpha, sigma = sigma, base_seed = base_seed,
                    sims_null = sims_null, sims_alternative = sims_alternative,
                    workers = workers),
    session = utils::sessionInfo(),
    git_commit = Sys.getenv("SOLOMONR_COMMIT", NA_character_),
    started = started,
    finished = Sys.time()
  ),
  output
)

cat("Scenarios:", nrow(scenarios),
    " rows:", nrow(performance),
    " failures:", sum(performance$failures),
    " minutes:", round(as.numeric(difftime(Sys.time(), started, units = "mins")), 1), "\n")
