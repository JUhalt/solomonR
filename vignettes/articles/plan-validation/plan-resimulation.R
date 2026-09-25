# One-off verification for issue #24: plan each design analytically with
# plan_solomon(), then re-simulate the returned design with power_solomon() and
# compare the GLM rejection rate with the target of 0.80.
#
# Usage:
#   Rscript plan-resimulation.R <repo> <output.rds> [workers]

args <- commandArgs(trailingOnly = TRUE)
repo <- if (length(args) >= 1) args[[1]] else "."
output <- if (length(args) >= 2) args[[2]] else "plan-resimulation.rds"
workers <- if (length(args) >= 3) as.integer(args[[3]]) else max(1L, parallel::detectCores() - 1L)

target <- 0.80
sims <- 2000L
base_seed <- 20260925L

layouts <- list("equal" = c(1, 1, 1, 1), "1:1:2:2" = c(1, 1, 2, 2))
designs <- expand.grid(
  delta = c(0.3, 0.6),
  sens = c(0, 0.3),
  rho = c(0, 0.5),
  layout = names(layouts),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

suppressMessages(pkgload::load_all(repo, quiet = TRUE))

plans <- do.call(rbind, lapply(seq_len(nrow(designs)), function(i) {
  d <- designs[i, ]
  p <- plan_solomon(power = target, delta = d$delta, sens = d$sens, rho = d$rho,
                    allocation = layouts[[d$layout]])
  p <- p[!is.na(p$total_n), ]
  cbind(design = i, d[rep(1, nrow(p)), ], p)
}))
rownames(plans) <- NULL

check_one <- function(j) {
  row <- plans[j, ]
  res <- power_solomon(
    n = c(n1 = row$n1, n2 = row$n2, n3 = row$n3, n4 = row$n4),
    delta = row$delta, rho = row$rho, sens = row$sens, sims = sims,
    stouffer = FALSE, seed = base_seed + j
  )
  glm <- res[res$test == "GLM (HC3, t)" & res$estimand == row$estimand, ]
  data.frame(row = j, simulated_power = glm$power, simulated_mcse = glm$mcse,
             failures = glm$failures)
}

started <- Sys.time()
cl <- parallel::makeCluster(workers)
on.exit(parallel::stopCluster(cl), add = TRUE)
parallel::clusterExport(cl, c("repo", "plans", "sims", "base_seed", "check_one"))
invisible(parallel::clusterEvalQ(cl, {
  suppressMessages(pkgload::load_all(repo, quiet = TRUE))
  NULL
}))
checks <- do.call(rbind, parallel::parLapplyLB(cl, seq_len(nrow(plans)), check_one))

result <- cbind(plans, checks[order(checks$row), c("simulated_power", "simulated_mcse", "failures")])
result$difference <- result$simulated_power - target
result$smallest_cell <- pmin(result$n1, result$n2, result$n3, result$n4)

saveRDS(list(result = result, target = target, sims = sims, base_seed = base_seed,
             session = utils::sessionInfo(),
             git_commit = Sys.getenv("SOLOMONR_COMMIT", NA_character_),
             started = started, finished = Sys.time()), output)

cat("plans:", nrow(result), " minutes:",
    round(as.numeric(difftime(Sys.time(), started, units = "mins")), 1), "\n")
