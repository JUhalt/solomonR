# The check behind the warning of perm_solomon(statistic = "difference")
# (issue #113). Run from the package root:
#
#   Rscript tools/perm-difference-check.R
#
# Design: 8 treated and 24 control participants in each pretest condition,
# outcomes normal with standard deviation 2 in the treated arms and 1 in the
# control arms, and no effect, so only the average effect is zero (the sharp
# null hypothesis is false). Each replication tests the average treatment
# effect with both statistics, 199 permutations each, and rejects when the
# permutation p-value is at most .05. Result on October 3, 2026 (R 4.6.1,
# 1,000 replications, six workers, about 7 minutes):
#
#   difference   0.144 (Monte Carlo standard error 0.011)
#   studentized  0.064 (Monte Carlo standard error 0.008)

n_rep <- as.integer(Sys.getenv("N_REP", "1000"))
pkg <- normalizePath(".")

cl <- parallel::makeCluster(6)
parallel::clusterExport(cl, "pkg")
invisible(parallel::clusterEvalQ(cl, suppressMessages(pkgload::load_all(pkg, quiet = TRUE))))

one <- function(i) {
  set.seed(113000 + i)
  n_t <- 8
  n_c <- 24
  pretested <- rep(c(1, 0), each = n_t + n_c)
  treat <- rep(c(rep(1, n_t), rep(0, n_c)), 2)
  y <- stats::rnorm(length(treat), 0, ifelse(treat == 1, 2, 1))
  fit <- fit_solomon_glm(y, treat, pretested, robust = "none")
  # The warning this check motivated is expected here.
  c(
    difference = suppressWarnings(
      perm_solomon(fit, reps = 199, seed = i, statistic = "difference")
    )$p_perm,
    studentized = perm_solomon(fit, reps = 199, seed = i, statistic = "studentized")$p_perm
  )
}

p <- do.call(rbind, parallel::parLapply(cl, seq_len(n_rep), one))
parallel::stopCluster(cl)

rejection <- colMeans(p <= 0.05)
print(data.frame(
  statistic = names(rejection),
  type1 = rejection,
  mcse = sqrt(rejection * (1 - rejection) / n_rep),
  row.names = NULL
))
