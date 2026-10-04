# Monte Carlo check of the pretest-effect standard errors after adding the
# variance of the estimated mean pretest (issue #104 review).
#
# Rscript sim_pretest_se.R <out.csv> [reps_scale]

args <- commandArgs(trailingOnly = TRUE)
out_file <- args[1]
scale_reps <- if (length(args) > 1) as.numeric(args[2]) else 1
wt <- "C:/Users/TheGr/OneDrive/Documents/Articles/R Packages/solomonR-wt26"

library(parallel)
n_workers <- 5L
cl <- makeCluster(n_workers)
clusterExport(cl, "wt")
invisible(clusterEvalQ(cl, {
  suppressMessages(pkgload::load_all(wt, quiet = TRUE, export_all = TRUE))
  NULL
}))

pretest <- c("Pretest effect | control", "Pretest effect | treated", "Pretest main effect")

# ---- one replication per scenario type ----------------------------------------

rep_fourgroup <- function(n, rho, seed) {
  d <- simulate_solomon(n = n, delta = 0.5, sens = 0.3, pretest_effect = 0.3, rho = rho,
                        seed = seed)
  truth <- attr(d, "truth")
  tv <- setNames(truth$true_value, truth$estimand)
  rows <- list()
  add <- function(method, e) {
    e <- e[e$contrast %in% pretest, ]
    rows[[length(rows) + 1]] <<- data.frame(method = method, contrast = e$contrast,
                                            estimate = e$estimate, se = e$std.error,
                                            lo = e$conf.low, hi = e$conf.high,
                                            truth = tv[e$contrast])
  }
  for (rob in c("HC3", "none")) {
    fit <- fit_solomon_glm(y_post, treat, pretested, y_pre, robust = rob, data = d)
    add(paste0("glm ", rob), fit$effects)
    # The same contrasts with the mean pretest treated as fixed (before the fix).
    cf <- coef(fit$model)
    Ls <- list(c(pretested = 1), c(pretested = 1, "treat:pretested" = 1),
               c(pretested = 1, "treat:pretested" = 0.5))
    old <- do.call(rbind, lapply(seq_along(Ls), function(j) {
      L <- setNames(numeric(length(cf)), names(cf))
      L[names(Ls[[j]])] <- Ls[[j]]
      est <- sum(L * cf)
      se <- sqrt(drop(t(L) %*% fit$vcov %*% L))
      q <- qt(0.975, df.residual(fit$model))
      data.frame(contrast = pretest[j], estimate = est, std.error = se,
                 conf.low = est - q * se, conf.high = est + q * se)
    }))
    add(paste0("glm ", rob, " (mean fixed)"), old)
  }
  for (inf in c("wald", "satterthwaite")) {
    ml <- fit_solomon_ml(y_post, treat, pretested, y_pre, inference = inf, data = d)
    add(paste0("ml ", inf), ml$effects)
  }
  do.call(rbind, rows)
}

# Heterogeneous slopes: the pretest-posttest slope is 0.8 among controls and
# 0.3 among treated participants, so the common-slope model is misspecified.
rep_slopes <- function(n, seed) {
  set.seed(seed)
  treat <- rep(c(1, 0, 1, 0), each = n)
  pretested <- rep(c(1, 1, 0, 0), each = n)
  A <- rnorm(4 * n)
  slope <- ifelse(treat == 1, 0.3, 0.8)
  y_post <- 0.5 * treat + 0.3 * pretested + slope * A + rnorm(4 * n, sd = 0.6)
  y_pre <- ifelse(pretested == 1, A, NA)
  fit <- fit_solomon_glm(y_post, treat, pretested, y_pre)
  e <- fit$effects[fit$effects$contrast %in% pretest, ]
  data.frame(method = "glm HC3", contrast = e$contrast, estimate = e$estimate, se = e$std.error,
             lo = e$conf.low, hi = e$conf.high, truth = 0.3)
}

rep_ngroup <- function(n, rho, seed) {
  d <- simulate_solomon(n = n, delta = c(0.5, 0.3), sens = c(0.3, 0), pretest_effect = 0.3,
                        rho = rho, seed = seed)
  truth <- attr(d, "truth")
  fit <- fit_solomon_glm(y_post, treat, pretested, y_pre, control = "Control", data = d)
  e <- fit$effects[fit$effects$contrast %in% pretest, ]
  key <- paste(e$comparison, e$contrast)
  tk <- paste(truth$comparison, truth$contrast)
  data.frame(method = "glm HC3 (two treatments)", contrast = key, estimate = e$estimate,
             se = e$std.error, lo = e$conf.low, hi = e$conf.high,
             truth = truth$true_value[match(key, tk)])
}

# Cluster-randomized cells: 10 clusters of 8 participants per cell, a
# cluster effect on the posttest (ICC about .12) and on the pretest.
rep_cr2 <- function(seed) {
  set.seed(seed)
  k <- 10; m <- 8
  cells <- expand.grid(treat = c(1, 0), pretested = c(1, 0))
  d <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i) {
    data.frame(treat = cells$treat[i], pretested = cells$pretested[i],
               site = paste(i, rep(seq_len(k), each = m)))
  }))
  sites <- unique(d$site)
  u <- setNames(rnorm(length(sites), sd = 0.35), sites)
  w <- setNames(rnorm(length(sites), sd = 0.3), sites)
  A <- u[d$site] + rnorm(nrow(d), sd = 0.94)
  d$y_post <- 0.5 * d$treat + 0.3 * d$pretested + 0.8 * A + w[d$site] + rnorm(nrow(d), sd = 0.5)
  d$y_pre <- ifelse(d$pretested == 1, A, NA)
  fit <- suppressWarnings(fit_solomon_glm(y_post, treat, pretested, y_pre, robust = "CR2",
                                          cluster = site, data = d))
  e <- fit$effects[fit$effects$contrast %in% pretest, ]
  data.frame(method = "glm CR2", contrast = e$contrast, estimate = e$estimate, se = e$std.error,
             lo = e$conf.low, hi = e$conf.high, truth = 0.3)
}

rep_mmrm <- function(n, seed) {
  set.seed(seed)
  treat <- rep(c(1, 0, 1, 0), each = n)
  pretested <- rep(c(1, 1, 0, 0), each = n)
  N <- 4 * n
  A <- rnorm(N)
  long <- do.call(rbind, lapply(1:3, function(t) {
    data.frame(id = seq_len(N), occasion = t, treat = treat, pretested = pretested,
               y_pre = ifelse(pretested == 1, A, NA),
               y_post = 0.5 * treat + 0.3 * pretested + 0.8 * A + rnorm(N, sd = 0.6))
  }))
  fit <- suppressWarnings(fit_solomon_mmrm(y_post, treat, pretested, id, occasion, y_pre,
                                           data = long))
  e <- fit$effects[fit$effects$contrast %in% pretest, ]
  data.frame(method = "mmrm KR", contrast = paste(e$occasion, e$contrast), estimate = e$estimate,
             se = e$std.error, lo = e$conf.low, hi = e$conf.high, truth = 0.3)
}

# Marginal contrasts, with no pretest effect: logit slope 2 on the pretest
# (binary), or log-rate slope 0.5 (counts).
rep_marginal <- function(n, outcome, seed) {
  set.seed(seed)
  treat <- rep(c(1, 0, 1, 0), each = n)
  pretested <- rep(c(1, 1, 0, 0), each = n)
  A <- rnorm(4 * n)
  y <- if (outcome == "binary") {
    rbinom(4 * n, 1, plogis(-0.3 + 0.5 * treat + 2 * A))
  } else {
    rpois(4 * n, exp(0.5 + 0.3 * treat + 0.5 * A))
  }
  y_pre <- ifelse(pretested == 1, A, NA)
  fam <- if (outcome == "binary") binomial() else poisson()
  fit <- suppressWarnings(fit_solomon_glm(y, treat, pretested, y_pre, family = fam))
  m <- marginal_solomon(fit, method = "delta")$effects
  m <- m[m$contrast %in% pretest, ]
  ratio <- m$scale != m$scale[1]
  data.frame(method = paste("marginal delta", outcome), contrast = paste(m$scale, m$contrast),
             estimate = m$estimate, se = m$std.error, lo = m$conf.low, hi = m$conf.high,
             truth = ifelse(grepl("ratio", m$scale), 1, 0))
}

clusterExport(cl, c("pretest", "rep_fourgroup", "rep_slopes", "rep_ngroup", "rep_cr2",
                    "rep_mmrm", "rep_marginal"))

run <- function(label, reps, f) {
  reps <- ceiling(reps * scale_reps)
  t0 <- Sys.time()
  res <- parLapply(cl, seq_len(reps), function(r, g) {
    tryCatch(g(r), error = function(e) NULL)
  }, g = f)
  failed <- sum(vapply(res, is.null, logical(1)))
  res <- do.call(rbind, res)
  res$scenario <- label
  message(sprintf("%s: %d reps, %d failed, %.1f min", label, reps, failed,
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  attr(res, "failed") <- failed
  res
}

results <- list()
for (n in c(30, 100)) for (rho in c(0.5, 0.8)) {
  label <- sprintf("four-group n=%d rho=%.1f", n, rho)
  results[[label]] <- run(label, 4000, eval(bquote(function(r) rep_fourgroup(.(n), .(rho), 1e5 + r))))
}
results[["slopes"]] <- run("heterogeneous slopes n=100", 3000, function(r) rep_slopes(100, 2e5 + r))
results[["ngroup"]] <- run("two treatments n=30 rho=0.8", 3000,
                           function(r) rep_ngroup(30, 0.8, 3e5 + r))
results[["cr2"]] <- run("CR2 10 clusters of 8 per cell", 2000, function(r) rep_cr2(4e5 + r))
results[["binary"]] <- run("marginal binary n=150", 2000,
                           function(r) rep_marginal(150, "binary", 5e5 + r))
results[["count"]] <- run("marginal count n=100", 2000,
                          function(r) rep_marginal(100, "count", 6e5 + r))
results[["mmrm"]] <- run("mmrm n=60, 3 occasions", 1000, function(r) rep_mmrm(60, 7e5 + r))

stopCluster(cl)

all <- do.call(rbind, results)

summ <- do.call(rbind, lapply(split(all, list(all$scenario, all$method, all$contrast), drop = TRUE),
                              function(x) {
  R <- nrow(x)
  cover <- mean(x$lo <= x$truth & x$truth <= x$hi)
  is_ratio <- grepl("ratio", x$contrast)
  est <- if (is_ratio[1]) log(x$estimate) else x$estimate
  tr <- if (is_ratio[1]) 0 else x$truth[1]
  data.frame(scenario = x$scenario[1], method = x$method[1], contrast = x$contrast[1], reps = R,
             bias = mean(est) - tr, emp_sd = sd(est), mean_se = mean(x$se),
             se_ratio = mean(x$se) / sd(est), coverage = cover,
             mcse = sqrt(0.95 * 0.05 / R))
}))
summ <- summ[order(summ$scenario, summ$method, summ$contrast), ]
rownames(summ) <- NULL

# Coverage within the Bonferroni-adjusted Monte Carlo tolerance.
corrected <- !grepl("mean fixed", summ$method)
z <- stats::qnorm(1 - 0.05 / (2 * sum(corrected)))
summ$tolerance <- z * summ$mcse
summ$within_tolerance <- abs(summ$coverage - 0.95) <= summ$tolerance
write.csv(summ, out_file, row.names = FALSE)
print(summ, digits = 3)
cat(sprintf("\nCorrected methods: %d of %d cells within tolerance; coverage %.3f to %.3f.\n",
            sum(summ$within_tolerance[corrected]), sum(corrected),
            min(summ$coverage[corrected]), max(summ$coverage[corrected])))
