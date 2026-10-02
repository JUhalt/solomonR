# Simulation study for issue #45: Solomon designs with several treatments.
#
# Protocol (ADEMP): posted on issue #45 before any run,
# https://github.com/JUhalt/solomonR/issues/45#issuecomment-5950791880
# Amendment 1 (Scheffe post hoc tests in M5), posted before the study was run:
# https://github.com/JUhalt/solomonR/issues/45#issuecomment-5956070600
#
# Run from the package root of a detached snapshot of the repository:
#
#   Rscript vignettes/articles/ngroup-validation/ngroup-simulation.R
#
# Writes performance.csv, agreement.csv, and run-information.csv next to this
# script. Each task (one scenario, one substream of 500 replications) is
# cached in the checkpoint directory, so an interrupted run resumes.

library(parallel)

out_dir <- file.path("vignettes", "articles", "ngroup-validation")
cache_dir <- file.path(out_dir, "checkpoints")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

reps <- 5000L
chunk <- 500L
alpha <- 0.05
agree_glm <- 100L   # replications per scenario checked against fit_solomon_glm()
agree_steyn <- 20L  # replications per scenario checked against fit_solomon_steyn()

# ---- Scenarios -----------------------------------------------------------------

effects <- list(
  S0 = list(pi = 0, delta = function(k) rep(0, k), sens = function(k) rep(0, k)),
  S1 = list(pi = 5, delta = function(k) rep(0, k), sens = function(k) rep(0, k)),
  S2 = list(pi = 0, delta = function(k) rep(5, k), sens = function(k) rep(0, k)),
  S3 = list(pi = 0, delta = function(k) c(5, rep(0, k - 1)), sens = function(k) c(5, rep(0, k - 1)))
)

scenarios <- do.call(rbind, lapply(2:3, function(k) {
  sizes <- c("10", "25", "50", if (k == 2) "mai")
  expand.grid(k = k, n = sizes, rho = c(0.3, 0.7), spread = c("equal", "unequal"),
              effects = names(effects), stringsAsFactors = FALSE)
}))
scenarios$scenario <- seq_len(nrow(scenarios))

# Group sizes in the order pretested treatments, pretested control,
# unpretested treatments, unpretested control (.solomon_cells()).
cell_sizes <- function(k, n) {
  if (n == "mai") return(c(24, 23, 27, 22, 15, 22))  # Mai et al. (2020) posttests
  rep(as.integer(n), 2 * (k + 1))
}

# ---- Data ------------------------------------------------------------------------

# Condition 0 is the control; 1..k the treatments.
design_of <- function(k, n) {
  sizes <- cell_sizes(k, n)
  cond <- rep(c(1:k, 0L, 1:k, 0L), times = sizes)
  pretested <- rep(rep(c(1L, 0L), each = k + 1), times = sizes)
  list(cond = cond, pretested = pretested, k = k)
}

simulate_data <- function(des, s) {
  k <- des$k
  e <- effects[[s$effects]]
  delta <- c(0, e$delta(k))
  sens <- c(0, e$sens(k))
  h <- c(1, if (s$spread == "unequal") 1.5 else 1, rep(1, k - 1))
  N <- length(des$cond)
  x <- stats::rnorm(N)
  u <- stats::rnorm(N)
  idx <- des$cond + 1L
  y <- 50 + e$pi * des$pretested + delta[idx] + sens[idx] * des$pretested +
    10 * (s$rho * x + sqrt(1 - s$rho^2) * h[idx] * u)
  list(y = y, pre = 50 + 10 * x)
}

# True values of the four contrasts for condition a minus condition b.
types <- c("ATE (avg over pretest)", "Pretest x Treatment",
           "Treatment | pretested", "Treatment | unpretested")
truth_pair <- function(s, k, a, b) {
  e <- effects[[s$effects]]
  delta <- c(0, e$delta(k))
  sens <- c(0, e$sens(k))
  d <- delta[a + 1] - delta[b + 1]
  sd <- sens[a + 1] - sens[b + 1]
  c(d + sd / 2, sd, d + sd, d)
}

# Pairs of conditions in the package's order: each treatment against the
# control, then each pair of treatments.
pairs_of <- function(k) {
  vc <- lapply(1:k, function(j) c(j, 0L))
  tt <- if (k > 1) utils::combn(1:k, 2, simplify = FALSE) else list()
  c(vc, tt)
}
pair_names <- function(k) {
  lab <- function(c) if (c == 0) "Control" else paste0("T", c)
  vapply(pairs_of(k), function(p) paste(lab(p[1]), "vs", lab(p[2])), "")
}

# ---- The joint model, computed directly (M1, M2, M3) ---------------------------

joint_fit <- function(des, y, pre, adjust_pre) {
  k <- des$k
  T <- sapply(1:k, function(j) as.numeric(des$cond == j))
  P <- des$pretested
  X <- cbind(1, T, P, if (adjust_pre) P * pre, T * P)
  p <- ncol(X)
  main_col <- 1 + (1:k)
  inter_col <- p - k + (1:k)
  XtXi <- solve(crossprod(X))
  b <- drop(XtXi %*% crossprod(X, y))
  r <- y - drop(X %*% b)
  df <- length(y) - p
  hat <- rowSums((X %*% XtXi) * X)
  w <- r^2 / (1 - hat)^2

  pairs <- pairs_of(k)
  L <- do.call(rbind, lapply(types, function(type) {
    do.call(rbind, lapply(pairs, function(pr) {
      main <- numeric(p); inter <- numeric(p)
      if (pr[1] > 0) { main[main_col[pr[1]]] <- 1; inter[inter_col[pr[1]]] <- 1 }
      if (pr[2] > 0) { main[main_col[pr[2]]] <- main[main_col[pr[2]]] - 1
                       inter[inter_col[pr[2]]] <- inter[inter_col[pr[2]]] - 1 }
      switch(type,
             "ATE (avg over pretest)" = main + 0.5 * inter,
             "Pretest x Treatment" = inter,
             "Treatment | pretested" = main + inter,
             "Treatment | unpretested" = main)
    }))
  }))
  est <- drop(L %*% b)
  A <- L %*% XtXi %*% t(X)
  V_hc3 <- A %*% (w * t(A))
  V_conv <- sum(r^2) / df * (L %*% XtXi %*% t(L))
  list(est = est, V_hc3 = V_hc3, V_conv = V_conv, df = df, k = k,
       n_pairs = length(pairs))
}

# Tests from a joint fit for one covariance: per-contrast results for all
# pairs and, for the comparisons against the control, the omnibus tests.
joint_tests <- function(f, V) {
  k <- f$k
  np <- f$n_pairs
  se <- sqrt(diag(V))
  t <- f$est / se
  p <- 2 * stats::pt(-abs(t), f$df)
  q <- stats::qt(1 - alpha / 2, f$df)
  omni <- vapply(seq_along(types), function(ti) {
    rows <- (ti - 1) * np + seq_len(k)
    bb <- f$est[rows]
    W <- drop(t(bb) %*% solve(V[rows, rows], bb))
    stats::pf(W / k, k, f$df, lower.tail = FALSE)
  }, numeric(1))
  list(est = f$est, se = se, p = p, lo = f$est - q * se, hi = f$est + q * se,
       omni = omni)
}

# Holm-adjusted p-values within each contrast type, for the comparisons
# against the control (the first k pairs) or for all pairs.
holm_by_type <- function(p, k, np, all_pairs) {
  out <- rep(NA_real_, length(p))
  for (ti in seq_along(types)) {
    rows <- (ti - 1) * np + if (all_pairs) seq_len(np) else seq_len(k)
    out[rows] <- stats::p.adjust(p[rows], "holm")
  }
  out
}

# ---- Overlapping four-group analyses (M4) -------------------------------------

# For each pair of conditions, the interaction F test of the 2 x 2 ANOVA of
# its four posttest groups (equal to the t test of the interaction
# coefficient of lm(y ~ treat * pretested) on those groups).
overlapping <- function(des, y) {
  vapply(pairs_of(des$k), function(pr) {
    keep <- des$cond %in% pr
    g <- interaction(des$cond[keep] == pr[1], des$pretested[keep])
    m <- tapply(y[keep], g, mean)
    n <- tapply(y[keep], g, length)
    sse <- sum((y[keep] - m[g])^2)
    dfe <- sum(keep) - 4
    # cells: FALSE.0, TRUE.0, FALSE.1, TRUE.1
    est <- (m["TRUE.1"] - m["FALSE.1"]) - (m["TRUE.0"] - m["FALSE.0"])
    se <- sqrt(sse / dfe * sum(1 / n))
    2 * stats::pt(-abs(est / se), dfe)
  }, numeric(1))
}

# ---- Steyn's (2009) sequence (M5) -----------------------------------------------

# Post hoc tests of every pair of groups, with a pooled SD: Scheffe tests
# (Steyn, 2005), or pairwise t tests with Holm's adjustment, as
# stats::pairwise.t.test().
pairwise_pooled <- function(y, g, posthoc) {
  m <- tapply(y, g, mean)
  n <- tapply(y, g, length)
  G <- length(m)
  dfe <- length(y) - G
  s2 <- sum((y - m[as.character(g)])^2) / dfe
  cmb <- utils::combn(names(m), 2)
  t <- apply(cmb, 2, function(ab) {
    (m[ab[1]] - m[ab[2]]) / sqrt(s2 * (1 / n[ab[1]] + 1 / n[ab[2]]))
  })
  p <- if (posthoc == "scheffe") {
    stats::pf(t^2 / (G - 1), G - 1, dfe, lower.tail = FALSE)
  } else {
    stats::p.adjust(2 * stats::pt(-abs(t), dfe), "holm")
  }
  data.frame(a = cmb[1, ], b = cmb[2, ], p = p, stringsAsFactors = FALSE)
}

oneway_p <- function(y, g) {
  m <- tapply(y, g, mean)
  n <- tapply(y, g, length)
  G <- length(m)
  ssb <- sum(n * (m - mean(y))^2)
  ssw <- sum((y - m[as.character(g)])^2)
  stats::pf((ssb / (G - 1)) / (ssw / (length(y) - G)), G - 1, length(y) - G,
            lower.tail = FALSE)
}

steyn <- function(des, y, posthoc) {
  k <- des$k
  # Groups: Ob_j / Oe_j for treatment j pretested / unpretested; Od, Of.
  g <- ifelse(des$cond == 0,
              ifelse(des$pretested == 1, "Od", "Of"),
              paste0(ifelse(des$pretested == 1, "Ob", "Oe"), des$cond))
  e1_p <- oneway_p(y, g)
  e1 <- e1_p < alpha
  pw <- pairwise_pooled(y, g, posthoc)
  differs <- function(a, b) {
    hit <- (pw$a == a & pw$b == b) | (pw$a == b & pw$b == a)
    pw$p[hit] < alpha
  }
  effective <- vapply(1:k, function(j) {
    e1 && ((differs(paste0("Ob", j), "Od") && differs(paste0("Ob", j), "Of")) ||
           (differs(paste0("Oe", j), "Od") && differs(paste0("Oe", j), "Of")))
  }, logical(1))
  e4_p <- vapply(1:k, function(j) {
    a <- y[g == paste0("Ob", j)]; b <- y[g == paste0("Oe", j)]
    stats::t.test(a, b, var.equal = TRUE)$p.value
  }, numeric(1))
  flag <- any(e4_p < alpha)
  differ <- FALSE
  e5_p <- NA_real_
  if (!flag) {
    gt <- paste0("T", des$cond)[des$cond > 0]
    yt <- y[des$cond > 0]
    pw5 <- pairwise_pooled(yt, gt, posthoc)
    e5_p <- min(pw5$p)
    differ <- any(pw5$p < alpha)
  }
  list(e1_p = e1_p, effective = effective, e4_p = e4_p, flag = flag,
       differ = differ, e5_min_p = e5_p, pairs = pw)
}

# ---- One task: a scenario and a substream ----------------------------------------

run_task <- function(task) {
  s <- scenarios[task$scenario, ]
  file <- file.path(cache_dir, sprintf("s%03d_c%02d.rds", task$scenario, task$chunk))
  if (file.exists(file)) return(readRDS(file))

  assign(".Random.seed", task$seed, envir = .GlobalEnv)
  k <- s$k
  des <- design_of(k, s$n)
  pairs <- pairs_of(k)
  np <- length(pairs)
  truth <- unlist(lapply(seq_along(types), function(ti) {
    vapply(pairs, function(pr) truth_pair(s, k, pr[1], pr[2])[ti], numeric(1))
  }))
  null_pair <- abs(truth) < 1e-12
  vc <- rep(seq_len(np) <= k, length(types))           # rows against the control
  type_of <- rep(seq_along(types), each = np)
  e <- effects[[s$effects]]
  delta <- e$delta(k); sens <- e$sens(k)

  acc <- list()
  add <- function(name, value) acc[[name]] <<- (if (is.null(acc[[name]])) 0 else acc[[name]]) + value
  agree <- list()

  for (r in seq_len(chunk)) {
    d <- simulate_data(des, s)

    fits <- list(adj = joint_fit(des, d$y, d$pre, TRUE),
                 post = joint_fit(des, d$y, d$pre, FALSE))
    res <- list(
      M1 = joint_tests(fits$adj, fits$adj$V_hc3),
      M2 = joint_tests(fits$post, fits$post$V_hc3),
      M3 = joint_tests(fits$post, fits$post$V_conv)
    )

    for (m in names(res)) {
      x <- res[[m]]
      padj_vc <- holm_by_type(x$p, k, np, all_pairs = FALSE)
      padj_all <- holm_by_type(x$p, k, np, all_pairs = TRUE)
      dev <- x$est - truth
      add(paste(m, "bias_sum", sep = "|"), dev)
      add(paste(m, "bias_sq", sep = "|"), dev^2)
      add(paste(m, "se_sum", sep = "|"), x$se)
      add(paste(m, "cover", sep = "|"), as.numeric(x$lo <= truth & truth <= x$hi))
      add(paste(m, "reject", sep = "|"), as.numeric(x$p < alpha))
      add(paste(m, "omni_reject", sep = "|"), as.numeric(x$omni < alpha))
      fw_vc <- vapply(seq_along(types), function(ti) {
        rows <- which(type_of == ti & vc & null_pair)
        if (length(rows)) as.numeric(any(padj_vc[rows] < alpha)) else NA_real_
      }, numeric(1))
      fw_all <- vapply(seq_along(types), function(ti) {
        rows <- which(type_of == ti & null_pair)
        if (length(rows)) as.numeric(any(padj_all[rows] < alpha)) else NA_real_
      }, numeric(1))
      fw_unadj <- vapply(seq_along(types), function(ti) {
        rows <- which(type_of == ti & vc & null_pair)
        if (length(rows)) as.numeric(any(x$p[rows] < alpha)) else NA_real_
      }, numeric(1))
      add(paste(m, "fwer_control", sep = "|"), fw_vc)
      add(paste(m, "fwer_pairwise", sep = "|"), fw_all)
      add(paste(m, "fwer_unadjusted", sep = "|"), fw_unadj)

      if (task$chunk == 1L && r <= agree_glm) {
        agree[[length(agree) + 1L]] <- list(rep = r, method = m, d = d, x = x,
                                            padj_vc = padj_vc, padj_all = padj_all)
      }
    }

    p4 <- overlapping(des, d$y)
    sens_null <- null_pair[type_of == 2]
    add("M4|any_control", as.numeric(any(p4[seq_len(k)][sens_null[seq_len(k)]] < alpha)))
    add("M4|any_pair", as.numeric(any(p4[sens_null] < alpha)))
    add("M4|n_null_control", as.numeric(any(sens_null[seq_len(k)])))
    add("M4|n_null_pair", as.numeric(any(sens_null)))

    # M5: Scheffe post hoc tests (amendment 1); M5h: Holm-adjusted pairwise
    # t tests, as first specified.
    null_treat <- delta == 0 & sens == 0
    for (m5 in c("M5", "M5h")) {
      st <- steyn(des, d$y, if (m5 == "M5") "scheffe" else "holm")
      add(paste0(m5, "|e1"), as.numeric(st$e1_p < alpha))
      add(paste0(m5, "|null_effective"),
          if (any(null_treat)) as.numeric(any(st$effective[null_treat])) else 0)
      add(paste0(m5, "|flag"), as.numeric(st$flag))
      add(paste0(m5, "|differ"), as.numeric(st$differ))
      if (task$chunk == 1L && r <= agree_steyn) {
        agree[[length(agree) + 1L]] <- list(rep = r, method = m5, d = d, steyn = st)
      }
    }
  }

  out <- list(scenario = task$scenario, chunk = task$chunk, sums = acc,
              agree = agree, n = chunk)
  saveRDS(out, file)
  out
}

# ---- Run ---------------------------------------------------------------------------

if (sys.nframe() == 0L) {
  commit <- tryCatch(system("git rev-parse HEAD", intern = TRUE), error = function(e) NA)
  started <- Sys.time()

  RNGkind("L'Ecuyer-CMRG")
  set.seed(45045)
  tasks <- list()
  stream <- .Random.seed
  for (sc in scenarios$scenario) {
    stream <- parallel::nextRNGStream(stream)
    sub <- stream
    for (ch in seq_len(reps / chunk)) {
      tasks[[length(tasks) + 1L]] <- list(scenario = sc, chunk = ch, seed = sub)
      sub <- parallel::nextRNGSubStream(sub)
    }
  }

  workers <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
  cl <- parallel::makeCluster(workers)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterEvalQ(cl, RNGkind("L'Ecuyer-CMRG"))
  parallel::clusterExport(cl, c(
    "scenarios", "effects", "types", "cell_sizes", "design_of", "simulate_data",
    "truth_pair", "pairs_of", "pair_names", "joint_fit", "joint_tests",
    "holm_by_type", "overlapping", "pairwise_pooled", "oneway_p", "steyn",
    "run_task", "cache_dir", "chunk", "alpha", "agree_glm", "agree_steyn"
  ))
  results <- parallel::parLapplyLB(cl, tasks, run_task, chunk.size = 1)

  source(file.path(out_dir, "ngroup-summarize.R"))
  summarize_run(results, commit = commit, started = started, workers = workers)
}
