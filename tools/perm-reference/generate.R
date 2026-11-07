# Generate the reference set for issue #134 (design: p134/design.md).
#
#   Rscript tools/perm-reference/generate.R [full|test] [out_dir]
#
# Writes tools/perm-reference/perm-reference-full.rds and generate.log,
# tests/testthat/fixtures/perm-reference.rds (the compact set), and, when
# out_dir is given, reference-full.rds and reference.csv there.
#
# Rules (design Section 11): this script does not load the package under
# test, fits no reference statistic with glm(), and decides no equality with
# a tolerance; every seed is recorded.

args <- commandArgs(TRUE)
tier <- if (length(args) >= 1) args[1] else "full"
out_dir <- if (length(args) >= 2) args[2] else NULL
REF_DIR <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))),
                         winslash = "/")
ROOT <- normalizePath(file.path(REF_DIR, "..", ".."), winslash = "/")
REF_PYTHON <- Sys.getenv("PYTHON", Sys.which("python"))
stopifnot(nzchar(REF_PYTHON))
stopifnot(identical(RNGkind(), c("Mersenne-Twister", "Inversion", "Rejection")))
source(file.path(REF_DIR, "lib.R"))
source(file.path(REF_DIR, "reference.R"))
source(file.path(REF_DIR, "worlds.R"))
if ("solomonR" %in% loadedNamespaces()) stop("the generator must not load solomonR")

LOG <- file.path(REF_DIR, "generate.log")
cat("", file = LOG)
logf <- function(...) {
  msg <- sprintf(...)
  cat(msg, "\n", sep = "", file = LOG, append = TRUE)
  cat(msg, "\n", sep = "")
}

# ---------------------------------------------------------------------------
# Self-checks (design Section 11, rule 3): not used for any decision.

self_check <- function(w, nc, tabs, theta = NULL) {
  nP <- length(nc$P$keys)
  out <- NULL
  set <- tabs$tables
  for (cn in ALLC) for (s in STATS) {
    tb <- set[[cn]][[s]]
    cand <- which(tb$code == "finite" & is.finite(tb$approx))
    if (!length(cand)) next
    pick <- unique(cand[round(seq(1, length(cand), length.out = min(2, length(cand))))])
    for (i in pick) {
      iU <- (i - 1) %/% nP + 1
      iP <- (i - 1) %% nP + 1
      rows <- class_rows(w, nc$U$comps[iU, ], nc$P$comps[iP, ])
      ds <- design_of(rows)
      rep_i <- rep(seq_len(nrow(ds$X)), ds$f)
      X <- ds$X[rep_i, , drop = FALSE]
      y <- ds$y[rep_i]
      off <- log(ds$e[rep_i])
      q <- qr(X)
      X <- X[, q$pivot[seq_len(q$rank)], drop = FALSE]
      L <- c(cell_L(cn), rep(0, ncol(ds$X) - 4))[q$pivot[seq_len(q$rank)]]
      fam <- switch(w$family, gaussian = stats::gaussian(), binomial = stats::binomial(),
                    poisson = stats::poisson(), negbin = MASS::negative.binomial(theta))
      one <- function(ctrl) {
        g <- suppressWarnings(stats::glm(y ~ 0 + X, family = fam, offset = off, control = ctrl))
        b <- stats::coef(g)
        if (anyNA(b)) return(NA_real_)
        V <- sandwich::vcovHC(g, type = "HC3")
        est <- sum(L * b)
        if (s == "difference") est else est / sqrt(drop(t(L) %*% V %*% L))
      }
      val <- one(stats::glm.control())
      val_tight <- one(stats::glm.control(epsilon = 1e-14, maxit = 200))
      if (is.na(val)) next
      rd <- abs(val - tb$approx[i]) / max(1, abs(tb$approx[i]))
      rdt <- abs(val_tight - tb$approx[i]) / max(1, abs(tb$approx[i]))
      out <- rbind(out, data.frame(world = w$id, family = w$family, contrast = cn, statistic = s, class = i,
                                   reference = tb$approx[i], glm_sandwich = val, rel_diff = rd,
                                   rel_diff_tight = rdt))
    }
  }
  out
}

# ---------------------------------------------------------------------------
# Participant worlds

ref_participant <- function(w) {
  t0 <- proc.time()[["elapsed"]]
  w <- prep_world(w)
  spaces <- build_spaces(w)
  nc <- needed_comps(w, spaces)
  nP <- length(nc$P$keys)
  thetas <- NULL
  tabs <- list()
  if (w$family == "negbin") {
    thetas <- list()
    for (o in names(w$obs)) {
      rows <- class_rows(w, labeling_comp(w, w$obs[[o]], "U"), labeling_comp(w, w$obs[[o]], "P"))
      ds <- design_of(rows)
      th <- nb_theta(ds$X, ds$y, ds$f, log(ds$e))
      thetas[[o]] <- th
      if (th$boundary || th$theta > 1e3) {
        tabs[[o]] <- NULL
        next
      }
      tabs[[o]] <- lapply(c(lo = 1 - 1e-6, mid = 1, hi = 1 + 1e-6), function(m)
        ref_world_tables(w, nc, th$theta * m))
    }
  } else {
    tb <- ref_world_tables(w, nc)
    for (o in names(w$obs)) tabs[[o]] <- list(mid = tb)
  }
  checks <- NULL
  cases <- list()
  for (o in names(w$obs)) {
    th <- if (!is.null(thetas)) thetas[[o]]$theta else NA_real_
    for (cn in ALLC) for (s in STATS) {
      base <- data.frame(world = w$id, block = w$block, items = w$items, family = w$family,
                         path = if (w$prep$sep) "separable" else "shared covariate",
                         hc3_path = if (any(vapply(w$covs, function(v) is.numeric(w$data[[v]]), logical(1))))
                           ".glm_influence" else "sandwich",
                         obs = o, contrast = cn, statistic = s, theta = th, stringsAsFactors = FALSE)
      if (is.null(tabs[[o]])) {
        cases[[length(cases) + 1]] <- cbind(base, space = "all", outcome = "no-reference:NB-BND",
                                            stringsAsFactors = FALSE)
        next
      }
      rule <- data_rule(w, cn)
      for (sp in c("all", if (!is.null(spaces$model)) "model")) {
        wU <- if (sp == "all") space_weights(spaces$all$U, nc$U$keys) else space_weights(spaces$model$U[[o]], nc$U$keys)
        wP <- if (sp == "all") space_weights(spaces$all$P, nc$P$keys) else space_weights(spaces$model$P[[o]], nc$P$keys)
        res <- lapply(tabs[[o]], function(tt) {
          tab <- tt$tables[[cn]][[s]]
          rl <- tab$rel[[o]]
          list(tab = tab, rl = rl, oi = tt$obs_idx[[o]],
               cp = case_counts(tab, rl$rel, wU, wP, nP, s, "proposal"),
               cd = case_counts(tab, rl$rel, wU, wP, nP, s, "drop"),
               ci = case_counts(tab, rl$rel, wU, wP, nP, s, "inf"))
        })
        m <- res$mid
        oc <- m$tab$code[m$oi]
        outcome <- "p"
        if (rule == "ERR") outcome <- "error:no-variation"
        if (!is.null(m$rl$obs_undefined)) outcome <- paste0("error:observed-", m$rl$obs_undefined)
        else if (!is.null(m$cp$error) && outcome == "p") outcome <- paste0("error:", m$cp$error)
        if (is.null(m$rl$obs_undefined) && (isTRUE(m$cp$undecided) || !is.null(m$rl$unexplained)))
          outcome <- "excluded:undecided"
        stable <- TRUE
        if (length(res) == 3) {
          stable <- identical(res$lo$cp$b, m$cp$b) && identical(res$hi$cp$b, m$cp$b) &&
            identical(res$lo$cd$b, m$cd$b) && identical(res$hi$cd$b, m$cd$b)
          if (!stable) outcome <- "excluded:theta-unstable"
        }
        num <- function(x) if (is.null(x)) NA_real_ else as.numeric(x)
        cases[[length(cases) + 1]] <- cbind(base, data.frame(
          space = sp, outcome = outcome, data_rule = rule, obs_code = oc,
          T_obs = m$tab$approx[m$oi], b = m$cp$b, M = m$cp$M,
          p = num(m$cp$b) / num(m$cp$M), b_tied = m$cp$b_tied, b_strict = m$cp$b_strict,
          b_drop = m$cd$b, M_drop = m$cd$M,
          p_drop = if (oc %in% c("ZZ", "NE", "H1") || (oc == "SEP-L" && s == "studentized")) NA_real_ else
            num(m$cd$b) / num(m$cd$M),
          b_inf = m$ci$b, M_inf = m$ci$M, p_inf = num(m$ci$b) / num(m$ci$M),
          n_inf = m$cp$n_inf, n_sepl = m$cp$n_sepl, n_h1 = m$cp$n_h1, n_ne = m$cp$n_ne,
          n_zz = m$cp$n_zz, n_lim = m$cp$n_lim, n_alias = m$cp$n_alias,
          band_n = if (is.null(m$rl$band_n)) 0 else m$rl$band_n,
          how = if (is.null(m$rl$how)) NA_character_ else m$rl$how[1],
          theta_stable = stable, stringsAsFactors = FALSE))
      }
    }
    if (!is.null(tabs[[o]])) checks <- rbind(checks, self_check(w, nc, tabs[[o]]$mid, th))
  }
  cases <- do.call(rbind, lapply(cases, function(z) {
    for (nm in setdiff(CASE_COLS, names(z))) z[[nm]] <- NA
    z[CASE_COLS]
  }))
  # compact per-class store: codes, signed values, relation to each observed
  # class (rel: -1, 0, 1, NA) for the proposal space
  store <- lapply(names(w$obs), function(o) {
    if (is.null(tabs[[o]])) return(NULL)
    tt <- tabs[[o]]$mid
    lapply(stats::setNames(ALLC, ALLC), function(cn) lapply(stats::setNames(STATS, STATS), function(s) {
      tab <- tt$tables[[cn]][[s]]
      list(code = tab$code, approx = tab$approx, rel = as.integer(tab$rel[[o]]$rel))
    }))
  })
  names(store) <- names(w$obs)
  list(world = w[c("id", "block", "items", "family", "data", "pre", "covs", "exposure", "obs", "seed", "mode",
                   "reps", "base", "shift")],
       types = lapply(w$prep$cond, function(ck) list(idx = ck$idx, row_type = ck$row_type, m = ck$m,
                                                     n_out = ck$n_out)),
       keys = list(U = nc$U$keys, P = nc$P$keys),
       weights = list(all = list(U = space_weights(spaces$all$U, nc$U$keys), P = space_weights(spaces$all$P, nc$P$keys))),
       store = store, cases = cases, thetas = thetas, checks = checks,
       secs = proc.time()[["elapsed"]] - t0)
}

CASE_COLS <- c("world", "block", "items", "family", "path", "hc3_path", "obs", "contrast", "statistic", "theta",
               "space", "outcome", "data_rule", "obs_code", "T_obs", "b", "M", "p", "b_tied", "b_strict",
               "b_drop", "M_drop", "p_drop", "b_inf", "M_inf", "p_inf", "n_inf", "n_sepl", "n_h1", "n_ne",
               "n_zz", "n_lim", "n_alias", "band_n", "how", "theta_stable")

# Location worlds: the reference is the base world's, after asserting that
# the stored data are an exact shift of the base data (design Section 7).
ref_shifted <- function(w, base_ref, base_w) {
  wh <- w$shift$what
  col <- switch(wh, y = "y", pre = "pre", x = "x", year = "year", xmat = "xm")
  a <- as.numeric(as.matrix(w$data[[col]]))
  b <- as.numeric(as.matrix(base_w$data[[col]]))
  ok <- isTRUE(all.equal(is.na(a), is.na(b))) && all((a - w$shift$c)[!is.na(a)] == b[!is.na(b)])
  grid <- all(vapply(a[!is.na(a)], function(v) {
    k <- 0
    while (v * 2^k != round(v * 2^k) && k < 60) k <- k + 1
    abs(v * 2^k) < 2^53
  }, logical(1)))
  if (!ok || !grid) stop("shift is not exact: ", w$id)
  r <- base_ref
  r$world <- w[c("id", "block", "items", "family", "data", "pre", "covs", "exposure", "obs", "seed", "mode",
                 "reps", "base", "shift")]
  r$cases$world <- w$id
  r$cases$block <- "L"
  r$cases$items <- w$items
  r$cases$how <- paste0("invariance: reference of ", base_w$id, " (exact shift by ", format(w$shift$c, digits = 17),
                        " of ", wh, ")")
  r$shift_check <- list(exact = ok, dyadic_below_2_53 = grid)
  r$band_rows <- NULL
  r$band_fail <- NULL
  r$checks <- NULL
  r
}

# ---------------------------------------------------------------------------
# Clustered worlds (exact Stage 1 scores in exact.py)

ref_cluster <- function(w) {
  t0 <- proc.time()[["elapsed"]]
  d <- w$data
  inm <- !is.na(d$y)
  d <- d[inm, ]
  fam <- w$family
  stage1 <- if (fam == "gaussian") "ols" else "means"
  if (fam != "gaussian" && (isTRUE(w$pre) || length(w$covs))) stop("cluster world needs an exact Stage 1: ", w$id)
  z <- cbind(if (isTRUE(w$pre)) ifelse(d$pretested == 1, d$pre, 0), if (length(w$covs)) as.matrix(d[w$covs]))
  e <- if (is.null(w$exposure)) rep(1, nrow(d)) else d[[w$exposure]]
  rows <- lapply(seq_len(nrow(d)), function(i) list(cl = as.character(d$cl[i]), y = dyad(d$y[i]), p = d$pretested[i],
                                                    z = if (is.null(z)) list() else as.list(dyad(z[i, ])),
                                                    e = dyad(e[i])))
  cases <- list()
  cl_treat <- tapply(d$treat, d$cl, `[`, 1)
  cl_pre <- tapply(d$pretested, d$cl, `[`, 1)
  for (cn in ALLC) {
    wt <- c(P = CW[[cn]][["P"]], U = CW[[cn]][["U"]])
    if (w$design == "A") {
      strata <- list()
      for (k in c("P", "U")) if (wt[[k]] != 0) {
        ids <- names(cl_pre)[cl_pre == (k == "P")]
        strata[[length(strata) + 1]] <- list(ids = as.list(ids), t = as.list(unname(cl_treat[ids])), w = dyad(wt[[k]]))
      }
      inp <- list(rows = rows, stage1 = stage1, design = "A", strata = strata)
    } else {
      ids <- names(cl_treat)
      inp <- list(rows = rows, stage1 = stage1, design = "B",
                  strata = list(list(ids = as.list(ids), t = as.list(unname(cl_treat[ids])), w = "1")),
                  wB = list(dyad(wt[["P"]]), dyad(wt[["U"]])))
    }
    r <- py_result("cluster", inp)
    o <- r$obs + 1
    for (s in STATS) {
      rr <- r[[if (s == "studentized") "stud" else "diff"]]
      rank <- vapply(rr$rank, function(v) if (is.null(v)) NA_real_ else as.numeric(v), 0)
      code <- unlist(rr$code)
      M <- length(rank)
      if (is.na(rank[o])) {
        outcome <- "error:observed-undefined"
        b <- NA
      } else {
        b <- sum(rank >= rank[o], na.rm = TRUE)
        outcome <- if (code[o] == "INF-C") "p (observed +Inf)" else "p"
      }
      cases[[length(cases) + 1]] <- data.frame(
        world = w$id, block = w$block, items = w$items, family = fam, path = paste("cluster, design", w$design),
        hc3_path = NA, obs = "gen", contrast = cn, statistic = s, theta = NA, space = "allocations",
        outcome = outcome, data_rule = NA, obs_code = code[o],
        T_obs = vapply(rr$approx, function(v) if (is.null(v)) NA_real_ else as.numeric(v), 0)[o],
        b = b, M = M, p = b / M, b_tied = sum(rank == rank[o], na.rm = TRUE),
        b_strict = sum(rank > rank[o], na.rm = TRUE), b_drop = NA, M_drop = NA, p_drop = NA, b_inf = NA,
        M_inf = NA, p_inf = NA, n_inf = sum(code == "INF-C"), n_sepl = 0, n_h1 = 0, n_ne = 0,
        n_zz = sum(code == "ZZ"), n_lim = 0, n_alias = 0, band_n = 0,
        how = "exact rationals (Stage 1 and allocations)",
        theta_stable = NA, stringsAsFactors = FALSE)
    }
  }
  list(world = w[c("id", "block", "items", "family", "data", "pre", "covs", "exposure", "obs", "seed", "mode",
                   "cluster", "design", "base", "shift")],
       cases = do.call(rbind, cases), secs = proc.time()[["elapsed"]] - t0)
}

# ---------------------------------------------------------------------------
# Run

t_all <- proc.time()[["elapsed"]]
W <- build_grid(tier)
if (nzchar(Sys.getenv("REF_ONLY"))) W <- W[grepl(Sys.getenv("REF_ONLY"), names(W))]
logf("perm-reference generator, tier %s, %d worlds; R %s; Python %s", tier, length(W), getRversion(),
     paste(system2(REF_PYTHON, "--version", stdout = TRUE), collapse = ""))
REFS <- list()
# Optional checkpoints (REF_CKPT, a folder): a participant world's result is
# reused only when the computing code and the world are byte-for-byte those
# that produced it, so a resumed run gives what a fresh run gives.
CKPT <- Sys.getenv("REF_CKPT")
code_hash <- digest::digest(list(
  lapply(c("lib.R", "reference.R", "exact.py"), function(f) readLines(file.path(REF_DIR, f))),
  deparse(ref_participant), deparse(self_check), CASE_COLS), algo = "sha1")
for (w in W) {
  nb0 <- length(BAND_LOG$rows)
  nf0 <- length(BAND_LOG$fail)
  ck <- if (nzchar(CKPT) && is.null(w$base) && !identical(w$mode, "cluster"))
    file.path(CKPT, paste0(w$id, "-", digest::digest(list(code_hash, w), algo = "sha1"), ".rds")) else ""
  r <- tryCatch({
    if (nzchar(ck) && file.exists(ck)) {
      readRDS(ck)
    } else if (identical(w$mode, "cluster")) {
      if (!is.null(w$base)) ref_shifted(w, REFS[[w$base]], W[[w$base]]) else ref_cluster(w)
    } else if (!is.null(w$base)) {
      ref_shifted(w, REFS[[w$base]], W[[w$base]])
    } else {
      ref_participant(w)
    }
  }, error = function(e) {
    logf("ERROR in %s: %s", w$id, conditionMessage(e))
    NULL
  })
  if (is.null(r)) next
  if (is.null(r$band_rows)) {
    r$band_rows <- if (length(BAND_LOG$rows) > nb0) do.call(rbind, BAND_LOG$rows[(nb0 + 1):length(BAND_LOG$rows)])
    r$band_fail <- if (length(BAND_LOG$fail) > nf0) do.call(rbind, BAND_LOG$fail[(nf0 + 1):length(BAND_LOG$fail)])
    if (!is.null(r$band_rows)) r$band_rows$world <- w$id
  }
  if (nzchar(ck) && !file.exists(ck)) saveRDS(r, ck)
  REFS[[w$id]] <- r
  cs <- r$cases
  logf("%-32s %6.1f s  cases %3d  p-cases %3d  errors %2d  excluded %2d  band %d", w$id, r$secs, nrow(cs),
       sum(cs$outcome == "p" | grepl("^p ", cs$outcome)), sum(grepl("^error", cs$outcome)),
       sum(grepl("^excluded", cs$outcome)), sum(cs$band_n, na.rm = TRUE))
}
cases <- do.call(rbind, lapply(REFS, function(r) {
  z <- r$cases
  for (nm in setdiff(CASE_COLS, names(z))) z[[nm]] <- NA
  z[CASE_COLS]
}))
rownames(cases) <- NULL
cases$case_id <- sprintf("%s|%s|%s|%s|%s", cases$world, cases$obs, CSHORT[cases$contrast], cases$statistic, cases$space)
checks <- do.call(rbind, lapply(REFS, `[[`, "checks"))
seeds <- do.call(rbind, SEED_LOG$rows)
band <- do.call(rbind, lapply(REFS, `[[`, "band_rows"))
band_fail <- do.call(rbind, lapply(REFS, `[[`, "band_fail"))
if (!is.null(band_fail)) logf("70-digit fits that failed (their classes stay undecided): %d", nrow(band_fail))
if (!is.null(band)) logf("band: %d comparisons at 50 digits; %d unexplained (gap below 1e-35); %d failed the eps check",
                         nrow(band), sum(band$rel_gap50 < 1e-35), sum(!band$eps_check))
logf("worlds %d, cases %d; python calls %d (%.0f s); total %.0f s", length(REFS), nrow(cases), PY_CALLS$n,
     PY_CALLS$secs, proc.time()[["elapsed"]] - t_all)
if (!is.null(checks)) logf(paste("self-checks: %d classes; largest relative difference from glm() + sandwich:",
                                 "Gaussian %.2e; iterated %.2e at glm()'s default convergence, %.2e at epsilon 1e-14"),
                           nrow(checks), max(c(0, checks$rel_diff[checks$family == "gaussian"])),
                           max(c(0, checks$rel_diff[checks$family != "gaussian"])),
                           max(c(0, checks$rel_diff_tight[checks$family != "gaussian"]), na.rm = TRUE))
full <- list(created = format(Sys.time(), "%Y-%m-%d %H:%M:%S"), tier = tier, R = R.version.string,
             packages = vapply(c("MASS", "sandwich", "lpSolve", "jsonlite"), function(p) as.character(utils::packageVersion(p)), ""),
             python = system2(REF_PYTHON, "--version", stdout = TRUE), rng = RNGkind(),
             refs = REFS, cases = cases, checks = checks, seeds = seeds, band = band, band_fail = band_fail)
saveRDS(full, file.path(REF_DIR, "perm-reference-full.rds"))

# Compact fixture: the test tier (worlds with at most 400 labelings and the
# blocks D, L, G), with per-class codes and relations, and the case table.
compact <- lapply(REFS, function(r) {
  M <- suppressWarnings(max(as.numeric(r$cases$M), na.rm = TRUE))
  keep <- r$world$block %in% c("D", "L", "G") || (is.finite(M) && M <= 400)
  if (!keep) return(NULL)
  out <- r[c("world", "types", "keys", "store", "cases", "weights", "thetas")]
  # worlds with more than 5,000 classes keep their cases but not the
  # per-class tables (in the full file), to keep the fixture small
  if (length(r$keys$U) * length(r$keys$P) > 5000) out$store <- NULL
  out
})
compact <- compact[!vapply(compact, is.null, logical(1))]
fx <- file.path(ROOT, "tests", "testthat", "fixtures")
dir.create(fx, showWarnings = FALSE, recursive = TRUE)
saveRDS(list(created = full$created, cases = cases, worlds = compact), file.path(fx, "perm-reference.rds"),
        compress = "xz", version = 2)
logf("compact fixture: %d worlds, %.0f KB", length(compact), file.size(file.path(fx, "perm-reference.rds")) / 1024)
if (!is.null(out_dir)) {
  saveRDS(full, file.path(out_dir, "reference-full.rds"))
  utils::write.csv(cases, file.path(out_dir, "reference.csv"), row.names = FALSE)
  if (!is.null(seeds)) utils::write.csv(seeds, file.path(out_dir, "reference-seeds.csv"), row.names = FALSE)
  if (!is.null(band)) utils::write.csv(band, file.path(out_dir, "reference-band.csv"), row.names = FALSE)
}
