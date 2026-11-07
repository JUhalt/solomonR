# The grid of worlds (design Section 9). Each world is a data set as passed
# to fit_solomon_glm(), the model, and one or more observed labelings. Data
# are integers or dyadic rationals, so that every value is exact in binary
# and the generator can compute with them exactly. Every random element uses
# a recorded seed: data worlds 1347001-1347999 (a block of 10 per world: the
# first is used, the rest replace it in order when it fails a screen), and
# near-tie constructions 1348001-1348499.

SEED_LOG <- new.env()
SEED_LOG$rows <- list()
log_seed <- function(world, seed, purpose, used) {
  SEED_LOG$rows[[length(SEED_LOG$rows) + 1]] <- data.frame(world = world, seed = seed, purpose = purpose,
                                                           used = used, stringsAsFactors = FALSE)
}

groups4 <- function(n) {
  if (length(n) == 1) n <- rep(n, 4)
  g <- rep(1:4, times = n)  # 1 PT, 2 PC, 3 UT, 4 UC
  data.frame(treat = c(1, 0, 1, 0)[g], pretested = c(1, 1, 0, 0)[g], g = g)
}

dyadic_pre <- function(n) round(stats::rnorm(n, 0, 2) * 4) / 4

# Draw the outcome for a family.
draw_y <- function(fam, eta, nb_size = 1) {
  n <- length(eta)
  switch(fam,
         gaussian = round(4 * (eta + stats::rnorm(n))),
         binomial = stats::rbinom(n, 1, stats::plogis(eta)),
         poisson = stats::rpois(n, exp(eta)),
         negbin = stats::rnbinom(n, size = nb_size, mu = exp(eta)))
}

emptied_obs <- function(fam, d, treat) {
  if (fam == "gaussian") return(FALSE)
  for (p in 0:1) for (a in 0:1) {
    sel <- d$pretested == p & treat == a & !is.na(d$y)
    S <- sum(d$y[sel])
    if (S == 0 || (fam == "binomial" && S == sum(sel))) return(TRUE)
  }
  FALSE
}

# Negative binomial worlds are screened for an interior theta (design 6 and
# T5): theta* from nb_theta() (lib.R; the highest local maximum of the
# profile, above the Poisson limit) on the observed labeling's data, and
# glm.nb, fitted as the package fits it, converging without the warnings
# behind solomonR_theta_boundary_warning and without an error (glmnb_check()).
nb_interior <- function(dat, treat, pre, covs, exposure) {
  ex <- if (isTRUE(exposure)) "ex" else NULL
  cd <- cell_design(dat, treat, pre, covs, ex)
  th <- nb_theta(cd$X, cd$y, rep(1, length(cd$y)), cd$off)
  if (th$boundary || !is.finite(th$theta) || th$theta > 1e3) return(FALSE)
  gn <- glmnb_check(dat, treat, pre, covs, ex)
  is.null(gn$error) && !isTRUE(gn$boundary)
}

random_labeling <- function(d) {
  t <- d$treat
  for (p in 0:1) {
    i <- which(d$pretested == p)
    t[i] <- sample(d$treat[i])
  }
  t
}

# A standard world: four groups of n, a pretest covariate or not, an
# optional shared covariate (numeric or two-level factor), optional exposure.
# screen_finite (blocks M and L; design T4): every observed labeling has a
# maximum-likelihood estimate (no separation and no recession), certified by
# mle_exists() (lib.R), so that its statistic is finite.
std_world <- function(id, block, fam, n, pre, cov = c("none", "num", "fac"), exposure = FALSE,
                      seed, items = "", nb_size = 1, eta0 = NULL, extra_obs = TRUE, screen_theta = TRUE,
                      screen_finite = FALSE) {
  cov <- match.arg(cov)
  finite_ok <- function(dat, t, covs) {
    !screen_finite || mle_exists(dat, t, fam, pre, covs, if (exposure) "ex" else NULL)
  }
  for (s in seed + 0:9) {
    set.seed(s)
    d <- groups4(n)
    N <- nrow(d)
    # distinct dyadic pretests (multiples of 1/4)
    prev <- rep(NA_real_, N)
    prev[d$pretested == 1] <- sample(seq(-6, 6, by = 0.25), sum(d$pretested == 1))
    for (attempt in 1:50) {
    x <- round(stats::rnorm(N) * 4) / 4
    fx <- factor(ifelse(stats::runif(N) < 0.5, "a", "b"))
    ex <- sample(1:3, N, replace = TRUE)
    base <- if (is.null(eta0)) switch(fam, gaussian = 0, binomial = -0.2, poisson = 0.8, negbin = 1) else eta0
    eta <- base + 0.6 * d$treat + ifelse(is.na(prev), 0, 0.4 * prev) +
      (if (cov == "num") 0.3 * x else 0) + (if (cov == "fac") 0.3 * (fx == "b") else 0) +
      (if (exposure) log(ex) else 0)
    y <- draw_y(fam, eta, nb_size)
    dat <- data.frame(y = y, treat = d$treat, pretested = d$pretested)
    if (pre) dat$pre <- prev
    covs <- character(0)
    if (cov == "num") {
      dat$x <- x
      covs <- "x"
    }
    if (cov == "fac") {
      dat$f <- fx
      covs <- "f"
    }
    if (exposure) dat$ex <- ex
    obs <- list(gen = dat$treat)
    if (extra_obs) {
      # the first of up to 200 random labelings that empties no arm
      for (r in 1:200) {
        alt <- random_labeling(dat)
        if (!emptied_obs(fam, dat, alt) && !identical(alt, dat$treat) && finite_ok(dat, alt, covs)) break
      }
      obs$alt <- alt
    }
    bad <- any(vapply(obs, function(t) emptied_obs(fam, dat, t), logical(1)))
    for (p in 0:1) if (length(unique(dat$y[dat$pretested == p])) < 2) bad <- TRUE
    if (cov == "fac" && length(unique(fx)) < 2) bad <- TRUE
    if (!bad && screen_finite) bad <- !all(vapply(obs, function(t) finite_ok(dat, t, covs), logical(1)))
    if (!bad && fam == "negbin" && screen_theta) {
      bad <- !all(vapply(obs, function(t) nb_interior(dat, t, pre, covs, exposure), logical(1)))
    }
    if (bad) next
    log_seed(id, s, sprintf(paste0("data (outcome draw %d in the seed's stream; earlier draws emptied an arm, left a ",
                                   "condition constant%s)"), attempt,
                            paste0(if (screen_finite) ", had no estimate for an observed labeling" else "",
                                   if (fam == "negbin" && screen_theta) ", or had theta at the boundary" else "")), TRUE)
    return(list(id = id, block = block, items = items, family = fam, data = dat, pre = pre,
                covs = covs, exposure = if (exposure) "ex" else NULL, obs = obs, seed = s,
                mode = "A"))
    }
    log_seed(id, s, "data: 50 outcome draws failed the screens, replaced", FALSE)
  }
  stop("no seed in the block passed the screens: ", id)
}

# Rows given literally, four groups in the order PT, PC, UT, UC.
lit_world <- function(id, block, fam, ylist, pre = NULL, items = "", mode = "C", seed = NA) {
  n <- lengths(ylist)
  d <- groups4(n)
  dat <- data.frame(y = unlist(ylist), treat = d$treat, pretested = d$pretested)
  if (!is.null(pre)) dat$pre <- pre
  list(id = id, block = block, items = items, family = fam, data = dat, pre = !is.null(pre),
       covs = character(0), exposure = NULL, obs = list(gen = dat$treat), seed = seed, mode = mode)
}

build_grid <- function(tier = c("full", "test")) {
  tier <- match.arg(tier)
  W <- list()
  add <- function(w) W[[w$id]] <<- w
  s <- 1347001
  nxt <- function() {
    v <- s
    s <<- s + 10
    v
  }
  fams <- c("gaussian", "binomial", "poisson", "negbin")
  # --- Block M: main
  for (fam in fams) {
    for (sv in c("S0", "S1", "S2", "S3")) {
      for (rep in c("a", "b")) {
        add(std_world(sprintf("M-%s-%s-%s", sv, fam, rep), "M", fam, 3, pre = sv != "S0",
                      cov = switch(sv, S0 = "none", S1 = "none", S2 = "num", S3 = "fac"), seed = nxt(),
                      items = "1", screen_finite = TRUE))
      }
    }
    add(std_world(sprintf("M-S1-%s-n4", fam), "M", fam, 4, pre = TRUE, seed = nxt(), items = "1", screen_finite = TRUE))
  }
  add(std_world("M-Cexp-poisson-a", "M", "poisson", 3, pre = TRUE, exposure = TRUE, seed = nxt(), screen_finite = TRUE))
  add(std_world("M-Cexp-negbin-a", "M", "negbin", 3, pre = TRUE, exposure = TRUE, seed = nxt(), screen_finite = TRUE))

  # --- Block D: the issue's cases
  # D1a: binary outcome, binary pretest, 20 per group (item 1), compositions.
  sd <- nxt()
  set.seed(sd)
  d <- groups4(20)
  pre <- ifelse(d$pretested == 1, stats::rbinom(80, 1, 0.5), NA)
  y <- stats::rbinom(80, 1, stats::plogis(-0.3 + 0.5 * d$treat + ifelse(is.na(pre), 0, 0.8 * pre)))
  log_seed("D1a-binomial", sd, "data", TRUE)
  W[["D1a-binomial"]] <- list(id = "D1a-binomial", block = "D", items = "1", family = "binomial",
                              data = data.frame(y = y, treat = d$treat, pretested = d$pretested, pre = pre),
                              pre = TRUE, covs = character(0), exposure = NULL, obs = list(gen = d$treat),
                              seed = sd, mode = "C", reps = 999)
  # D1b: integer pretest, 3 per group (binomial).
  w <- std_world("D1b-binomial", "D", "binomial", 3, pre = TRUE, seed = nxt(), items = "1")
  w$data$pre <- ifelse(w$data$pretested == 1, round(w$data$pre), NA)
  add(w)
  # D2: logistic, 10 per group, successes 6, 4, 4, 0 (item 2).
  bin <- function(n, k) lapply(k, function(kk) rep(c(1, 0), c(kk, n - kk)))
  add(lit_world("D2-binomial", "D", "binomial", bin(10, c(6, 4, 4, 0)), items = "2"))
  # D3: Poisson, the issue's unpretested counts (item 3); pretested counts drawn.
  sd <- nxt()
  set.seed(sd)
  log_seed("D3-poisson", sd, "data (pretested counts only)", TRUE)
  add(lit_world("D3-poisson", "D", "poisson",
                list(stats::rpois(10, 1.2), stats::rpois(10, 1), c(1, 0, 0, 2, 0, 0, 1, 0, 0, 0),
                     c(0, 0, 1, 0, 0, 0, 0, 0, 0, 0)), items = "3", seed = sd))
  # D4: Gaussian binary scores, groups 5, 5, 4, 4, 1 of 4 against 3 of 4
  # unpretested successes (item 4); pretested scores drawn.
  sd <- nxt()
  set.seed(sd)
  log_seed("D4-gaussian", sd, "data (pretested scores only)", TRUE)
  add(lit_world("D4-gaussian", "D", "gaussian",
                list(stats::rbinom(5, 1, 0.6), stats::rbinom(5, 1, 0.4), c(1, 0, 0, 0), c(1, 1, 1, 0)),
                items = "4", mode = "A", seed = sd))
  # D5: 5 per group, every unpretested score 0.3 (item 5); pretested drawn.
  sd <- nxt()
  set.seed(sd)
  log_seed("D5-gaussian", sd, "data (pretested scores only)", TRUE)
  add(lit_world("D5-gaussian", "D", "gaussian",
                list(round(stats::rnorm(5) * 4) / 4, round(stats::rnorm(5) * 4) / 4, rep(0.3, 5), rep(0.3, 5)),
                items = "5", mode = "C", seed = sd))
  # D131: the deterministic cases of test-permutation.R (binary_data()).
  add(lit_world("D131a-gaussian", "D", "gaussian", bin(10, c(6, 5, 7, 3)), items = "131"))
  add(lit_world("D131b-binomial", "D", "binomial", bin(10, c(6, 4, 4, 1)), items = "131"))
  add(lit_world("D131c-gaussian", "D", "gaussian", bin(5, c(3, 1, 2, 2)), items = "131,obs-zero"))
  # D8: Poisson, difference statistic, tied scores, 299 draws (item 8). Four
  # per group without a pretest covariate, so that every cell is saturated
  # and the classes stay few (at 10 per group they number about 9.5 million).
  w <- std_world("D8-poisson", "D", "poisson", 4, pre = FALSE, seed = nxt(), items = "8", extra_obs = FALSE)
  w$mode <- "C"
  w$reps <- 299
  add(w)

  # --- Block F: other degenerate categories
  w <- std_world("F-sep-binomial-n4", "F", "binomial", 4, pre = TRUE, seed = nxt(), items = "sep")
  add(w)
  add(std_world("F-rec-poisson", "F", "poisson", 3, pre = TRUE, seed = nxt(), items = "recession", eta0 = -0.9))
  add(std_world("F-lim-binomial", "F", "binomial", 3, pre = TRUE, cov = "num", seed = nxt(),
                items = "LIM", eta0 = -1.2))
  add(lit_world("F-infinf-binomial", "F", "binomial", list(c(1, 0, 0), c(0, 1, 0), c(1, 0, 0), c(0, 1, 0)),
                items = "INF-X", mode = "A"))
  add(lit_world("F-00-gaussian", "F", "gaussian", list(c(1, 1, 2), c(2, 2, 1), c(2, 1, 1), c(1, 2, 2)),
                items = "ZZ,INF-C", mode = "A"))
  add(lit_world("F-const-poisson", "F", "poisson", list(c(1, 3, 1), c(3, 1, 3), c(1, 1, 3), c(3, 3, 1)),
                items = "INF-C", mode = "A"))
  w <- lit_world("F-lev1-gaussian", "F", "gaussian", list(c(3, 5, 2), c(4, 1, 6), c(2, 7, 4), c(5, 1, 3)),
                 pre = c(1, 1, 5, 2, 2, 2, NA, NA, NA, NA, NA, NA), items = "H1", mode = "A")
  add(w)
  # F-miss: 16 rows, 2 missing posttests in each unpretested arm.
  for (cv in c("none", "num")) {
    w <- std_world(paste0("F-miss-gaussian-", cv), "F", "gaussian", 4, pre = TRUE, cov = cv, seed = nxt(),
                   items = "4,missing", extra_obs = FALSE)
    w$data$y[w$data$pretested == 0][c(1, 2, 5, 6)] <- NA
    add(w)
  }
  # F-unbal: 5 and 3 treated in the pretested condition.
  sd <- nxt()
  set.seed(sd)
  log_seed("F-unbal-gaussian", sd, "data", TRUE)
  w <- lit_world("F-unbal-gaussian", "F", "gaussian",
                 lapply(c(5, 3, 4, 4), function(k) round(stats::rnorm(k) * 4)), items = "unbalanced", mode = "A",
                 seed = sd)
  add(w)
  add(lit_world("F-obsinf-binomial", "F", "binomial", list(c(1, 1, 0), c(1, 0, 0), c(1, 1, 1), c(0, 1, 0)),
                items = "obs-inf", mode = "A"))
  add(lit_world("F-obszero-poisson", "F", "poisson", list(c(2, 0, 1), c(1, 3, 0), c(1, 2, 0), c(2, 0, 1)),
                items = "obs-zero", mode = "A"))
  w <- std_world("F-thbound-negbin", "F", "negbin", 3, pre = TRUE, seed = nxt(), items = "NB-BND",
                 nb_size = 50, extra_obs = FALSE, screen_theta = FALSE)
  add(w)

  # --- Block L: location (references by invariance from the base world)
  shift_world <- function(base, id, what, c, items) {
    w <- base
    w$id <- id
    w$block <- "L"
    w$items <- items
    w$base <- base$id
    w$shift <- list(what = what, c = c)
    if (what == "y") w$data$y <- w$data$y + c
    if (what == "pre") w$data$pre <- w$data$pre + c
    if (what %in% c("x", "year")) w$data[[what]] <- w$data[[what]] + c
    if (what == "xmat") w$data$xm <- I(matrix(w$data$xm + c, ncol = 1))
    w
  }
  gb <- std_world("L-base-gaussian", "L", "gaussian", 3, pre = TRUE, seed = nxt(), items = "6", screen_finite = TRUE)
  add(gb)
  add(shift_world(gb, "L-y1e9-gaussian", "y", 1e9, "6"))
  add(shift_world(gb, "L-y1e12-gaussian", "y", 1e12, "6"))
  gb8 <- gb
  gb8$id <- "L-base8-gaussian"
  gb8$data$y <- gb8$data$y / 8
  add(gb8)
  add(shift_world(gb8, "L-y1e12-gaussian-grid8", "y", 1e12, "6"))
  add(shift_world(gb, "L-pre2e40-gaussian", "pre", 2^40, "9 (pretest)"))
  for (fam in fams) {
    b <- std_world(paste0("L-base-x-", fam), "L", fam, 3, pre = TRUE, cov = "num", seed = nxt(), items = "9",
                   screen_finite = TRUE)
    add(b)
    sdx <- stats::sd(b$data$x)
    for (k in c(200, 1e4, 1e6)) {
      cc <- 2^ceiling(log2(k * sdx))
      add(shift_world(b, sprintf("L-x%g-%s", k, fam), "x", cc, "9"))
    }
    # calendar year 2019 to 2024, 3 per group
    yb <- std_world(paste0("L-base-year-", fam), "L", fam, 3, pre = TRUE, cov = "none", seed = nxt(), items = "9",
                    screen_finite = TRUE)
    sy0 <- nxt()
    ok <- FALSE
    for (sy in sy0 + 0:9) {
      set.seed(sy)
      # up to 20 draws of the year values in the seed's stream: each observed
      # labeling must keep an estimate with the year in the model (and, for
      # the negative binomial, an interior theta)
      for (k in 1:20) {
        yr <- sample(2019:2024, nrow(yb$data), replace = TRUE) - 2021
        dd <- yb$data
        dd$year <- yr
        ok <- all(vapply(yb$obs, function(t) mle_exists(dd, t, fam, TRUE, "year", NULL), logical(1))) &&
          (fam != "negbin" || all(vapply(yb$obs, function(t) nb_interior(dd, t, TRUE, "year", FALSE), logical(1))))
        if (ok) break
      }
      log_seed(yb$id, sy, if (ok) sprintf("year values (draw %d in the seed's stream)", k) else
        "year values: 20 draws failed the screens, replaced", ok)
      if (ok) break
    }
    if (!ok) stop("no year values passed the screens: ", yb$id)
    yb$data$year <- yr
    yb$covs <- "year"
    add(yb)
    add(shift_world(yb, paste0("L-year-", fam), "year", 2021, "9"))
  }
  # two-valued year with a binary outcome, 6 and 10 per group (no pretest)
  for (nn in c(6, 10)) {
    sdv <- nxt()
    set.seed(sdv)
    log_seed(paste0("L-base-year2-binomial-n", nn), sdv, "data", TRUE)
    d <- groups4(nn)
    yr <- sample(0:1, nrow(d), replace = TRUE)
    yy <- stats::rbinom(nrow(d), 1, stats::plogis(-0.2 + 0.6 * d$treat + 0.5 * yr))
    yb <- list(id = paste0("L-base-year2-binomial-n", nn), block = "L", items = "9", family = "binomial",
               data = data.frame(y = yy, treat = d$treat, pretested = d$pretested, year = yr), pre = FALSE,
               covs = "year", exposure = NULL, obs = list(gen = d$treat), seed = sdv, mode = "C", reps = 499)
    add(yb)
    add(shift_world(yb, paste0("L-year2-binomial-n", nn), "year", 2023, "9"))
  }
  # one-column numeric matrix covariate far from zero (Gaussian)
  mb <- std_world("L-base-xmat-gaussian", "L", "gaussian", 3, pre = TRUE, cov = "num", seed = nxt(), items = "9",
                  screen_finite = TRUE)
  mb$data$xm <- I(matrix(mb$data$x, ncol = 1))
  mb$data$x <- NULL
  mb$covs <- "xm"
  add(mb)
  add(shift_world(mb, "L-xmat2e20-gaussian", "xmat", 2^20, "9 (matrix)"))

  # --- Block G: clustered
  W <- c(W, cluster_worlds(nxt))

  # --- Added after the review of the first build (appended, so that the
  # replay seeds of the worlds above keep their positions)
  # F-sepobs: the data of the first build's M-S1-binomial-b (seed 1347121,
  # outcome draw 1, not screened), whose observed labelings are separated by
  # the pretest with a Treatment | pretested contrast that has no limit
  add(std_world("F-sepobs-binomial", "F", "binomial", 3, pre = TRUE, seed = 1347121, items = "SEP-I (observed)"))
  # F-exactfit: pretested counts on an exact log-linear curve (every
  # pretested residual 0; the studentized Treatment | pretested is +Inf)
  add(lit_world("F-exactfit-poisson", "F", "poisson", list(c(1, 4, 2), c(2, 8, 4), c(2, 0, 1), c(1, 3, 0)),
                pre = c(0, 2, 1, 0, 2, 1, NA, NA, NA, NA, NA, NA), items = "INF-C (covariate), exact fit", mode = "A"))
  # F-arm1: an unpretested treated arm of one row (groups 3, 3, 1, 3)
  add(lit_world("F-arm1-poisson", "F", "poisson", list(c(2, 0, 3), c(1, 4, 2), c(3), c(0, 1, 2)),
                items = "H1 (saturated arm of one)", mode = "A"))
  # F-miss for counts: 16 rows, 2 missing posttests in each unpretested arm
  w <- std_world("F-miss-poisson-none", "F", "poisson", 4, pre = TRUE, seed = nxt(), items = "missing (counts)",
                 extra_obs = FALSE)
  w$data$y[w$data$pretested == 0][c(1, 2, 5, 6)] <- NA
  add(w)
  # F-const-negbin: constant unpretested arms in some labelings (1, 1, 1
  # against 3, 3, 3), overdispersed pretested counts so that theta* is
  # interior; no covariate, so theta cancels from the statistics
  add(lit_world("F-const-negbin", "F", "negbin", list(c(0, 9, 1), c(12, 0, 2), c(1, 3, 1), c(3, 1, 3)),
                items = "INF-C (negative binomial)", mode = "A"))
  # D3-negbin: the issue's unpretested counts (item 3) under the negative
  # binomial; pretested counts drawn overdispersed
  sd3 <- nxt()
  okd3 <- FALSE
  for (s3 in sd3 + 0:9) {
    set.seed(s3)
    for (k in 1:20) {
      w <- lit_world("D3-negbin", "D", "negbin",
                     list(stats::rnbinom(10, size = 1, mu = 2.5), stats::rnbinom(10, size = 1, mu = 1.5),
                          c(1, 0, 0, 2, 0, 0, 1, 0, 0, 0), c(0, 0, 1, 0, 0, 0, 0, 0, 0, 0)),
                     items = "3 (negative binomial)", seed = s3)
      okd3 <- nb_interior(w$data, w$obs$gen, FALSE, character(0), FALSE)
      if (okd3) break
    }
    log_seed("D3-negbin", s3, if (okd3) sprintf("data (pretested counts only; draw %d in the seed's stream)", k) else
      "data: 20 draws had theta at the boundary, replaced", okd3)
    if (okd3) break
  }
  if (okd3) add(w)
  # G-conststrat: every unpretested score equal, so every cluster score of
  # that stratum is equal (the cluster analogue of item 5)
  w <- cluster_const_world(nxt())
  add(w)
  # Block E: constructed near-ties (item 7), seeds 1348001-1348499
  W <- c(W, near_tie_worlds())

  if (tier == "test") W <- W[vapply(W, function(w) w$block %in% c("D", "L", "G", "E") || grepl("-a$", w$id), logical(1))]
  W
}

# Clustered worlds. Design A: whole clusters assigned to the four conditions;
# Design B: treatment by cluster, pretesting within clusters.
cluster_worlds <- function(nxt) {
  W <- list()
  mk_A <- function(id, fam, k, m, pre, exposure = FALSE, cov = 0, items = "", seed) {
    set.seed(seed)
    log_seed(id, seed, "data", TRUE)
    cl <- rep(seq_len(4 * k), each = m)
    grp <- rep(rep(1:4, each = k), each = m)
    treat <- c(1, 0, 1, 0)[grp]
    pretested <- c(1, 1, 0, 0)[grp]
    N <- length(cl)
    prev <- ifelse(pretested == 1, round(stats::rnorm(N) * 2), NA)
    ce <- stats::rnorm(4 * k, 0, 0.5)[cl]
    ex <- sample(1:3, N, replace = TRUE)
    x <- round(stats::rnorm(N) * 2)
    eta <- switch(fam, gaussian = 0, binomial = -0.2, poisson = 0.5, negbin = 0.8) + 0.5 * treat + ce +
      (if (pre) 0.4 * ifelse(is.na(prev), 0, prev) else 0) + (if (exposure) log(ex) else 0) +
      (if (cov != 0) 0.3 * x else 0)
    y <- draw_y(fam, eta)
    dat <- data.frame(y = y, treat = treat, pretested = pretested, cl = cl)
    if (pre) dat$pre <- prev
    if (exposure) dat$ex <- ex
    covs <- character(0)
    if (cov != 0) {
      dat$x <- x
      covs <- "x"
    }
    list(id = id, block = "G", items = items, family = fam, data = dat, pre = pre, covs = covs,
         exposure = if (exposure) "ex" else NULL, cluster = "cl", design = "A", obs = list(gen = treat),
         seed = seed, mode = "cluster")
  }
  add <- function(w) W[[w$id]] <<- w
  add(mk_A("G-A-gaussian-k3", "gaussian", 3, 3, TRUE, seed = nxt()))
  add(mk_A("G-A-gaussian-k4", "gaussian", 4, 2, TRUE, seed = nxt()))
  add(mk_A("G-A-binomial-k3", "binomial", 3, 4, FALSE, seed = nxt()))
  add(mk_A("G-A-poisson-exp-k3", "poisson", 3, 3, FALSE, exposure = TRUE, seed = nxt()))
  add(mk_A("G-A-N-negbin-k3", "negbin", 3, 3, FALSE, seed = nxt(), items = "D7"))
  b10 <- mk_A("G-10-base-gaussian", "gaussian", 3, 3, TRUE, cov = 1, seed = nxt(), items = "10")
  add(b10)
  for (k in c(1e6, 1e8)) {
    w <- b10
    w$id <- sprintf("G-10-x%g-gaussian", k)
    cc <- 2^ceiling(log2(k * stats::sd(w$data$x)))
    w$data$x <- w$data$x + cc
    w$base <- b10$id
    w$shift <- list(what = "x", c = cc)
    add(w)
  }
  # G-zero: observed statistic exactly 0 (3 treated and 3 control clusters,
  # unpretested cluster means 1, 2, 3 in both arms).
  w <- mk_A("G-zero-gaussian", "gaussian", 3, 2, FALSE, seed = nxt(), items = "10 (Oct), obs-zero")
  u <- w$data$pretested == 0
  w$data$y[u] <- rep(c(1, 2, 3, 1, 2, 3), each = 2)
  add(w)
  # G-const: arms constant in the unpretested stratum.
  w <- mk_A("G-const-gaussian", "gaussian", 3, 2, FALSE, seed = nxt(), items = "INF-C")
  w$data$y[u] <- rep(c(1, 1, 1, 2, 2, 2), each = 2)
  add(w)
  # Design B: 6 and 10 clusters, every cluster with both conditions.
  for (k in c(6, 10)) {
    sd <- nxt()
    set.seed(sd)
    id <- sprintf("G-B-gaussian-k%d", k)
    log_seed(id, sd, "data", TRUE)
    cl <- rep(seq_len(k), each = 4)
    treat <- rep(rep(c(1, 0), each = k / 2), each = 4)
    pretested <- rep(c(1, 1, 0, 0), k)
    N <- length(cl)
    prev <- ifelse(pretested == 1, round(stats::rnorm(N) * 2), NA)
    y <- round(2 * (0.5 * treat + stats::rnorm(k, 0, 0.5)[cl] + ifelse(is.na(prev), 0, 0.4 * prev) + stats::rnorm(N)))
    add(list(id = id, block = "G", items = "", family = "gaussian",
             data = data.frame(y = y, treat = treat, pretested = pretested, cl = cl, pre = prev), pre = TRUE,
             covs = character(0), exposure = NULL, cluster = "cl", design = "B", obs = list(gen = treat),
             seed = sd, mode = "cluster"))
  }
  W
}

# G-conststrat: Design A, 3 clusters of 2 per condition, every unpretested
# score equal (2), so that the cluster scores of the unpretested stratum are
# all equal: Treatment | unpretested is 0 for every allocation (an error),
# and ATE and Pretest x Treatment vary with the pretested stratum only.
cluster_const_world <- function(seed) {
  set.seed(seed)
  log_seed("G-conststrat-gaussian", seed, "data", TRUE)
  k <- 3
  m <- 2
  cl <- rep(seq_len(4 * k), each = m)
  grp <- rep(rep(1:4, each = k), each = m)
  treat <- c(1, 0, 1, 0)[grp]
  pretested <- c(1, 1, 0, 0)[grp]
  N <- length(cl)
  ce <- stats::rnorm(4 * k, 0, 0.5)[cl]
  y <- round(4 * (0.5 * treat + ce + stats::rnorm(N)))
  y[pretested == 0] <- 2
  list(id = "G-conststrat-gaussian", block = "G", items = "5 (cluster level)", family = "gaussian",
       data = data.frame(y = y, treat = treat, pretested = pretested, cl = cl), pre = FALSE,
       covs = character(0), exposure = NULL, cluster = "cl", design = "A", obs = list(gen = treat),
       seed = seed, mode = "cluster")
}

# ---------------------------------------------------------------------------
# Block E: constructed near-ties (design Section 9, item 7).
#
# A world with a pretest and no other covariate, 4 per group, so that
# Treatment | pretested depends on the pretested rows only. The observed
# labeling is chosen among the 70 labelings of the pretested condition for a
# small studentized statistic (|T_obs| between 0.005 and 0.05, where the
# package's tolerance of 1e-6 is absolute), and two pretest values are moved
# so that one other labeling sits 5e-7 below |T_obs| and another 5e-7 above,
# inside that tolerance; the moved values are rounded to multiples of 2^-40.
# Kind "rel": |T_obs| between 1.2 and 3, one labeling at a relative 1e-5
# below it. Kind "gauss": the Gaussian analogue, a relative 5e-11 below and
# above an |T_obs| between 1.2 and 3 (the Gaussian tolerance is 1e-10). The
# construction uses doubles; the generator then decides every order exactly
# (Gaussian) or by the rules of design Section 5.

# Studentized Treatment | pretested of the pretested rows (doubles).
e_tp_stat <- function(dat, treat, fam, theta) {
  P <- which(dat$pretested == 1)
  X <- cbind(treat[P] == 0, treat[P] == 1, dat$pre[P]) * 1
  f <- rep(1, length(P))
  fit <- newton_fit(X, dat$y[P], f, rep(0, length(P)), fam, theta)
  if (!fit$converged) return(NA_real_)
  hc <- hc3_contrast(X, f, fit, c(-1, 1, 0))
  if (!is.finite(hc$V) || hc$V <= 0) return(NA_real_)
  hc$est / sqrt(hc$V)
}

e_theta <- function(dat, treat, fam) {
  if (fam != "negbin") return(NULL)
  cd <- cell_design(dat, treat, TRUE, character(0), NULL)
  th <- nb_theta(cd$X, cd$y, rep(1, length(cd$y)), cd$off)
  if (th$boundary || !is.finite(th$theta) || th$theta > 1e3) return(NA_real_)
  th$theta
}

# Move pretest values r of `dat` so that the gaps g_j = |T_aj| - |T_obs| equal
# targets[j] (one or two targets), theta fixed; Newton's method with a
# numerical Jacobian, from the nearest root of the first gap along row r[1].
e_move <- function(dat, r, to, ta, targets, fam, theta) {
  gaps <- function(p) {
    d2 <- dat
    d2$pre[r] <- p
    To <- e_tp_stat(d2, to, fam, theta)
    vapply(seq_along(ta), function(j) abs(e_tp_stat(d2, ta[[j]], fam, theta)) - abs(To) - targets[j], 0)
  }
  p <- dat$pre[r]
  # start: the nearest root of the first gap along the first row
  g1 <- function(dl) {
    q <- p
    q[1] <- q[1] + dl
    gaps(q)[1]
  }
  grid <- seq(-2, 2, by = 0.02)
  v <- vapply(grid, g1, 0)
  ok <- is.finite(v)
  idx <- which(ok[-length(v)] & ok[-1] & sign(v[-length(v)]) != sign(v[-1]))
  if (!length(idx)) return(NULL)
  j <- idx[which.min(abs(grid[idx]))]
  u <- tryCatch(stats::uniroot(g1, grid[c(j, j + 1)], tol = 1e-15, maxiter = 200), error = function(e) NULL)
  if (is.null(u)) return(NULL)
  p[1] <- p[1] + u$root
  for (it in 1:60) {
    g <- gaps(p)
    if (!all(is.finite(g))) return(NULL)
    if (max(abs(g)) < 1e-14) break
    h <- 1e-6
    J <- vapply(seq_along(p), function(k) {
      q <- p
      q[k] <- q[k] + h
      (gaps(q) - g) / h
    }, numeric(length(g)))
    J <- matrix(J, nrow = length(g))
    st <- tryCatch(if (length(p) == length(g)) solve(J, -g) else NULL, error = function(e) NULL)
    if (is.null(st) || !all(is.finite(st))) return(NULL)
    if (max(abs(st)) > 0.25) st <- st * 0.25 / max(abs(st))
    p <- p + st
  }
  if (max(abs(gaps(p))) > 1e-12) return(NULL)
  dat$pre[r] <- p
  dat
}

e_construct <- function(dat, fam, kind) {
  P <- which(dat$pretested == 1)
  comps <- utils::combn(length(P), length(P) / 2)
  lab <- function(j) {
    t <- dat$treat
    t[P] <- 0
    t[P[comps[, j]]] <- 1
    t
  }
  labs <- lapply(seq_len(ncol(comps)), lab)
  rng <- if (kind == "pair") c(0.005, 0.05) else c(1.2, 3)
  th0 <- e_theta(dat, dat$treat, fam)
  if (fam == "negbin" && !is.finite(th0)) return(NULL)
  T0 <- vapply(labs, function(t) e_tp_stat(dat, t, fam, th0), 0)
  cand_o <- which(is.finite(T0) & abs(T0) >= rng[1] & abs(T0) <= rng[2])
  for (io in cand_o) {
    to <- labs[[io]]
    th <- e_theta(dat, to, fam)
    if (fam == "negbin" && !is.finite(th)) next
    Tv <- vapply(labs, function(t) e_tp_stat(dat, t, fam, th), 0)
    To <- Tv[io]
    if (!is.finite(To) || abs(To) < rng[1] || abs(To) > rng[2]) next
    swap <- vapply(labs, function(t) all(t[P] == 1 - to[P]), logical(1))
    other <- is.finite(Tv) & seq_along(Tv) != io & !swap
    below <- which(other & abs(Tv) < abs(To))
    above <- which(other & abs(Tv) > abs(To))
    if (!length(below) || (kind != "rel" && !length(above))) next
    ia <- below[which.max(abs(Tv[below]))]
    ta <- list(labs[[ia]])
    targets <- if (kind == "pair") -5e-7 else if (kind == "rel") -1e-5 * abs(To) else -5e-11 * abs(To)
    if (kind != "rel") {
      ib <- above[which.min(abs(Tv[above]))]
      ta <- c(ta, list(labs[[ib]]))
      targets <- c(targets, if (kind == "pair") 5e-7 else 5e-11 * abs(To))
    }
    nr <- length(targets)
    rows <- if (nr == 1) as.list(P) else utils::combn(P, 2, simplify = FALSE)
    for (r in rows) {
      for (outer in 1:4) {
        d2 <- e_move(dat, r, to, ta, targets, fam, th)
        if (is.null(d2)) break
        d2$pre[r] <- round(d2$pre[r] * 2^40) / 2^40
        th2 <- e_theta(d2, to, fam)
        if (fam == "negbin" && !is.finite(th2)) {
          d2 <- NULL
          break
        }
        To2 <- e_tp_stat(d2, to, fam, th2)
        g <- vapply(ta, function(t) abs(e_tp_stat(d2, t, fam, th2)) - abs(To2), 0)
        tol <- if (kind == "gauss") 0.4 * abs(targets) else 0.2 * abs(targets)
        good <- all(is.finite(g)) && all(abs(g - targets) <= tol) && abs(To2) >= rng[1] && abs(To2) <= rng[2]
        if (good) {
          d2$treat <- to
          return(list(data = d2, gaps = g, T_obs = To2, theta = th2, rows = r))
        }
        # the negative binomial's theta moved with the pretest: solve again
        if (fam != "negbin") break
        dat <- d2
        th <- th2
      }
    }
  }
  NULL
}

near_tie_worlds <- function() {
  W <- list()
  specs <- list(list(id = "E-binomial-pair", fam = "binomial", kind = "pair", seed = 1348001),
                list(id = "E-poisson-pair", fam = "poisson", kind = "pair", seed = 1348011),
                list(id = "E-negbin-pair", fam = "negbin", kind = "pair", seed = 1348021),
                list(id = "E-poisson-rel", fam = "poisson", kind = "rel", seed = 1348031),
                list(id = "E-gaussian-pair", fam = "gaussian", kind = "gauss", seed = 1348041))
  for (sp in specs) {
    done <- FALSE
    for (s in sp$seed + 0:9) {
      set.seed(s)
      for (attempt in 1:10) {
        d <- groups4(4)
        N <- nrow(d)
        prev <- rep(NA_real_, N)
        prev[d$pretested == 1] <- sample(seq(-3, 3, by = 0.25), sum(d$pretested == 1))
        base <- switch(sp$fam, gaussian = 0, binomial = 0, poisson = 0.8, negbin = 1)
        eta <- base + 0.3 * d$treat + ifelse(is.na(prev), 0, 0.15 * prev)
        y <- draw_y(sp$fam, eta, nb_size = 1)
        dat <- data.frame(y = y, treat = d$treat, pretested = d$pretested, pre = prev)
        if (any(vapply(0:1, function(p) length(unique(dat$y[dat$pretested == p])) < 2, logical(1)))) next
        res <- e_construct(dat, sp$fam, sp$kind)
        if (is.null(res)) next
        ok <- sp$fam == "gaussian" || mle_exists(res$data, res$data$treat, sp$fam, TRUE, character(0), NULL)
        if (ok && sp$fam == "negbin") ok <- nb_interior(res$data, res$data$treat, TRUE, character(0), FALSE)
        if (!ok) next
        log_seed(sp$id, s, sprintf(paste0("near-tie construction (draw %d in the seed's stream): pretest rows %s moved; ",
                                          "gaps %s at |T_obs| = %.6g (doubles)"), attempt,
                                   paste(res$rows, collapse = ", "), paste(format(res$gaps, digits = 4), collapse = ", "),
                                   abs(res$T_obs)), TRUE)
        W[[sp$id]] <- list(id = sp$id, block = "E", items = "7", family = sp$fam, data = res$data, pre = TRUE,
                           covs = character(0), exposure = NULL, obs = list(gen = res$data$treat), seed = s,
                           mode = "A", construction = list(kind = sp$kind, rows = res$rows, gaps = res$gaps,
                                                           T_obs = res$T_obs))
        done <- TRUE
        break
      }
      if (done) break
      log_seed(sp$id, s, "near-tie construction failed in 10 draws, replaced", FALSE)
    }
    if (!done) stop("no near-tie construction succeeded: ", sp$id)
  }
  W
}
