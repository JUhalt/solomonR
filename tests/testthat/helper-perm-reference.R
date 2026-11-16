# The reference set for perm_solomon() and plot_perm() (issue #134):
# reading the fixture, and running the package on exactly the labelings the
# reference counts, by the rules of tools/perm-reference/compare.R (the
# reference-set comment on #134, Section 6) with the two harness changes and
# the conventions D1 to D16 posted on #134 before the remedy. No tolerance
# decides a count: the package's integer count b and its number of labelings
# used are compared with the reference's exactly.

perm_reference_cache <- new.env(parent = emptyenv())

# The fixture as posted on #134: SHA-256
# e236850dbee4b6be930fd0dc90b2c0beb74d76271bcc81f24c05fd85ddf59419. R before
# 4.5 has no tools::sha256sum(), so there its MD5 is checked instead.
perm_reference <- function() {
  if (is.null(perm_reference_cache$ref)) {
    path <- testthat::test_path("fixtures", "perm-reference.rds")
    ok <- if (exists("sha256sum", envir = asNamespace("tools"), inherits = FALSE)) {
      identical(unname(getExportedValue("tools", "sha256sum")(path)),
                "e236850dbee4b6be930fd0dc90b2c0beb74d76271bcc81f24c05fd85ddf59419")
    } else {
      identical(unname(tools::md5sum(path)), "9a3704e1fa6f0cbd6a5e83ad61beccbb")
    }
    if (!ok) stop("tests/testthat/fixtures/perm-reference.rds is not the reference set posted on #134.")
    perm_reference_cache$ref <- readRDS(path)
  }
  perm_reference_cache$ref
}

# "default" (CI) or "full" (SOLOMONR_PERM_REFERENCE=full: every case of the
# fixture's worlds with a class table, every cluster world).
perm_reference_tier <- function() {
  tier <- Sys.getenv("SOLOMONR_PERM_REFERENCE", "default")
  if (!tier %in% c("default", "full")) stop("SOLOMONR_PERM_REFERENCE must be unset, \"default\", or \"full\".")
  tier
}

# The Mode C replay seed of a world: 1348501 plus its position in the grid
# (compare.R, REPLAY_SEED0). The cases list every world of the full
# reference in grid order.
perm_replay_seed <- function(ref, world) 1348501 + match(world, unique(ref$cases$world)) - 1

perm_family_arg <- function(f) {
  switch(f, gaussian = stats::gaussian(), binomial = stats::binomial(), poisson = stats::poisson(),
         negbin = "negative_binomial")
}

# The fit of a world under a labeling, as compare.R's fit_world() makes it;
# the fit's own warnings are muffled, and an error is returned.
perm_fit_world <- function(w, treat, cluster = FALSE) {
  d <- w$data
  args <- list(y_post = d$y, treat = treat, pretested = d$pretested, family = perm_family_arg(w$family))
  if (isTRUE(w$pre)) args$y_pre <- d$pre
  if (length(w$covs)) args$covariates <- d[w$covs]
  if (!is.null(w$exposure)) args$exposure <- d[[w$exposure]]
  if (cluster) {
    args$cluster <- d$cl
    args$robust <- if (w$family == "negbin") "HC3" else "CR2"
  }
  tryCatch(suppressWarnings(do.call(fit_solomon_glm, args)), error = function(e) e)
}

# The class (composition) of a labeling: how many rows of each row type it
# treats in each condition (compare.R, class_of()).
perm_class_of <- function(r, treat) {
  key <- function(k) {
    ty <- r$types[[k]]
    tr <- treat[ty$idx] == 1
    paste(tabulate(ty$row_type[tr & !is.na(ty$row_type)], nbins = length(ty$m)), collapse = ",")
  }
  (match(key("U"), r$keys$U) - 1) * length(r$keys$P) + match(key("P"), r$keys$P)
}

# The reference's count on given classes under the adopted conventions
# (compare.R, ref_count(), its proposal columns): an emptied arm, separation,
# and constant arms counted as +Inf (D4, D10); leverage 1 an error with the
# studentized statistic (D3); a labeling under which the contrast is not
# estimable left out (D16). `counted` and `kept` are per labeling.
perm_ref_count <- function(st, cls, statistic) {
  code <- st$code[cls]
  rel <- st$rel[cls]
  counted <- !is.na(rel) & rel >= 0
  counted[code %in% c("SEP-L", "SEP-I")] <- TRUE
  kept <- code != "ALIAS"
  list(b = sum(counted[kept]), m = sum(kept),
       error = any(code == "NE") || (statistic == "studentized" && any(code == "H1")),
       undecided = any(code == "UNDECIDED") ||
         any(is.na(rel) & !(code %in% c("NE", "H1", "SEP-L", "SEP-I", "ALIAS"))),
       counted = counted & kept, kept = kept)
}

# Every labeling of the rows in the model (decision D2; compare.R with
# CMP_SPACE=model): in each pretest condition, every choice of treated rows
# among that condition's rows in the model, with as many treated as the
# observed labeling has there; one column per labeling, 0 outside the model.
perm_labelings <- function(d, obs, used) {
  combos <- lapply(c(U = 0, P = 1), function(p) {
    idx <- which(used & d$pretested == p)
    list(idx = idx, cb = utils::combn(length(idx), sum(obs[idx])))
  })
  grid <- expand.grid(u = seq_len(ncol(combos$U$cb)), p = seq_len(ncol(combos$P$cb)))
  vapply(seq_len(nrow(grid)), function(j) {
    t <- numeric(nrow(d))
    t[combos$U$idx[combos$U$cb[, grid$u[j]]]] <- 1
    t[combos$P$idx[combos$P$cb[, grid$p[j]]]] <- 1
    t
  }, numeric(nrow(d)))
}

# The value of an expression, or its error, with the classes of the
# warnings it gave (muffled).
perm_capture <- function(expr) {
  warns <- character(0)
  value <- withCallingHandlers(tryCatch(expr, error = function(e) e), warning = function(w) {
    warns <<- c(warns, class(w)[1])
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = unique(warns))
}

# plot_perm()'s p-value as printed, against a target: either neighbour at a
# half-way point of the third decimal (compare.R, plot_ok()).
perm_plot_ok <- function(p, target) {
  s <- plot_perm(p)$labels$subtitle
  got <- sub(".*permutation p = ([0-9.]+).*", "\\1", s)
  x <- target * 1000
  ok <- sprintf("%.3f", target)
  if (abs(x - floor(x) - 0.5) < 1e-6) ok <- sprintf("%.3f", c(floor(x), ceiling(x)) / 1000)
  got %in% ok
}

perm_without_p <- function(p) {
  p$p_perm <- NULL
  p
}

# Whether a case expects a p-value: the stored outcome, or D5's proposal
# (alt_outcome), which the remedy adopts for the ERR-S cases.
perm_wants_p <- function(case) {
  startsWith(case$outcome, "p") || (!is.na(case$alt_outcome) && startsWith(case$alt_outcome, "p"))
}

# The error class that meets a reference error raised by perm_solomon().
perm_error_class <- function(outcome) {
  switch(sub("^error:", "", outcome), "no-variation" = "solomonR_perm_no_variation_error",
         "H1" = , "observed-H1" = "solomonR_perm_leverage_error", "solomonR_perm_error")
}

# One participant case: Mode A (every labeling of the rows in the model, fed
# through the sampler seam .perm_treat()) or Mode C (the package's own draws
# under the recorded seed, replayed over the rows in the model and classed
# by the fixture's table).
perm_expect_case <- function(r, case, seed) {
  label <- case$case_id
  # no verdict: theta at the boundary (D14), or a case the set excludes
  if (grepl("^(no-reference|excluded)", case$outcome)) return(invisible(NULL))
  w <- r$world
  d <- w$data
  obs <- w$obs[[case$obs]]
  wants_p <- perm_wants_p(case)
  fit <- perm_fit_world(w, obs)
  if (inherits(fit, "error")) {
    # A reference error is met by a classed error of the fit. An arm of one
    # under HC3 is refused by fit_solomon_glm() (#173), a refusal kept for
    # 1.0 (decision M5 of #134), so F-arm1-poisson's cases with a p-value
    # stay unmet, by design.
    expect_true(any(startsWith(class(fit), "solomonR_")), label = paste(label, "(fit error classed)"))
    if (wants_p) {
      expect_true(case$world == "F-arm1-poisson" && inherits(fit, "solomonR_hc3_singleton_error"),
                  label = paste(label, "(fit error only where decision M5 keeps it):", conditionMessage(fit)))
    }
    return(invisible(NULL))
  }
  used <- rep(TRUE, nrow(d))
  if (!is.null(fit$model$na.action)) used[fit$model$na.action] <- FALSE
  # the rows in the model, by the package and by the reference
  expect_identical(as.numeric(sum(used)), sum(vapply(r$types, function(ty) sum(!is.na(ty$row_type)), 0)),
                   label = paste(label, "(rows in the model)"))
  st <- r$store[[case$obs]][[case$contrast]][[case$statistic]]
  M <- prod(vapply(0:1, function(p) choose(sum(used & d$pretested == p), sum(obs[used & d$pretested == p])), 0))
  if (identical(w$mode, "A") && M <= 4900) {
    labs <- perm_labelings(d, obs, used)
    cls <- apply(labs, 2, function(t) perm_class_of(r, t))
    k <- 0L
    res <- perm_capture(testthat::with_mocked_bindings(
      perm_solomon(fit, case$contrast, reps = ncol(labs), statistic = case$statistic, return_dist = TRUE),
      .perm_treat = function(treat, pretested) {
        k <<- k + 1L
        labs[used, k]
      }, .package = "solomonR"))
  } else {
    reps <- if (is.null(w$reps)) 999L else as.integer(w$reps)
    draws <- withr::with_seed(seed, lapply(seq_len(reps), function(i) {
      t <- numeric(nrow(d))
      t[used] <- .perm_treat(fit$data$treat[used], fit$data$pretested[used])
      t
    }))
    cls <- vapply(draws, function(t) perm_class_of(r, t), 0)
    res <- perm_capture(perm_solomon(fit, case$contrast, reps = reps, seed = seed, statistic = case$statistic,
                                     return_dist = TRUE))
  }
  v <- res$value
  rc <- perm_ref_count(st, cls, case$statistic)
  if (!wants_p) {
    expect_true(rc$error || startsWith(case$outcome, "error"), label = paste(label, "(reference error)"))
    expect_true(inherits(v, perm_error_class(case$outcome)),
                label = paste0(label, " (error of class ", perm_error_class(case$outcome), ")"))
    return(invisible(NULL))
  }
  if (inherits(v, "error")) {
    expect_false(TRUE, label = paste(label, "stopped:", conditionMessage(v)))
    return(invisible(NULL))
  }
  expect_false(rc$error || rc$undecided, label = paste(label, "(the reference decides every labeling)"))
  b <- round(v$p_perm * (v$valid_reps + 1)) - 1
  expect_identical(c(b = b, m = as.numeric(v$valid_reps)), c(b = as.numeric(rc$b), m = as.numeric(rc$m)),
                   label = label)
  # No draw is classified differently from the reference: each is left out,
  # or counted as at least as extreme, by both or by neither (compare.R,
  # misclassified(), with the window and unit floor the result records).
  kept <- !is.na(v$z_perm)
  counted <- kept
  counted[kept] <- .perm_at_least(v$z_perm[kept], v$z_obs,
                                  .perm_unit(v$z_perm[kept], case$statistic, v$tie_scale), v$tie_tolerance)
  expect_identical(kept, rc$kept, label = paste(label, "(draws left out)"))
  expect_identical(counted, rc$counted, label = paste(label, "(draws counted)"))
  # D5's warning in exactly the ERR-S cases
  expect_identical("solomonR_perm_no_variation_warning" %in% res$warnings, !startsWith(case$outcome, "p"),
                   label = paste(label, "(no-variation warning)"))
  # plot_perm() prints the reference's p, with and without p_perm
  target <- (rc$b + 1) / (rc$m + 1)
  expect_true(perm_plot_ok(v, target), label = paste(label, "(plot_perm() with p_perm)"))
  expect_true(perm_plot_ok(perm_without_p(v), target), label = paste(label, "(plot_perm() without p_perm)"))
  invisible(NULL)
}

# One cluster case: every allocation (reps above their number), so the
# package enumerates and p = b / M.
perm_expect_cluster_case <- function(r, case) {
  label <- case$case_id
  fit <- perm_fit_world(r$world, r$world$obs$gen, cluster = TRUE)
  res <- perm_capture(perm_solomon(fit, case$contrast, reps = 100000L, statistic = case$statistic,
                                   return_dist = TRUE))
  v <- res$value
  if (!perm_wants_p(case)) {
    expect_true(inherits(v, perm_error_class(case$outcome)),
                label = paste0(label, " (error of class ", perm_error_class(case$outcome), ")"))
    return(invisible(NULL))
  }
  if (inherits(v, "error")) {
    expect_false(TRUE, label = paste(label, "stopped:", conditionMessage(v)))
    return(invisible(NULL))
  }
  B <- as.numeric(case$b)
  M <- as.numeric(case$M)
  expect_true(isTRUE(v$exact), label = paste(label, "(enumerated)"))
  expect_identical(c(b = round(v$p_perm * v$valid_reps), m = as.numeric(v$valid_reps)), c(b = B, m = M),
                   label = label)
  expect_identical("solomonR_perm_no_variation_warning" %in% res$warnings, !startsWith(case$outcome, "p"),
                   label = paste(label, "(no-variation warning)"))
  expect_true(perm_plot_ok(v, B / M), label = paste(label, "(plot_perm() with p_perm)"))
  expect_true(perm_plot_ok(perm_without_p(v), B / M), label = paste(label, "(plot_perm() without p_perm)"))
  invisible(NULL)
}
