# Compare the package with the reference set (issue #134; design Section 8).
#
#   Rscript tools/perm-reference/compare.R <reference-full.rds> <out.csv> [workers]
#
# Loads the package from this worktree with pkgload::load_all(). The package
# is run on exactly the labelings the reference counts:
#   Mode A  every labeling, fed through the package's own refit by mocking
#           its sampler (.perm_treat) with testthat::with_mocked_bindings();
#   Mode C  the package's own draws under a recorded seed, replayed with its
#           own sampler and classified by the reference's class table;
#   cluster every allocation (reps >= the number of allocations), so the
#           package enumerates.
# The integer count b is compared, not p (design Section 8). The plot_perm()
# subtitle is checked against the reference p for every object. Also per
# case (design Section 8): the rows in the model against the reference's,
# the package's negative binomial theta against theta* and the band over
# which the reference certified its decisions, and the classed warnings of
# the fit; and every draw that the package classifies differently from the
# reference is listed (misclassified-draws.csv beside the output), in every
# case with a class table whether or not its counts agree, including the
# draws the package drops (a statistic that is not finite) and those it uses
# where the reference leaves the labeling out.

args <- commandArgs(TRUE)
ref_file <- args[1]
out_csv <- args[2]
workers <- if (length(args) >= 3) as.integer(args[3]) else 1L
ROOT <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))),
                                "..", ".."), winslash = "/")
REPLAY_SEED0 <- 1348501  # Mode C seeds 1348501-1348799, one per world, in grid order

ref <- readRDS(ref_file)

fam_arg <- function(f) switch(f, gaussian = stats::gaussian(), binomial = stats::binomial(),
                              poisson = stats::poisson(), negbin = "negative_binomial")

fit_world <- function(wd, treat, cluster = FALSE) {
  d <- wd$data
  args <- list(y_post = d$y, treat = treat, pretested = d$pretested, family = fam_arg(wd$family))
  if (isTRUE(wd$pre)) args$y_pre <- d$pre
  if (length(wd$covs)) args$covariates <- d[wd$covs]
  if (!is.null(wd$exposure)) args$exposure <- d[[wd$exposure]]
  if (cluster) {
    args$cluster <- d$cl
    args$robust <- if (wd$family == "negbin") "HC3" else "CR2"
  }
  # the classed warnings of the fit are kept (theta at the boundary, ...)
  warns <- character(0)
  fit <- withCallingHandlers(do.call(solomonR::fit_solomon_glm, args), warning = function(w) {
    cl <- setdiff(class(w), c("warning", "condition", "simpleWarning"))
    warns <<- c(warns, if (length(cl)) paste(cl, collapse = "/") else "warning")
    invokeRestart("muffleWarning")
  })
  attr(fit, "fit_warnings") <- paste(unique(warns), collapse = ";")
  fit
}

# Rows in the model, by the reference (row types not "OUT") and by the package.
ref_rows <- function(r) sum(vapply(r$types, function(ty) sum(!is.na(ty$row_type)), 0))
pkg_rows <- function(fit) tryCatch(as.numeric(stats::nobs(fit$model)), error = function(e) NA_real_)

# The package's theta (that of the model it refits) against theta* and the
# certified band (design Section 8: reported as a theta mismatch, not a p
# mismatch).
theta_cols <- function(fit, cr) {
  if (is.null(cr$theta) || !is.finite(cr$theta)) {
    return(data.frame(theta_ref = NA, theta_pkg = NA, theta_rel_diff = NA, theta_band = NA, theta_ok = NA))
  }
  tp <- tryCatch(fit$model$theta, error = function(e) NULL)
  if (is.null(tp)) tp <- NA_real_
  rd <- tp / cr$theta - 1
  data.frame(theta_ref = cr$theta, theta_pkg = tp, theta_rel_diff = rd, theta_band = cr$theta_band,
             theta_ok = is.finite(rd) && abs(rd) <= cr$theta_band)
}

class_of <- function(r, treat) {
  key <- function(k) {
    ty <- r$types[[k]]
    tr <- treat[ty$idx] == 1
    v <- tabulate(ty$row_type[tr & !is.na(ty$row_type)], nbins = length(ty$m))
    paste(v, collapse = ",")
  }
  iU <- match(key("U"), r$keys$U)
  iP <- match(key("P"), r$keys$P)
  (iU - 1) * length(r$keys$P) + iP
}

# Reference count on given labelings (classes), proposal conventions.
ref_count <- function(st, cls, s) {
  code <- st$code[cls]
  rel <- st$rel[cls]
  under <- if (is.null(st$under)) code else st$under[cls]
  cnt <- as.numeric(!is.na(rel) & rel >= 0)
  cnt[code %in% c("SEP-L", "SEP-I")] <- 1
  keep <- code != "ALIAS"  # left out (see reference.R, case_counts())
  err <- if (any(code == "NE")) "NE" else if (any(code == "H1" & s == "studentized")) "H1" else NULL
  und <- any(code == "UNDECIDED") || any(is.na(rel) & !(code %in% c("NE", "H1", "SEP-L", "SEP-I", "ALIAS")))
  # the alternative that leaves out labelings without a defined statistic
  kd <- keep & !(code %in% c("NE", "H1", "ZZ", "SEP-I")) & !(code == "SEP-L" & s == "studentized")
  # the alternative that gives a leverage-1 labeling its infinite value
  hinf <- code == "H1" & under %in% c("INF-E", "INF-X", "INF-C", "SEP-L", "SEP-I")
  kh <- keep & !(code == "NE") & !(code == "H1" & !hinf)
  cnt_h <- cnt
  cnt_h[hinf] <- 1
  list(b = sum(cnt[keep]), m = sum(keep), error = err, undecided = und,
       b_drop = sum(cnt[kd]), m_drop = sum(kd), b_h1inf = sum(cnt_h[kh]), m_h1inf = sum(kh),
       err_h1inf = any(code == "NE") || any(code == "H1" & !hinf), cnt = cnt, keep = keep)
}

# The package returns only its finite permutation statistics. To list the
# draws it drops, the full vector is read from perm_solomon()'s frame when it
# calls .perm_at_least(), which is mocked by a wrapper that returns what the
# package's own function returns (no change to R/).
ZCAP <- new.env()
capture_pal <- function(z, z0, unit = 1, tol = 1e-10) {
  pf <- parent.frame()
  if (exists("z_perm", envir = pf, inherits = FALSE)) assign("z_full", get("z_perm", envir = pf), envir = ZCAP)
  ZCAP$orig(z, z0, unit, tol)
}
run_captured <- function(expr_fun) {
  ZCAP$z_full <- NULL
  ZCAP$orig <- get(".perm_at_least", envir = asNamespace("solomonR"))
  expr_fun()
}

# Draws that the package classifies differently from the reference (design
# Section 8): index, class, reference code and |T|, package z, whether each
# keeps the draw (the package drops a statistic that is not finite; the
# reference leaves out ALIAS labelings), and whether each counts it as at
# least as extreme. Listed for every case, whether or not the counts agree.
MISCLASS <- list()
misclassified <- function(case_id, cls, st, z_full, z_obs, family_obj, statistic) {
  if (is.null(z_full) || length(z_full) != length(cls)) return(NULL)
  rc <- ref_count(st, cls, statistic)
  pkg_kept <- is.finite(z_full)
  pk <- rep(NA, length(z_full))
  if (any(pkg_kept)) {
    zf <- z_full[pkg_kept]
    pk[pkg_kept] <- ZCAP$orig(zf, z_obs, solomonR:::.perm_unit(zf, statistic),
                              solomonR:::.perm_tie_tolerance(family_obj))
  }
  ref_kept <- rc$keep
  rf <- rc$cnt == 1 & ref_kept
  bad <- which(pkg_kept != ref_kept | (pkg_kept & ref_kept & pk != rf))
  if (!length(bad)) return(NULL)
  data.frame(case_id = case_id, draw = bad, class = cls[bad], ref_code = st$code[cls[bad]],
             ref_T = st$approx[cls[bad]], pkg_z = z_full[bad], pkg_z_obs = z_obs, ref_kept = ref_kept[bad],
             pkg_kept = pkg_kept[bad], ref_counted = rf[bad], pkg_counted = pk[bad], stringsAsFactors = FALSE)
}

safe_bind <- function(rows) {
  cols <- unique(unlist(lapply(rows, names)))
  do.call(rbind, lapply(rows, function(z) {
    for (nm in setdiff(cols, names(z))) z[[nm]] <- NA
    z[cols]
  }))
}

run_pkg <- function(expr) {
  warns <- character(0)
  val <- withCallingHandlers(
    tryCatch(expr, error = function(e) structure(list(msg = conditionMessage(e), class = class(e)), class = "pkg_error")),
    warning = function(w) {
      warns <<- c(warns, paste(setdiff(class(w), c("warning", "condition", "simpleWarning")), collapse = "/"))
      invokeRestart("muffleWarning")
    })
  list(val = val, warns = unique(warns))
}

subtitle_p <- function(p) {
  s <- tryCatch(solomonR::plot_perm(p)$labels$subtitle, error = function(e) NA_character_)
  if (is.na(s)) return(NA_character_)
  sub(".*permutation p = ([0-9.]+).*", "\\1", s)
}
plot_ok <- function(got, target) {
  if (is.na(got)) return(NA)
  tgt <- sprintf("%.3f", target)
  if (identical(got, tgt)) return(TRUE)
  # either neighbour at a half-way point of the third decimal (design 8)
  x <- target * 1000
  if (abs(x - floor(x) - 0.5) < 1e-6) return(got %in% c(sprintf("%.3f", floor(x) / 1000), sprintf("%.3f", ceiling(x) / 1000)))
  FALSE
}

compare_world <- function(r, replay_seed) {
  wd <- r$world
  rows <- list()
  cs <- r$cases
  CSHORT <- c("ATE (avg over pretest)" = "ATE", "Pretest x Treatment" = "PxT",
              "Treatment | pretested" = "T|P", "Treatment | unpretested" = "T|U")
  cs$case_id <- sprintf("%s|%s|%s|%s|%s", cs$world, cs$obs, CSHORT[cs$contrast], cs$statistic, cs$space)
  mode <- wd$mode
  if (identical(mode, "cluster")) {
    fit <- fit_world(wd, wd$obs$gen, cluster = TRUE)
    for (i in seq_len(nrow(cs))) {
      cr <- cs[i, ]
      res <- run_pkg(solomonR::perm_solomon(fit, cr$contrast, reps = 100000L, statistic = cr$statistic,
                                            return_dist = TRUE))
      v <- res$val
      out <- data.frame(case_id = cr$case_id, mode = "cluster (enumerated)", seed = NA, reps = NA,
                        ref_outcome = cr$outcome, ref_b = cr$b, ref_M = cr$M, ref_p = cr$p,
                        ref_b_same = cr$b, ref_m_same = cr$M, stringsAsFactors = FALSE)
      if (inherits(v, "pkg_error")) {
        out <- cbind(out, pkg_outcome = paste("error:", v$msg), pkg_b = NA, pkg_valid = NA, pkg_p = NA,
                     pkg_p_same = NA, agree = grepl("^error", cr$outcome) && any(grepl("^solomonR_", v$class)),
                     agree_alt = if (is.na(cr$alt_outcome)) NA else FALSE,
                     z_maxrel = NA, plot_with_p = NA, plot_without_p = NA)
      } else {
        bp <- round(v$p_perm * v$valid_reps)
        q <- v
        q$p_perm <- NULL
        same <- isTRUE(v$exact) && v$valid_reps == cr$M && isTRUE(bp == cr$b)
        out <- cbind(out, pkg_outcome = "p", pkg_b = bp, pkg_valid = v$valid_reps, pkg_p = v$p_perm,
                     pkg_p_same = bp / v$valid_reps,
                     agree = grepl("^p", cr$outcome) && same,
                     agree_alt = if (is.na(cr$alt_outcome)) NA else same,
                     z_maxrel = NA, plot_with_p = plot_ok(subtitle_p(v), cr$p),
                     plot_without_p = plot_ok(subtitle_p(q), cr$p))
      }
      out$warnings <- paste(res$warns, collapse = ";")
      rows[[length(rows) + 1]] <- out
    }
    return(safe_bind(rows))
  }
  cs <- cs[cs$space == "all", ]
  d <- wd$data
  # every labeling of the space over all rows (the reference's M can be
  # smaller, when a convention leaves labelings out)
  M_all <- prod(vapply(0:1, function(p) choose(sum(d$pretested == p), sum(wd$obs$gen[d$pretested == p])), 0))
  use_A <- identical(mode, "A") && is.finite(M_all) && M_all <= 4900
  if (use_A) {
    combos <- lapply(c(U = 0, P = 1), function(p) {
      idx <- which(d$pretested == p)
      cb <- utils::combn(length(idx), sum(wd$obs$gen[idx]))
      list(idx = idx, cb = cb)
    })
    grid <- expand.grid(u = seq_len(ncol(combos$U$cb)), p = seq_len(ncol(combos$P$cb)))
    labs <- vapply(seq_len(nrow(grid)), function(j) {
      t <- numeric(nrow(d))
      t[combos$U$idx[combos$U$cb[, grid$u[j]]]] <- 1
      t[combos$P$idx[combos$P$cb[, grid$p[j]]]] <- 1
      t
    }, numeric(nrow(d)))
    if (ncol(labs) != M_all) stop("labelings do not match M for ", wd$id)
  }
  for (o in names(wd$obs)) {
    fit <- tryCatch(fit_world(wd, wd$obs[[o]]), error = function(e) e)
    fit_warn <- if (inherits(fit, "error")) "" else attr(fit, "fit_warnings")
    rows_ok <- if (inherits(fit, "error")) NA else isTRUE(pkg_rows(fit) == ref_rows(r))
    if (use_A) {
      cls <- apply(labs, 2, function(t) class_of(r, t))
    } else {
      reps <- if (is.null(wd$reps)) 999L else as.integer(wd$reps)
      draws <- if (inherits(fit, "error")) NULL else
        withr::with_seed(replay_seed, lapply(seq_len(reps), function(i) solomonR:::.perm_treat(fit$data$treat, fit$data$pretested)))
      cls <- if (is.null(draws)) NULL else vapply(draws, function(t) class_of(r, t), 0)
    }
    for (i in which(cs$obs == o)) {
      cr <- cs[i, ]
      st <- r$store[[o]][[cr$contrast]][[cr$statistic]]
      base <- data.frame(case_id = cr$case_id, mode = if (use_A) "A (all labelings)" else "C (replay)",
                         seed = if (use_A) NA else replay_seed, reps = if (use_A) M_all else length(cls),
                         ref_outcome = cr$outcome, ref_b = cr$b, ref_M = cr$M, ref_p = cr$p,
                         stringsAsFactors = FALSE)
      rc <- if (!is.null(st) && !is.null(cls)) ref_count(st, cls, cr$statistic) else list(b = NA, m = NA)
      base$ref_b_same <- rc$b
      base$ref_m_same <- rc$m
      base$fit_warnings <- fit_warn
      base$rows_ok <- rows_ok
      if (inherits(fit, "error")) {
        rows[[length(rows) + 1]] <- cbind(base, pkg_outcome = paste("fit error:", conditionMessage(fit)), pkg_b = NA,
                                          pkg_valid = NA, pkg_p = NA, pkg_p_same = NA,
                                          # a reference error is met by a classed solomonR error,
                                          # whether perm_solomon() or the fit raises it
                                          agree = if (grepl("^(no-reference|excluded)", cr$outcome)) NA else
                                            grepl("^error", cr$outcome) && any(grepl("^solomonR_", class(fit))),
                                          z_maxrel = NA, plot_with_p = NA, plot_without_p = NA, warnings = "")
        next
      }
      base <- cbind(base, theta_cols(fit, cr))
      if (use_A) {
        env <- new.env()
        env$k <- 0L
        feeder <- function(treat, pretested) {
          env$k <- env$k + 1L
          labs[, env$k]
        }
        res <- run_captured(function() run_pkg(testthat::with_mocked_bindings(
          solomonR::perm_solomon(fit, cr$contrast, reps = M_all, statistic = cr$statistic, return_dist = TRUE),
          .perm_treat = feeder, .perm_at_least = capture_pal, .package = "solomonR")))
      } else {
        res <- run_captured(function() run_pkg(testthat::with_mocked_bindings(
          solomonR::perm_solomon(fit, cr$contrast, reps = length(cls), seed = replay_seed,
                                 statistic = cr$statistic, return_dist = TRUE),
          .perm_at_least = capture_pal, .package = "solomonR")))
      }
      v <- res$val
      if (inherits(v, "pkg_error")) {
        out <- cbind(base, pkg_outcome = paste("error:", v$msg), pkg_b = NA, pkg_valid = NA, pkg_p = NA,
                     pkg_p_same = NA, agree = grepl("^error", cr$outcome) && any(grepl("^solomonR_", v$class)),
                     agree_alt = if (is.na(cr$alt_outcome)) NA else FALSE,
                     agree_h1inf = if (identical(cr$outcome_h1inf, "p") && !(cr$n_h1_inf %in% c("0", 0))) FALSE else NA,
                     z_maxrel = NA, plot_with_p = NA, plot_without_p = NA)
      } else {
        bp <- round(v$p_perm * (v$valid_reps + 1)) - 1
        zr <- NA
        if (v$valid_reps == length(cls) && !is.null(st)) {
          a <- st$approx[cls]
          fin <- is.finite(a) & st$code[cls] %in% c("finite", "LIM", "ZERO")
          if (any(fin)) zr <- max(abs(v$z_perm[fin] - a[fin]) / pmax(1, abs(a[fin])))
        }
        q <- v
        q$p_perm <- NULL
        target <- (rc$b + 1) / (rc$m + 1)
        same <- is.null(rc$error) && !isTRUE(rc$undecided) && v$valid_reps == rc$m && bp == rc$b
        ok <- grepl("^p", cr$outcome) && same
        if (!is.null(st) && !is.null(cls)) {
          mc <- misclassified(cr$case_id, cls, st, ZCAP$z_full, v$z_obs, stats::family(fit$model), cr$statistic)
          if (!is.null(mc)) MISCLASS[[length(MISCLASS) + 1]] <<- mc
        }
        out <- cbind(base, pkg_outcome = "p", pkg_b = bp, pkg_valid = v$valid_reps, pkg_p = v$p_perm,
                     pkg_p_same = bp / v$valid_reps, agree = ok,
                     agree_alt = if (is.na(cr$alt_outcome)) NA else same,
                     agree_h1inf = if (identical(cr$outcome_h1inf, "p") && cr$n_h1_inf %in% c("0", 0)) NA else
                       if (identical(cr$outcome_h1inf, "p")) !isTRUE(rc$undecided) && !isTRUE(rc$err_h1inf) &&
                         v$valid_reps == rc$m_h1inf && bp == rc$b_h1inf else NA,
                     agree_drop = !is.null(rc$b_drop) && !isTRUE(rc$undecided) && v$valid_reps == rc$m_drop &&
                       bp == rc$b_drop,
                     z_maxrel = zr,
                     plot_with_p = if (grepl("^p", cr$outcome)) plot_ok(subtitle_p(v), target) else NA,
                     plot_without_p = if (grepl("^p", cr$outcome)) plot_ok(subtitle_p(q), target) else NA)
      }
      out$warnings <- paste(res$warns, collapse = ";")
      # no reference p (theta at the boundary, or a case excluded): no verdict
      if (grepl("^(no-reference|excluded)", cr$outcome)) out$agree <- NA
      rows[[length(rows) + 1]] <- out
    }
  }
  out <- safe_bind(rows)
  attr(out, "misclassified") <- if (length(MISCLASS)) do.call(rbind, MISCLASS)
  out
}

# Scope block S: refused spellings must be refused by perm_solomon() with a
# classed error; accepted spellings must give a p-value (design Section 9).
scope_block <- function() {
  set.seed(1349501)  # reserve block; recorded in build.md
  g <- rep(1:4, each = 6)
  treat <- c(1, 0, 1, 0)[g]
  pretested <- c(1, 1, 0, 0)[g]
  yb <- stats::rbinom(24, 1, 0.5)
  yc <- stats::rpois(24, 3) + 1
  yg <- round(stats::rnorm(24, 10, 2)) + 0.5
  specs <- list(
    list("gaussian(link = \"log\")", yg, stats::gaussian(link = "log"), "refuse"),
    list("Gamma(link = \"log\")", yg, stats::Gamma(link = "log"), "refuse"),
    list("inverse.gaussian(link = \"log\")", yg, stats::inverse.gaussian(link = "log"), "refuse"),
    list("binomial(link = \"probit\")", yb, stats::binomial(link = "probit"), "refuse"),
    list("binomial(link = \"cloglog\")", yb, stats::binomial(link = "cloglog"), "refuse"),
    list("binomial(link = \"cauchit\")", yb, stats::binomial(link = "cauchit"), "refuse"),
    list("quasibinomial()", yb, stats::quasibinomial(), "refuse"),
    list("poisson(link = \"identity\")", yc, stats::poisson(link = "identity"), "refuse"),
    list("poisson(link = \"sqrt\")", yc, stats::poisson(link = "sqrt"), "refuse"),
    list("quasipoisson()", yc, stats::quasipoisson(), "refuse"),
    list("MASS::negative.binomial(2)", yc, MASS::negative.binomial(2), "refuse (pending D1)"),
    list("MASS::negative.binomial(2, link = \"sqrt\")", yc, MASS::negative.binomial(2, link = "sqrt"), "refuse"),
    list("\"poisson\"", yc, "poisson", "accept"),
    list("poisson", yc, stats::poisson, "accept"),
    list("poisson()", yc, stats::poisson(), "accept"),
    list("\"binomial\"", yb, "binomial", "accept"),
    list("binomial(link = \"logit\")", yb, stats::binomial(link = "logit"), "accept"),
    list("gaussian", yg, stats::gaussian, "accept"),
    list("\"negative_binomial\"", yc, "negative_binomial", "accept"))
  rows <- list()
  for (sp in specs) for (cl in c(FALSE, TRUE)) {
    fit <- tryCatch(suppressWarnings(solomonR::fit_solomon_glm(
      sp[[2]], treat, pretested, family = sp[[3]],
      cluster = if (cl) rep(1:12, each = 2) else NULL,
      robust = if (cl && !identical(sp[[3]], "negative_binomial")) "CR2" else "HC3")), error = function(e) e)
    if (inherits(fit, "error")) {
      res <- list(val = structure(list(msg = paste("fit_solomon_glm:", conditionMessage(fit)), class = "fit"),
                                  class = "pkg_error"), warns = character(0))
    } else {
      res <- run_pkg(solomonR::perm_solomon(fit, "Treatment | unpretested", reps = 99L, seed = 1349502L))
    }
    v <- res$val
    got <- if (inherits(v, "pkg_error")) paste0("error: ", v$msg) else sprintf("p = %.3f", v$p_perm)
    classed <- inherits(v, "pkg_error") && any(grepl("^solomonR_", v$class))
    agree <- if (startsWith(sp[[4]], "refuse")) classed else !inherits(v, "pkg_error")
    rows[[length(rows) + 1]] <- data.frame(
      case_id = sprintf("S|%s|%s", sp[[1]], if (cl) "cluster" else "participant"), mode = "S (scope)",
      seed = 1349502, reps = 99, ref_outcome = sp[[4]], ref_b = NA, ref_M = NA, ref_p = NA, ref_b_same = NA,
      ref_m_same = NA, pkg_outcome = got, pkg_b = NA, pkg_valid = NA,
      pkg_p = if (inherits(v, "pkg_error")) NA else v$p_perm, pkg_p_same = NA, agree = agree, z_maxrel = NA,
      plot_with_p = NA, plot_without_p = NA, warnings = paste(res$warns, collapse = ";"), stringsAsFactors = FALSE)
  }
  safe_bind(rows)
}

t0 <- proc.time()[["elapsed"]]
refs <- ref$refs
seeds <- stats::setNames(REPLAY_SEED0 + seq_along(refs) - 1L, names(refs))
# CMP_ONLY (a regular expression) runs some worlds only, without the scope
# block; their seeds stay those of the full grid
only <- Sys.getenv("CMP_ONLY")
if (nzchar(only)) {
  refs <- refs[grepl(only, names(refs))]
  seeds <- seeds[names(refs)]
}
work <- function(nm) {
  t1 <- proc.time()[["elapsed"]]
  MISCLASS <<- list()
  out <- tryCatch(compare_world(refs[[nm]], seeds[[nm]]), error = function(e) {
    data.frame(case_id = paste0(nm, "|harness error"), mode = NA, seed = NA, reps = NA, ref_outcome = NA, ref_b = NA,
               ref_M = NA, ref_p = NA, ref_b_same = NA, ref_m_same = NA, pkg_outcome = conditionMessage(e), pkg_b = NA,
               pkg_valid = NA, pkg_p = NA, pkg_p_same = NA, agree = NA, z_maxrel = NA, plot_with_p = NA,
               plot_without_p = NA, warnings = "", stringsAsFactors = FALSE)
  })
  attr(out, "secs") <- proc.time()[["elapsed"]] - t1
  out
}
if (workers > 1) {
  cl <- parallel::makeCluster(workers)
  parallel::clusterExport(cl, c("ROOT", "refs", "seeds", "fam_arg", "fit_world", "class_of", "ref_count", "run_pkg", "safe_bind",
                                "subtitle_p", "plot_ok", "compare_world", "work", "ref_rows", "pkg_rows", "theta_cols",
                                "misclassified", "MISCLASS", "ZCAP", "capture_pal", "run_captured"))
  invisible(parallel::clusterEvalQ(cl, suppressMessages(pkgload::load_all(ROOT, quiet = TRUE))))
  ord <- names(refs)[order(-vapply(refs, function(r) suppressWarnings(max(as.numeric(r$cases$M), na.rm = TRUE)), 0))]
  res <- parallel::parLapplyLB(cl, ord, work)
  parallel::stopCluster(cl)
  names(res) <- ord
  res <- res[names(refs)]
} else {
  suppressMessages(pkgload::load_all(ROOT, quiet = TRUE))
  res <- lapply(names(refs), work)
  names(res) <- names(refs)
}
suppressMessages(pkgload::load_all(ROOT, quiet = TRUE))
secs <- vapply(res, function(x) attr(x, "secs"), 0)
mis <- do.call(rbind, lapply(res, function(x) attr(x, "misclassified")))
cmp <- safe_bind(c(unname(res), if (!nzchar(only)) list(scope_block())))
key <- ref$cases[, c("case_id", "world", "block", "items", "family", "obs", "contrast", "statistic", "space")]
cmp <- merge(key, cmp, by = "case_id", all.y = TRUE, sort = FALSE)
utils::write.csv(cmp, out_csv, row.names = FALSE)
if (!is.null(mis)) utils::write.csv(mis, file.path(dirname(out_csv), "misclassified-draws.csv"), row.names = FALSE)
saveRDS(list(cmp = cmp, secs = secs, total = proc.time()[["elapsed"]] - t0, seeds = seeds, misclassified = mis),
        sub("\\.csv$", ".rds", out_csv))
cat(sprintf("cases compared %d; agree %d; disagree %d; NA %d; total %.0f s\n", nrow(cmp), sum(cmp$agree %in% TRUE),
            sum(cmp$agree %in% FALSE), sum(is.na(cmp$agree)), proc.time()[["elapsed"]] - t0))
