# Reference set for the participant- and cluster-level permutation tests
# (issue #134): enumeration, exact keys, fits at the exact estimate, and the
# bridge to exact.py. This file must not load the package under test or name
# its internal functions (design, Section 11, rule 1); it fits nothing with
# glm() for a reference statistic and decides no equality with a tolerance.
#
# The statistic (design, Section 3.2) is computed in cell-means form: one
# indicator for each of the four cells (UC, UT, PC, PT), the pretest (0 for
# the unpretested), the covariates, and the offset. HC3 (MacKinnon & White,
# 1985) as sandwich computes it for a glm (Zeileis, 2006): bread (X'WX)^-1
# with Fisher weights, scores x_i (y_i - mu_i) (dmu/deta) / Var(mu), hat
# values of W^(1/2) X, each score divided by 1 - h_i.

CONTRASTS <- c("ATE (avg over pretest)", "Pretest x Treatment",
               "Treatment | pretested", "Treatment | unpretested")
CSHORT <- c("ATE (avg over pretest)" = "ATE", "Pretest x Treatment" = "PxT",
            "Treatment | pretested" = "T|P", "Treatment | unpretested" = "T|U")
# Condition weights (unpretested, pretested): T = w_U e_U + w_P e_P.
CW <- list("ATE (avg over pretest)" = c(U = 0.5, P = 0.5),
           "Pretest x Treatment" = c(U = -1, P = 1),
           "Treatment | pretested" = c(U = 0, P = 1),
           "Treatment | unpretested" = c(U = 1, P = 0))
STATS <- c("studentized", "difference")

cell_L <- function(cn) {
  w <- CW[[cn]]
  c(UC = -w[["U"]], UT = w[["U"]], PC = -w[["P"]], PT = w[["P"]])
}

# ---------------------------------------------------------------------------
# Exact numbers and the Python bridge

dyad1 <- function(x) {
  if (!is.finite(x)) stop("dyad1: not finite")
  if (x == round(x) && abs(x) < 2^53) return(sprintf("%.0f", x))
  k <- 0
  m <- x
  while (m != round(m)) {
    m <- m * 2
    k <- k + 1
  }
  sprintf("%.0f/%.0f", m, 2^k)
}
dyad <- function(x) vapply(as.numeric(x), dyad1, character(1))
rows_json <- function(X) lapply(seq_len(nrow(X)), function(i) as.list(dyad(X[i, ])))
frac_inv <- function(s) {
  p <- strsplit(s, "/", fixed = TRUE)[[1]]
  if (length(p) == 1) p <- c(p, "1")
  paste0(p[2], "/", p[1])
}

PY_CALLS <- new.env()
PY_CALLS$n <- 0
PY_CALLS$secs <- 0
# One Python process serves every request (exact.py serve), so that the
# interpreter starts once.
py_server <- function() {
  if (is.null(PY_CALLS$proc) || !PY_CALLS$proc$is_alive()) {
    PY_CALLS$proc <- processx::process$new(REF_PYTHON, c("-I", file.path(REF_DIR, "exact.py"), "serve"),
                                           stdin = "|", stdout = "|", stderr = "|")
  }
  PY_CALLS$proc
}
py_run <- function(mode, input) {
  fin <- tempfile(fileext = ".json")
  fout <- tempfile(fileext = ".json")
  on.exit(unlink(c(fin, fout)))
  jsonlite::write_json(input, fin, auto_unbox = TRUE, digits = NA, null = "null")
  t0 <- proc.time()[["elapsed"]]
  pr <- py_server()
  pr$write_input(paste0(mode, "\t", fin, "\t", fout, "\n"))
  reply <- character(0)
  while (!length(reply)) {
    pr$poll_io(-1)
    reply <- pr$read_output_lines()
    if (!pr$is_alive() && !length(reply)) stop("exact.py server died:\n", paste(pr$read_all_error_lines(), collapse = "\n"))
  }
  PY_CALLS$n <- PY_CALLS$n + 1
  PY_CALLS$secs <- PY_CALLS$secs + proc.time()[["elapsed"]] - t0
  if (!identical(reply[1], "done")) stop("exact.py ", mode, " failed: ", reply[1])
  jsonlite::read_json(fout, simplifyVector = FALSE)
}
py_result <- function(mode, input) {
  out <- py_run(mode, input)
  if (!isTRUE(out$ok)) stop("UNEXPLAINED: ", out$unexplained)
  out$result
}

# ---------------------------------------------------------------------------
# Compositions

enum_comps <- function(m, t) {
  J <- length(m)
  out <- list()
  rec <- function(j, left, cur) {
    if (j > J) {
      if (left == 0) out[[length(out) + 1L]] <<- cur
      return(invisible())
    }
    rest <- if (j < J) sum(m[(j + 1):J]) else 0
    lo <- max(0, left - rest)
    hi <- min(m[j], left)
    if (lo > hi) return(invisible())
    for (v in lo:hi) rec(j + 1L, left - v, c(cur, v))
  }
  rec(1L, t, integer(0))
  do.call(rbind, out)
}
comp_key <- function(v) paste(v, collapse = ",")
comp_mult <- function(m, v) prod(choose(m, v))

# ---------------------------------------------------------------------------
# World preparation: row types, labeling spaces, classes

type_string <- function(d, vars) {
  parts <- lapply(vars, function(v) {
    x <- d[[v]]
    if (is.factor(x) || is.character(x) || is.logical(x)) as.character(x) else sprintf("%a", as.numeric(x))
  })
  do.call(paste, c(parts, sep = "|"))
}

prep_world <- function(w) {
  d <- w$data
  n <- nrow(d)
  inmod <- !is.na(d$y)
  if (isTRUE(w$pre)) inmod <- inmod & (d$pretested == 0 | !is.na(d$pre))
  for (v in w$covs) inmod <- inmod & !is.na(d[[v]])
  if (!is.null(w$exposure)) inmod <- inmod & !is.na(d[[w$exposure]])
  dd <- d
  if (isTRUE(w$pre)) dd$pre[dd$pretested == 0] <- 0
  vars <- c("y", if (isTRUE(w$pre)) "pre", w$covs, w$exposure)
  key <- type_string(dd, vars)
  key[!inmod] <- "OUT"
  # design columns other than the cells
  zl <- list()
  if (isTRUE(w$pre)) zl$pre <- ifelse(d$pretested == 1, d$pre, 0)
  for (v in w$covs) {
    x <- d[[v]]
    if (is.factor(x) || is.character(x)) {
      lev <- levels(factor(x[inmod]))
      for (l in lev[-1]) zl[[paste0(v, l)]] <- as.numeric(as.character(x) == l)
    } else {
      zl[[v]] <- as.numeric(as.matrix(x)[, 1])
    }
  }
  Z <- if (length(zl)) do.call(cbind, zl) else matrix(numeric(0), n, 0)
  Z[is.na(Z)] <- 0
  Z[!inmod, ] <- 0
  e <- if (is.null(w$exposure)) rep(1, n) else d[[w$exposure]]
  cond <- list()
  for (k in c("U", "P")) {
    idx <- which(d$pretested == (k == "P"))
    keys_k <- key[idx]
    ut <- unique(keys_k[keys_k != "OUT"])
    ut <- ut[order(ut)]
    rep_row <- idx[match(ut, keys_k)]
    m_in <- as.integer(table(factor(keys_k[keys_k != "OUT"], levels = ut)))
    n_out <- sum(keys_k == "OUT")
    cond[[k]] <- list(idx = idx, types = ut, rep = rep_row, m = m_in, n_out = n_out,
                      y = d$y[rep_row], Z = Z[rep_row, , drop = FALSE], e = e[rep_row],
                      row_type = match(keys_k, ut))
  }
  w$prep <- list(inmod = inmod, Z = Z, cond = cond, sep = length(w$covs) == 0)
  w
}

# In-model composition of a labeling in condition k: treated rows per type.
labeling_comp <- function(w, treat, k) {
  ck <- w$prep$cond[[k]]
  tr <- treat[ck$idx] == 1
  rt <- ck$row_type
  v <- integer(length(ck$types))
  tab <- table(factor(rt[tr & !is.na(rt)], levels = seq_along(ck$types)))
  v[] <- as.integer(tab)
  v
}

# Labeling spaces. Lambda_all: labels permuted over every row of the
# condition (rows outside the model hold labels too). Lambda_model: over the
# rows in the model only.
build_spaces <- function(w) {
  out <- list()
  t_all <- c(U = sum(w$obs[[1]][w$prep$cond$U$idx]), P = sum(w$obs[[1]][w$prep$cond$P$idx]))
  for (o in names(w$obs)) {
    for (k in c("U", "P")) {
      s <- sum(w$obs[[o]][w$prep$cond[[k]]$idx])
      if (s != t_all[[k]]) stop(w$id, ": observed labeling ", o, " is not in the labeling space")
    }
  }
  for (k in c("U", "P")) {
    ck <- w$prep$cond[[k]]
    m_all <- c(ck$m, if (ck$n_out > 0) ck$n_out)
    comps <- enum_comps(m_all, t_all[[k]])
    mult <- apply(comps, 1, function(v) comp_mult(m_all, v))
    inm <- comps[, seq_along(ck$m), drop = FALSE]
    out$all[[k]] <- list(keys = apply(inm, 1, comp_key), mult = mult, comps = inm)
    if (w$prep$cond$U$n_out + w$prep$cond$P$n_out > 0) {
      out$model[[k]] <- list()
      for (o in names(w$obs)) {
        tm <- sum(labeling_comp(w, w$obs[[o]], k))
        cm <- enum_comps(ck$m, tm)
        out$model[[k]][[o]] <- list(keys = apply(cm, 1, comp_key),
                                    mult = apply(cm, 1, function(v) comp_mult(ck$m, v)), comps = cm)
      }
    }
  }
  out
}

# Unique in-model compositions needed per condition, with the weight of each
# in each space.
needed_comps <- function(w, spaces) {
  res <- list()
  for (k in c("U", "P")) {
    allk <- spaces$all[[k]]
    keys <- allk$keys
    comps <- allk$comps
    if (!is.null(spaces$model)) {
      for (o in names(spaces$model[[k]])) {
        keys <- c(keys, spaces$model[[k]][[o]]$keys)
        comps <- rbind(comps, spaces$model[[k]][[o]]$comps)
      }
    }
    u <- !duplicated(keys)
    res[[k]] <- list(keys = keys[u], comps = comps[u, , drop = FALSE])
  }
  res
}

space_weights <- function(space_k, keys) {
  wt <- tapply(space_k$mult, factor(space_k$keys, levels = keys), sum)
  wt[is.na(wt)] <- 0
  as.numeric(wt)
}

# Aggregated rows of a class: cell (0 UC, 1 UT, 2 PC, 3 PT), z, y, e, f.
class_rows <- function(w, cU, cP, conds = c("U", "P")) {
  rows <- list()
  for (k in conds) {
    ck <- w$prep$cond[[k]]
    v <- if (k == "U") cU else cP
    base <- if (k == "U") 0L else 2L
    for (j in seq_along(ck$m)) {
      for (arm in c(1L, 0L)) {
        f <- if (arm == 1L) v[j] else ck$m[j] - v[j]
        if (f > 0) rows[[length(rows) + 1L]] <- list(cell = base + arm, z = ck$Z[j, ], y = ck$y[j],
                                                     e = ck$e[j], f = f, type = j, cond = k)
      }
    }
  }
  rows
}

design_of <- function(rows, conds = c("U", "P"), zcols = TRUE) {
  cells <- if (identical(conds, "P")) 2:3 else if (identical(conds, "U")) 0:1 else 0:3
  X <- t(vapply(rows, function(r) {
    cx <- as.numeric(cells == r$cell)
    c(cx, if (zcols) r$z)
  }, numeric(length(cells) + if (zcols) length(rows[[1]]$z) else 0)))
  if (length(rows) == 1) X <- matrix(X, nrow = 1)
  list(X = X, y = vapply(rows, `[[`, 0, "y"), e = vapply(rows, `[[`, 0, "e"),
       f = vapply(rows, `[[`, 0, "f"), cell = vapply(rows, `[[`, 0L, "cell"))
}

# ---------------------------------------------------------------------------
# Fits at the exact estimate (iterated models)

fam_parts <- function(fam, eta, y, theta) {
  if (fam == "binomial") {
    mu <- stats::plogis(eta)
    w <- mu * (1 - mu)
    s <- y - mu
    hw <- w
    ll <- y * eta - log1p(exp(-abs(eta))) - pmax(eta, 0)
  } else if (fam == "poisson") {
    mu <- exp(eta)
    w <- mu
    s <- y - mu
    hw <- w
    ll <- y * eta - mu
  } else if (fam == "negbin") {
    mu <- exp(eta)
    w <- mu * theta / (theta + mu)
    s <- (y - mu) * theta / (theta + mu)
    hw <- (y + theta) * theta * mu / (theta + mu)^2
    ll <- y * eta - (y + theta) * log(theta + mu)
  } else {
    mu <- eta
    w <- rep(1, length(eta))
    s <- y - mu
    hw <- w
    ll <- -(y - mu)^2 / 2
  }
  list(mu = mu, w = w, s = s, hw = hw, ll = ll)
}

# Damped Newton on standardized columns. Converged when the relative step is
# below 1e-12 on two successive iterations and the scaled score below 1e-10.
newton_fit <- function(X, y, f, off, fam, theta = NULL, maxit = 400, reverse = FALSE) {
  ord <- if (reverse) rev(seq_len(nrow(X))) else seq_len(nrow(X))
  X <- X[ord, , drop = FALSE]
  y <- y[ord]
  f <- f[ord]
  off <- off[ord]
  sc <- apply(abs(X), 2, max)
  sc[sc == 0] <- 1
  Xs <- sweep(X, 2, sc, "/")
  p <- ncol(Xs)
  b <- numeric(p)
  if (fam == "gaussian") {
    A <- crossprod(Xs * f, Xs)
    b <- solve(A, crossprod(Xs * f, y - off))[, 1]
    conv <- TRUE
    it <- 1
  } else {
    conv <- FALSE
    small <- 0
    eta <- drop(Xs %*% b) + off
    fp <- fam_parts(fam, eta, y, theta)
    ll <- sum(f * fp$ll)
    for (it in seq_len(maxit)) {
      g <- drop(crossprod(Xs, f * fp$s))
      H <- crossprod(Xs * (f * fp$hw), Xs)
      step <- tryCatch(solve(H, g), error = function(e) NULL)
      if (is.null(step)) break
      lam <- 1
      repeat {
        bn <- b + lam * step
        etan <- drop(Xs %*% bn) + off
        fpn <- fam_parts(fam, etan, y, theta)
        lln <- sum(f * fpn$ll)
        if (is.finite(lln) && lln >= ll - 1e-12 * abs(ll)) break
        lam <- lam / 2
        if (lam < 1e-10) break
      }
      b <- bn
      fp <- fpn
      ll <- lln
      rel <- max(abs(lam * step)) / max(1, max(abs(b)))
      gs <- drop(crossprod(Xs, f * fp$s)) / sqrt(pmax(diag(crossprod(Xs * (f * fp$w), Xs)), 1e-300))
      small <- if (rel < 1e-12) small + 1 else 0
      if (small >= 2 && max(abs(gs)) < 1e-10) {
        conv <- TRUE
        break
      }
      if (fam == "binomial" && max(abs(etan - off)) > 60) break
      if (fam != "binomial" && min(fp$mu) < 1e-300) break
    }
  }
  eta <- drop(Xs %*% b) + off
  fp <- fam_parts(fam, eta, y, theta)
  back <- order(ord)
  list(b_std = b, sc = sc, b = b / sc, eta = eta[back], mu = fp$mu[back], w = fp$w[back],
       s = fp$s[back], converged = conv, iter = it)
}

# Rows with leverage exactly 1 (`skip`, certified by exact.py struct) on
# which the contrast does not depend add exactly 0 to the HC3 sum (their
# c_i is 0 and so is 1 - h_i); they are left out rather than computed as
# rounding noise divided by rounding noise.
hc3_contrast <- function(X, f, fit, L, skip = integer(0)) {
  A <- crossprod(X * (f * fit$w), X)
  Ai <- tryCatch(solve(A), error = function(e) NULL)
  if (is.null(Ai)) return(list(est = NA_real_, V = NA_real_, hmax = NA_real_, kappa = Inf))
  est <- sum(L * fit$b)
  g <- drop(Ai %*% L)
  ci <- drop(X %*% g)
  h <- fit$w * rowSums((X %*% Ai) * X)
  use <- setdiff(seq_len(nrow(X)), skip)
  V <- sum((f * (ci * fit$s / (1 - h))^2)[use])
  h <- h[use]
  sc <- sqrt(diag(A))
  kappa <- kappa(A / outer(sc, sc), exact = TRUE)
  list(est = est, V = V, hmax = max(h), kappa = kappa)
}

# Recession directions by a linear program (a screen only; every direction is
# verified exactly by exact.py, and the absence of one is certified only by
# mlecert on the reduced model).
lp_recession <- function(X, y, f, fam) {
  sc <- apply(abs(X), 2, max)
  sc[sc == 0] <- 1
  Xs <- sweep(X, 2, sc, "/")
  n <- nrow(Xs)
  p <- ncol(Xs)
  B <- 1e4
  obj <- c(rep(0, 2 * p), rep(1, n))
  con <- NULL
  dir <- NULL
  rhs <- NULL
  for (i in seq_len(n)) {
    if (fam == "binomial") {
      s <- if (y[i] == 1) 1 else -1
      con <- rbind(con, c(s * Xs[i, ], -s * Xs[i, ], -as.numeric(seq_len(n) == i)))
      dir <- c(dir, ">=")
      rhs <- c(rhs, 0)
    } else if (y[i] > 0) {
      con <- rbind(con, c(Xs[i, ], -Xs[i, ], rep(0, n)))
      dir <- c(dir, "=")
      rhs <- c(rhs, 0)
      con <- rbind(con, c(rep(0, 2 * p), as.numeric(seq_len(n) == i)))
      dir <- c(dir, "<=")
      rhs <- c(rhs, 0)
    } else {
      con <- rbind(con, c(-Xs[i, ], Xs[i, ], -as.numeric(seq_len(n) == i)))
      dir <- c(dir, ">=")
      rhs <- c(rhs, 0)
    }
  }
  con <- rbind(con, cbind(diag(2 * p), matrix(0, 2 * p, n)), cbind(matrix(0, n, 2 * p), diag(n)))
  dir <- c(dir, rep("<=", 2 * p + n))
  rhs <- c(rhs, rep(B, 2 * p), rep(1, n))
  sol <- lpSolve::lp("max", obj, con, dir, rhs)
  if (sol$status != 0) return(NULL)
  x <- sol$solution
  d <- (x[seq_len(p)] - x[p + seq_len(p)]) / sc
  t <- x[2 * p + seq_len(n)]
  list(I0 = which(t > 1e-7), d = d)
}

# Rank facts depend on the design only, not on the outcome or theta: cached.
STRUCT_CACHE <- new.env()
cached_struct <- function(units) {
  keys <- vapply(units, function(u) digest::digest(u, algo = "sha1"), "")
  miss <- which(!vapply(keys, exists, logical(1), envir = STRUCT_CACHE, inherits = FALSE))
  miss <- miss[!duplicated(keys[miss])]
  if (length(miss)) {
    res <- py_result("struct", list(units = units[miss]))
    for (j in seq_along(miss)) assign(keys[miss[j]], res[[j]], envir = STRUCT_CACHE)
  }
  lapply(keys, get, envir = STRUCT_CACHE, inherits = FALSE)
}

# Evaluate a batch of iterated fits. Each item: X, y, f, e (exposure), L (a
# named list of contrast vectors on X's columns). Returns, per item, the
# status of the estimate ("mle" or "limit" with the rows kept, or
# "undecided") and per contrast: est, V, eps (measured error), code
# ("finite", "H1", "NE", "DIV" = not estimable in the limiting model), and
# what band() needs.
eval_iter_batch <- function(items, fam, theta = NULL) {
  nI <- length(items)
  if (nI == 0) return(list())
  res <- vector("list", nI)
  keep <- lapply(items, function(it) seq_len(nrow(it$X)))
  status <- rep("todo", nI)
  # structure of the full designs (estimability, leverage 1)
  st_full <- cached_struct(lapply(items, function(it)
    list(X = rows_json(it$X), f = as.list(it$f), L = lapply(it$L, function(l) as.list(dyad(l))))))
  fits <- vector("list", nI)
  stt <- st_full
  for (round in 1:4) {
    todo <- which(status == "todo")
    if (!length(todo)) break
    cert_items <- list()
    cert_idx <- integer(0)
    for (i in todo) {
      it <- items[[i]]
      r <- keep[[i]]
      X <- it$X[r, , drop = FALSE]
      piv <- unlist(stt[[i]]$piv) + 1
      if (!length(piv)) piv <- seq_len(ncol(X))
      fit <- newton_fit(X[, piv, drop = FALSE], it$y[r], it$f[r], log(it$e[r]), fam, theta)
      fits[[i]] <- list(fit = fit, piv = piv)
      # The negative binomial (log link, theta fixed) has the directions of
      # recession of the Poisson model, so the Poisson fit's means certify
      # that its estimate exists.
      cfam <- if (fam == "negbin") "poisson" else fam
      cfit <- if (fam == "negbin") newton_fit(X[, piv, drop = FALSE], it$y[r], it$f[r], log(it$e[r]), "poisson") else fit
      if (cfit$converged && fit$converged) {
        cert_items[[length(cert_items) + 1]] <- list(fam = cfam, X = rows_json(X[, piv, drop = FALSE]),
                                                     y = as.list(dyad(it$y[r])), f = as.list(it$f[r]),
                                                     mu = as.list(dyad(cfit$mu)))
        cert_idx <- c(cert_idx, i)
      } else if (cfit$converged && !fit$converged) {
        # the estimate exists but the solver did not converge: the class is
        # left undecided and reported (design Section 5)
        warning("solver failure on a class with an estimate")
      }
    }
    ok <- rep(FALSE, nI)
    if (length(cert_items)) {
      cr <- unlist(py_result("mlecert", list(items = cert_items)))
      ok[cert_idx] <- cr
    }
    for (i in todo[ok[todo]]) status[i] <- if (length(keep[[i]]) == nrow(items[[i]]$X)) "mle" else "limit"
    # the rest: directions of recession
    rest <- todo[!ok[todo]]
    if (!length(rest)) break
    dir_items <- list()
    dir_idx <- integer(0)
    cand <- list()
    for (i in rest) {
      it <- items[[i]]
      r <- keep[[i]]
      lp <- lp_recession(it$X[r, , drop = FALSE], it$y[r], it$f[r], fam)
      if (is.null(lp) || !length(lp$I0)) {
        status[i] <- "undecided"
        next
      }
      cand[[as.character(i)]] <- lp
      dir_items[[length(dir_items) + 1]] <- list(fam = fam, X = rows_json(it$X[r, , drop = FALSE]),
                                                 y = as.list(dyad(it$y[r])), I0 = as.list(lp$I0 - 1L),
                                                 d = as.list(unname(lp$d)))
      dir_idx <- c(dir_idx, i)
    }
    if (length(dir_items)) {
      dr <- unlist(py_result("dircert", list(items = dir_items)))
      for (j in seq_along(dir_idx)) {
        i <- dir_idx[j]
        if (!dr[j]) {
          status[i] <- "undecided"
          next
        }
        keep[[i]] <- keep[[i]][-cand[[as.character(i)]]$I0]
        # every row diverges: no contrast is estimable in the limiting model
        if (!length(keep[[i]])) status[i] <- "limit_noest"
      }
    }
    # structure of the reduced designs
    red <- which(status == "todo")
    if (length(red)) {
      sr <- cached_struct(lapply(red, function(i) {
        it <- items[[i]]
        list(X = rows_json(it$X[keep[[i]], , drop = FALSE]), f = as.list(it$f[keep[[i]]]),
             L = lapply(it$L, function(l) as.list(dyad(l))))
      }))
      for (j in seq_along(red)) stt[[red[j]]] <- sr[[j]]
      # a reduced model in which every contrast is not estimable needs no fit
      for (i in red) {
        if (!any(vapply(stt[[i]]$c, function(z) isTRUE(z$est), logical(1)))) status[i] <- "limit_noest"
      }
    }
  }
  status[status == "todo"] <- "undecided"
  for (i in seq_len(nI)) {
    it <- items[[i]]
    out <- list(status = status[i], keep = keep[[i]], c = list())
    for (cn in names(it$L)) {
      cf <- st_full[[i]]$c[[cn]]
      if (!isTRUE(cf$est)) {
        # NE: an arm the contrast uses has no rows; ALIAS: a covariate is
        # collinear with the cells under this labeling
        out$c[[cn]] <- list(code = if (isTRUE(cf$alias)) "ALIAS" else "NE")
        next
      }
      if (status[i] == "undecided") {
        out$c[[cn]] <- list(code = "UNDECIDED", h1 = isTRUE(cf$h1))
        next
      }
      cs <- stt[[i]]$c[[cn]]
      if (!isTRUE(cs$est) || status[i] == "limit_noest") {
        out$c[[cn]] <- list(code = "DIV", h1 = isTRUE(cf$h1))
        next
      }
      r <- keep[[i]]
      X <- it$X[r, , drop = FALSE]
      piv <- fits[[i]]$piv
      fit <- fits[[i]]$fit
      L <- it$L[[cn]][piv]
      lev1 <- unlist(stt[[i]]$lev1_rows) + 1L
      hc <- hc3_contrast(X[, piv, drop = FALSE], it$f[r], fit, L, lev1)
      fit2 <- newton_fit(X[, piv, drop = FALSE], it$y[r], it$f[r], log(it$e[r]), fam, theta, reverse = TRUE)
      hc2 <- hc3_contrast(X[, piv, drop = FALSE], it$f[r], fit2, L, lev1)
      out$c[[cn]] <- list(code = if (isTRUE(cf$h1) || isTRUE(cs$h1)) "H1" else "finite",
                          h1 = isTRUE(cf$h1) || isTRUE(cs$h1),
                          est = hc$est, V = hc$V, est2 = hc2$est, V2 = hc2$V,
                          ill = !is.finite(hc$kappa) || hc$kappa > 1e10 || 1 - hc$hmax < 1e-6,
                          piv = piv, b = fit$b, skip = lev1)
    }
    res[[i]] <- out
  }
  res
}

# ---------------------------------------------------------------------------
# Negative binomial: theta by maximum likelihood (profile score solved
# tightly), on the observed data.

nb_theta <- function(X, y, f, off) {
  score <- function(lt) {
    th <- exp(lt)
    fit <- newton_fit(X, y, f, off, "negbin", th)
    if (!fit$converged) return(NA_real_)
    mu <- fit$mu
    sum(f * (digamma(y + th) - digamma(th) + log(th) + 1 - log(th + mu) - (y + th) / (th + mu)))
  }
  grid <- log(c(1e-3, 1e-2, 0.1, 0.5, 1, 2, 5, 10, 30, 100, 300, 1e3))
  sv <- vapply(grid, score, 0)
  j <- which(diff(sign(sv)) != 0)
  if (!length(j) || anyNA(sv[c(j[1], j[1] + 1)])) return(list(theta = NA_real_, boundary = TRUE, sv = sv))
  j <- j[1]
  r <- stats::uniroot(score, c(grid[j], grid[j + 1]), tol = 1e-14, maxiter = 200)
  list(theta = exp(r$root), boundary = FALSE, sv = sv)
}

# ---------------------------------------------------------------------------
# Data-level rules (design, Section 3.4, rules 3 and 4)

data_rule <- function(w, cn) {
  wts <- CW[[cn]]
  used <- names(wts)[wts != 0]
  constant <- character(0)
  boundary <- FALSE
  for (k in used) {
    ck <- w$prep$cond[[k]]
    yv <- unique(ck$y)
    if (length(yv) == 1) {
      constant <- c(constant, k)
      if ((w$family == "binomial" && yv %in% c(0, 1)) || (w$family %in% c("poisson", "negbin") && yv == 0))
        boundary <- TRUE
    }
  }
  if (!length(constant)) return("ok")
  if (boundary) return("ERR")
  if (w$prep$sep && length(used) == 1) return("ERR")
  "ERR-S"
}
