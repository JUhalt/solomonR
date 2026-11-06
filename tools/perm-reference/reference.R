# The reference for one world: every class of labelings, its codes, its
# value, and its relation (-1, 0, 1) to each observed class, decided exactly
# or by the rule of design Section 5. Sourced by generate.R after lib.R.

ALLC <- CONTRASTS

# ---------------------------------------------------------------------------
# Gaussian identity: every class in exact rationals (exact.py gauss).

ref_gauss <- function(w, nc) {
  nU <- length(nc$U$keys)
  nP <- length(nc$P$keys)
  k <- ncol(w$prep$Z)
  classes <- vector("list", nU * nP)
  for (iU in seq_len(nU)) for (iP in seq_len(nP)) {
    rows <- class_rows(w, nc$U$comps[iU, ], nc$P$comps[iP, ])
    classes[[(iU - 1) * nP + iP]] <- list(
      cell = lapply(rows, `[[`, "cell"), z = lapply(rows, function(r) as.list(dyad(r$z))),
      y = as.list(dyad(vapply(rows, `[[`, 0, "y"))), f = as.list(vapply(rows, `[[`, 0, "f")))
  }
  inp <- list(k = k, classes = classes,
              contrasts = lapply(stats::setNames(ALLC, ALLC), function(cn) as.list(dyad(c(cell_L(cn), rep(0, k))))))
  r <- py_result("gauss", inp)
  out <- list()
  for (cn in ALLC) for (s in STATS) {
    rr <- r[[cn]][[if (s == "studentized") "stud" else "diff"]]
    code <- unlist(rr$code)
    rank <- vapply(rr$rank, function(v) if (is.null(v)) NA_real_ else as.numeric(v), 0)
    approx <- vapply(rr$approx, function(v) if (is.null(v)) NA_real_ else as.numeric(v), 0)
    out[[cn]][[s]] <- list(code = code, approx = approx, rank = rank, how = "exact rationals")
  }
  out
}

rel_from_rank <- function(rank, o) {
  r <- sign(rank - rank[o])
  r
}

# ---------------------------------------------------------------------------
# Exact keys of iterated fits (design E1 to E3): equal keys give equal |T| in
# exact arithmetic. Rows are grouped by (cell, covariates, exposure); `perm`
# relabels the cells (a swap of the arms). For the difference statistic in a
# canonical-link model (logit, Poisson log) the likelihood, and so the
# estimate, depends on the outcome only through the design and X'y; with the
# negative binomial (theta fixed) the score depends on the group sums. The
# HC3 variance needs, per group, n, the sum, and the sum of squares.

suff_key <- function(rows, fam, stat, perm = NULL) {
  cell <- vapply(rows, function(r) as.integer(r$cell), 0L)
  if (!is.null(perm)) cell <- perm[cell + 1L]
  zs <- vapply(rows, function(r) paste(sprintf("%a", c(r$z, r$e)), collapse = ","), "")
  grp <- paste(cell, zs, sep = "|")
  f <- vapply(rows, `[[`, 0, "f")
  y <- vapply(rows, `[[`, 0, "y")
  n_g <- tapply(f, grp, sum)
  design <- paste(names(n_g), n_g, sep = ":", collapse = ";")
  if (stat == "difference" && fam %in% c("binomial", "poisson")) {
    Z <- do.call(rbind, lapply(rows, function(r) r$z))
    X <- cbind(outer(cell, 0:3, "==") * 1, if (length(Z)) Z)
    xty <- colSums(X * (f * y))
    return(paste(design, paste(sprintf("%a", xty), collapse = ","), sep = "#"))
  }
  S_g <- tapply(f * y, grp, sum)
  if (stat == "difference") return(paste(design, paste(sprintf("%a", S_g), collapse = ","), sep = "#"))
  Q_g <- tapply(f * y^2, grp, sum)
  paste(design, paste(sprintf("%a", S_g), collapse = ","), paste(sprintf("%a", Q_g), collapse = ","), sep = "#")
}

PERMS <- list(id = c(0L, 1L, 2L, 3L), sU = c(1L, 0L, 2L, 3L), sP = c(0L, 1L, 3L, 2L), sUP = c(1L, 0L, 3L, 2L))
# swaps that keep |T| for each contrast (design E2), and those that reverse its sign
SWAPS <- list("Treatment | unpretested" = c("id", "sU", "sP", "sUP"),
              "Treatment | pretested" = c("id", "sU", "sP", "sUP"),
              "ATE (avg over pretest)" = c("id", "sUP"), "Pretest x Treatment" = c("id", "sUP"))
REVERSE <- list("Treatment | unpretested" = c("sU", "sUP"), "Treatment | pretested" = c("sP", "sUP"),
                "ATE (avg over pretest)" = "sUP", "Pretest x Treatment" = "sUP")

# ---------------------------------------------------------------------------
# Parts of a separable iterated model: one condition at a time.

part_eval <- function(w, k, comps, theta, saturated) {
  ck <- w$prep$cond[[k]]
  fam <- w$family
  n <- nrow(comps)
  parts <- vector("list", n)
  sat_in <- list()
  sat_idx <- integer(0)
  num_items <- list()
  num_idx <- integer(0)
  for (i in seq_len(n)) {
    v <- comps[i, ]
    mc <- ck$m - v
    nt <- sum(v)
    ncn <- sum(mc)
    St <- sum(v * ck$y)
    Sc <- sum(mc * ck$y)
    p <- list(nt = nt, nc = ncn, St = St, Sc = Sc, ne = nt == 0 || ncn == 0,
              key_o = paste0("t:", comp_key(v), "|c:", comp_key(mc)),
              key_s = paste0("t:", comp_key(mc), "|c:", comp_key(v)),
              zero_self = all(v == mc), sigma = 0, sat = saturated)
    if (!p$ne) {
      emp <- function(S, nn) if (fam == "binomial") S == 0 || S == nn else S == 0
      top <- function(S, nn) fam == "binomial" && S == nn
      et <- emp(St, nt)
      ec <- emp(Sc, ncn)
      if (et || ec) {
        st <- if (et) (if (top(St, nt)) 1 else -1) else 0
        sc <- if (ec) (if (top(Sc, ncn)) -1 else 1) else 0
        p$sigma <- if (st != 0) st else sc
        p$emptied <- TRUE
      } else {
        p$emptied <- FALSE
      }
    }
    if (saturated) {
      p$h1 <- !p$ne && (nt == 1 || ncn == 1)
      if (!p$ne && !p$emptied) {
        arm <- function(vv) lapply(which(vv > 0), function(j) list(dyad(ck$y[j]), dyad(ck$e[j]), vv[j]))
        sat_in[[length(sat_in) + 1]] <- list(t = arm(v), c = arm(mc))
        sat_idx <- c(sat_idx, i)
      }
    } else if (!p$ne && !p$emptied) {
      rows <- list()
      for (j in seq_along(ck$m)) for (arm in c(1L, 0L)) {
        f <- if (arm == 1L) v[j] else mc[j]
        if (f > 0) rows[[length(rows) + 1]] <- list(cell = arm, z = ck$Z[j, ], y = ck$y[j], e = ck$e[j], f = f)
      }
      X <- t(vapply(rows, function(r) c(as.numeric(c(0, 1) == r$cell), r$z), numeric(2 + ncol(ck$Z))))
      num_items[[length(num_items) + 1]] <- list(
        X = X, y = vapply(rows, `[[`, 0, "y"), e = vapply(rows, `[[`, 0, "e"),
        f = vapply(rows, `[[`, 0, "f"), L = list(own = c(-1, 1, rep(0, ncol(ck$Z)))), rows = rows)
      num_idx <- c(num_idx, i)
      # exact keys of the part, ordered (control, treated) and swapped
      p$key_o <- suff_key(rows, fam, "studentized")
      p$key_s <- suff_key(rows, fam, "studentized", perm = c(1L, 0L))
      p$dkey_o <- suff_key(rows, fam, "difference")
      p$dkey_s <- suff_key(rows, fam, "difference", perm = c(1L, 0L))
    }
    parts[[i]] <- p
  }
  if (length(sat_in)) {
    sp <- py_result("satparts", list(fam = fam, parts = sat_in))
    for (j in seq_along(sat_idx)) {
      i <- sat_idx[j]
      R <- sp[[j]]$R
      V <- sp[[j]]$V
      parts[[i]]$R <- R
      parts[[i]]$V <- V
      parts[[i]]$e <- log(eval(parse(text = R)))
      parts[[i]]$Vd <- eval(parse(text = V))
      parts[[i]]$zero <- R == "1/1"
      parts[[i]]$Vzero <- V == "0/1"
      parts[[i]]$key_o <- paste0("R:", R, "|V:", V)
      parts[[i]]$key_s <- paste0("R:", frac_inv(R), "|V:", V)
      parts[[i]]$dkey_o <- paste0("R:", R)
      parts[[i]]$dkey_s <- paste0("R:", frac_inv(R))
    }
  }
  if (length(num_items)) {
    ev <- eval_iter_batch(num_items, fam, theta)
    for (j in seq_along(num_idx)) {
      i <- num_idx[j]
      cc <- ev[[j]]$c$own
      parts[[i]]$num <- cc
      parts[[i]]$status <- ev[[j]]$status
      parts[[i]]$h1 <- isTRUE(cc$h1)
      parts[[i]]$band <- list(item = num_items[[j]], keep = ev[[j]]$keep, piv = cc$piv, b = cc$b, skip = cc$skip)
      if (identical(ev[[j]]$status, "limit")) {
        rk <- num_items[[j]]$rows[ev[[j]]$keep]
        parts[[i]]$key_o <- paste0("LIM:", suff_key(rk, fam, "studentized"))
        parts[[i]]$key_s <- paste0("LIM:", suff_key(rk, fam, "studentized", perm = c(1L, 0L)))
        parts[[i]]$dkey_o <- paste0("LIM:", suff_key(rk, fam, "difference"))
        parts[[i]]$dkey_s <- paste0("LIM:", suff_key(rk, fam, "difference", perm = c(1L, 0L)))
      }
      # exactly zero when the part's keys are invariant under the arm swap
      parts[[i]]$zero <- identical(parts[[i]]$dkey_o, parts[[i]]$dkey_s)
    }
  }
  parts
}

# ---------------------------------------------------------------------------
# Combination of parts into class codes, values, and exact keys.

combine_parts <- function(w, pU, pP, cn, s) {
  wt <- CW[[cn]]
  used <- names(wt)[wt != 0]
  P <- list(U = pU, P = pP)
  stud <- s == "studentized"
  for (k in used) if (isTRUE(P[[k]]$ne)) return(list(code = "NE"))
  for (k in used) if (!P[[k]]$sat && identical(P[[k]]$num$code, "ALIAS")) return(list(code = "ALIAS"))
  if (stud) for (k in used) if (isTRUE(P[[k]]$h1)) return(list(code = "H1"))
  emp <- used[vapply(used, function(k) isTRUE(P[[k]]$emptied), logical(1))]
  if (length(emp)) {
    if (length(emp) == 2) {
      sg <- sign(wt[emp]) * c(P$U$sigma, P$P$sigma)[match(emp, c("U", "P"))]
      if (sg[1] != sg[2]) return(list(code = "INF-X"))
    }
    return(list(code = "INF-E"))
  }
  for (k in used) {
    if (!P[[k]]$sat) {
      cd <- P[[k]]$num$code
      if (identical(cd, "DIV")) return(list(code = "SEP-L"))
      if (identical(cd, "UNDECIDED")) return(list(code = "UNDECIDED"))
    }
  }
  num <- 0
  V <- 0
  num2 <- 0
  V2 <- 0
  ill <- FALSE
  exact_zero <- TRUE
  vzero <- TRUE
  limit <- FALSE
  for (k in used) {
    p <- P[[k]]
    if (p$sat) {
      num <- num + wt[[k]] * p$e
      num2 <- num2 + wt[[k]] * p$e
      V <- V + wt[[k]]^2 * p$Vd
      V2 <- V2 + wt[[k]]^2 * p$Vd
      exact_zero <- exact_zero && p$zero
      vzero <- vzero && p$Vzero
    } else {
      num <- num + wt[[k]] * p$num$est
      num2 <- num2 + wt[[k]] * p$num$est2
      V <- V + wt[[k]]^2 * p$num$V
      V2 <- V2 + wt[[k]]^2 * p$num$V2
      exact_zero <- exact_zero && p$zero
      vzero <- FALSE
      ill <- ill || isTRUE(p$num$ill)
      limit <- limit || identical(p$status, "limit")
    }
  }
  # keys: equal keys give equal |T| exactly (design E1 to E3a)
  ko <- function(k) if (stud) P[[k]]$key_o else P[[k]]$dkey_o
  ks <- function(k) if (stud) P[[k]]$key_s else P[[k]]$dkey_s
  if (length(used) == 1) {
    key <- min(ko(used), ks(used))
  } else {
    v <- c(paste(ko("U"), ko("P")), paste(ks("U"), ks("P")))
    if (isTRUE(P$U$zero) || isTRUE(P$P$zero)) v <- c(v, paste(ks("U"), ko("P")), paste(ko("U"), ks("P")))
    key <- min(v)
  }
  if (!stud) {
    if (exact_zero) return(list(code = "ZERO", key = key, val = 0, signed = 0, eps = 0))
    val <- abs(num)
    return(list(code = if (limit) "LIM" else "finite", key = key, val = val, signed = num,
                eps = abs(num - num2) + 1e-15 * val, ill = ill))
  }
  if (vzero) {
    if (exact_zero) return(list(code = "ZZ", key = key, val = 0, signed = 0, eps = 0))
    return(list(code = "INF-C", key = key))
  }
  if (exact_zero) return(list(code = "ZERO", key = key, val = 0, signed = 0, eps = 0))
  t1 <- num / sqrt(V)
  t2 <- num2 / sqrt(V2)
  list(code = if (limit) "LIM" else "finite", key = key, val = abs(t1), signed = t1,
       eps = abs(abs(t1) - abs(t2)) + 1e-15 * abs(t1), ill = ill)
}

INF_CODES <- c("INF-E", "INF-X", "INF-C")
CONV_CODES <- c("NE", "H1", "SEP-L", "UNDECIDED")

# ---------------------------------------------------------------------------
# Relations of numeric classes to an observed class (design Section 5).

numeric_relations <- function(tab, o, s, band_fun) {
  n <- length(tab$code)
  rel <- rep(NA_real_, n)
  co <- tab$code[o]
  # SEP-L is +Inf for the difference statistic (fixed) and, as proposed, for
  # the studentized statistic (D10); case_counts() also gives the counts with
  # those labelings left out
  isinf <- tab$code %in% c(INF_CODES, "SEP-L")
  if (co %in% c("NE", "H1", "UNDECIDED", "ALIAS")) return(list(rel = rel, obs_undefined = co))
  fin <- tab$code %in% c("finite", "LIM", "ZERO", "ZZ")
  if (isinf[o]) {
    rel[isinf] <- 0
    rel[fin] <- -1
    return(list(rel = rel))
  }
  rel[isinf] <- 1
  vo <- tab$val[o]
  if (vo == 0) {
    rel[fin] <- ifelse(tab$val[fin] == 0 & tab$code[fin] %in% c("ZERO", "ZZ"), 0, 1)
    # a finite value that is not certified zero but tiny is decided below
    unc <- which(fin & !(tab$code %in% c("ZERO", "ZZ")))
    rel[unc] <- 1
    return(list(rel = rel, band = integer(0)))
  }
  same <- fin & tab$key == tab$key[o]
  rel[same] <- 0
  rest <- which(fin & !same)
  g <- abs(tab$val[rest] - vo) / vo
  epsm <- pmax(tab$eps[rest], tab$eps[o])
  cert <- g >= 1e-9 & g * vo >= 1e4 * epsm & !tab$ill[rest] & !isTRUE(tab$ill[o])
  cert[tab$code[rest] %in% c("ZERO", "ZZ")] <- TRUE
  rel[rest[cert]] <- sign(tab$val[rest[cert]] - vo)
  bandi <- rest[!cert]
  list(rel = rel, band = bandi)
}

# ---------------------------------------------------------------------------
# Band: 70-digit Newton for the classes doubles could not order.

band_values <- function(w, specs, theta) {
  if (!length(specs)) return(list())
  items <- lapply(specs, function(sp) {
    X <- sp$X[sp$keep, sp$piv, drop = FALSE]
    list(fam = w$family, X = rows_json(X), y = as.list(dyad(sp$y[sp$keep])), f = as.list(sp$f[sp$keep]),
         e = if (any(sp$e != 1)) as.list(dyad(sp$e[sp$keep])) else NULL,
         theta = if (!is.null(theta)) dyad(theta) else NULL,
         b0 = as.list(dyad(sp$b)), parts = sp$parts, skip = as.list(as.integer(sp$skip) - 1L))
  })
  py_result("band", list(items = items))
}

# ---------------------------------------------------------------------------
# The world's class tables.

ref_world_tables <- function(w, nc, theta = NULL) {
  fam <- w$family
  nU <- length(nc$U$keys)
  nP <- length(nc$P$keys)
  obs_idx <- vapply(w$obs, function(t) {
    iU <- match(comp_key(labeling_comp(w, t, "U")), nc$U$keys)
    iP <- match(comp_key(labeling_comp(w, t, "P")), nc$P$keys)
    (iU - 1) * nP + iP
  }, 0)
  if (fam == "gaussian") {
    tabs <- ref_gauss(w, nc)
    out <- list()
    for (cn in ALLC) for (s in STATS) {
      tb <- tabs[[cn]][[s]]
      rels <- lapply(seq_along(obs_idx), function(j) {
        o <- obs_idx[j]
        if (tb$code[o] %in% c("NE", "H1", "ALIAS")) return(list(rel = rep(NA_real_, length(tb$code)), obs_undefined = tb$code[o]))
        list(rel = rel_from_rank(tb$rank, o), how = rep("exact", length(tb$code)))
      })
      names(rels) <- names(w$obs)
      out[[cn]][[s]] <- c(tb, list(rel = rels))
    }
    return(list(tables = out, obs_idx = obs_idx, nU = nU, nP = nP))
  }
  sepU_sat <- w$prep$sep && (is.null(w$exposure) || fam == "poisson")
  if (sepU_sat) return(ref_iter_sep(w, nc, theta, obs_idx))
  ref_iter_full(w, nc, theta, obs_idx)
}

ref_iter_sep <- function(w, nc, theta, obs_idx) {
  fam <- w$family
  nU <- length(nc$U$keys)
  nP <- length(nc$P$keys)
  satP <- !isTRUE(w$pre) && (is.null(w$exposure) || fam == "poisson")
  pU <- part_eval(w, "U", nc$U$comps, theta, TRUE)
  pP <- part_eval(w, "P", nc$P$comps, theta, satP)
  out <- list()
  log <- character(0)
  for (cn in ALLC) for (s in STATS) {
    wt <- CW[[cn]]
    used <- names(wt)[wt != 0]
    allsat <- all(c(U = TRUE, P = satP)[used])
    n <- nU * nP
    code <- character(n)
    val <- signed <- eps <- rep(NA_real_, n)
    key <- character(n)
    ill <- logical(n)
    for (iU in seq_len(nU)) for (iP in seq_len(nP)) {
      i <- (iU - 1) * nP + iP
      cb <- combine_parts(w, pU[[iU]], pP[[iP]], cn, s)
      code[i] <- cb$code
      if (!is.null(cb$val)) {
        val[i] <- cb$val
        signed[i] <- cb$signed
        eps[i] <- cb$eps
        key[i] <- cb$key
        ill[i] <- isTRUE(cb$ill)
      }
    }
    tab <- list(code = code, val = val, approx = signed, eps = eps, key = key, ill = ill)
    rels <- list()
    if (allsat) {
      # exact relations from exact.py satcmp, for the classes that are finite
      # in exact arithmetic (codes decided above: NE, H1, INF-E, INF-X)
      fin <- which(code %in% c("finite", "ZERO", "ZZ", "INF-C"))
      parts_in <- list()
      pidx <- list(U = integer(nU), P = integer(nP))
      for (k in used) {
        pl <- if (k == "U") pU else pP
        ck <- w$prep$cond[[k]]
        comps <- if (k == "U") nc$U$comps else nc$P$comps
        for (j in seq_along(pl)) {
          if (!is.null(pl[[j]]$R)) {
            v <- comps[j, ]
            mc <- ck$m - v
            arm <- function(vv) lapply(which(vv > 0), function(q) list(dyad(ck$y[q]), dyad(ck$e[q]), vv[q]))
            parts_in[[length(parts_in) + 1]] <- list(t = arm(v), c = arm(mc))
            pidx[[k]][j] <- length(parts_in)
          }
        }
      }
      cls <- lapply(fin, function(i) {
        iU <- (i - 1) %/% nP + 1
        iP <- (i - 1) %% nP + 1
        list(if (wt[["U"]] != 0) pidx$U[iU] - 1 else 0, if (wt[["P"]] != 0) pidx$P[iP] - 1 else 0)
      })
      sgn <- list(c(sign(wt[["U"]]), sign(wt[["P"]])))
      names(sgn) <- cn
      obs_fin <- match(obs_idx, fin)
      ok_obs <- which(!is.na(obs_fin))
      sc <- if (length(fin)) py_result("satcmp", list(fam = fam, parts = parts_in, contrasts = sgn,
                                                    classes = cls, obs = as.list(obs_fin[ok_obs] - 1)))[[cn]] else NULL
      ss <- if (s == "studentized") "stud" else "diff"
      if (length(fin)) {
        ap <- sc[[ss]]$approx
        # satcmp's value is sum(sign * log R); the difference statistic
        # carries the condition weight (1/2 for the ATE)
        tab$approx[fin] <- vapply(ap, function(v) if (is.null(v)) NA_real_ else as.numeric(v), 0) *
          (if (s == "difference") abs(wt[used])[[1]] else 1)
        tab$code[fin] <- unlist(sc[[ss]]$code)
      }
      code <- tab$code
      for (j in seq_along(obs_idx)) {
        o <- obs_idx[j]
        rel <- rep(NA_real_, n)
        isinf <- code %in% INF_CODES
        if (code[o] %in% c("NE", "H1", "ALIAS")) {
          rels[[j]] <- list(rel = rel, obs_undefined = code[o])
          next
        }
        if (code[o] %in% INF_CODES) {
          rel[isinf] <- 0
          rel[code %in% c("finite", "ZERO", "ZZ")] <- -1
        } else {
          q <- match(j, ok_obs)
          rel[fin] <- unlist(sc[[ss]]$rel[[q]])
          rel[isinf] <- 1
        }
        rels[[j]] <- list(rel = rel, how = "exact (closed form, integer relations)")
      }
    } else {
      for (j in seq_along(obs_idx)) {
        nr <- numeric_relations(tab, obs_idx[j], s)
        nr$how <- "doubles with certified gaps, 70-digit band"
        rels[[j]] <- nr
      }
    }
    names(rels) <- names(w$obs)
    out[[cn]][[s]] <- c(tab, list(rel = rels))
  }
  # one 70-digit fit per pretested part in the band, for every contrast and
  # unpretested part it is needed with
  spec_of <- function(iP, needs) {
    bi <- pP[[iP]]$band
    parts <- list()
    for (nd in needs) {
      wt <- CW[[nd$cn]]
      fixed <- list()
      if (wt[["U"]] != 0) fixed <- list(list(w = dyad(wt[["U"]]), R = pU[[nd$iU]]$R, V = pU[[nd$iU]]$V))
      parts[[paste(nd$cn, nd$iU, sep = "#")]] <- list(L = as.list(dyad(bi$item$L$own[bi$piv])), w = dyad(wt[["P"]]),
                                                      fixed = fixed)
    }
    list(X = bi$item$X, y = bi$item$y, f = bi$item$f, e = bi$item$e, keep = bi$keep, piv = bi$piv, b = bi$b,
         parts = parts, skip = bi$skip)
  }
  out <- resolve_band(w, out, obs_idx, nP, theta, spec_of,
                      unit_of = function(i) (i - 1) %% nP + 1,
                      part_of = function(i, cn) paste(cn, (i - 1) %/% nP + 1, sep = "#"))
  list(tables = out, obs_idx = obs_idx, nU = nU, nP = nP)
}

# Resolve every band entry with one call: the classes in the band, and their
# observed classes, are grouped by fitting unit (a class, or a pretested
# part), and each unit is fitted once at 70 digits.
resolve_band <- function(w, out, obs_idx, nP, theta, spec_of, unit_of, part_of) {
  needs <- list()
  for (cn in ALLC) for (s in STATS) {
    rels <- out[[cn]][[s]]$rel
    for (j in seq_along(rels)) {
      bi <- rels[[j]]$band
      if (!length(bi)) next
      for (i in c(bi, obs_idx[j])) {
        u <- as.character(unit_of(i))
        needs[[u]][[part_of(i, cn)]] <- list(cn = cn, iU = (i - 1) %/% nP + 1, i = i)
      }
    }
  }
  if (!length(needs)) return(out)
  units <- names(needs)
  specs <- lapply(units, function(u) spec_of(as.integer(u), needs[[u]]))
  bv <- band_values(w, specs, theta)
  names(bv) <- units
  for (u in units) if (!is.null(bv[[u]]$error)) {
    BAND_LOG$fail[[length(BAND_LOG$fail) + 1]] <- data.frame(world = w$id, unit = u, error = bv[[u]]$error,
                                                             stringsAsFactors = FALSE)
  }
  val50 <- function(i, cn, s) {
    v <- bv[[as.character(unit_of(i))]][[part_of(i, cn)]][[if (s == "studentized") "stud" else "diff"]]
    if (is.null(v)) NA_character_ else v
  }
  for (cn in ALLC) for (s in STATS) {
    tab <- out[[cn]][[s]]
    for (j in seq_along(tab$rel)) {
      bi <- tab$rel[[j]]$band
      if (!length(bi)) next
      o <- obs_idx[j]
      s50 <- vapply(bi, function(i) val50(i, cn, s), "")
      so <- val50(o, cn, s)
      if (is.na(so)) {
        out[[cn]][[s]]$rel[[j]]$rel[bi] <- NA
        out[[cn]][[s]]$rel[[j]]$unexplained <- length(bi)
        out[[cn]][[s]]$rel[[j]]$band_n <- length(bi)
        next
      }
      ok <- !is.na(s50)
      r <- list(rel = rep(NA_real_, length(bi)))
      if (any(ok)) r$rel[ok] <- band_decide(s50[ok], so, tab$val[bi][ok], tab$val[o],
                                           pmax(tab$eps[bi][ok], tab$eps[o]))$rel
      out[[cn]][[s]]$rel[[j]]$rel[bi] <- r$rel
      out[[cn]][[s]]$rel[[j]]$band_n <- length(bi)
      if (any(is.na(r$rel))) out[[cn]][[s]]$rel[[j]]$unexplained <- sum(is.na(r$rel))
    }
  }
  out
}

# 50-digit decision: the order of |T| from the 50-digit values when their
# relative gap is at least 1e-35; otherwise unexplained (NA), and the case is
# excluded (design Section 5, rules 2 and 3). The double value must agree
# with the 50-digit value within 1e4 eps, which checks the measured error.
BAND_LOG <- new.env()
BAND_LOG$rows <- list()
band_decide <- function(s50, s50o, vd, vdo, eps) {
  g <- py_result("decgap", list(pairs = lapply(seq_along(s50), function(q) list(s50[[q]], s50o))))
  rel <- rep(NA_real_, length(s50))
  for (q in seq_along(s50)) {
    gr <- as.numeric(g[[q]]$rel)
    if (gr >= 1e-35) rel[q] <- as.numeric(g[[q]]$sign)
    agree <- abs(vd[q] - abs(as.numeric(s50[[q]]))) <= 1e4 * eps[q] + 1e-300
    BAND_LOG$rows[[length(BAND_LOG$rows) + 1]] <- data.frame(
      double = vd[q], digits50 = s50[[q]], observed50 = s50o, rel_gap50 = gr,
      sign = as.numeric(g[[q]]$sign), eps = eps[q], eps_check = agree, stringsAsFactors = FALSE)
  }
  list(rel = rel, log = character(0))
}

ref_iter_full <- function(w, nc, theta, obs_idx) {
  fam <- w$family
  nU <- length(nc$U$keys)
  nP <- length(nc$P$keys)
  n <- nU * nP
  mU <- w$prep$cond$U$m
  mP <- w$prep$cond$P$m
  items <- vector("list", n)
  info <- vector("list", n)
  k <- ncol(w$prep$Z)
  for (iU in seq_len(nU)) for (iP in seq_len(nP)) {
    i <- (iU - 1) * nP + iP
    cU <- nc$U$comps[iU, ]
    cP <- nc$P$comps[iP, ]
    rows <- class_rows(w, cU, cP)
    ds <- design_of(rows)
    items[[i]] <- list(X = ds$X, y = ds$y, f = ds$f, e = ds$e,
                       L = lapply(stats::setNames(ALLC, ALLC), function(cn) c(cell_L(cn), rep(0, k))))
    arm <- function(cell) {
      sel <- ds$cell == cell
      c(n = sum(ds$f[sel]), S = sum(ds$f[sel] * ds$y[sel]))
    }
    a <- lapply(0:3, arm)
    info[[i]] <- list(a = a, cU = cU, cP = cP, rows = rows)
  }
  ev <- eval_iter_batch(items, fam, theta)
  out <- list()
  for (cn in ALLC) for (s in STATS) {
    wt <- CW[[cn]]
    used <- names(wt)[wt != 0]
    stud <- s == "studentized"
    code <- character(n)
    val <- signed <- eps <- rep(NA_real_, n)
    key <- character(n)
    ill <- logical(n)
    for (i in seq_len(n)) {
      a <- info[[i]]$a
      emp_sig <- function(cond) {
        cells <- if (cond == "U") c(2, 1) else c(4, 3)  # treated, control (1-based: UC=1, UT=2, PC=3, PT=4)
        at <- a[[cells[1]]]
        ac <- a[[cells[2]]]
        emp <- function(x) if (fam == "binomial") x[["S"]] == 0 || x[["S"]] == x[["n"]] else x[["S"]] == 0
        if (at[["n"]] == 0 || ac[["n"]] == 0) return(NA)
        if (emp(at)) return(if (fam == "binomial" && at[["S"]] == at[["n"]]) 1 else -1)
        if (emp(ac)) return(if (fam == "binomial" && ac[["S"]] == ac[["n"]]) -1 else 1)
        0
      }
      cc <- ev[[i]]$c[[cn]]
      if (cc$code %in% c("NE", "ALIAS")) {
        code[i] <- cc$code
        next
      }
      if (stud && isTRUE(cc$h1)) {
        code[i] <- "H1"
        next
      }
      sg <- vapply(used, emp_sig, 0)
      sg[is.na(sg)] <- 0
      if (any(sg != 0)) {
        if (length(used) == 2 && all(sg != 0) && sign(wt[used[1]]) * sg[1] != sign(wt[used[2]]) * sg[2]) {
          code[i] <- "INF-X"
        } else {
          code[i] <- "INF-E"
        }
        next
      }
      if (cc$code == "DIV") {
        code[i] <- "SEP-L"
        next
      }
      if (cc$code == "UNDECIDED") {
        code[i] <- "UNDECIDED"
        next
      }
      lim <- identical(ev[[i]]$status, "limit")
      # in a limiting model the statistic is that of the rows kept (the rows
      # removed are fitted exactly and add nothing), so the key is theirs
      rk <- if (lim) info[[i]]$rows[ev[[i]]$keep] else info[[i]]$rows
      sws <- SWAPS[[cn]]
      if (length(used) == 2 && !stud) {
        # design E3a: when one condition's contrast is exactly 0, the swap of
        # the other condition alone keeps |T| of the difference statistic
        kd0 <- suff_key(rk, fam, "difference")
        if (identical(suff_key(rk, fam, "difference", PERMS$sP), kd0) ||
            identical(suff_key(rk, fam, "difference", PERMS$sU), kd0)) sws <- names(PERMS)
      }
      kv <- vapply(sws, function(sw) suff_key(rk, fam, s, PERMS[[sw]]), "")
      key[i] <- paste0(if (lim) "LIM:", min(kv))
      # exactly zero when a sign-reversing swap leaves the exact key unchanged
      kd <- suff_key(rk, fam, "difference")
      zero <- any(vapply(REVERSE[[cn]], function(sw) identical(suff_key(rk, fam, "difference", PERMS[[sw]]), kd),
                         logical(1)))
      if (zero) {
        code[i] <- "ZERO"
        val[i] <- 0
        signed[i] <- 0
        eps[i] <- 0
        next
      }
      if (!stud) {
        val[i] <- abs(cc$est)
        signed[i] <- cc$est
        eps[i] <- abs(cc$est - cc$est2) + 1e-15 * abs(cc$est)
      } else {
        t1 <- cc$est / sqrt(cc$V)
        t2 <- cc$est2 / sqrt(cc$V2)
        val[i] <- abs(t1)
        signed[i] <- t1
        eps[i] <- abs(abs(t1) - abs(t2)) + 1e-15 * abs(t1)
      }
      ill[i] <- isTRUE(cc$ill)
      code[i] <- if (lim) "LIM" else "finite"
    }
    tab <- list(code = code, val = val, approx = signed, eps = eps, key = key, ill = ill)
    rels <- list()
    for (j in seq_along(obs_idx)) {
      nr <- numeric_relations(tab, obs_idx[j], s)
      nr$how <- "doubles with certified gaps, 70-digit band"
      rels[[j]] <- nr
    }
    names(rels) <- names(w$obs)
    out[[cn]][[s]] <- c(tab, list(rel = rels))
  }
  spec_of <- function(i, needs) {
    it <- items[[i]]
    parts <- list()
    for (nd in needs) {
      cc <- ev[[i]]$c[[nd$cn]]
      parts[[nd$cn]] <- list(L = as.list(dyad(it$L[[nd$cn]][cc$piv])), w = "1")
    }
    cc <- ev[[i]]$c[[needs[[1]]$cn]]
    list(X = it$X, y = it$y, f = it$f, e = it$e, keep = ev[[i]]$keep, piv = cc$piv, b = cc$b, parts = parts,
         skip = cc$skip)
  }
  out <- resolve_band(w, out, obs_idx, nP, theta, spec_of, unit_of = function(i) i, part_of = function(i, cn) cn)
  list(tables = out, obs_idx = obs_idx, nU = nU, nP = nP, status = vapply(ev, `[[`, "", "status"))
}

# ---------------------------------------------------------------------------
# Counts: b and M under the conventions.

exact_sum <- function(wU, inner) {
  # sum of wU[i] * inner[i], exactly
  if (sum(wU) * max(1, max(inner)) < 2^53) return(sprintf("%.0f", sum(wU * inner)))
  r <- py_result("bigsum", list(items = list(list(a = as.list(sprintf("%.0f", wU)),
                                                  b = as.list(sprintf("%.0f", inner))))))
  r[[1]]
}

case_counts <- function(tab, rel, wU, wP, nP, s, rule) {
  code <- tab$code
  n <- length(code)
  # value per class under the rule: count (1), not count (0), drop
  cnt <- ifelse(is.na(rel), NA_real_, as.numeric(rel >= 0))
  keep <- rep(TRUE, n)
  conv <- code %in% c("NE", "H1") | (code == "SEP-L" & s == "studentized")
  if (rule == "drop") {
    # labelings without a defined statistic (H1, NE, SEP-L studentized) and
    # 0/0 (ZZ) are left out of both counts
    keep[conv] <- FALSE
    keep[code == "ZZ"] <- FALSE
  } else if (rule == "inf") {
    # H1 and SEP-L counted as +Inf; NE left out
    cnt[code == "H1" | (code == "SEP-L" & s == "studentized")] <- 1
    keep[code == "NE"] <- FALSE
  } else {
    # proposal: SEP-L counted as +Inf (relations already say so); H1 and NE
    # are errors
    keep[code %in% c("NE", "H1")] <- FALSE
  }
  # ALIAS (a covariate collinear with the cells under the labeling; not in
  # the design's table, a decision for the lead): left out under every rule
  # here, and counted in n_alias so that the error alternative can be read
  keep[code == "ALIAS"] <- FALSE
  mat_w <- function(sel) {
    m <- matrix(as.numeric(sel), ncol = nP, byrow = TRUE)
    inner <- drop(m %*% wP)
    exact_sum(wU, inner)
  }
  wpos <- outer(wU, wP) > 0
  wposv <- as.vector(t(wpos))
  err <- NULL
  if (rule == "proposal") {
    if (any(code[wposv] == "H1")) err <- "H1"
    if (any(code[wposv] == "NE")) err <- "NE"
  }
  und <- any(code[wposv] == "UNDECIDED") || any(is.na(cnt[keep & wposv]))
  if (rule == "drop" && any(code[wposv] %in% c("NE", "H1"))) err <- NULL
  b <- mat_w(keep & !is.na(cnt) & cnt == 1)
  M <- mat_w(keep)
  list(b = b, M = M, error = err, undecided = und,
       b_tied = mat_w(keep & !is.na(rel) & rel == 0), b_strict = mat_w(keep & !is.na(rel) & rel == 1),
       n_inf = mat_w(keep & code %in% INF_CODES), n_sepl = mat_w(code == "SEP-L"),
       n_h1 = mat_w(code == "H1"), n_ne = mat_w(code == "NE"), n_zz = mat_w(code == "ZZ"),
       n_alias = mat_w(code == "ALIAS"),
       n_lim = mat_w(code == "LIM"))
}
