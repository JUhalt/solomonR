# Classes of labelings of perm_solomon() (issue #134).
#
# A class is the set of labelings that give each arm the same rows, where two
# rows count as the same when everything the model sees in them is equal, bit
# for bit (the posted reference set, Section 3). Each class is computed once,
# from a representative labeling built from the class alone, so labelings of
# one class get the same double. Swapping the arms of a condition is a
# reparameterization that reverses or keeps the contrast (the set's table of
# swaps), so swapped classes share one computation with the sign kept.
# Grouping labelings by what the statistic depends on is the device of the
# network algorithms for exact tests (Mehta & Patel, 1983).

# Row types of the rows in the model, by pretest condition. Two rows are of
# one type when everything the model sees in them is equal, bit for bit
# (sprintf("%a")): the outcome as glm() fits it (`y`, from .perm_outcome()),
# the covariate columns of the design `X` (every column but the intercept and
# the treatment and pretest terms), the offset `off` (NULL without an
# exposure), and the pretest condition `pretested`. `separable` says that no
# covariate is shared by the two conditions, so that a contrast on one
# condition depends only on that condition's labeling.
.perm_class_setup <- function(X, y, off, pretested, contrast, separable) {
  Z <- X[, !colnames(X) %in% c("(Intercept)", "treat", "pretested", "treat:pretested"), drop = FALSE]
  parts <- list(as.matrix(y), Z, off, pretested)
  parts <- parts[!vapply(parts, is.null, logical(1))]
  M <- do.call(cbind, lapply(parts, function(p) as.matrix(p) + 0))
  key <- do.call(paste, c(lapply(seq_len(ncol(M)), function(j) sprintf("%a", M[, j])), sep = "|"))
  per <- lapply(c(U = 0, P = 1), function(p) {
    idx <- which(pretested == p)
    type <- match(key[idx], unique(key[idx]))
    # each row's position among the rows of its type, in row order, so that a
    # class's representative labeling (the first rows of each type treated)
    # is built in one step
    rank <- stats::ave(seq_along(type), type, FUN = seq_along)
    m <- tabulate(type, nbins = max(c(0L, type)))
    list(idx = idx, type = type, m = m, rank = rank, single = all(m == 1L))
  })
  w_c <- .perm_contrast_weights(contrast)  # (pretested, unpretested)
  one <- separable && sum(w_c != 0) == 1L
  list(per = per,
       # conditions whose labeling the statistic depends on
       conds = if (one) c("U", "P")[c(w_c[2] != 0, w_c[1] != 0)] else c("U", "P"),
       # how swaps act: "each" (T|U, T|P: each condition's pair unordered, the
       # sign from the contrast's own condition) or "joint" (PxT, ATE: both
       # conditions swapped together)
       swap = if (contrast %in% c("Treatment | unpretested", "Treatment | pretested")) "each" else "joint",
       sign_cond = if (identical(contrast, "Treatment | unpretested")) "U" else
         if (identical(contrast, "Treatment | pretested")) "P" else NA_character_)
}

# A composition's key: single-row types as packed bits, the others as counts.
.perm_comp_key <- function(v, m) {
  s <- m == 1L
  bits <- if (any(s)) {
    b <- v[s] == 1L
    paste(as.character(packBits(c(b, logical((8L - length(b) %% 8L) %% 8L)))), collapse = "")
  } else ""
  paste0(bits, "/", paste(v[!s], collapse = ","))
}

# Orientation of a composition against its swap m - v: 1 when v is the
# lexicographically larger, -1 when smaller, 0 when equal.
.perm_orient <- function(v, m) {
  dlt <- v - (m - v)
  k <- which(dlt != 0)
  if (!length(k)) 0L else if (dlt[k[1]] > 0) 1L else -1L
}

# The composition of a condition: how many rows of each type are treated.
# When every row is its own type, it is the labeling itself.
.perm_comp_of <- function(p, tr) {
  if (p$single) as.integer(tr[p$idx] == 1) else tabulate(p$type[tr[p$idx] == 1], nbins = length(p$m))
}

# The class of a labeling `tr` (treatment of the rows in the model): its key,
# the canonical compositions, the sign of the labeling's statistic relative to
# the canonical one, and whether the swap maps the class to itself (contrast
# exactly 0).
.perm_class_of <- function(cc, tr) {
  pU <- cc$per$U
  pP <- cc$per$P
  comp <- list(U = .perm_comp_of(pU, tr), P = .perm_comp_of(pP, tr))
  orient <- c(U = .perm_orient(comp$U, pU$m), P = .perm_orient(comp$P, pP$m))
  conds <- cc$conds
  if (identical(cc$swap, "joint")) {
    # the first condition (in U, P order) that is not its own swap decides
    o <- orient[conds]
    lead <- o[o != 0]
    sgn <- if (length(lead)) lead[[1]] else 0L
    flip <- c(U = sgn < 0, P = sgn < 0)
  } else {
    flip <- orient < 0
    sgn <- orient[[cc$sign_cond]]
  }
  self <- sgn == 0L
  canon <- list(U = if (flip[["U"]]) pU$m - comp$U else comp$U,
                P = if (flip[["P"]]) pP$m - comp$P else comp$P)
  key <- if (length(conds) == 2L) {
    paste0(.perm_comp_key(canon$U, pU$m), "#", .perm_comp_key(canon$P, pP$m))
  } else {
    .perm_comp_key(canon[[conds]], cc$per[[conds]]$m)
  }
  list(key = key, canon = canon, sign = if (self) 1L else as.integer(sgn), self = self)
}

# The canonical compositions of a class, decoded from its key, so that the
# cache need not hold them (memory stays one short string per class).
.perm_key_canon <- function(cc, key) {
  parts <- strsplit(key, "#", fixed = TRUE)[[1]]
  canon <- list(U = NULL, P = NULL)
  for (j in seq_along(cc$conds)) {
    k <- cc$conds[j]
    m <- cc$per[[k]]$m
    s <- m == 1L
    bc <- strsplit(parts[j], "/", fixed = TRUE)[[1]]
    v <- integer(length(m))
    if (any(s)) {
      hex <- bc[1]
      bytes <- as.raw(strtoi(substring(hex, seq(1, nchar(hex), 2), seq(2, nchar(hex), 2)), 16L))
      v[s] <- as.integer(rawToBits(bytes))[seq_len(sum(s))]
    }
    if (any(!s)) v[!s] <- as.integer(strsplit(if (length(bc) > 1) bc[2] else "", ",", fixed = TRUE)[[1]])
    canon[[k]] <- v
  }
  canon
}

# The representative labeling of a class: within each row type, its first
# rows (in row order) treated. Conditions the statistic does not depend on
# keep the labeling `base`.
.perm_rep_labeling <- function(cc, canon, base) {
  tr <- base
  for (k in cc$conds) {
    p <- cc$per[[k]]
    tr[p$idx] <- as.numeric(p$rank <= canon[[k]][p$type])
  }
  tr
}
