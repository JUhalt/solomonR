# The small-sample omnibus test of a design with several treatments (issue
# #45, amendment 2).
#
# With `robust = "HC3"`, a Gaussian family, and the identity link,
# fit_solomon_glm() tests each omnibus hypothesis L b = 0, with one row of L
# for each treatment against the control (q rows), by the approximate
# Hotelling's T-squared (AHT) test of Pustejovsky and Tipton (2018), with
# each participant as a cluster and the identity working model: independent
# errors with one variance. The CR2 covariance of clusters of one
# participant is the HC2 covariance of MacKinnon and White (1985), and the
# test is
#   Q = (L b)' (L V L')^{-1} (L b), with L V L' = sum_i c_i c_i' e_i^2 / (1 - h_i),
#   F = Q (eta - q + 1) / (eta q), referred to F(q, eta - q + 1),
#   eta = q (q + 1) / (2 sum_i sum_j (a_i' a_j)^2 m_ij^2).
# Here X is the design matrix, H = X (X'X)^{-1} X' with diagonal h_i, e holds
# the residuals, c_i = L (X'X)^{-1} x_i, G = L (X'X)^{-1} L',
# a_i = G^{-1/2} c_i / sqrt(1 - h_i), and m_ij are the entries of I - H.
# These are Equations 10, 12, and 13 of Pustejovsky and Tipton (2018; the
# numbers are those of arXiv:1601.01981v2) for clusters of one participant
# and ordinary least squares. eta depends on the design matrix alone, not on
# the outcomes. With one constraint, the test is the t test with the degrees
# of freedom of Bell and McCaffrey (2002).
#
# The former test of these fits referred the HC3 Wald statistic over q to
# F(q, N - p). In the package's first simulation study it rejected a true
# null hypothesis too often for Condition | pretested and Condition |
# unpretested with three treatments and 10 participants per group. Every
# other family and link keeps that test, `robust = "none"` the model-based
# one, and `robust = "CR2"` the AHT test of clubSandwich (.cr2_wald()).

# The label of each omnibus test in the `method` column of `omnibus`.
.omnibus_methods <- c(
  aht_hc2 = "AHT, HC2 (Pustejovsky & Tipton, 2018)",
  aht_cr2 = "AHT, CR2 (Pustejovsky & Tipton, 2018)",
  wald_f_hc3 = "Wald F, HC3",
  wald_chisq_hc3 = "Wald chi-square, HC3",
  f_model = "F, model-based",
  chisq_model = "chi-square, model-based"
)

# Whether the omnibus tests of an HC3 fit of this family are AHT tests: a
# Gaussian family with the identity link.
.aht_family <- function(family) {
  identical(family$family, "gaussian") && identical(family$link, "identity")
}

# Whether every test of the omnibus table `omnibus` is an AHT test, and
# whether every one is the AHT test on HC2, by the `method` column. A fit
# made before the column was added has none, and holds the former tests, so
# both are FALSE for it: its print and its report describe what it holds.
.omnibus_aht <- function(omnibus) {
  !is.null(omnibus$method) && all(omnibus$method %in% .omnibus_methods[c("aht_hc2", "aht_cr2")])
}
.omnibus_aht_hc2 <- function(omnibus) {
  !is.null(omnibus$method) && all(omnibus$method == .omnibus_methods[["aht_hc2"]])
}

# The AHT test of the constraints L b = 0 (one per row of `L`) on the
# coefficients of the least-squares fit `model`, a glm with a Gaussian
# family and the identity link, with each participant as a cluster.
#
# It is computed from the QR decomposition X = QR that the fit holds, as the
# HC3 covariance is (.glm_influence()), without the inverse of X'X and
# without an N x N matrix. With q_i the ith row of Q, x_i = R' q_i, so that
# the hat values are h_ij = q_i' q_j and c_i = B' q_i with B = R^{-T} L',
# which backsolve() gives. With B = Q_B R_B a second decomposition,
# G = B'B = R_B' R_B, and
#   g_i = Q_B' q_i = R_B^{-T} c_i
# has g_i' g_j = c_i' G^{-1} c_j, which is all that the test needs of
# G^{-1/2}: a_i' a_j = g_i' g_j / sqrt((1 - h_i)(1 - h_j)). Then, with
# a_i = g_i / sqrt(1 - h_i) and t = R_B^{-T} L b,
#   Q = t' (sum_i a_i a_i' e_i^2)^{-1} t,
# and, because m_ij is -h_ij off the diagonal and 1 - h_i on it,
#   sum_i sum_j (a_i' a_j)^2 m_ij^2
#     = sum_i (a_i' a_i)^2 (1 - 2 h_i) + sum_i sum_j (a_i' a_j)^2 (q_i' q_j)^2.
# The last sum is the sum, over the pairs (s, u) of constraints, of the
# squared entries of the p x p matrix sum_i a_is a_iu q_i q_i'. The cost is
# of the order N p^2 q^2, and no step grows with N^2: clubSandwich's
# Wald_test() builds an array of q x q x N x N numbers for the same test.
#
# A participant at leverage 1 (a hat value above .hc3_leverage_limit) adds
# the limit of its terms, 0, as in the HC3 covariance (.hc3_influence()):
# the test is then that of the other participants, and it is defined only
# when no constraint depends on that participant, which the caller checks
# (.hc3_undefined(); issue #173).
#
# Returns `statistic` (F), `df1` (q), `df2` (eta - q + 1), and `p.value`,
# with `eta` and `Q`. As in clubSandwich, the factor (eta - q + 1) / eta is
# taken as 0 when it is negative, and the p-value is NA when `df2` is not
# positive. When the constraints are linearly dependent, every value but
# `df1` is NA.
.aht_test <- function(model, L) {
  cf <- stats::coef(model)
  L <- matrix(L, ncol = length(cf))
  q <- nrow(L)
  w <- model$weights
  good <- w > 0
  qx <- model$qr
  if (is.null(qx) || nrow(qx$qr) != sum(good)) {
    X <- stats::model.matrix(model)[good, , drop = FALSE]
    qx <- qr(sqrt(w[good]) * X, tol = if (is.null(qx$tol)) 1e-11 else qx$tol)
  }
  rank <- seq_len(qx$rank)
  pivot <- qx$pivot[rank]
  Qx <- qr.Q(qx)[, rank, drop = FALSE]
  hat <- rowSums(Qx^2)
  B <- backsolve(qr.R(qx)[rank, rank, drop = FALSE], t(L[, pivot, drop = FALSE]),
                 transpose = TRUE)
  qb <- qr(B)
  if (qb$rank < q) {
    return(list(statistic = NA_real_, df1 = q, df2 = NA_real_, p.value = NA_real_,
                eta = NA_real_, Q = NA_real_))
  }
  leverage <- hat > .hc3_leverage_limit
  A <- (Qx %*% qr.Q(qb)) / sqrt(ifelse(leverage, 1, 1 - hat))
  A[leverage, ] <- 0
  e <- (sqrt(w) * model$residuals)[good]
  tb <- backsolve(qr.R(qb), drop(L %*% cf)[qb$pivot], transpose = TRUE)
  Q <- sum(tb * solve(crossprod(A * e), tb))

  total <- sum(rowSums(A^2)^2 * (1 - 2 * hat))
  for (s in seq_len(q)) {
    for (u in s:q) {
      block <- crossprod(Qx, Qx * (A[, s] * A[, u]))
      total <- total + (if (s == u) 1 else 2) * sum(block^2)
    }
  }
  eta <- q * (q + 1) / (2 * total)
  df2 <- eta - q + 1
  statistic <- max(df2 / eta, 0) * Q / q
  list(
    statistic = statistic, df1 = q, df2 = df2,
    p.value = if (df2 > 0) stats::pf(statistic, q, df2, lower.tail = FALSE) else NA_real_,
    eta = eta, Q = Q
  )
}
