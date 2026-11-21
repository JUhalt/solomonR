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

# Which rows of the omnibus table `omnibus` are AHT tests on HC2 that have
# no statistic because L V L' is singular to working precision
# (.aht_singular()). Unlike a test that
# is undefined at leverage 1 (issue #173), whose `df2` is NA too, such a
# test keeps its `df2`, which depends on the design alone.
.omnibus_singular <- function(omnibus) {
  if (is.null(omnibus$method)) return(rep(FALSE, NROW(omnibus)))
  omnibus$method == .omnibus_methods[["aht_hc2"]] & is.na(omnibus$statistic) &
    !is.na(omnibus$df2)
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
# Q is undefined when L V L' is singular, that is, when the residuals that
# enter it do not span the q constraints (.aht_singular()). Q, the
# statistic, and the p-value are then NA, and `singular` is TRUE; eta and
# `df2`, which depend on the design alone, are kept. clubSandwich does not
# test for this: on such a matrix it stops, as not positive definite, or
# gives a value, depending on the rounding error of the residuals.
#
# Returns `statistic` (F), `df1` (q), `df2` (eta - q + 1), and `p.value`,
# with `eta`, `Q`, and `singular`. As in clubSandwich, the factor
# (eta - q + 1) / eta is taken as 0 when it is negative, and the p-value is
# NA when `df2` is not positive. When the constraints are linearly
# dependent, every value but `df1` is NA.
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
                eta = NA_real_, Q = NA_real_, singular = FALSE))
  }
  leverage <- hat > .hc3_leverage_limit
  A <- (Qx %*% qr.Q(qb)) / sqrt(ifelse(leverage, 1, 1 - hat))
  A[leverage, ] <- 0
  e <- (sqrt(w) * model$residuals)[good]
  tb <- backsolve(qr.R(qb), drop(L %*% cf)[qb$pivot], transpose = TRUE)
  M <- crossprod(A * e)
  singular <- .aht_singular(M, (sqrt(w) * (model$fitted.values + model$residuals))[good])
  Q <- if (singular) NA_real_ else sum(tb * solve(M, tb))

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
    eta = eta, Q = Q, singular = singular
  )
}

# Whether L V L', held as M = sum_i a_i a_i' e_i^2 (.aht_test()), is
# singular to working precision. `outcome` holds the outcome of the N
# participants, as the fitted values plus the residuals, times the root of
# their weights.
#
# The matrix is singular when the residuals that enter it do not span the
# q constraints. Without pretest scores or covariates, the residuals of a
# group are 0 when its outcome does not vary. Condition | pretested and
# Condition | unpretested are then singular when two of the groups that
# they compare have residuals of 0, and the other two tests when the
# pretested and the unpretested groups of two conditions do: the
# contrasts of the other groups still span the constraints otherwise.
# With pretest scores or covariates, the residuals of a group are 0 when
# its outcome is an exact linear function of them, with the slopes of the
# fit: when every pretested participant's score changes by the same
# amount, or none changes, Condition | pretested is singular. The
# unpretested groups have no pretest scores.
#
# Computed residuals are 0 only up to rounding error, so the matrix is
# taken as singular when either
# - its reciprocal condition number is below the machine epsilon, the rule
#   by which base R's solve() refuses it; or
# - its smallest eigenvalue is at most tau^2, with tau = 10 N eps times the
#   root mean square of `outcome`.
# In a direction v of unit length, v' M v is a weighted mean of
# e_i^2 / (1 - h_i), with weights (g_i' v)^2 that sum to 1, as the g_i of
# .aht_test() sum to the identity in their outer products. So the second
# rule asks whether the residuals that carry some direction of the
# contrasts are 0 up to rounding error. The first misses that when the
# residuals that enter the test are rounding error in every direction:
# when the outcome does not vary within any of the groups that a test of
# group means compares, as when no unpretested participant of a linear
# probability model has an event, M is made of rounding error alone, is
# not singular by its condition number, and gave an arbitrary statistic. It
# misses it too for two groups with residuals of 0 among groups that vary,
# once the outcome lies so far from zero (10^8 with a spread of 0.5) that
# the rounding error of those residuals is no longer negligible against
# the others.
#
# The rounding error of a residual is of the order of eps times the size
# of the outcome values that the decomposition combines to compute it,
# those of the other groups included, and it grows with the size of the
# fit. So tau is set by the outcome, not by the fitted values: with the
# groups that do not vary at 0 and the other groups' means near 0 against
# their spread, the fitted values are near 0 but the rounding error is
# not, and a bound on the fitted values missed such tests. Of 9,478 tests
# (two to five treatments, 12 to 14,740 participants, binary, integer,
# and normal outcomes 0 to 10^10 from zero, with and without the
# pretest), the 2,673 whose residuals are 0 in exact arithmetic were all
# found singular: 499 by the first rule, and the other 2,174 by the
# second, with a root of the smallest eigenvalue of at most 0.41 N eps
# times the root mean square of the outcome, 24 times below tau. Of the
# 6,805 others, the second rule found 30 singular, and the first none.
# Outcomes that vary give residuals of about their spread, which is within
# tau only when the spread of the groups that a test compares is below
# about 10 N eps of the outcome's root mean square: 16 of the 30 tests had
# an outcome 10^8 from zero with a spread of 10^-4 (2,432 to 9,820
# participants), and 14 had groups with a spread of 10^-9 or less against
# 10^3 in the other groups. The eigenvalues of M are found only to within
# about eps times the largest, so the second rule does not replace the
# first.
.aht_singular <- function(M, outcome) {
  if (rcond(M) < .Machine$double.eps) return(TRUE)
  tau <- 10 * length(outcome) * .Machine$double.eps * sqrt(mean(outcome^2))
  min(eigen(M, symmetric = TRUE, only.values = TRUE)$values) <= tau^2
}
