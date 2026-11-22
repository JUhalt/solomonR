# The small-sample omnibus tests of designs with several treatments (issue
# #45, amendment 2): the reference computations, by clubSandwich and from
# the test's equations, and the lines and sentences that the print and the
# report give for an omnibus table. test-ngroup-aht.R uses them, and so do
# the two tests that compare new fits with results stored from master
# (test-covariate-location.R and test-pretest-center.R), where only these
# tests' numbers and wording may differ from master's.

# The approximate Hotelling's T-squared test of Pustejovsky and Tipton
# (2018) for each omnibus hypothesis, from
# clubSandwich::Wald_test(test = "HTZ") on the CR2 covariance with each
# participant as a cluster, of a model fitted by lm() with one indicator for
# each treatment. `y_pre` is centered at 0 among the pretested participants
# in the model and each numeric covariate at its mean there, which changes
# no test and keeps clubSandwich accurate: it inverts X'X, and lost accuracy
# on a covariate far from zero (#137). Returns the tests in the order and
# under the labels of the `omnibus` table. With `engine = "literal"`, the
# tests are computed by aht_literal() instead, which shares no code with
# clubSandwich or the package. `cluster`, if given, groups the participants
# into clusters (CR2), in place of one cluster for each participant. `only`,
# if given, names the tests to compute, as when one of them has no
# statistic.
aht_reference <- function(y, condition, pretested, y_pre = NULL, covariates = NULL, control,
                          cluster = NULL, engine = c("clubSandwich", "literal"), only = NULL) {
  engine <- match.arg(engine)
  treatments <- setdiff(unique(as.character(stats::na.omit(condition))), control)
  d <- data.frame(y = y, pretested = pretested)
  d$cluster_ <- if (is.null(cluster)) seq_along(y) else cluster
  terms <- paste0("t", seq_along(treatments))
  for (j in seq_along(treatments)) d[[terms[j]]] <- as.numeric(condition == treatments[j])
  rhs <- sprintf("(%s) * pretested", paste(terms, collapse = " + "))
  if (!is.null(y_pre)) {
    d$pre0 <- ifelse(pretested == 1, y_pre, 0)
    rhs <- c(rhs, "pre0")
  }
  if (!is.null(covariates)) {
    names(covariates) <- paste0("z", seq_along(covariates))
    d <- cbind(d, covariates)
    rhs <- c(rhs, names(covariates))
  }
  d <- d[stats::complete.cases(d), , drop = FALSE]
  if (!is.null(y_pre)) {
    d$pre0 <- ifelse(d$pretested == 1, d$pre0 - mean(d$pre0[d$pretested == 1]), 0)
  }
  for (name in names(covariates)) {
    if (is.numeric(d[[name]])) d[[name]] <- d[[name]] - mean(d[[name]])
  }
  model <- stats::lm(stats::reformulate(rhs, "y"), data = d)
  V <- clubSandwich::vcovCR(model, cluster = d$cluster_, type = "CR2")
  cn <- names(stats::coef(model))
  constraints <- function(main, interaction) {
    L <- matrix(0, length(terms), length(cn), dimnames = list(NULL, cn))
    for (j in seq_along(terms)) {
      L[j, terms[j]] <- main
      L[j, paste0(terms[j], ":pretested")] <- interaction
    }
    L
  }
  tests <- list(
    "Condition (avg over pretest)" = constraints(1, 0.5),
    "Pretest x Condition" = constraints(0, 1),
    "Condition | pretested" = constraints(1, 1),
    "Condition | unpretested" = constraints(1, 0)
  )
  if (!is.null(only)) tests <- tests[only]
  out <- do.call(rbind, lapply(names(tests), function(test) {
    if (engine == "literal") {
      lt <- aht_literal(model, tests[[test]], d$cluster_)
      return(data.frame(test = test, statistic = lt[["statistic"]], df1 = length(terms),
                        df2 = lt[["df2"]], p.value = lt[["p.value"]], stringsAsFactors = FALSE))
    }
    wt <- as.data.frame(clubSandwich::Wald_test(model, constraints = tests[[test]], vcov = V,
                                                test = "HTZ"))
    data.frame(test = test, statistic = wt$Fstat, df1 = wt$df_num, df2 = wt$df_denom,
               p.value = wt$p_val, stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  out
}

# The AHT test of the constraints L b = 0 on the coefficients of the
# least-squares fit `model`, computed literally from the definitions of
# Pustejovsky and Tipton (2018; Equations 10, 12, and 13 of
# arXiv:1601.01981v2) under the identity working model, with explicit
# inverses and N x N matrices, for small data sets only:
# - the CR2 adjustment of cluster j, A_j = (I - H_jj)^(-1/2), the
#   symmetric inverse square root, and V = (X'X)^-1 (sum_j X_j' A_j e_j
#   e_j' A_j X_j) (X'X)^-1;
# - Q = (L b)' (L V L')^-1 (L b);
# - with Omega = L (X'X)^-1 L', the normalized covariance estimate
#   D = Omega^(-1/2) L V L' Omega^(-1/2) has entries
#   d_st = sum_j (g_sj' e*)(g_tj' e*), where e* is the error vector,
#   g_sj = (I - H)_j' A_j X_j (X'X)^-1 L' Omega^(-1/2) u_s (u_s the sth unit
#   vector), and (I - H)_j the rows of I - H for cluster j. Under the
#   working model (independent errors with one variance, taken as 1),
#   Var(d_st) = sum_i sum_j [(g_si' g_sj)(g_ti' g_tj) + (g_si' g_tj)(g_ti' g_sj)],
#   and eta = q (q + 1) / sum_s sum_t Var(d_st);
# - F = Q (eta - q + 1) / (eta q), on F(q, eta - q + 1).
# With one participant in each cluster, this is the AHT test on HC2 that
# fit_solomon_glm() computes with `robust = "HC3"`; with other clusters, the
# test of `robust = "CR2"`. It forms (X'X)^-1, so it is given covariates
# near zero.
aht_literal <- function(model, L, cluster = seq_len(stats::nobs(model))) {
  X <- stats::model.matrix(model)
  e <- stats::residuals(model)
  n <- nrow(X)
  q <- nrow(L)
  XtXinv <- solve(crossprod(X))
  IH <- diag(n) - X %*% XtXinv %*% t(X)
  inv_sqrt <- function(S) {
    ev <- eigen(S, symmetric = TRUE)
    ev$vectors %*% diag(1 / sqrt(ev$values), length(ev$values)) %*% t(ev$vectors)
  }
  groups <- unique(cluster)
  adj <- matrix(0, n, n)
  meat <- matrix(0, ncol(X), ncol(X))
  for (j in groups) {
    at <- which(cluster == j)
    adj[at, at] <- inv_sqrt(IH[at, at, drop = FALSE])
    u <- t(X[at, , drop = FALSE]) %*% adj[at, at, drop = FALSE] %*% e[at]
    meat <- meat + u %*% t(u)
  }
  V <- XtXinv %*% meat %*% XtXinv
  b <- L %*% stats::coef(model)
  Q <- drop(t(b) %*% solve(L %*% V %*% t(L), b))
  P <- adj %*% X %*% XtXinv %*% t(L) %*% inv_sqrt(L %*% XtXinv %*% t(L))
  # G[[s]] has the vector g_sj in column j.
  G <- lapply(seq_len(q), function(s) {
    vapply(groups, function(j) {
      at <- which(cluster == j)
      drop(crossprod(IH[at, , drop = FALSE], P[at, s]))
    }, numeric(n))
  })
  total <- 0
  for (s in seq_len(q)) {
    for (t in seq_len(q)) {
      st <- crossprod(G[[s]], G[[t]])
      total <- total + sum(crossprod(G[[s]]) * crossprod(G[[t]])) + sum(st * t(st))
    }
  }
  eta <- q * (q + 1) / total
  df2 <- eta - q + 1
  statistic <- Q * df2 / (eta * q)
  c(statistic = statistic, df2 = df2, p.value = stats::pf(statistic, q, df2, lower.tail = FALSE),
    eta = eta, Q = Q)
}

# The columns of an omnibus table that the AHT test changed.
aht_columns <- c("statistic", "df2", "p.value")

# The label of the AHT test on HC2 in the `method` column.
aht_hc2_label <- "AHT, HC2 (Pustejovsky & Tipton, 2018)"

# The heading of the omnibus tests in the print of a fit whose tests are AHT
# tests, and the phrase that names the test in the method of its report in
# place of master's "Omnibus Wald F tests".
aht_heading <- "Omnibus tests (small-sample AHT tests; Pustejovsky & Tipton, 2018)"
aht_method_phrase <- paste(
  "Omnibus Wald tests, using the small-sample test of Pustejovsky and Tipton (2018) with",
  "the HC2 covariance matrix (MacKinnon & White, 1985) and each participant as a cluster,"
)

# The rows of the omnibus table `om` as the print gives them, with runs of
# spaces reduced to one: the test, its F statistic with `df2` to one decimal
# unless it is whole, and the p-value.
aht_print_rows <- function(om) {
  whole <- abs(om$df2 - round(om$df2)) < 1e-8
  df2 <- ifelse(whole, sprintf("%.0f", om$df2), sprintf("%.1f", om$df2))
  p <- ifelse(om$p.value < .001, "<.001", sprintf("%.3f", om$p.value))
  sprintf("%s F(%d, %s) = %.2f %s", om$test, as.integer(om$df1), df2, om$statistic, p)
}

# The lines of the printed fit `text`, a new fit whose omnibus table is
# `om`, put back in the form master printed them in, `master`: the heading
# of the omnibus tests, the table under it, and the note that ends the
# print are the only lines that the AHT test changed. They are checked here,
# and master's lines are returned in their place, so that a caller compares
# every other line with master's as before. Returns NULL when a block is
# not found.
aht_text_as_master <- function(text, master, om, info = NULL) {
  squeeze <- function(x) gsub(" +", " ", trimws(x))
  at <- match(aht_heading, text)
  ref_at <- match("Omnibus tests", master)
  testthat::expect_false(is.na(at), info = info)
  testthat::expect_false(is.na(ref_at), info = info)
  if (is.na(at) || is.na(ref_at)) return(NULL)
  block <- seq_len(nrow(om) + 2L) - 1L
  testthat::expect_identical(squeeze(text[at + 1L]), "Test Statistic p", info = info)
  testthat::expect_identical(squeeze(text[at + 1L + seq_len(nrow(om))]), aht_print_rows(om),
                             info = info)
  text <- append(text[-(at + block)], master[ref_at + block], after = at - 1L)
  # The note that ends the print, after its empty line.
  note <- utils::tail(which(startsWith(text, "Experimental: ")), 1L)
  ref_note <- utils::tail(which(startsWith(master, "Experimental: ")), 1L)
  testthat::expect_length(note, 1L)
  testthat::expect_length(ref_note, 1L)
  if (!length(note) || !length(ref_note)) return(NULL)
  new_note <- text[note:length(text)]
  testthat::expect_identical(new_note, .ngroup_lifecycle_note(list(omnibus = om))[-1L], info = info)
  testthat::expect_identical(master[ref_note:length(master)], .ngroup_lifecycle_note()[-1L],
                             info = info)
  c(text[seq_len(note - 1L)], master[ref_note:length(master)])
}

# The sentence of the report on the two omnibus tests it gives, for the
# omnibus table `om`.
aht_report_sentence <- function(om) {
  gave <- function(test) {
    r <- om[om$test == test, ]
    sprintf("%s, %s", .apa_omnibus(r$statistic, r$df1, r$df2, r$reference, FALSE, 2),
            .apa_p(r$p.value, FALSE))
  }
  sprintf(paste(
    "The omnibus test of the Pretest x Condition interaction (pretest sensitization) gave",
    "%s, and the omnibus test of the conditions, averaged over pretest conditions, gave %s."
  ), gave("Pretest x Condition"), gave("Condition (avg over pretest)"))
}
