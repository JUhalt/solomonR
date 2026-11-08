# HC3 standard errors when a participant's leverage is 1 (issue #173).
#
# HC3 divides each participant's squared residual by (1 - h)^2, with h the
# participant's hat value, or leverage (MacKinnon & White, 1985). At h = 1
# the model fits that participant exactly, so the term is 0/0, and the HC3
# standard error of every estimate that depends on that participant is
# undefined. This happens when the pretest varies within the pretested
# groups only through one participant, or when a covariate varies only
# through one participant, such as a factor level that participant alone
# holds. A group with a single participant also has h = 1; for it the fits
# refuse HC3 (.stop_hc3_singleton()).
#
# The rule. A participant is at leverage 1 when the hat value is above
# 1 - 1e-8 (.hc3_leverage_limit). An estimate is a combination l'b of the
# coefficients, and participant i's weight in it is a_i = q_i' R^{-T} l,
# with W^(1/2) X = QR the fit's decomposition (.glm_influence()). When
# |a_i| is at most 1e-8 times the largest |a_j| over the participants
# (.hc3_weight_tolerance), the estimate does not depend on participant i,
# whose term takes its limit, 0. Otherwise the estimate's HC3 standard
# error is undefined, and its standard error, statistic, p-value, interval,
# and Wald R-squared are reported as NA. sandwich::vcovHC() returns NaN
# throughout in that case, so the HC3 covariance then comes from the QR
# decomposition (.hc3_vcov()). The rule is that of the reference set posted
# on issue #134 (H1).
.hc3_leverage_limit <- 1 - 1e-8
.hc3_weight_tolerance <- 1e-8

# The influences of the participants of the glm fit `model` on its
# coefficients, or on the combinations `map` (.glm_influence()), each
# divided by 1 - h_i as HC3 divides them, with the rows of the participants
# at leverage 1 set to their limit, 0. Returns them as `G`, so that the HC3
# covariance is G'G; `undefined`, which flags the columns of `G` whose HC3
# variance is undefined; `leverage`, which flags the rows at leverage 1; and
# `weight`, the participants' weights in each column (NULL when no row is at
# leverage 1).
.hc3_influence <- function(model, map = NULL) {
  infl <- .glm_influence(model, map)
  G <- infl$influence / (1 - infl$hat)
  leverage <- infl$hat > .hc3_leverage_limit
  undefined <- rep(FALSE, ncol(G))
  if (any(leverage)) {
    G[leverage, ] <- 0
    undefined <- .hc3_depends(infl$weight, leverage)
  }
  list(G = G, undefined = undefined, leverage = leverage, weight = infl$weight)
}

# For the estimates whose weights are the columns of `weight` (one row per
# participant), whether each participant at leverage 1 (the rows flagged by
# `leverage`) has a weight in each estimate above 1e-8 times the largest
# weight in it: a matrix with a row for each such participant and a column
# for each estimate.
.hc3_weighted <- function(weight, leverage) {
  weight <- abs(as.matrix(weight))
  largest <- apply(weight, 2, max)
  weight[leverage, , drop = FALSE] >
    .hc3_weight_tolerance * rep(largest, each = sum(leverage))
}

# Whether the HC3 variance of each estimate (column of `weight`) is
# undefined: a participant at leverage 1 has weight in it.
.hc3_depends <- function(weight, leverage) {
  unname(colSums(.hc3_weighted(weight, leverage)) > 0)
}

# The entries of the HC3 covariance of the estimates (columns of `weight`)
# that are undefined: those of two estimates, or of one with itself, in
# which the same participant at leverage 1 has weight.
.hc3_undefined_entries <- function(weight, leverage) {
  big <- .hc3_weighted(weight, leverage)
  crossprod(big + 0) > 0
}

# The participants at leverage 1 of the model `fit` of fit_solomon_glm(),
# from which its HC3 standard errors come: NULL when there are none, and
# otherwise `leverage` (the rows of the model at leverage 1), `weight` (the
# participants' weights in its coefficients), and `rows` (their row numbers
# in the data as given; `in_model` flags the rows the model used).
.hc3_leverage <- function(fit, in_model) {
  infl <- .glm_influence(fit)
  leverage <- infl$hat > .hc3_leverage_limit
  if (!any(leverage)) return(NULL)
  list(leverage = leverage, weight = infl$weight, rows = which(in_model)[leverage])
}

# Whether the HC3 standard error of each combination (each row of `L`, or
# `L` itself as a vector) of the coefficients of the fit that `hc3`
# (.hc3_leverage()) describes is undefined. FALSE for each when `hc3` is
# NULL.
.hc3_undefined <- function(hc3, L) {
  n <- if (is.matrix(L)) nrow(L) else 1L
  if (is.null(hc3)) return(rep(FALSE, n))
  .hc3_depends(hc3$weight %*% t(matrix(L, nrow = n)), hc3$leverage)
}

# The covariance `V` of the coefficients of the fit that `hc3` describes
# (.hc3_leverage()), with NA in its undefined entries: the HC3 covariance
# that a fit returns. `V` itself when `hc3` is NULL.
.hc3_mark <- function(V, hc3) {
  if (is.null(hc3)) return(V)
  V[.hc3_undefined_entries(hc3$weight, hc3$leverage)] <- NA
  V
}

# The classed warning of a fit with participants at leverage 1 under HC3,
# once for each fit. `rows` are their row numbers in the data as given;
# `pretest` and `covariates` say whether the model has a pretest and
# covariates, which name the likely cause.
.warn_hc3_leverage <- function(rows, pretest, covariates) {
  shown <- if (length(rows) > 10L) {
    c(rows[1:10], sprintf("and %d more", length(rows) - 10L))
  } else {
    rows
  }
  one <- length(rows) == 1L
  causes <- c(
    if (pretest) {
      paste0("the pretest varies within the pretested groups only through ",
             if (one) "this participant" else "them")
    },
    if (covariates) {
      paste0("a covariate varies only through ", if (one) "this participant" else "them",
             ", such as a factor level that ", if (one) "the participant holds alone" else
               "they alone hold")
    }
  )
  warning(structure(
    class = c("solomonR_hc3_leverage_warning", "warning", "condition"),
    list(message = paste0(
      "HC3 standard errors are undefined for the estimates that depend on ",
      if (one) "participant " else "participants ", paste(shown, collapse = ", "),
      if (one) " (row number of the data)" else " (row numbers of the data)",
      ", whose leverage (hat value) is 1: the model fits ",
      if (one) "this participant" else "each of them", " exactly. ",
      if (length(causes)) paste0("Likely cause: ", paste(causes, collapse = ", or "), ". "),
      "Those estimates have NA standard errors, tests, and intervals; the others ",
      "keep theirs. Use `robust = \"none\"` for model-based standard errors, which ",
      "are defined."
    ), call = NULL)
  ))
}

# "row 3 of the data", or "among rows 1, 25 of the data", from the row
# numbers `rows` of the participants at leverage 1, for the messages on a
# contrast that depends on one of them.
.hc3_rows_phrase <- function(rows) {
  paste0(if (length(rows) == 1L) "row " else "among rows ", paste(rows, collapse = ", "),
         " of the data")
}

# The classed refusal of HC3 when a group has a single analyzable
# participant, whose hat value is 1, so that HC3 is undefined for most
# contrasts. `groups` names the groups (.group_names()); with several
# treatments the message suggests leaving the condition out.
.stop_hc3_singleton <- function(groups, ngroup = FALSE) {
  stop(structure(
    class = c("solomonR_hc3_singleton_error", "error", "condition"),
    list(message = paste0(
      "HC3 standard errors are undefined because these groups have a single ",
      "analyzable participant: ", paste(groups, collapse = "; "), ". ",
      if (ngroup) "Use `robust = \"none\"` or leave the condition out; " else
        "Use `robust = \"none\"`; ",
      "see validate_solomon()."
    ), call = NULL)
  ))
}
