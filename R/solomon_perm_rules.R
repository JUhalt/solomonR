# Rules of perm_solomon() for supported models, degenerate labelings, and
# the statistic at the exact estimate (issue #134). The conventions are
# those of the reference set posted on #134 (Section 5, decisions D1-D16).

# Classed errors of perm_solomon(). Every one also has the class
# solomonR_perm_error, so one handler can catch them all.
.perm_stop <- function(class, ...) {
  stop(structure(
    class = c(class, "solomonR_perm_error", "error", "condition"),
    list(message = paste0(...), call = NULL)
  ))
}

# The four models the package documents and validates (decision D1): the
# Gaussian (identity), binomial (logit), and Poisson (log) models of glm(),
# and the NB2 model of family = "negative_binomial" (MASS::glm.nb(), theta
# by maximum likelihood). Checked before the participant and cluster paths
# split, so that both refuse the same fits.
.perm_check_family <- function(fit) {
  fam <- stats::family(fit$model)
  ok <- (identical(fam$family, "gaussian") && identical(fam$link, "identity")) ||
    (identical(fam$family, "binomial") && identical(fam$link, "logit")) ||
    (identical(fam$family, "poisson") && identical(fam$link, "log")) ||
    (inherits(fit$model, "negbin") && identical(fam$link, "log"))
  if (ok) return(invisible(NULL))
  hint <- if (fam$family %in% c("quasipoisson", "quasibinomial")) {
    sprintf(" Refit with family = %s(), which gives the same coefficients and HC3 standard errors.",
            sub("^quasi", "", fam$family))
  } else if (.is_negbin(fam)) {
    " For a negative binomial model, fit with family = \"negative_binomial\" (theta estimated by MASS::glm.nb())."
  } else {
    ""
  }
  .perm_stop(
    "solomonR_perm_family_error",
    "perm_solomon() supports four models: Gaussian (identity link), binomial (logit link), ",
    "Poisson (log link), and family = \"negative_binomial\". This fit uses ", fam$family,
    " with the ", fam$link, " link, for which the permutation p-value has not been checked ",
    "(issue #134).", hint
  )
}

# The outcome as glm() fits it: a factor counts its first level as 0 and the
# others as 1, and a logical outcome is 0 or 1 (binomial()$initialize), so
# that a binary outcome meets the same rules however it is coded.
.perm_outcome <- function(y) {
  if (is.factor(y)) return(as.numeric(y != levels(y)[1L]))
  as.numeric(y)
}

# No variation in a pretest condition the contrast uses (the acceptance
# criterion of #134 and decision D5). `y` (as .perm_outcome() gives it),
# `pretested`, `shared` (the model-matrix columns of the covariates the two
# conditions share, NULL or none without them), and `offset` (NULL without
# an exposure) are those of the rows in the model.
#
# D5: an error when the contrast is 0 or undefined under every labeling
# (ERR). With the outcome constant in a condition the contrast uses, that
# holds when
# - the constant is at a bound (0 or 1 for a binary outcome, 0 for a count)
#   in a condition the contrast uses: the contrast is undefined;
# - every condition the contrast uses is fitted exactly whatever the labels,
#   with a treatment effect of 0: its outcome and offset are constant, and no
#   shared covariate links it to the other condition (there is none, or each
#   is constant within the condition, where its own intercept absorbs it);
# - the outcome and the offset are constant in both conditions: the model
#   without treatment fits every row, so every treatment coefficient is 0.
# Otherwise the contrast still varies, through the other condition, a shared
# covariate, or the exposure (ERR-S): the p-value, with a classed warning
# (D5's proposal).
.perm_check_variation <- function(y, pretested, contrast, family, shared = NULL, offset = NULL) {
  uses <- .perm_contrast_weights(contrast) != 0
  rows <- list(pretested = pretested == 1, unpretested = pretested == 0)
  one_value <- function(v) length(unique(v)) == 1L
  constant <- vapply(rows, function(r) one_value(y[r]), logical(1))
  if (!any(constant[uses])) return(invisible(NULL))
  flat <- constant & vapply(rows, function(r) is.null(offset) || one_value(offset[r]), logical(1))
  apart <- vapply(rows, function(r) {
    is.null(shared) || ncol(shared) == 0L ||
      all(vapply(seq_len(ncol(shared)), function(j) one_value(shared[r, j]), logical(1)))
  }, logical(1))
  at_bound <- constant & vapply(rows, function(r) {
    v <- y[r][1L]
    switch(family$family, gaussian = FALSE, binomial = v %in% c(0, 1), isTRUE(v == 0))
  }, logical(1))
  and_list <- function(x) {
    if (length(x) <= 2L) return(paste(x, collapse = " and "))
    paste0(paste(x[-length(x)], collapse = ", "), ", and ", x[length(x)])
  }
  named <- function(keep) {
    x <- names(rows)[keep & uses]
    paste0(and_list(x), " condition", if (length(x) > 1L) "s" else "")
  }
  undefined <- any(at_bound & uses)
  if (undefined || all((flat & apart)[uses]) || all(flat)) {
    .perm_stop(
      "solomonR_perm_no_variation_error",
      if (all(constant) && !undefined) {
        paste0("The outcome does not vary in either pretest condition, so the contrast \"",
               contrast, "\" is 0 under every labeling and the test has no information.")
      } else {
        paste0("The outcome does not vary in the ", named(if (undefined) at_bound else constant),
               ", which the contrast \"", contrast, "\" uses, so the contrast is ",
               if (undefined) "undefined" else "0",
               " under every labeling and the test has no information.",
               if (!all(constant)) " Choose a contrast on a condition whose outcome varies.")
      }
    )
  }
  # The contrast still varies (decision D5's proposal): a p-value, with a
  # classed warning that names the constant condition and what the contrast
  # varies through.
  through <- c(
    if (any(!constant[uses])) "the other condition",
    if (any((constant & !apart)[uses])) "the covariates the conditions share",
    if (any((constant & !flat)[uses])) "the exposure"
  )
  warning(structure(
    class = c("solomonR_perm_no_variation_warning", "warning", "condition"),
    list(message = paste0(
      "The outcome does not vary in the ", named(constant), "; the contrast \"", contrast,
      "\" varies only through ", and_list(through), "."), call = NULL)
  ))
  invisible(NULL)
}
