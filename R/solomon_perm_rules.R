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

# No variation in a pretest condition the contrast uses (the acceptance
# criterion of #134 and decision D5). An error when the statistic is 0, or
# infinite, for every labeling: the outcome is constant at 0 (or 1 for a
# binary outcome) in such a condition, or constant in the only condition a
# contrast of a separable model uses. When the statistic still varies (the
# ATE or Pretest x Treatment with the other condition varying, or a
# covariate shared by the conditions), the p-value is computed (D5's
# proposal).
.perm_check_variation <- function(y, pretested, contrast, family, separable) {
  w <- .perm_contrast_weights(contrast)
  conds <- c(pretested = 1, unpretested = 0)[w != 0]
  values <- lapply(conds, function(p) unique(y[pretested == p]))
  constant <- vapply(values, function(v) length(v) == 1L, logical(1))
  if (!any(constant)) return(invisible(NULL))
  v <- unlist(values[constant])
  boundary <- if (identical(family$family, "binomial")) any(v %in% c(0, 1)) else
    if (identical(family$family, "gaussian")) FALSE else any(v == 0)
  if (boundary || (separable && length(conds) == 1L)) {
    .perm_stop(
      "solomonR_perm_no_variation_error",
      "The outcome does not vary in the ", paste(names(conds)[constant], collapse = " and "),
      " condition", if (sum(constant) > 1) "s" else "", ", which the contrast \"", contrast,
      "\" uses, so every permutation gives the same statistic and the test has no information. ",
      "Choose a contrast on a condition whose outcome varies."
    )
  }
  # The statistic still varies (decision D5's proposal): a p-value, with a
  # classed warning so that the constant condition is seen.
  warning(structure(
    class = c("solomonR_perm_no_variation_warning", "warning", "condition"),
    list(message = paste0(
      "The outcome does not vary in the ", paste(names(conds)[constant], collapse = " and "),
      " condition; the contrast \"", contrast, "\" varies only through the other condition",
      if (!separable) " and the covariates" else "", "."), call = NULL)
  ))
  invisible(NULL)
}
