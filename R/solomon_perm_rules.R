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
