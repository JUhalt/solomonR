# Negative-binomial (NB2) fits for fit_solomon_glm(family = "negative_binomial")
# (issue #62), following Cameron and Trivedi (2013, sec. 3.3).

# Families whose dispersion is fixed rather than estimated, so that tests use
# the normal distribution, as summary.glm() and MASS's summary for glm.nb()
# fits do.
.fixed_dispersion <- function(family) {
  family$family %in% c("binomial", "poisson") ||
    startsWith(family$family, "Negative Binomial")
}

.is_negbin <- function(family) {
  startsWith(family$family, "Negative Binomial")
}

# NB2 maximum-likelihood fit. When the counts show little overdispersion, the
# estimate of theta grows without bound and glm.nb() warns that its iteration
# limit was reached; that is reported once, as a classed warning.
.fit_negbin <- function(fml, df) {
  boundary <- FALSE
  fit <- withCallingHandlers(
    MASS::glm.nb(fml, data = df, na.action = stats::na.exclude),
    warning = function(w) {
      if (grepl("iteration limit reached|alternation limit reached", conditionMessage(w))) {
        boundary <<- TRUE
        invokeRestart("muffleWarning")
      }
    }
  )

  if (boundary || !is.null(fit$th.warn)) {
    warning(structure(
      class = c("solomonR_theta_boundary_warning", "warning", "condition"),
      list(
        message = sprintf(paste0(
          "The negative-binomial dispersion parameter did not converge ",
          "(theta = %.3g), which happens when the counts show little ",
          "overdispersion. The coefficients are then close to those of ",
          "family = poisson(), which may be preferred."
        ), fit$theta),
        call = NULL
      )
    ))
  }

  fit
}
