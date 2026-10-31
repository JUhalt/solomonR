# Numeric covariates and the pretest standardized for the inference of
# fit_solomon_glm() (issue #137).
#
# A numeric covariate far from zero relative to its spread, such as a
# calendar year, makes the design of the model ill-conditioned: the robust
# covariances computed from it lose accuracy, as do the delta method and the
# bootstrap of marginal_solomon() and the cell means of
# plot_sensitization(), and the cross-products that the influence functions
# of the pretest effects were computed from became singular (see
# .glm_influence()). A covariate or a pretest on a scale far from 1, such as
# 1e-8 or 1e8, does the same to the CR2 covariance and to the numerical
# derivatives of the delta method. The inference is therefore computed from
# the model refitted with each numeric covariate standardized, centered at
# its mean and divided by its standard deviation over the rows the model
# used, and with the pretest, already centered, divided by the power of two
# nearest its standard deviation. That is a reparameterization, X_s = X A
# for a fixed matrix A, so the fitted values, the model-based, HC3, and CR2
# covariances of the contrasts, and their Satterthwaite degrees of freedom
# are the same in exact arithmetic. The contrasts give no weight to the
# intercept or the covariates; those that weight the pretest's coefficient
# (the pretest effects at the other center) weight it in the units of the
# scaled pretest (.pretest_scale()). The model, the coefficients, and their
# covariance are returned for the covariates and the pretest as given
# (.covariate_map()), and the standardized coefficients and covariance in
# `covariate_scaling`, for the functions that take the fit.

# The numeric covariates of `df` standardized over the rows the model used
# (`in_model`). `covariates` names the covariates among the columns of
# `df`. Every numeric vector is standardized, whatever its class (such as a
# labeled column read from SPSS or Stata, or one wrapped in I()), as the
# numbers it enters the model as, and so are dates and times, which enter
# it as days or seconds since 1970, or in the units of a difference
# (.covariate_values()); a factor, a logical, or a character column enters
# as indicators and is left as given, as is a matrix. The pretest,
# `pre_obs`, is centered already, at the mean at which the pretest effects
# are evaluated, so it is only divided, by the power of two nearest its
# standard deviation over the pretested participants in the model; dividing
# by a power of two is exact. Returns the data with them standardized,
# `data`, and their `means` and `scales` (standard deviations, and that
# power of two), named by covariate, the pretest first, with a mean of 0;
# `means` is empty when there is neither. Without a pretest, a covariate may
# be named `pre_obs`; it is standardized as a covariate.
.covariate_scaling <- function(df, covariates, in_model) {
  means <- numeric()
  scales <- numeric()
  if ("pre_obs" %in% setdiff(names(df), covariates)) {
    s <- stats::sd(df$pre_obs[in_model & df$pretested %in% 1])
    s <- if (is.finite(s) && s > 0) 2^round(log2(s)) else 1
    df$pre_obs <- df$pre_obs / s
    means[["pre_obs"]] <- 0
    scales[["pre_obs"]] <- s
  }
  for (name in covariates) {
    x <- .covariate_values(df[[name]])
    if (is.null(x)) next
    m <- mean(x[in_model])
    s <- stats::sd(x[in_model])
    if (!is.finite(m)) next
    if (!is.finite(s) || s == 0) s <- 1
    df[[name]] <- (x - m) / s
    means[[name]] <- m
    scales[[name]] <- s
  }
  list(data = df, means = means, scales = scales)
}

# The number the pretest was divided by for the inference: a pretest
# coefficient of the standardized model is in its units, and a shift of the
# pretest in the units of the pretest as given is divided by it (`scaling`
# from .covariate_scaling()). 1 without a pretest.
.pretest_scale <- function(scaling) {
  if ("pre_obs" %in% names(scaling$scales)) scaling$scales[["pre_obs"]] else 1
}

# The numbers a covariate enters the model as, without its class: those of
# a numeric vector, and those of a date, a time, or a time difference. NULL
# for any other column.
.covariate_values <- function(x) {
  if (inherits(x, c("Date", "POSIXct", "difftime"))) return(as.numeric(x))
  if (is.numeric(x) && is.null(dim(x))) return(as.vector(unclass(x)))
  NULL
}

# The columns of a data frame (`frame`) or of a design matrix (`X`)
# standardized as `scaling` (.covariate_scaling()) standardized the data.
.standardized_frame <- function(frame, scaling) {
  for (name in names(scaling$means)) {
    frame[[name]] <- (.covariate_values(frame[[name]]) - scaling$means[[name]]) /
      scaling$scales[[name]]
  }
  frame
}

.standardized_columns <- function(X, scaling) {
  columns <- .formula_names(names(scaling$means))
  for (j in seq_along(columns)) {
    X[, columns[j]] <- (X[, columns[j]] - scaling$means[[j]]) / scaling$scales[[j]]
  }
  X
}

# The map from the coefficients of the model with standardized covariates
# to those of the covariates as given, b = M b_s, so that V = M V_s M'. A
# covariate's coefficient is divided by its standard deviation, and the
# intercept loses each covariate's mean times its coefficient; the
# pretest's coefficient is divided by its scale, and its mean is 0. `terms`
# are the names of the coefficients.
.covariate_map <- function(terms, scaling) {
  M <- diag(length(terms))
  dimnames(M) <- list(terms, terms)
  columns <- .formula_names(names(scaling$means))
  for (j in seq_along(columns)) {
    M[columns[j], columns[j]] <- 1 / scaling$scales[[j]]
    M["(Intercept)", columns[j]] <- -scaling$means[[j]] / scaling$scales[[j]]
  }
  M
}

# Whether .covariate_scaling() standardizes anything in `df`: a pretest, or
# a numeric covariate among `covariates`.
.standardizes <- function(df, covariates) {
  "pre_obs" %in% setdiff(names(df), covariates) ||
    any(vapply(covariates, function(name) !is.null(.covariate_values(df[[name]])), logical(1)))
}

# The model of fit_solomon_glm() fitted to `data`, a negative-binomial one
# when `negbin` is TRUE. With `quiet`, without its warnings.
.fit_model <- function(fml, data, family, negbin, quiet = FALSE) {
  fit <- function() {
    if (negbin) {
      .fit_negbin(fml, data)
    } else {
      stats::glm(fml, data = data, family = family, na.action = stats::na.exclude)
    }
  }
  if (quiet) suppressWarnings(fit()) else fit()
}

# The standardized model of a fit_solomon_glm() fit with `covariate_scaling`,
# refitted. The fit gave its warnings when it was made, so this one gives
# none.
.standardized_model <- function(fit) {
  .fit_model(stats::formula(fit$model), .standardized_frame(fit$data, fit$covariate_scaling),
             fit$family, .is_negbin(fit$family), quiet = TRUE)
}

# The coefficient table of fit_solomon_glm(), in the columns of
# .solomon_coefficient_columns, from the model `fit`, its covariance `V`,
# and, for CR2, the complete-case refit `fit_cr` that `V` was computed on.
# With standardized covariates (`scaling` from .covariate_scaling(), with
# `fit` and `fit_cr` fitted to them), the coefficients and their covariance
# are mapped to the covariates as given; the Satterthwaite degrees of
# freedom of a CR2 intercept are then those of the contrast that maps it,
# and a covariate's are those of its standardized coefficient. Returns the
# table, `table`, and the covariance of the coefficients as given, `vcov`:
# `V` itself without standardized covariates, and otherwise a matrix.
.coefficient_table <- function(fit, fit_cr, V, scaling, robust, df_model, conf_level) {
  tidy <- broom::tidy(fit)
  standardized <- length(scaling$means) > 0L
  V_given <- V
  if (standardized) {
    cf <- stats::coef(fit)
    M <- .covariate_map(names(cf), scaling)
    # broom names the estimates of a negative-binomial fit by term; they
    # keep their names.
    estimate <- drop(M %*% cf)[tidy$term]
    tidy$estimate <- if (is.null(names(tidy$estimate))) unname(estimate) else estimate
    V_given <- M %*% as.matrix(V) %*% t(M)
  }
  se_vec <- sqrt(diag(V_given))
  tidy$std.error <- unname(se_vec[match(tidy$term, names(se_vec))])
  tidy$statistic <- tidy$estimate / tidy$std.error
  if (robust == "CR2") {
    ct <- clubSandwich::coef_test(fit_cr, vcov = V, test = "Satterthwaite")
    tidy$df <- ct$df_Satt[match(tidy$term, ct$Coef)]
    if (standardized) {
      tidy$df[tidy$term == "(Intercept)"] <- as.data.frame(clubSandwich::linear_contrast(
        fit_cr, vcov = V, contrasts = M["(Intercept)", , drop = FALSE], test = "Satterthwaite"
      ))$df
    }
  } else {
    tidy$df <- df_model
  }
  tidy$p.value <- 2 * stats::pt(-abs(tidy$statistic), df = tidy$df)
  coef_ci <- .wald_ci(tidy$estimate, tidy$std.error, tidy$df, conf_level)
  tidy$conf.low <- unname(coef_ci[, "conf.low"])
  tidy$conf.high <- unname(coef_ci[, "conf.high"])
  list(table = tidy[, .solomon_coefficient_columns], vcov = V_given)
}

# The `covariate_scaling` element of a fit: the `means` and `scales` of the
# standardized covariates, and the `coefficients` and the covariance `vcov`
# of the standardized model `fit` (for CR2, the clubSandwich object); NULL
# without numeric covariates.
.covariate_scaling_element <- function(scaling, fit, V) {
  if (!length(scaling$means)) return(NULL)
  list(means = scaling$means, scales = scaling$scales, coefficients = stats::coef(fit),
       vcov = V)
}
