# Numeric covariates standardized for the inference of fit_solomon_glm()
# (issue #137).
#
# A numeric covariate far from zero relative to its spread, such as a
# calendar year, makes the design of the model ill-conditioned: the robust
# covariances computed from it lose accuracy, as do the delta method and the
# bootstrap of marginal_solomon() and the cell means of
# plot_sensitization(), and the cross-products that the influence functions
# of the pretest effects were computed from became singular (see
# .glm_influence()). The inference is therefore computed from the model
# refitted with each numeric covariate standardized: centered at its mean
# and divided by its standard deviation over the rows the model used. That
# is a reparameterization, X_s = X A for a fixed matrix A, so the fitted
# values, the model-based, HC3, and CR2 covariances of the contrasts, and
# their Satterthwaite degrees of freedom are the same in exact arithmetic;
# the contrasts give no weight to the intercept or the covariates. The
# model, the coefficients, and their covariance are returned for the
# covariates as given (.covariate_map()), and the standardized coefficients
# and covariance in `covariate_scaling`, for the functions that take the
# fit.

# The numeric covariates of `df` standardized over the rows the model used
# (`in_model`). `covariates` names the covariates among the columns of
# `df`. Every numeric vector is standardized, whatever its class (such as a
# labelled column read from SPSS or Stata, or one wrapped in I()), as the
# numbers it enters the model as, and so are dates and times, which enter
# it as days or seconds since 1970, or in the units of a difference
# (.covariate_values()); a factor, a logical, or a character column enters
# as indicators and is left as given, as is a matrix. Returns the data with
# them standardized, `data`, and their `means` and `scales` (standard
# deviations), named by covariate; `means` is empty when there is none.
.covariate_scaling <- function(df, covariates, in_model) {
  means <- numeric()
  scales <- numeric()
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
# intercept loses each covariate's mean times its coefficient. `terms` are
# the names of the coefficients.
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

# The model refitted to the standardized covariates, as fit_solomon_glm()
# fitted it to the covariates as given. The fit as given has already given
# any warning, such as one that theta did not converge, so the refit gives
# none.
.fit_standardized <- function(fml, data, family, negbin) {
  suppressWarnings(if (negbin) {
    .fit_negbin(fml, data)
  } else {
    stats::glm(fml, data = data, family = family, na.action = stats::na.exclude)
  })
}

# The standardized model of a fit_solomon_glm() fit with `covariate_scaling`.
.standardized_model <- function(fit) {
  .fit_standardized(stats::formula(fit$model),
                    .standardized_frame(fit$data, fit$covariate_scaling), fit$family,
                    .is_negbin(fit$family))
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
    tidy$estimate <- unname(drop(M %*% cf)[tidy$term])
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
