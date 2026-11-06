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
# (.covariate_map(), .given_vcov()), and the standardized coefficients and
# covariance in `covariate_scaling`, for the functions that take the fit.

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
# The columns of `X` are those of the model matrix of `model`, and are
# found by position (.scaled_columns()).
.standardized_frame <- function(frame, scaling) {
  for (name in names(scaling$means)) {
    frame[[name]] <- (.covariate_values(frame[[name]]) - scaling$means[[name]]) /
      scaling$scales[[name]]
  }
  frame
}

.standardized_columns <- function(X, scaling, model) {
  columns <- .scaled_columns(model, scaling)
  for (j in seq_along(columns)) {
    X[, columns[[j]]] <- (X[, columns[[j]]] - scaling$means[[j]]) / scaling$scales[[j]]
  }
  X
}

# The positions of the standardized covariates and of the scaled pretest of
# `scaling` among the columns of the model matrix of `model`, which are
# those of its coefficients, named as in `scaling`. Each is found from the
# term it enters the model as, not from the name of its column, which
# another column can share: a numeric covariate named `wave2` and the
# indicator of the second level of a factor `wave` both give a column
# `wave2`, so that a map indexed by name rescaled the indicator (issue
# #137). A covariate whose term is not found is looked up by name.
.scaled_columns <- function(model, scaling) {
  X <- stats::model.matrix(model)
  labels <- attr(stats::terms(model), "term.labels")
  assign <- attr(X, "assign")
  columns <- .formula_names(names(scaling$means))
  index <- vapply(columns, function(column) {
    j <- which(assign == match(column, labels))
    if (length(j) == 1L) j else match(column, colnames(X))
  }, integer(1), USE.NAMES = FALSE)
  if (anyNA(index)) {
    stop("The model has no column for ", paste(columns[is.na(index)], collapse = ", "), ".",
         call. = FALSE)
  }
  stats::setNames(index, names(scaling$means))
}

# The map from the coefficients of the model with standardized covariates
# to those of the covariates as given, b = M b_s, so that V = M V_s M'. A
# covariate's coefficient is divided by its standard deviation, and the
# intercept loses each covariate's mean times its coefficient; the
# pretest's coefficient is divided by its scale, and its mean is 0. `model`
# is either fit, as given or standardized: they have the same terms, whose
# columns are found by position (.scaled_columns()).
.covariate_map <- function(model, scaling) {
  terms <- names(stats::coef(model))
  M <- diag(length(terms))
  dimnames(M) <- list(terms, terms)
  columns <- .scaled_columns(model, scaling)
  intercept <- match("(Intercept)", terms)
  for (j in seq_along(columns)) {
    M[columns[[j]], columns[[j]]] <- 1 / scaling$scales[[j]]
    M[intercept, columns[[j]]] <- -scaling$means[[j]] / scaling$scales[[j]]
  }
  M
}

# The covariance of the coefficients as given, M V_s M', from the
# covariance `V` of the standardized model `fit` and the map `M`
# (.covariate_map()). Formed as M V_s M', it can cancel: nearly collinear
# covariates have standardized coefficients with large variances whose
# combination, the variance of the intercept as given, is small, so that
# their rounding error swamps it. With a covariate and its copy in single
# precision, the HC3 and CR2 standard errors of the intercept came out NaN,
# or off by a fifth of their value for covariates 10^-6 standard deviations
# apart, where the fit to the covariates as given was accurate. It is
# therefore computed through U = R^{-T} M', with R from the QR
# decomposition W^(1/2) X = QR of the standardized model, which
# backsolve() computes without that cancellation:
# - model-based ("none"): sigma^2 U'U;
# - HC3: G'G, a row of G being a participant's influence on M b divided by
#   1 - h_i (.glm_influence());
# - CR2: U' V_o U, with V_o the CR2 covariance of the standardized model in
#   an orthonormal basis and R the decomposition that gives the basis
#   (`basis`, from .cr2_basis()).
# A map that is diagonal only divides, which cancels nothing: the pretest
# alone, divided by a power of two, which is exact, so that such fits keep
# the covariance they had, or covariates whose mean is zero. Without a
# basis it is applied as M V_s M'. Returns a matrix named by the terms.
.given_vcov <- function(fit, V, M, robust, basis = NULL) {
  mapped <- function() M %*% as.matrix(V) %*% t(M)
  if (robust == "CR2") {
    if (is.null(basis)) return(mapped())
    U <- backsolve(basis$R, t(M), transpose = TRUE)
    out <- crossprod(U, .plain_matrix(basis$vcov) %*% U)
    out <- (out + t(out)) / 2
    dimnames(out) <- dimnames(M)
    return(out)
  }
  if (.diagonal_map(M)) return(mapped())
  q <- fit$qr
  p <- ncol(M)
  if (is.null(q) || q$rank < p || !identical(q$pivot, seq_len(p))) return(mapped())
  U <- backsolve(qr.R(q), t(M), transpose = TRUE)
  if (robust == "none") {
    G <- sqrt(summary(fit)$dispersion) * U
  } else {
    infl <- .glm_influence(fit, map = M)
    G <- infl$influence / (1 - infl$hat)
  }
  out <- crossprod(G)
  dimnames(out) <- dimnames(M)
  out
}

# Whether the map `M` (.covariate_map()) is diagonal: it only divides.
.diagonal_map <- function(M) all(M[row(M) != col(M)] == 0)

# The robust covariances of the standardized model are computed from its QR
# decomposition W^(1/2) X = QR when it has a numeric covariate, rather than
# by sandwich::vcovHC() and clubSandwich::vcovCR() (issue #137). Those form
# them from the inverse of X'WX, whose condition number is the square of
# the design's, so that with nearly collinear covariates their entries lost
# their accuracy, though the coefficient table, computed through R, kept
# it. With covariates 10^-8 standard deviations apart near zero, the HC3
# standard error of the average treatment effect came out 0.117 for 0.151,
# without a warning; at 10^-10 those of other contrasts were up to 54 times
# too large. CR2 fits stopped inside clubSandwich ("the leading minor ...
# is not positive"), or its adjustment matrices, computed from that
# inverse with its eigenvalues below 1e-12 set to zero, put the
# coefficients' standard errors and degrees of freedom off by percents.
# A fit without a numeric covariate, whose design is well conditioned,
# keeps the computations of sandwich and clubSandwich, and its results.
# `scaling` comes from .covariate_scaling().
.numeric_covariates <- function(scaling) any(names(scaling$means) != "pre_obs")

# A covariance as a plain matrix, without the attributes of a clubSandwich
# object.
.plain_matrix <- function(V) matrix(as.numeric(V), nrow(V), ncol(V))

# The HC3 covariance of a glm fit (MacKinnon & White, 1985): the sum of the
# cross-products of the participants' influences, each divided by 1 - h_i,
# from the fit's QR decomposition (.glm_influence()). It is the covariance
# that sandwich::vcovHC() computes, to about 1e-13 on well-conditioned
# designs, without forming the inverse of X'WX. Named by the coefficients.
.hc3_vcov <- function(model) {
  infl <- .glm_influence(model)
  crossprod(infl$influence / (1 - infl$hat))
}

# The variances of the combinations `L` (one per row) of the coefficients
# of the glm fit `model`, model-based ("none") or HC3, through its QR
# decomposition W^(1/2) X = QR: sigma^2 ||R^{-T} l||^2, or the sum of the
# squared influences on l'b, each divided by 1 - h_i (.glm_influence()).
# The covariance of nearly collinear covariates has large entries, which
# cancel in a combination that weights them nearly equally, as the
# gradients of the delta method do; these do not form it.
.combination_variance <- function(model, robust, L) {
  if (robust == "HC3") {
    infl <- .glm_influence(model, map = L)
    return(colSums((infl$influence / (1 - infl$hat))^2))
  }
  q <- model$qr
  rank <- seq_len(q$rank)
  U <- backsolve(qr.R(q)[rank, rank, drop = FALSE], t(L[, q$pivot[rank], drop = FALSE]),
                 transpose = TRUE)
  summary(model)$dispersion * colSums(U^2)
}

# The CR2 covariance (Bell & McCaffrey, 2002) of the glm fit `fit_cr`, with
# `cluster` the clusters of its rows, computed by clubSandwich in an
# orthonormal basis of the same columns. With W^(1/2) X = QR the fit's
# decomposition, the model is refitted to X_o = X R^{-1}, whose weighted
# columns are those of Q and so orthonormal. That is a reparameterization,
# b_o = R b, with the same fitted values, iterations, and weights, and the
# CR2 adjustments with clubSandwich's working model, the covariance, and
# the Satterthwaite degrees of freedom (Pustejovsky & Tipton, 2018) do not
# depend on the parameterization: a combination L b is the combination
# L R^{-1} b_o, which backsolve() computes. In the basis, the inverse of
# X'WX is the identity, and clubSandwich's computations are accurate. On
# well-conditioned designs the covariance equals that of clubSandwich on
# `fit_cr` to about 1e-12, and the degrees of freedom to about 1e-13.
# Returns `R`, the refit `fit`, and its CR2 covariance `vcov`, the
# clubSandwich object; NULL when the fit holds no decomposition of full
# rank in the order of its terms.
.cr2_basis <- function(fit_cr, cluster) {
  q <- fit_cr$qr
  p <- length(stats::coef(fit_cr))
  if (is.null(q) || q$rank < p || !identical(q$pivot, seq_len(p))) return(NULL)
  R <- qr.R(q)
  frame <- list(y = fit_cr$y,
                basis = t(backsolve(R, t(stats::model.matrix(fit_cr)), transpose = TRUE)))
  fit <- suppressWarnings(stats::glm(y ~ 0 + basis, family = stats::family(fit_cr), data = frame,
                                     offset = fit_cr$offset, control = fit_cr$control))
  list(R = R, fit = fit, vcov = clubSandwich::vcovCR(fit, cluster = cluster, type = "CR2"))
}

# The covariance of the standardized coefficients from the basis of
# .cr2_basis(), R^{-1} V_o R^{-T}, through backsolve(), named by `terms`.
.cr2_standardized_vcov <- function(basis, terms) {
  V <- backsolve(basis$R, t(backsolve(basis$R, .plain_matrix(basis$vcov))))
  V <- (V + t(V)) / 2
  dimnames(V) <- list(terms, terms)
  V
}

# The Satterthwaite degrees of freedom (Pustejovsky & Tipton, 2018) of the
# CR2 tests of the combinations `L` (one per row) of the coefficients of
# `fit_cr`, whose CR2 covariance is `V`, from clubSandwich: in the
# orthonormal basis of .cr2_basis() when `basis` is given, and with `fit_cr`
# and `V` otherwise.
.cr2_df <- function(L, fit_cr, V, basis = NULL) {
  L <- matrix(L, ncol = length(stats::coef(fit_cr)))
  if (!is.null(basis)) {
    fit_cr <- basis$fit
    V <- basis$vcov
    L <- t(backsolve(basis$R, t(L), transpose = TRUE))
  }
  as.data.frame(clubSandwich::linear_contrast(fit_cr, vcov = V, contrasts = L,
                                              test = "Satterthwaite"))$df
}

# The HTZ test (Pustejovsky & Tipton, 2018) of the constraints L b = 0 with
# the CR2 covariance `V` of `fit_cr`, from clubSandwich, in the basis of
# .cr2_basis() when `basis` is given, as for .cr2_df().
.cr2_wald <- function(L, fit_cr, V, basis = NULL) {
  if (!is.null(basis)) {
    fit_cr <- basis$fit
    V <- basis$vcov
    L <- t(backsolve(basis$R, t(L), transpose = TRUE))
  }
  clubSandwich::Wald_test(fit_cr, constraints = L, vcov = V, test = "HTZ")
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
# are mapped to the covariates as given (.given_vcov()). Under CR2 with a
# numeric covariate, `basis` is the orthonormal basis of .cr2_basis(), and
# the Satterthwaite degrees of freedom of each coefficient as given are
# those of the combination of the standardized coefficients that its row
# of the map gives (.cr2_df()). Returns the table, `table`, and the
# covariance of the coefficients as given, `vcov`: `V` itself without
# standardized covariates, and otherwise a matrix (.given_vcov()). The rows
# of the table are the coefficients, in their order, and are matched to the
# coefficients and the covariance by position, since two terms can share a
# name (.scaled_columns()).
.coefficient_table <- function(fit, fit_cr, V, scaling, robust, df_model, conf_level,
                               basis = NULL) {
  tidy <- .tidy_coefficients(fit)
  standardized <- length(scaling$means) > 0L
  cf <- stats::coef(fit)
  rows <- if (identical(tidy$term, names(cf))) seq_along(cf) else match(tidy$term, names(cf))
  V_given <- V
  if (standardized) {
    M <- .covariate_map(fit, scaling)
    # .tidy_coefficients() names the estimates of a negative-binomial fit by
    # term, as broom did; they keep their names.
    estimate <- drop(M %*% cf)[rows]
    tidy$estimate <- if (is.null(names(tidy$estimate))) unname(estimate) else estimate
    V_given <- .given_vcov(fit, V, M, robust, if (robust == "CR2") basis)
  }
  se_vec <- sqrt(diag(V_given))
  tidy$std.error <- unname(se_vec[rows])
  tidy$statistic <- tidy$estimate / tidy$std.error
  if (robust == "CR2" && standardized && !is.null(basis)) {
    tidy$df <- .cr2_df(M, fit_cr, V, basis)[rows]
  } else if (robust == "CR2") {
    # Only the degrees of freedom are taken from clubSandwich. With
    # standardized covariates the standard errors it computes as well, of
    # the standardized coefficients, are not used.
    # clubSandwich names the rows of its table by the coefficients, which
    # must then be unique; the rows are the coefficients in their order.
    df_only <- function(expr) if (standardized) suppressWarnings(expr) else expr
    fit_named <- fit_cr
    names(fit_named$coefficients) <- make.unique(names(fit_named$coefficients))
    ct <- df_only(clubSandwich::coef_test(fit_named, vcov = V, test = "Satterthwaite"))
    tidy$df <- if (nrow(ct) == length(cf)) ct$df_Satt[rows] else
      ct$df_Satt[match(tidy$term, ct$Coef)]
    if (standardized) {
      # The intercept as given is the combination of the standardized
      # coefficients that its row of the map gives: the standardized
      # intercept for a diagonal map.
      tidy$df[tidy$term == "(Intercept)"] <-
        df_only(.cr2_df(M["(Intercept)", , drop = FALSE], fit_cr, V))
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
# standardized covariates and the scaled pretest, and the `coefficients`
# and the covariance `vcov` of the standardized model `fit` (for CR2, the
# clubSandwich object, or with a numeric covariate the matrix of
# .cr2_standardized_vcov()); NULL with neither a pretest nor a numeric
# covariate.
.covariate_scaling_element <- function(scaling, fit, V) {
  if (!length(scaling$means)) return(NULL)
  list(means = scaling$means, scales = scaling$scales, coefficients = stats::coef(fit),
       vcov = V)
}
