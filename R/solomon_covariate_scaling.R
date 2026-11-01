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
# covariance `V` of the standardized model and the map `M`
# (.covariate_map()). Formed as M V_s M', it can cancel: nearly collinear
# covariates have standardized coefficients with large variances whose
# combination, the variance of the intercept as given, is small, so that
# their rounding error swamps it. With a covariate and its copy in single
# precision, the HC3 and CR2 standard errors of the intercept came out NaN,
# or off by a fifth of their value for covariates 10^-6 standard deviations
# apart, where the fit to the covariates as given was accurate. It is
# therefore computed from a factor of V_s = S S', as G'G with G = (M S)',
# and G from the QR decomposition W^(1/2) X = QR of the standardized model
# through R^{-T} M', which backsolve() computes without that cancellation:
# - model-based ("none"): S = sigma R^{-1}, so that G = sigma R^{-T} M';
# - HC3: a row of G is a participant's influence on M b divided by 1 - h_i
#   (.glm_influence());
# - CR2: a column of S is a cluster's adjusted estimating function,
#   (X'WX)^{-1} X_j' W_j A_j e_j (Bell & McCaffrey, 2002), from the matrices
#   X_j' W_j and A_j and the clusters that clubSandwich keeps with `V`, and
#   the residuals of the complete-case refit `fit_cr` it was computed on.
#   They are first checked to give `V` as clubSandwich computes it, and its
#   bread over its scale to be the inverse of R'R; should clubSandwich keep
#   other matrices, the check fails and the covariance is M V_s M'.
# A map that is diagonal only divides, which cancels nothing: the pretest
# alone, divided by a power of two, which is exact, so that such fits keep
# the covariance they had, or covariates whose mean is zero. It is applied
# as M V_s M'. `parts`, for CR2, are those of .cr2_parts(), computed here
# when not given. Returns a matrix named by the terms.
.given_vcov <- function(fit, fit_cr, V, M, robust, parts = NULL) {
  model <- if (robust == "CR2") fit_cr else fit
  mapped <- function() M %*% as.matrix(V) %*% t(M)
  if (.diagonal_map(M)) return(mapped())
  q <- model$qr
  p <- ncol(M)
  if (is.null(q) || q$rank < p || !identical(q$pivot, seq_len(p))) return(mapped())
  R <- qr.R(q)
  U <- backsolve(R, t(M), transpose = TRUE)
  if (robust == "none") {
    G <- sqrt(summary(model)$dispersion) * U
  } else if (robust == "HC3") {
    infl <- .glm_influence(model, map = M)
    G <- infl$influence / (1 - infl$hat)
  } else {
    if (is.null(parts)) parts <- .cr2_parts(model, V)
    if (is.null(parts)) return(mapped())
    G <- crossprod(backsolve(R, parts$comp, transpose = TRUE), U)
  }
  out <- crossprod(G)
  dimnames(out) <- dimnames(M)
  out
}

# Whether the map `M` (.covariate_map()) is diagonal: it only divides.
.diagonal_map <- function(M) all(M[row(M) != col(M)] == 0)

# The matrices of a CR2 covariance `V` of the glm fit `fit_cr` that
# .given_vcov() and .cr2_contrast_df() use, from those clubSandwich keeps
# with `V`: the R of the fit's decomposition W^{1/2} X = QR; the clusters'
# estimating matrices X_j' W_j (`est`), adjustment matrices A_j (`adj`),
# and working variances Phi_j (`target`); and their adjusted estimating
# functions X_j' W_j A_j e_j (`comp`, one column per cluster), with e_j the
# residuals. They are checked to give `V` as clubSandwich computes it, and
# its bread over its scale to be the inverse of R'R; should clubSandwich
# keep other matrices, the check fails and NULL is returned, as it is when
# the fit has no decomposition of full rank in the order of its terms.
.cr2_parts <- function(fit_cr, V) {
  q <- fit_cr$qr
  p <- nrow(as.matrix(V))
  if (is.null(q) || q$rank < p || !identical(q$pivot, seq_len(p))) return(NULL)
  R <- qr.R(q)
  est <- attr(V, "est_mats")
  adj <- attr(V, "adjustments")
  cluster <- attr(V, "cluster")
  bread <- attr(V, "bread")
  v_scale <- attr(V, "v_scale")
  e <- stats::residuals(fit_cr, type = "response")
  if (!is.list(est) || !is.list(adj) || length(est) != length(adj) || !is.factor(cluster) ||
      length(cluster) != length(e) || nlevels(cluster) != length(est) ||
      !is.matrix(bread) || !identical(dim(bread), c(p, p)) || !is.numeric(v_scale)) {
    return(NULL)
  }
  comp <- tryCatch(
    do.call(cbind, Map(function(x, a, r) x %*% a %*% r, est, adj, split(e, cluster))),
    error = function(err) NULL
  )
  if (is.null(comp) || !identical(dim(comp), c(p, length(est))) ||
      !isTRUE(all.equal(unname(bread %*% (tcrossprod(comp) / v_scale) %*% bread / v_scale),
                        unname(as.matrix(V)), tolerance = 1e-10)) ||
      !isTRUE(all.equal(unname(bread / v_scale), chol2inv(R), tolerance = 1e-10))) {
    return(NULL)
  }
  list(R = R, comp = comp, est = est, adj = adj, target = attr(V, "target"),
       inverse_var = isTRUE(attr(V, "inverse_var")))
}

# The Satterthwaite degrees of freedom of the CR2 test of the combination
# c'b of the coefficients (Pustejovsky & Tipton, 2018), `contrast` being c,
# as clubSandwich::linear_contrast() computes them, from the matrices of
# .cr2_parts(), without forming (X'WX)^{-1} c. For a glm, clubSandwich takes
# the working model of cluster j to be Phi_j = W_j^{-1}, the inverse of its
# weights, and then, with X_j the cluster's rows of the working design and
# A_j its adjustment,
#   f_j = A_j' W_j X_j (X'WX)^{-1} c,
#   P_jk = [j = k] f_j' Phi_j f_j - f_j' X_j (X'WX)^{-1} X_k' f_k,
#   df = (sum_j P_jj)^2 / sum_jk P_jk^2.
# It computes W_j X_j (X'WX)^{-1} c from the inverse it forms, whose
# entries are large for nearly collinear covariates and cancel in such a
# combination as the intercept as given: for two covariates 10^-8 standard
# deviations apart at a mean of 5, the intercept's degrees of freedom came
# out 6.18 for 18.88, as the same model with the second covariate replaced
# by their difference, an exact reparameterization, gives. Here
# (X'WX)^{-1} = R^{-1} R^{-T}, so that with g = R^{-T} c and
# K_j = R^{-T} X_j' W_j, both from backsolve(), which form no large terms,
#   f_j = A_j' K_j' g,  R^{-T} X_j' f_j = K_j Phi_j f_j = h_j,
#   P_jk = [j = k] f_j' Phi_j f_j - h_j' h_k,
# the computation through R^{-T} that the standard error takes
# (.given_vcov()). On well-conditioned designs it gives clubSandwich's
# degrees of freedom to about 1e-12. NULL without the parts, or with
# another working model.
.cr2_contrast_df <- function(parts, contrast) {
  if (is.null(parts) || !parts$inverse_var || !is.list(parts$target) ||
      length(parts$target) != length(parts$est)) {
    return(NULL)
  }
  g <- backsolve(parts$R, contrast, transpose = TRUE)
  clusters <- length(parts$est)
  diagonal <- numeric(clusters)
  H <- matrix(0, length(g), clusters)
  for (j in seq_len(clusters)) {
    K <- backsolve(parts$R, parts$est[[j]], transpose = TRUE)
    f <- drop(crossprod(parts$adj[[j]], crossprod(K, g)))
    phi_f <- drop(parts$target[[j]] %*% f)
    diagonal[j] <- sum(f * phi_f)
    H[, j] <- K %*% phi_f
  }
  P <- -crossprod(H)
  diag(P) <- diag(P) + diagonal
  sum(diag(P))^2 / sum(P^2)
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
# freedom of a CR2 intercept are then those of the contrast that maps it
# (.cr2_contrast_df()), and a covariate's are those of its standardized
# coefficient. Returns the
# table, `table`, and the covariance of the coefficients as given, `vcov`:
# `V` itself without standardized covariates, and otherwise a matrix
# (.given_vcov()). The rows of the table are the coefficients, in their
# order, and are matched to the coefficients and the covariance by
# position, since two terms can share a name (.scaled_columns()).
.coefficient_table <- function(fit, fit_cr, V, scaling, robust, df_model, conf_level) {
  tidy <- broom::tidy(fit)
  standardized <- length(scaling$means) > 0L
  cf <- stats::coef(fit)
  rows <- if (identical(tidy$term, names(cf))) seq_along(cf) else match(tidy$term, names(cf))
  V_given <- V
  parts <- NULL
  if (standardized) {
    M <- .covariate_map(fit, scaling)
    # broom names the estimates of a negative-binomial fit by term; they
    # keep their names.
    estimate <- drop(M %*% cf)[rows]
    tidy$estimate <- if (is.null(names(tidy$estimate))) unname(estimate) else estimate
    if (robust == "CR2" && !.diagonal_map(M)) parts <- .cr2_parts(fit_cr, V)
    V_given <- .given_vcov(fit, fit_cr, V, M, robust, parts)
  }
  se_vec <- sqrt(diag(V_given))
  tidy$std.error <- unname(se_vec[rows])
  tidy$statistic <- tidy$estimate / tidy$std.error
  if (robust == "CR2") {
    # Only the degrees of freedom are taken from clubSandwich. With
    # standardized covariates the standard errors it computes as well, of
    # the standardized coefficients and of the intercept as given from
    # M V_s M', are not used; with nearly collinear covariates their
    # rounding error can make them NaN, with a warning (.given_vcov()).
    # clubSandwich names the rows of its table by the coefficients, which
    # must then be unique; the rows are the coefficients in their order.
    df_only <- function(expr) if (standardized) suppressWarnings(expr) else expr
    fit_named <- fit_cr
    names(fit_named$coefficients) <- make.unique(names(fit_named$coefficients))
    ct <- df_only(clubSandwich::coef_test(fit_named, vcov = V, test = "Satterthwaite"))
    tidy$df <- if (nrow(ct) == length(cf)) ct$df_Satt[rows] else
      ct$df_Satt[match(tidy$term, ct$Coef)]
    if (standardized) {
      # The intercept as given combines the standardized coefficients, whose
      # terms cancel for nearly collinear covariates, so its degrees of
      # freedom are computed through R^{-T} (.cr2_contrast_df()). A diagonal
      # map only divides: the contrast is the standardized intercept.
      intercept <- M["(Intercept)", , drop = FALSE]
      df <- .cr2_contrast_df(parts, drop(intercept))
      if (is.null(df)) {
        df <- suppressWarnings(as.data.frame(
          clubSandwich::linear_contrast(fit_cr, vcov = V, contrasts = intercept,
                                        test = "Satterthwaite")
        ))$df
      }
      tidy$df[tidy$term == "(Intercept)"] <- df
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
# clubSandwich object); NULL with neither a pretest nor a numeric covariate.
.covariate_scaling_element <- function(scaling, fit, V) {
  if (!length(scaling$means)) return(NULL)
  list(means = scaling$means, scales = scaling$scales, coefficients = stats::coef(fit),
       vcov = V)
}
