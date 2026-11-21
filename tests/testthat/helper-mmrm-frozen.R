# The internal fitting function of fit_solomon_mmrm() as it was before the
# arguments `structures` and `vcov` were added to it (issue #57). The third
# amendment of the protocol on that issue says of the two arguments that
# their "defaults reproduce the current behavior; a unit test checks that".
# Calling the function with the defaults named would only show that the
# arguments are passed through, so the test in test-mmrm.R compares the
# current function, at its defaults, with this frozen copy.
#
# mmrm_fit_contrasts_70022ef() is .mmrm_fit_contrasts() of commit
# 70022ef29915c0206ea2732d93567764a864454e (R/solomon_mmrm.R, lines 39 to
# 177), line for line, without its comments. It uses that commit's list of
# fallback structures, kept below, and the package's own .welch_df() and
# .effects_table(), which the two arguments do not touch. It is a record:
# it does not change when the package does. If the default fit is changed
# on purpose later, the test that fails is the place to say so.

mmrm_fit_contrasts_70022ef <- function(observed, occasions, pre, grouped, df, conf_level,
                                       center = NULL, other = NULL) {
  empty <- setdiff(occasions, as.character(observed$occ))
  if (length(empty)) {
    stop("A Solomon contrast cannot be estimated because no analyzed participant has a ",
         "posttest at these occasions: ", paste(empty, collapse = ", "), ". Remove their ",
         "rows, or check that every group is observed at every occasion.", call. = FALSE)
  }
  fixed <- if (pre) "y ~ occ * treat * pretested + occ:pre_obs" else "y ~ occ * treat * pretested"
  method <- if (df == "kenward-roger") "Kenward-Roger" else "Satterthwaite"
  fit_one <- function(type) {
    f <- stats::as.formula(sprintf("%s + %s(occ | %sid)", fixed, type,
                                   if (grouped) "pgrp / " else ""))
    fit <- tryCatch(
      if (df == "kenward-roger") {
        mmrm::mmrm(f, data = observed, reml = TRUE, method = method,
                   vcov = "Kenward-Roger-Linear")
      } else {
        mmrm::mmrm(f, data = observed, reml = TRUE, method = method)
      },
      error = function(e) NULL
    )
    if (!is.null(fit) && !identical(as.integer(mmrm::component(fit, "convergence")), 0L)) fit <- NULL
    fit
  }

  structure_used <- "unstructured"
  fit <- fit_one("us")
  if (is.null(fit)) {
    fits <- lapply(names(.mmrm_fallbacks), fit_one)
    ok <- !vapply(fits, is.null, logical(1))
    if (!any(ok)) {
      stop("The model did not converge with an unstructured covariance or any of the ",
           "fallback structures (", paste(.mmrm_fallbacks, collapse = ", "), ").", call. = FALSE)
    }
    aic <- vapply(fits[ok], stats::AIC, numeric(1))
    best <- which(ok)[which.min(aic)]
    fit <- fits[[best]]
    structure_used <- unname(.mmrm_fallbacks[best])
  }

  beta <- mmrm::component(fit, "beta_est")
  tt <- stats::delete.response(stats::terms(stats::as.formula(fixed)))
  contrasts <- mmrm::component(fit, "contrasts")
  point <- function(t, tr, pr, pre_obs = 0) {
    nd <- data.frame(occ = factor(occasions[t], levels = occasions), treat = tr, pretested = pr,
                     pre_obs = pre_obs)
    stats::model.matrix(tt, nd, contrasts.arg = contrasts)[1, ]
  }
  weights <- function(v) {
    missing_terms <- setdiff(names(v)[abs(v) > 0], names(beta))
    if (length(missing_terms)) {
      stop("A Solomon contrast cannot be estimated because these terms are not estimable: ",
           paste(missing_terms, collapse = ", "), ". Check that every group is observed at ",
           "every occasion.", call. = FALSE)
    }
    out <- stats::setNames(numeric(length(beta)), names(beta))
    out[intersect(names(v), names(beta))] <- v[intersect(names(v), names(beta))]
    out
  }
  one <- function(L, occasion, contrast, slope = 0, terms = center) {
    r <- mmrm::df_1d(fit, weights(L))
    est <- r$est
    se <- r$se
    dfc <- r$df
    stat <- r$t_stat
    p <- r$p_val
    if (!is.null(terms) && slope != 0) {
      v_center <- slope^2 * terms$variance
      dfc <- .welch_df(c(se^2, v_center), c(dfc, terms$df))
      se <- sqrt(se^2 + v_center)
      stat <- est / se
      p <- 2 * stats::pt(-abs(stat), dfc)
    }
    crit <- stats::qt(1 - (1 - conf_level) / 2, dfc)
    data.frame(occasion = occasion, contrast = contrast, estimate = est, std.error = se,
               statistic = stat, df = dfc, p.value = p,
               conf.low = est - crit * se, conf.high = est + crit * se,
               stringsAsFactors = FALSE)
  }
  sens <- list()
  rows <- list()
  other_rows <- list()
  slopes <- stats::setNames(numeric(length(occasions)), occasions)
  for (t in seq_along(occasions)) {
    pre_t <- point(t, 1, 1) - point(t, 0, 1)
    un_t <- point(t, 1, 0) - point(t, 0, 0)
    sens[[t]] <- pre_t - un_t
    control_t <- point(t, 0, 1) - point(t, 0, 0)
    treated_t <- point(t, 1, 1) - point(t, 1, 0)
    slope <- if (pre) sum(weights(point(t, 0, 1, pre_obs = 1) - point(t, 0, 1)) * beta) else 0
    slopes[t] <- slope
    rows <- c(rows, list(
      one((pre_t + un_t) / 2, occasions[t], "ATE (avg over pretest)"),
      one(sens[[t]], occasions[t], "Pretest x Treatment"),
      one(pre_t, occasions[t], "Treatment | pretested"),
      one(un_t, occasions[t], "Treatment | unpretested"),
      one(control_t, occasions[t], "Pretest effect | control", slope),
      one(treated_t, occasions[t], "Pretest effect | treated", slope),
      one((control_t + treated_t) / 2, occasions[t], "Pretest main effect", slope)
    ))
    if (!is.null(other)) {
      control_o <- point(t, 0, 1, pre_obs = other$shift) - point(t, 0, 0)
      treated_o <- point(t, 1, 1, pre_obs = other$shift) - point(t, 1, 0)
      other_rows <- c(other_rows, list(
        one(control_o, occasions[t], "Pretest effect | control", slope, other$center),
        one(treated_o, occasions[t], "Pretest effect | treated", slope, other$center),
        one((control_o + treated_o) / 2, occasions[t], "Pretest main effect", slope,
            other$center)
      ))
    }
  }
  last <- length(occasions)
  rows <- c(rows, list(one(sens[[last]] - sens[[1]],
                           sprintf("%s vs %s", occasions[last], occasions[1]),
                           "Change in Pretest x Treatment")))
  list(effects = .effects_table(do.call(rbind, rows), keys = "occasion"), model = fit,
       covariance = structure_used, slopes = if (pre) slopes,
       other = if (length(other_rows)) .effects_table(do.call(rbind, other_rows), keys = "occasion"))
}

# The copy sees the fallback structures of 70022ef, and everything else as
# the package's function does.
environment(mmrm_fit_contrasts_70022ef) <- local({
  .mmrm_fallbacks <- c(toeph = "heterogeneous Toeplitz", ar1h = "heterogeneous AR(1)",
                       csh = "heterogeneous compound symmetry")
  environment()
}, envir = new.env(parent = asNamespace("solomonR")))
