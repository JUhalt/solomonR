#' @export
summary.solomon_glm <- function(object, digits = 3, ...) {
  res <- list(
    call       = object$call,
    formula    = stats::formula(object$model),
    robust     = object$robust,
    cluster    = object$cluster,
    n_clusters = object$n_clusters,
    conf_level = object$conf_level,
    coefs      = object$coefficients,
    effects    = object$effects,
    pretest_mean = object$pretest_mean,
    family     = object$family
  )
  # The two pretest centers and the attrition, for the lines printed when a
  # posttest is missing (issue #135). Fits made by earlier versions have
  # none of them.
  if (!is.null(object$pretest_center)) {
    res$pretest_center <- object$pretest_center
    res$pretest_mean_analyzed <- object$pretest_mean_analyzed
    res$pretest_mean_observed <- object$pretest_mean_observed
    res$pretest_n <- object$pretest_n
    res$pretest_sensitivity <- object$pretest_sensitivity
    res$attrition <- .fit_attrition(object)
  }
  class(res) <- "summary.solomon_glm"
  res
}

#' @export
print.summary.solomon_glm <- function(x, digits = 3, ...) {
  level <- if (is.null(x$conf_level)) 0.95 else x$conf_level
  cat("Summary: Solomon GLM (unified model)\n")
  cat("Formula: ", .formula_line(x$formula), "\n", sep = "")
  cat("Covariance: ", .solomon_vcov_label(x), "\n\n", sep = "")
  .solomon_glm_tables(x$coefs, x$effects, digits, level)
  if (is.finite(.null_na(x$pretest_mean))) cat("\n")
  .glm_pretest_note(x, digits)
  .center_output(x, digits)
  invisible(x)
}

#' @export
summary.solomon_classic <- function(object, ...) {
  class(object) <- c("summary.solomon_classic", setdiff(class(object), "summary.solomon_classic"))
  object
}

#' @export
print.summary.solomon_classic <- function(x, digits = 3, ...) {
  print.solomon_classic(x, digits = digits, ...)
}
