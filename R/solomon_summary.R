#' @export
summary.solomon_glm <- function(object, digits = 3, ...) {
  res <- list(
    call    = object$call,
    formula = stats::formula(object$model),
    robust  = object$robust,
    cluster = object$cluster,
    coefs   = object$coefficients,
    effects = object$effects
  )
  class(res) <- "summary.solomon_glm"
  res
}

#' @export
print.summary.solomon_glm <- function(x, digits = 3, ...) {
  cat("Summary: Solomon GLM (unified model)\n")
  cat("Formula: ", paste(deparse(x$formula), collapse = " "), "\n", sep = "")
  cat("Covariance: ", .solomon_vcov_label(x), "\n\n", sep = "")
  .solomon_glm_tables(x$coefs, x$effects, x$robust, digits)
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
