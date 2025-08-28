#' @export
summary.solomon_glm <- function(object, digits = 3, ...) {
  res <- list(
    call   = object$model$call,
    formula= stats::formula(object$model),
    coefs  = object$coefficients[, c("term","estimate","std.error","statistic","p.value")],
    effects= object$effects[,  c("contrast","estimate","std.error","statistic","p.value","r2","r2_lo","r2_hi")]
  )
  class(res) <- "summary.solomon_glm"
  res
}

#' @export
print.summary.solomon_glm <- function(x, digits = 3, ...) {
  cat("Summary: Solomon GLM (unified model)\n")
  cat("Formula: "); print(x$formula)
  cat("\nCoefficients (robust SEs):\n")
  ct <- x$coefs
  ct$estimate <- round(ct$estimate, digits)
  ct$std.error <- round(ct$std.error, digits)
  ct$statistic <- round(ct$statistic, digits)
  ct$p.value <- signif(ct$p.value, 3)
  print.data.frame(ct, row.names = FALSE)
  cat("\nKey contrasts:\n")
  ef <- x$effects
  ef$estimate <- round(ef$estimate, digits)
  ef$std.error <- round(ef$std.error, digits)
  ef$statistic <- round(ef$statistic, digits)
  ef$p.value <- signif(ef$p.value, 3)
  if ("r2" %in% names(ef)) {
    ef$r2    <- round(ef$r2, digits)
    ef$r2_lo <- round(ef$r2_lo, digits)
    ef$r2_hi <- round(ef$r2_hi, digits)
  }
  print.data.frame(ef, row.names = FALSE)
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
