#' @export
print.solomon_glm <- function(x, digits = 3, ...) {
  cat("Solomon GLM (unified model)\n")
  cat("Formula: "); print(stats::formula(x$model))
  cat("\nCoefficients (robust SEs):\n")
  ct <- x$coefficients[, c("term","estimate","std.error","statistic","p.value")]
  print.data.frame(transform(ct,
                             estimate = round(estimate, digits),
                             std.error = round(std.error, digits),
                             statistic = round(statistic, digits),
                             p.value = signif(p.value, 3)), row.names = FALSE)
  cat("\nKey contrasts:\n")
  ef <- x$effects[, c("contrast","estimate","std.error","statistic","p.value")]
  print.data.frame(transform(ef,
                             estimate = round(estimate, digits),
                             std.error = round(std.error, digits),
                             statistic = round(statistic, digits),
                             p.value = signif(p.value, 3)), row.names = FALSE)
  invisible(x)
}

#' @export
print.solomon_classic <- function(x, digits = 3, ...) {
  cat("Classic Solomon analysis\n")
  cat("\n2x2 ANOVA on posttest (tests sensitization via interaction):\n")
  if (!is.null(x$aov)) print(x$aov)
  if (!is.null(x$ancova)) {
    cat("\nANCOVA on pretested cells (Groups 1 & 2):\n")
    print(x$ancova)
  }
  if (!is.null(x$t_unpretested)) {
    cat("\nPosttest-only Welch t (Groups 3 & 4):\n")
    print(x$t_unpretested)
  }
  if (!is.null(x$stouffer)) {
    cat("\nStouffer Z (optional combine of ANCOVA + posttest-only):\n")
    print(x$stouffer)
  }
  invisible(x)
}
