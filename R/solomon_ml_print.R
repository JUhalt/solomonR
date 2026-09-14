#' @export
print.solomon_ml <- function(x, digits = 3, ...) {

  p_fmt <- function(p) {
    ifelse(
      p < .001,
      "<.001",
      sprintf("%.3f", p)
    )
  }

  level <- if (is.null(x$conf_level)) 0.95 else x$conf_level

  cat("Solomon full-information maximum-likelihood model\n")
  cat("-------------------------------------------------\n")
  cat("Method: van Engelenburg (1999)\n\n")

  cat(
    sprintf(
      "Centered pretest mean: %.3f\n",
      x$pretest_mean
    )
  )

  cat(
    sprintf(
      "Residual SD, unpretested: %.3f\n",
      x$sigma["unpretested"]
    )
  )

  cat(
    sprintf(
      "Residual SD, pretested:   %.3f\n\n",
      x$sigma["pretested"]
    )
  )

  cat("Key Solomon estimands\n")
  cat("---------------------\n")

  for (i in seq_len(nrow(x$effects))) {

    e <- x$effects[i, ]

    ci <- if (!is.null(e$conf.low)) {
      sprintf(
        ", %s%% CI [%.*f, %.*f]",
        format(100 * level),
        digits,
        e$conf.low,
        digits,
        e$conf.high
      )
    } else {
      ""
    }

    cat(
      sprintf(
        "%-28s %.*f (SE = %.*f), z = %.2f, p = %s%s\n",
        e$contrast,
        digits,
        e$estimate,
        digits,
        e$std.error,
        e$statistic,
        p_fmt(e$p.value),
        ci
      )
    )
  }

  cat(
    sprintf(
      "\nlogLik = %.2f; optimizer convergence = %d\n",
      x$logLik,
      x$convergence
    )
  )

  cat("Tests and intervals use large-sample normal reference distributions.\n")

  invisible(x)
}
