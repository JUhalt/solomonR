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
  inference <- if (is.null(x$inference)) "wald" else x$inference

  cat("Solomon full-information maximum-likelihood model\n")
  cat("-------------------------------------------------\n")
  cat("Method: van Engelenburg (1999)\n")
  cat(
    if (identical(inference, "wald")) {
      "Inference: Wald (large-sample normal reference)\n\n"
    } else {
      "Inference: small-sample (t; Welch-Satterthwaite df for combined contrasts)\n\n"
    }
  )

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

    statistic <- if (!is.null(e$df) && is.finite(e$df)) {
      sprintf("t(%s) = %.2f", .df_fmt(e$df), e$statistic)
    } else {
      sprintf("z = %.2f", e$statistic)
    }

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
        "%-28s %.*f (SE = %.*f), %s, p = %s%s\n",
        e$contrast,
        digits,
        e$estimate,
        digits,
        e$std.error,
        statistic,
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

  if (identical(inference, "wald") && isTRUE(x$small_sample)) {
    cat(
      .wrap_lines(
        sprintf(
          paste(
            "Note: the smallest cell has %d participants. Wald intervals can be",
            "too narrow in small samples; see inference = \"satterthwaite\" in",
            "?fit_solomon_ml."
          ),
          x$min_cell_n
        )
      ),
      "\n",
      sep = ""
    )
  }

  invisible(x)
}
