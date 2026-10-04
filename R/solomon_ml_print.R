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
      "Inference: Wald (van Engelenburg, 1999; large-sample normal reference)"
    } else {
      .wrap_lines(
        paste(
          "Inference: Satterthwaite (Satterthwaite, 1946; Welch, 1947; t with",
          "residual df within a pretest condition, Welch-Satterthwaite df for",
          "contrasts that combine them)"
        ),
        exdent = 2
      )
    },
    "\n\n",
    sep = ""
  )

  cat(
    sprintf(
      "Centered pretest mean: %.3f\n",
      x$pretest_mean
    )
  )

  # The ML residual SDs (SSE / n) are the ones Wald standard errors use;
  # Satterthwaite standard errors use the unbiased residual variances
  # (SSE / residual df), so both are shown. Fits saved before
  # `sigma_unbiased` existed show the ML values only.
  sd_text <- function(group) {
    ml <- sprintf("%.3f (ML", x$sigma[[group]])
    if (identical(inference, "wald")) {
      paste0(ml, "; used for SEs)")
    } else if (!is.null(x$sigma_unbiased)) {
      sprintf("%s); %.3f (from the unbiased variance; used for SEs)", ml,
              x$sigma_unbiased[[group]])
    } else {
      paste0(ml, ")")
    }
  }

  cat("Residual SD, unpretested: ", sd_text("unpretested"), "\n", sep = "")
  cat("Residual SD, pretested:   ", sd_text("pretested"), "\n\n", sep = "")

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

  # The pretest effects' standard errors include the variance of the mean
  # pretest (issue #104); older fits have no pretest-effect rows.
  if (any(x$effects$contrast %in% .solomon_pretest_order)) {
    cat("\n", .wrap_lines(paste(
      "Pretest effects: at the centered pretest mean, with standard errors that",
      "include its sampling variance."
    )), "\n", sep = "")
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
            "Note: the smallest cell has %d participants. In the package's",
            "simulation validation, Wald intervals were too narrow with fewer",
            "than %d participants per cell; the default,",
            "inference = \"satterthwaite\", was calibrated. See ?fit_solomon_ml."
          ),
          x$min_cell_n,
          .solomon_ml_small_cell
        )
      ),
      "\n",
      sep = ""
    )
  }

  invisible(x)
}
