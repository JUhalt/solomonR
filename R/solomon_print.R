# ---- shared formatting helpers ----

p_fmt <- function(p) {
  ifelse(is.na(p), "NA",
         ifelse(p < .001, "<.001", sprintf("%.3f", p)))
}

estse_str <- function(est, se, digits = 3) {
  ifelse(
    is.na(se),
    sprintf("%.*f (NA)", digits, est),
    sprintf("%.*f (%.3f)", digits, est, se)
  )
}

# Describe the covariance estimator and reference distribution of a
# solomon_glm fit or its summary.
.solomon_vcov_label <- function(x) {
  switch(
    x$robust,
    none = "model-based (conventional); normal-reference tests",
    HC3 = "HC3 heteroskedasticity-consistent; normal-reference tests",
    CR2 = sprintf(
      "CR2 cluster-robust (%d clusters); Satterthwaite t tests",
      length(unique(stats::na.omit(x$cluster)))
    )
  )
}

# Print the coefficient and key-contrast tables of a solomon_glm fit.
.solomon_glm_tables <- function(coefs, effects, robust, digits = 3) {

  show_df <- identical(robust, "CR2")
  stat_label <- if (show_df) "t" else "z"

  print_table <- function(tab, first_col, first_label, r2 = FALSE) {

    cols <- list(
      as.character(tab[[first_col]]),
      estse_str(tab$estimate, tab$std.error, digits),
      sprintf("%.2f", tab$statistic)
    )
    heads <- c(first_label, "Est (SE)", stat_label)

    if (show_df && "df" %in% names(tab)) {
      cols <- c(cols, list(sprintf("%.1f", tab$df)))
      heads <- c(heads, "df")
    }

    cols <- c(cols, list(p_fmt(tab$p.value)))
    heads <- c(heads, "p")

    if (r2 && "r2" %in% names(tab)) {
      cols <- c(cols, list(ifelse(is.na(tab$r2), "", sprintf("%.3f", tab$r2))))
      heads <- c(heads, "Wald R2")
    }

    widths <- mapply(function(h, v) max(nchar(h), nchar(v)), heads, cols)
    left <- seq_along(heads) <= 2L

    print_row <- function(values) {
      cells <- ifelse(
        left,
        sprintf("%-*s", widths, values),
        sprintf("%*s", widths, values)
      )
      cat(paste(cells, collapse = "  "), "\n", sep = "")
    }

    print_row(heads)
    for (i in seq_len(nrow(tab))) {
      print_row(vapply(cols, `[`, character(1), i))
    }
  }

  print_table(coefs, "term", "Term")
  cat("\n")
  print_table(effects, "contrast", "Key contrasts", r2 = TRUE)
}

#' @export
print.solomon_glm <- function(x, digits = 3, ...) {
  # header + formula (without the environment garbage)
  f <- tryCatch(stats::formula(x$model), error = function(e) NULL)
  cat("Solomon GLM (unified model)\n")
  if (!is.null(f)) cat("Formula: ", paste(deparse(f), collapse = " "), "\n", sep = "")
  cat("Covariance: ", .solomon_vcov_label(x), "\n\n", sep = "")

  .solomon_glm_tables(x$coefficients, x$effects, x$robust, digits)

  cat(
    "\nWald R2: partial R-squared for conventional Gaussian OLS;\n",
    "a Wald-based descriptive approximation when robust covariance is used.\n",
    sep = ""
  )
  invisible(x)
}

#' @export
print.solomon_classic <- function(x, digits = 3, ...) {

  p_fmt <- function(p) {
    if (is.na(p)) {
      "NA"
    } else if (p < .001) {
      "<.001"
    } else {
      sprintf("%.3f", p)
    }
  }

  show_F <- function(letter) {

    z <- x$tests[[letter]]$result

    marker <- if (letter %in% x$path) {
      "[PATH]"
    } else {
      "      "
    }

    cat(
      sprintf(
        "%s Test %s: %-45s F(1, %.0f) = %.2f, p = %s\n",
        marker,
        letter,
        x$tests[[letter]]$label,
        z$df,
        z$F,
        p_fmt(z$p.value)
      )
    )
  }

  cat("Classic Solomon analysis (historical teaching workflow)\n")
  cat("-------------------------------------------------------\n")

  cat(
    sprintf(
      "Selected pretested-group method: Test %s (%s)\n",
      x$settings$selected_test,
      x$settings$pretested_test
    )
  )

  cat(
    sprintf(
      "Historical decision path: %s\n\n",
      x$path_string
    )
  )

  cat("Historical Tests A-I\n")
  cat("--------------------\n")

  cat(
    "All tests are shown below. Tests marked [PATH] were reached by\n",
    "the historical decision sequence for these data.\n\n",
    sep = ""
  )

  show_F("A")
  show_F("B")
  show_F("C")
  show_F("D")
  show_F("E")
  show_F("F")
  show_F("G")

  h <- x$tests$H$result

  h_marker <- if ("H" %in% x$path) {
    "[PATH]"
  } else {
    "      "
  }

  cat(
    sprintf(
      "%s Test H: %-45s t(%.0f) = %.2f, p = %s\n",
      h_marker,
      x$tests$H$label,
      h$df,
      h$statistic,
      p_fmt(h$p.value)
    )
  )

  if (!is.null(x$tests$I$result)) {

    ii <- x$tests$I$result

    i_marker <- if ("I" %in% x$path) {
      "[PATH]"
    } else {
      "      "
    }

    cat(
      sprintf(
        "%s Test I: %-45s Z = %.2f, p(one-tailed) = %s\n",
        i_marker,
        ii$source,
        ii$z,
        p_fmt(ii$p.value)
      )
    )
  }

  cat("\nHistorical interpretation\n")
  cat("-------------------------\n")
  cat(x$conclusion, "\n")

  if (!is.null(x$g_post)) {
    cat(
      sprintf(
        "\nGroups 3-4 effect size: Hedges g = %.3f, 95%% CI [%.3f, %.3f]\n",
        x$g_post["g"],
        x$g_post["lower"],
        x$g_post["upper"]
      )
    )
  }

  if (isTRUE(x$settings$combine_with_stouffer)) {
    cat(
      "\nCaution: Test I is reproduced for historical teaching and replication.\n",
      "Later simulation work raised concerns about Type I error for the\n",
      "conditional meta-analytic sequence; it is not the default modern\n",
      "inferential recommendation in solomonR.\n",
      sep = ""
    )
  }

  invisible(x)
}

#' @export
print.solomon_perm <- function(x, digits = 2, ...) {

  format_p <- function(p) {

    if (is.na(p)) {
      return("NA")
    }

    if (p < .001) {
      return("< .001")
    }

    sub(
      "^0",
      "",
      sprintf("%.3f", p)
    )
  }

  cat("Solomon randomization test\n")
  cat("--------------------------\n")

  cat(
    sprintf(
      "Contrast: %s\n",
      x$contrast
    )
  )

  cat(
    sprintf(
      "Observed studentized statistic: z = %.*f\n",
      digits,
      x$z_obs
    )
  )

  cat(
    sprintf(
      "Permutation p = %s\n",
      format_p(x$p_perm)
    )
  )

  cat(
    sprintf(
      "Valid permutations: %d of %d\n",
      x$valid_reps,
      x$reps
    )
  )

  if (!is.null(x$z_perm)) {
    cat(
      sprintf(
        "Permutation distribution retained (%d draws); use plot_perm() to visualize it.\n",
        length(x$z_perm)
      )
    )
  }

  invisible(x)
}
