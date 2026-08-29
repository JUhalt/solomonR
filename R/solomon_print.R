p_fmt <- function(p) {
  ifelse(is.na(p), "NA",
         ifelse(p < .001, "<.001", sprintf("%.3f", p)))
}
estse_str <- function(est, se, digits = 3) sprintf("%.*f (%.3f)", digits, est, se)

p_fmt <- function(p) {
  ifelse(is.na(p), "NA",
         ifelse(p < .001, "<.001", sprintf("%.3f", p)))
}
estse_str <- function(est, se, digits = 3) sprintf("%.*f (%.3f)", digits, est, se)

#' @export
print.solomon_glm <- function(x, digits = 3, ...) {
  # header + formula (without the environment garbage)
  f <- tryCatch(stats::formula(x$model), error = function(e) NULL)
  cat("Solomon GLM (unified model)\n")
  if (!is.null(f)) cat("Formula: ", paste(deparse(f), collapse = " "), "\n\n", sep = "")

  # --- Coefficients table ---
  ct <- x$coefficients[, c("term","estimate","std.error","statistic","p.value")]
  ct$`Est (SE)` <- estse_str(ct$estimate, ct$std.error, digits)
  ct$z <- sprintf("%.2f", ct$statistic)
  ct$p <- p_fmt(ct$p.value)

  # compute widths
  term_w  <- max(nchar("Term"),      nchar(ct$term))
  estse_w <- max(nchar("Est (SE)"),  nchar(ct$`Est (SE)`))
  z_w     <- max(nchar("z"),         nchar(ct$z))
  p_w     <- max(nchar("p"),         nchar(ct$p))

  cat(sprintf("%-*s  %-*s  %*s  %*s\n",
              term_w, "Term", estse_w, "Est (SE)", z_w, "z", p_w, "p"))
  for (i in seq_len(nrow(ct))) {
    cat(sprintf("%-*s  %-*s  %*s  %*s\n",
                term_w,  ct$term[i],
                estse_w, ct$`Est (SE)`[i],
                z_w,     ct$z[i],
                p_w,     ct$p[i]))
  }

  # --- Key contrasts table ---
  ef <- x$effects[, c(
    "contrast",
    "estimate",
    "std.error",
    "statistic",
    "p.value",
    "r2"
  )]
  ef$`Est (SE)` <- estse_str(ef$estimate, ef$std.error, digits)
  ef$z <- sprintf("%.2f", ef$statistic)
  ef$p <- p_fmt(ef$p.value)
  ef$R2 <- ifelse(is.na(ef$r2), "", sprintf("%.3f", ef$r2))

  con_w  <- max(nchar("Contrast"), nchar(ef$contrast))
  est2_w <- max(nchar("Est (SE)"), nchar(ef$`Est (SE)`))
  z2_w   <- max(nchar("z"),        nchar(ef$z))
  p2_w   <- max(nchar("p"),        nchar(ef$p))
  r2_w <- max(
    nchar("Wald R2"),
    nchar(ef$R2)
  )

  cat("\n")
  cat(
    sprintf(
      "%-*s  %-*s  %*s  %*s  %*s\n",
      con_w, "Key contrasts",
      est2_w, "Est (SE)",
      z2_w, "z",
      p2_w, "p",
      r2_w, "Wald R2"
    )
  )
  for (i in seq_len(nrow(ef))) {
    cat(sprintf("%-*s  %-*s  %*s  %*s  %*s\n",
                con_w,  ef$contrast[i],
                est2_w, ef$`Est (SE)`[i],
                z2_w,   ef$z[i],
                p2_w,   ef$p[i],
                r2_w,   ef$R2[i]))
  }
  invisible(x)
  cat(
    "\nWald R2: partial R-squared for conventional Gaussian OLS;\n",
    "a Wald-based descriptive approximation when robust covariance is used.\n",
    sep = ""
  )
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
print.summary.solomon_glm <- function(x, digits = 3, ...) {
  cat("Summary: Solomon GLM (unified model)\n\n")
  # mirror the same printing using x$coefs / x$effects
  ct <- x$coefs[, c("term","estimate","std.error","statistic","p.value")]
  ct$`Est (SE)` <- estse_str(ct$estimate, ct$std.error, digits)
  ct$z <- sprintf("%.2f", ct$statistic)
  ct$p <- p_fmt(ct$p.value)

  term_w  <- max(nchar("Term"),      nchar(ct$term))
  estse_w <- max(nchar("Est (SE)"),  nchar(ct$`Est (SE)`))
  z_w     <- max(nchar("z"),         nchar(ct$z))
  p_w     <- max(nchar("p"),         nchar(ct$p))

  cat(sprintf("%-*s  %-*s  %*s  %*s\n",
              term_w, "Term", estse_w, "Est (SE)", z_w, "z", p_w, "p"))
  for (i in seq_len(nrow(ct))) {
    cat(sprintf("%-*s  %-*s  %*s  %*s\n",
                term_w,  ct$term[i],
                estse_w, ct$`Est (SE)`[i],
                z_w,     ct$z[i],
                p_w,     ct$p[i]))
  }

  ef <- x$effects[, c("contrast","estimate","std.error","statistic","p.value","r2")]
  ef$`Est (SE)` <- estse_str(ef$estimate, ef$std.error, digits)
  ef$z <- sprintf("%.2f", ef$statistic)
  ef$p <- p_fmt(ef$p.value)
  ef$R2 <- ifelse(is.na(ef$r2), "", sprintf("%.3f", ef$r2))

  con_w  <- max(nchar("Contrast"), nchar(ef$contrast))
  est2_w <- max(nchar("Est (SE)"), nchar(ef$`Est (SE)`))
  z2_w   <- max(nchar("z"),        nchar(ef$z))
  p2_w   <- max(nchar("p"),        nchar(ef$p))
  r2_w   <- max(nchar("R2"),       nchar(ef$R2))

  cat("\n")
  cat(sprintf("%-*s  %-*s  %*s  %*s  %*s\n",
              con_w, "Key contrasts", est2_w, "Est (SE)", z2_w, "z", p2_w, "p", r2_w, "R2"))
  for (i in seq_len(nrow(ef))) {
    cat(sprintf("%-*s  %-*s  %*s  %*s  %*s\n",
                con_w,  ef$contrast[i],
                est2_w, ef$`Est (SE)`[i],
                z2_w,   ef$z[i],
                p2_w,   ef$p[i],
                r2_w,   ef$R2[i]))
  }
  invisible(x)
}

#' @export
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
