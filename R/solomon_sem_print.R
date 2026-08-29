#' @export
print.solomon_sem <- function(x, digits = 3, ...) {

  p_fmt <- function(p) {
    if (is.na(p)) {
      "NA"
    } else if (p < .001) {
      "<.001"
    } else {
      sprintf("%.3f", p)
    }
  }

  # ----------------------------------------------------------
  # Header
  # ----------------------------------------------------------

  if (identical(x$mode, "mean")) {
    cat("Solomon SEM (multi-group mean structure)\n")
  } else if (identical(x$mode, "ancova_pretested")) {
    cat("Solomon SEM (ANCOVA in pretested groups)\n")
  } else {
    cat("Solomon SEM\n")
  }

  # ----------------------------------------------------------
  # Model fit
  # ----------------------------------------------------------

  fm <- x$fitmeasures

  if (!is.null(fm)) {

    df_fit <- unname(fm["df"])

    if (
      identical(x$mode, "mean") &&
      length(df_fit) == 1L &&
      is.finite(df_fit) &&
      df_fit == 0
    ) {

      cat(
        "Fit: saturated four-group mean structure (df = 0)\n",
        "Global CFI/RMSEA/SRMR are not diagnostic for this model.\n\n",
        sep = ""
      )

    } else {

      cat(
        sprintf(
          "Fit: CFI=%.3f, RMSEA=%.3f, SRMR=%.3f; df=%.0f\n\n",
          unname(fm["cfi"]),
          unname(fm["rmsea"]),
          unname(fm["srmr"]),
          df_fit
        )
      )
    }
  }

  # ----------------------------------------------------------
  # Key contrasts
  # ----------------------------------------------------------

  ef <- x$effects

  ef$`Est (SE)` <- sprintf(
    "%.*f (%.*f)",
    digits,
    ef$estimate,
    digits,
    ef$std.error
  )

  ef$z <- sprintf(
    "%.2f",
    ef$statistic
  )

  ef$p <- vapply(
    ef$p.value,
    p_fmt,
    character(1)
  )

  contrast_w <- max(
    nchar("Key contrasts"),
    nchar(ef$contrast)
  )

  est_w <- max(
    nchar("Est (SE)"),
    nchar(ef$`Est (SE)`)
  )

  z_w <- max(
    nchar("z"),
    nchar(ef$z)
  )

  p_w <- max(
    nchar("p"),
    nchar(ef$p)
  )

  cat(
    sprintf(
      "%-*s  %-*s  %*s  %*s\n",
      contrast_w,
      "Key contrasts",
      est_w,
      "Est (SE)",
      z_w,
      "z",
      p_w,
      "p"
    )
  )

  for (i in seq_len(nrow(ef))) {

    cat(
      sprintf(
        "%-*s  %-*s  %*s  %*s\n",
        contrast_w,
        ef$contrast[i],
        est_w,
        ef$`Est (SE)`[i],
        z_w,
        ef$z[i],
        p_w,
        ef$p[i]
      )
    )
  }

  invisible(x)
}
