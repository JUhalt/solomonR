p_fmt <- function(p) ifelse(is.na(p), "NA", ifelse(p < .001, "<.001", sprintf("%.3f", p)))
estse_str <- function(est, se, digits = 3) sprintf("%.*f (%.3f)", digits, est, se)

#' @export
print.solomon_sem <- function(x, digits = 3, ...) {
  if (isTRUE(x$mode == "ancova_pretested")) {
    cat("Solomon SEM (ANCOVA in pretested groups)\n")
  } else {
    cat("Solomon SEM (multi-group mean structure)\n")
  }
  fm <- x$fitmeasures
  cat(sprintf("Fit: CFI=%.3f, RMSEA=%.3f, SRMR=%.3f; df=%d\n\n",
              fm["cfi"], fm["rmsea"], fm["srmr"], as.integer(fm["df"])))

  ef <- x$effects
  ef$`Est (SE)` <- estse_str(ef$estimate, ef$std.error, digits)
  ef$z <- sprintf("%.2f", ef$statistic)
  ef$p <- p_fmt(ef$p.value)

  con_w  <- max(nchar("Contrast"), nchar(ef$contrast))
  est_w  <- max(nchar("Est (SE)"), nchar(ef$`Est (SE)`))
  z_w    <- max(nchar("z"),        nchar(ef$z))
  p_w    <- max(nchar("p"),        nchar(ef$p))

  cat(sprintf("%-*s  %-*s  %*s  %*s\n",
              con_w, "Key contrasts", est_w, "Est (SE)", z_w, "z", p_w, "p"))
  for (i in seq_len(nrow(ef))) {
    cat(sprintf("%-*s  %-*s  %*s  %*s\n",
                con_w, ef$contrast[i], est_w, ef$`Est (SE)`[i],
                z_w, ef$z[i], p_w, ef$p[i]))
  }
  invisible(x)
}

#' @export
summary.solomon_sem <- function(object, digits = 3, ...) {
  print.solomon_sem(object, digits = digits, ...)
  invisible(object)
}
