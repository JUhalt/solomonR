#' @export
print.solomon_sem_latent <- function(x, digits = 3, ...) {
  level <- if (is.null(x$settings$conf_level)) 0.95 else x$settings$conf_level

  cat("Solomon SEM (latent)\n")
  cat(sprintf("POST measurement invariance: %s\n", x$settings$invariance_post))
  cat("Latent mean reference: U0 (unpretested control) fixed at 0\n")
  if (isTRUE(x$settings$ancova))
    cat(sprintf("Pretested latent ANCOVA: %s (reference: P0)\n", x$settings$invariance_pre))
  cat("\n")

  # POST (4-group)
  fm <- x$fitmeasures_post
  cat(sprintf("POST model fit (4 groups): CFI=%.3f, RMSEA=%.3f, SRMR=%.3f; df=%d\n\n",
              fm["cfi"], fm["rmsea"], fm["srmr"], as.integer(fm["df"])))
  ef <- x$effects_post
  ef$`Est (SE)` <- estse_str(ef$estimate, ef$std.error, digits)
  ef$z <- ifelse(is.na(ef$statistic), "", sprintf("%.2f", ef$statistic))
  ef$p <- ifelse(is.na(ef$p.value),   "", p_fmt(ef$p.value))

  has_ci <- all(c("conf.low", "conf.high") %in% names(ef))
  ci_label <- if (has_ci) sprintf("%s%% CI", format(100 * level)) else ""
  ef$ci <- if (has_ci) {
    sprintf("[%.*f, %.*f]", digits, ef$conf.low, digits, ef$conf.high)
  } else {
    rep("", nrow(ef))
  }

  con_w  <- max(nchar("Key contrasts (latent POST)"), nchar(ef$contrast))
  est_w  <- max(nchar("Est (SE)"), nchar(ef$`Est (SE)`))
  z_w    <- max(nchar("z"), nchar(ef$z))
  p_w    <- max(nchar("p"), nchar(ef$p))
  ci_w   <- max(nchar(ci_label), nchar(ef$ci))
  cat(sprintf("%-*s  %-*s  %*s  %*s  %*s\n",
              con_w, "Key contrasts (latent POST)", est_w, "Est (SE)", z_w, "z", p_w, "p",
              ci_w, ci_label))
  for (i in seq_len(nrow(ef))) {
    cat(sprintf("%-*s  %-*s  %*s  %*s  %*s\n",
                con_w, ef$contrast[i],
                est_w, ef$`Est (SE)`[i],
                z_w,   ef$z[i],
                p_w,   ef$p[i],
                ci_w,  ef$ci[i]))
  }

  # PRETESTED ANCOVA (optional)
  if (!is.null(x$fit_pre)) {
    cat("\n")
    fm2 <- x$fitmeasures_pre
    cat(sprintf("Pretested latent ANCOVA fit: CFI=%.3f, RMSEA=%.3f, SRMR=%.3f; df=%d\n",
                fm2["cfi"], fm2["rmsea"], fm2["srmr"], as.integer(fm2["df"])))
    ef2 <- x$effects_pre
    if (!is.null(ef2)) {
      ci2 <- if (!is.null(ef2$conf.low)) {
        sprintf(", %s [%.*f, %.*f]", ci_label, digits, ef2$conf.low[1], digits, ef2$conf.high[1])
      } else {
        ""
      }
      cat(sprintf("  %s: %.*f%s\n", ef2$contrast[1], digits, ef2$estimate[1], ci2))
    }
  }
  invisible(x)
}

#' @export
summary.solomon_sem_latent <- function(object, digits = 3, ...) {
  print.solomon_sem_latent(object, digits = digits, ...)
  invisible(object)
}
