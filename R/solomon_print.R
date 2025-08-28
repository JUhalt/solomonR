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
print.solomon_classic <- function(x, digits = 3, ...) {
  line <- function(...) cat(sprintf(...), "\n", sep = "")
  p_fmt <- function(p) ifelse(p < .001, "<.001", sprintf("%.3f", p))

  line("Classic Solomon analysis")

  # ANOVA interaction summary
  if (!is.null(x$aov)) {
    a <- x$aov
    ai <- a[a$term %in% c("treat:pretested","factor(treat):factor(preind)","treat:pr…", "treat:pretested"), , drop=FALSE]
    ar <- a[a$term %in% c("Residuals"), , drop=FALSE]
    if (nrow(ai) == 1 && nrow(ar) == 1) {
      df1 <- ai$df; df2 <- ar$df
      F   <- ai$statistic; p <- ai$p.value
      line("")
      line("Interaction (Pretest × Treatment): F(%d, %d) = %.2f, p = %s",
           df1, df2, F, p_fmt(p))
    }
  }

  # ANCOVA line (pretested)
  if (!is.null(x$ancova)) {
    an <- x$ancova
    tr <- an[grepl("^treat", an$term), , drop=FALSE]
    if (nrow(tr) == 1) {
      line("Pretested (ANCOVA): beta_treat = %.*f (SE = %.3f), t = %.2f, p = %s",
           digits, tr$estimate, tr$std.error, tr$statistic, p_fmt(tr$p.value))
    }
  }

  # Welch t (unpretested)
  if (!is.null(x$t_unpretested)) {
    tt <- x$t_unpretested
    line("Unpretested (Welch t): t(%0.1f) = %.2f, p = %s; Δ = %.*f, 95%% CI [%.*f, %.*f]",
         tt$parameter, tt$statistic, p_fmt(tt$p.value),
         digits, tt$estimate, digits, tt$conf.low, digits, tt$conf.high)
  }

  # Hedges g
  if (!is.null(x$g_post)) {
    line("Effect size (Groups 3–4): Hedges g = %.3f, 95%% CI [%.3f, %.3f]",
         x$g_post["g"], x$g_post["lower"], x$g_post["upper"])
  }

  # Optional Stouffer
  if (!is.null(x$stouffer)) {
    z <- x$stouffer$z_meta
    p1 <- x$stouffer$p_meta_one_tailed
    line("Stouffer Z (one-tailed): Z = %.2f, p = %s", z, p_fmt(p1))
  }

  invisible(x)
}
