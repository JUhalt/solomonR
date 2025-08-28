# Brown–Forsythe (median-based Levene) without new deps
bf_test <- function(y, group) {
  group <- factor(group)
  med <- tapply(y, group, stats::median, na.rm = TRUE)
  z <- abs(y - med[group])
  fit <- stats::lm(z ~ group)
  a <- stats::anova(fit)
  a$`Pr(>F)`[1]
}

#' Assumption checks for Solomon analyses
#' @param y_post numeric posttest
#' @param treat 0/1
#' @param pretested 0/1
#' @param y_pre numeric pretest (NA for unpretested)
#' @export
check_solomon_assumptions <- function(y_post, treat, pretested, y_pre) {
  df <- data.frame(y_post, treat=factor(treat), pretested=factor(pretested))
  df$cell <- interaction(df$pretested, df$treat, drop = TRUE)

  # (1) Heteroscedasticity across the four posttest cells (Brown–Forsythe)
  p_bf4 <- bf_test(df$y_post, df$cell)

  # (2) Heteroscedasticity in posttest-only cells (Groups 3 vs 4) for Welch t transparency
  p_bf_un <- bf_test(df$y_post[df$pretested==0], df$treat[df$pretested==0])

  # (3) Normality (Shapiro) within each posttest cell (warn if n<3)
  shaps <- tapply(df$y_post, df$cell, function(v) if (sum(is.finite(v)) >= 3) stats::shapiro.test(v)$p.value else NA_real_)

  # (4) Homogeneity of regression slopes for ANCOVA (pretested cells only):
  # y_post ~ treat * y_pre  (if interaction significant -> slope heterogeneity)
  p_slope <- NA_real_
  if (any(pretested==1) && all(is.finite(y_pre[pretested==1]))) {
    dd <- data.frame(y = y_post[pretested==1], treat = factor(treat[pretested==1]), pre = y_pre[pretested==1])
    fit <- stats::lm(y ~ treat * pre, data = dd)
    a <- stats::anova(fit)
    # last row is interaction (treat:pre)
    p_slope <- a$`Pr(>F)`[which(rownames(a) == "treat:pre")][1]
  }

  structure(list(
    brown_forsythe_4cell_p = p_bf4,
    brown_forsythe_unpre_p = p_bf_un,
    shapiro_p_by_cell = shaps,
    ancova_slope_homogeneity_p = p_slope
  ), class = "solomon_checks")
}

#' @export
print.solomon_checks <- function(x, ...) {
  line <- function(...) cat(sprintf(...), "\n", sep = "")
  p_fmt <- function(p) ifelse(is.na(p), "NA", ifelse(p < .001, "<.001", sprintf("%.3f", p)))
  ok <- function(flag) ifelse(flag, "OK", "FLAG")

  line("Assumption checks (α = .05)")
  bf4_ok  <- !is.na(x$brown_forsythe_4cell_p) && x$brown_forsythe_4cell_p >= .05
  bfu_ok  <- !is.na(x$brown_forsythe_unpre_p) && x$brown_forsythe_unpre_p >= .05
  shp_min <- suppressWarnings(min(x$shapiro_p_by_cell, na.rm = TRUE))
  shp_ok  <- is.finite(shp_min) && shp_min >= .05
  slp_ok  <- !is.na(x$ancova_slope_homogeneity_p) && x$ancova_slope_homogeneity_p >= .05

  line("  HoV across 4 posttest cells (Brown–Forsythe): p = %s  -> %s", p_fmt(x$brown_forsythe_4cell_p), ok(bf4_ok))
  line("  HoV in unpretested cells (Welch target):       p = %s  -> %s", p_fmt(x$brown_forsythe_unpre_p),  ok(bfu_ok))
  line("  Normality by cell (Shapiro, min p):            p = %s  -> %s", p_fmt(shp_min),                   ok(shp_ok))
  line("  ANCOVA slope homogeneity (pretested):          p = %s  -> %s", p_fmt(x$ancova_slope_homogeneity_p), ok(slp_ok))

  line("")
  line("Recommendations:")
  line("  • Robust SEs (HC3): %s", ifelse(bf4_ok, "fine", "recommended"))
  line("  • Welch t for groups 3–4: %s", ifelse(bfu_ok, "fine", "recommended (package uses Welch)"))
  line("  • Permutation p-values: %s", ifelse(all(bf4_ok, bfu_ok, shp_ok, slp_ok), "optional", "consider"))
  invisible(x)
}
