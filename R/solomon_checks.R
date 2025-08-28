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

  list(
    brown_forsythe_4cell_p = p_bf4,
    brown_forsythe_unpre_p = p_bf_un,
    shapiro_p_by_cell = shaps,
    ancova_slope_homogeneity_p = p_slope,
    notes = c(
      "Use Welch t for Groups 3 vs 4 if brown_forsythe_unpre_p < .05 (you already do).",
      "Prefer HC3 SEs (your default) if brown_forsythe_4cell_p < .05.",
      "If ancova_slope_homogeneity_p < .05, ANCOVA assumption is violated; prefer the unified GLM with interaction terms or a permutation p-value."
    )
  )
}
