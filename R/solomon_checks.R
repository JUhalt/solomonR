# Brown-Forsythe (median-based Levene) without new deps
bf_test <- function(y, group) {
  group <- factor(group)
  med <- tapply(y, group, stats::median, na.rm = TRUE)
  z <- abs(y - med[group])
  fit <- stats::lm(z ~ group)
  a <- stats::anova(fit)
  a$`Pr(>F)`[1]
}

#' Descriptive assumption diagnostics for Solomon analyses
#'
#' `r lifecycle::badge("stable")`
#' Reports Brown-Forsythe tests of equal posttest variance (Brown & Forsythe,
#' 1974), Shapiro-Wilk normality tests within cells, and a test of
#' homogeneous pretest-posttest slopes among pretested participants with
#' complete scores.
#'
#' These results are descriptive. Choosing an analysis according to whether
#' a preliminary assumption test is significant can distort Type I error
#' rates (Zimmerman, 2004), so solomonR does not use them as gates:
#' heteroskedasticity-consistent standard errors are a reasonable default for
#' the unified GLM regardless of these results (Long & Ervin, 2000). A small
#' slope-homogeneity p-value is substantively informative: it suggests that
#' the treatment effect among pretested participants depends on the pretest
#' score.
#'
#' Designs with several treatments: give `treat` as a factor or character
#' vector of conditions and name the control with `control`. The tests then
#' cover all 2(k + 1) cells of a design with k treatments (Steyn, 2009): the
#' Brown-Forsythe tests compare all the posttest cells and the k + 1
#' unpretested groups, and the slope-homogeneity test is a test on k degrees
#' of freedom that the pretest-posttest slope is the same in the k + 1
#' pretested groups.
#'
#' @param y_post numeric posttest
#' @param treat 0/1 (or logical) treatment indicator; or a factor or
#'   character vector of conditions, with the control named by `control`
#' @param pretested 0/1 (or logical) pretest indicator
#' @param y_pre numeric pretest (NA for unpretested)
#' @param control The control condition when `treat` is a factor or
#'   character vector with more than two conditions, a Solomon N-group
#'   design; see [fit_solomon_glm()]. With two conditions the result is the
#'   same as with a 0/1 `treat`.
#' @param data Optional data frame. When supplied, the other data arguments
#'   are looked up in it first, as bare column names (`y_post = post`) or as
#'   strings (`y_post = "post"`).
#' @return An object of class `solomon_checks` with the p-values
#'   `brown_forsythe_4cell_p` (all four posttest cells),
#'   `brown_forsythe_unpre_p` (the unpretested cells), `shapiro_p_by_cell`
#'   (one p-value for each cell, named by the pretest indicator and the
#'   treatment, such as `1.0`), and `ancova_slope_homogeneity_p`. For a
#'   design with several treatments, the test across all the posttest cells
#'   is `brown_forsythe_cells_p`, the cells of `shapiro_p_by_cell` are named
#'   by the pretest indicator and the condition, such as `1.RP`, and
#'   `conditions` names the control and the treatments.
#' @references
#' Brown, M. B., & Forsythe, A. B. (1974). Robust tests for the equality of
#' variances. *Journal of the American Statistical Association, 69*(346),
#' 364–367. https://doi.org/10.1080/01621459.1974.10482955
#'
#' Long, J. S., & Ervin, L. H. (2000). Using heteroscedasticity consistent
#' standard errors in the linear regression model. *The American Statistician,
#' 54*(3), 217–224. https://doi.org/10.1080/00031305.2000.10474549
#'
#' Mai, N. N., Takahashi, Y., & Oo, M. M. (2020). Testing the effectiveness of
#' transfer interventions using Solomon four-group designs. *Education
#' Sciences, 10*(4), Article 92. https://doi.org/10.3390/educsci10040092
#'
#' Steyn, R. (2009). Re-designing the Solomon four-group: Can we improve on
#' this exemplary model? *Design Principles and Practices: An International
#' Journal, 3*(1), 383–394. https://doi.org/10.18848/1833-1874/CGP/v03i01/37588
#'
#' Zimmerman, D. W. (2004). A note on preliminary tests of equality of
#' variances. *British Journal of Mathematical and Statistical Psychology,
#' 57*(1), 173–181. https://doi.org/10.1348/000711004849222
#' @examples
#' with(solomon_example, check_solomon_assumptions(y_post, treat, pretested, y_pre))
#'
#' # A six-group design: two treatments and a control (Mai et al., 2020).
#' check_solomon_assumptions(post_behavior, condition, pretested, pre_behavior,
#'                           control = "Control", data = mai2020)
#' @export
check_solomon_assumptions <- function(y_post, treat, pretested, y_pre,
                                      control = NULL, data = NULL) {
  .solomon_data_args(
    data, c("y_post", "treat", "pretested", "y_pre"),
    environment(), parent.frame()
  )
  design <- .solomon_conditions(treat, control)
  several <- design$k > 1L
  treat <- if (several) design$condition else design$treat
  pretested <- .solomon_indicator(pretested, "pretested")
  .solomon_check_lengths(
    y_post = y_post,
    treat = treat,
    pretested = pretested,
    y_pre = y_pre
  )

  df <- data.frame(y_post, treat=factor(treat), pretested=factor(pretested))
  df$cell <- interaction(df$pretested, df$treat, drop = TRUE)

  # (1) Heteroscedasticity across all the posttest cells (Brown-Forsythe):
  # four in the four-group design, 2(k + 1) with k treatments
  p_bf4 <- bf_test(df$y_post, df$cell)

  # (2) Heteroscedasticity in posttest-only cells (Groups 3 vs 4, or the
  # k + 1 unpretested groups)
  p_bf_un <- bf_test(df$y_post[df$pretested==0], df$treat[df$pretested==0])

  # (3) Normality (Shapiro) within each posttest cell (NA if n<3)
  shaps <- tapply(df$y_post, df$cell, function(v) if (sum(is.finite(v)) >= 3) stats::shapiro.test(v)$p.value else NA_real_)

  # (4) Homogeneity of regression slopes for ANCOVA among pretested
  # participants with complete scores: y_post ~ treat * y_pre
  # (if interaction significant -> slope heterogeneity). With k treatments
  # the interaction has k degrees of freedom, one slope difference for each
  # treatment; the test needs every one of the k + 1 pretested groups, and
  # is not reported on fewer degrees of freedom.
  p_slope <- NA_real_
  pre_rows <- which(pretested == 1L & is.finite(y_pre) & is.finite(y_post))
  if (several) {
    pre_rows <- pre_rows[!is.na(treat[pre_rows])]
  }
  n_groups <- design$k + 1L
  if (length(pre_rows) > 2L * n_groups && length(unique(treat[pre_rows])) == n_groups) {
    dd <- data.frame(y = y_post[pre_rows], treat = factor(treat[pre_rows]), pre = y_pre[pre_rows])
    fit <- stats::lm(y ~ treat * pre, data = dd)
    a <- stats::anova(fit)
    slope_row <- which(rownames(a) == "treat:pre")[1]
    if (!is.na(slope_row) && a$Df[slope_row] == design$k) {
      p_slope <- a$`Pr(>F)`[slope_row]
    }
  }

  if (several) {
    return(structure(list(
      brown_forsythe_cells_p = p_bf4,
      brown_forsythe_unpre_p = p_bf_un,
      shapiro_p_by_cell = shaps,
      ancova_slope_homogeneity_p = p_slope,
      conditions = .conditions_table(design)
    ), class = "solomon_checks"))
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
  row <- function(label, p) cat(sprintf("  %-52s p = %s\n", label, p_fmt(p)))

  shp_min <- suppressWarnings(min(x$shapiro_p_by_cell, na.rm = TRUE))
  if (!is.finite(shp_min)) shp_min <- NA_real_

  cat("Solomon assumption diagnostics (descriptive)\n")
  if (!is.null(x$conditions)) {
    cat(.ngroup_design_line(x$conditions), "\n", sep = "")
    row("Equal variance, all posttest cells (Brown-Forsythe)", x$brown_forsythe_cells_p)
  } else {
    row("Equal variance, four posttest cells (Brown-Forsythe)", x$brown_forsythe_4cell_p)
  }
  row("Equal variance, unpretested cells (Brown-Forsythe)", x$brown_forsythe_unpre_p)
  row("Normality within cells (Shapiro-Wilk, smallest p)", shp_min)
  row("Homogeneous slopes, pretested groups (Treat x Pre)", x$ancova_slope_homogeneity_p)

  cat(
    "\nThese p-values describe the data; they are not gates for choosing an\n",
    "analysis. Selecting a test because a preliminary assumption test was or\n",
    "was not significant can distort Type I error rates (Zimmerman, 2004).\n",
    "HC3 robust standard errors are a reasonable default for the unified GLM\n",
    "regardless of these results (Long & Ervin, 2000). A small slope p-value\n",
    "suggests the treatment effect in pretested groups depends on the pretest\n",
    "score, which is substantively informative.\n",
    sep = ""
  )
  invisible(x)
}
