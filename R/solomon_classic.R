# Internal helper: test a linear contrast from an lm object
.classic_contrast <- function(fit, L, test) {

  b <- stats::coef(fit)
  V <- stats::vcov(fit)

  L_full <- stats::setNames(
    numeric(length(b)),
    names(b)
  )

  common <- intersect(names(L), names(b))
  L_full[common] <- L[common]

  estimate <- sum(L_full * b)

  variance <- as.numeric(
    t(L_full) %*% V %*% L_full
  )

  se <- sqrt(variance)
  statistic <- estimate / se
  df <- stats::df.residual(fit)

  p <- 2 * stats::pt(
    -abs(statistic),
    df = df
  )

  data.frame(
    test = test,
    estimate = estimate,
    std.error = se,
    statistic = statistic,
    df = df,
    F = statistic^2,
    p.value = p,
    row.names = NULL
  )
}


# Internal helper: directional one-tailed p-value from a t statistic
.classic_one_tailed_p <- function(statistic, df, direction) {

  if (direction == "greater") {
    stats::pt(
      statistic,
      df = df,
      lower.tail = FALSE
    )
  } else {
    stats::pt(
      statistic,
      df = df,
      lower.tail = TRUE
    )
  }
}


# Internal helper: Braver & Braver Stouffer combination
.classic_stouffer_pair <- function(pre_result,
                                   unpre_result,
                                   direction,
                                   source) {

  p_pre <- .classic_one_tailed_p(
    pre_result$statistic,
    pre_result$df,
    direction
  )

  p_un <- .classic_one_tailed_p(
    unpre_result$statistic,
    unpre_result$df,
    direction
  )

  # Keep qnorm() away from exact 0 and 1.
  eps <- .Machine$double.eps

  p_one <- pmin(
    pmax(c(p_pre, p_un), eps),
    1 - eps
  )

  combined <- stouffer_solomon(p_one)

  data.frame(
    test = "I",
    procedure = "Braver & Braver (1988)",
    source = source,
    z = combined$z_meta,
    p.value = combined$p_meta_one_tailed,
    p_pretested_one_tailed = p_one[1],
    p_unpretested_one_tailed = p_one[2],
    row.names = NULL
  )
}


#' Historical Solomon Four-Group Analysis
#'
#' Implements the historical Test A-I framework associated with analysis of
#' the Solomon four-group design. All tests are calculated when possible,
#' while the historical decision pathway is stored separately.
#'
#' The historical sequence includes the four-group posttest factorial model,
#' simple treatment effects, ANCOVA, gain-score analysis, the equivalent
#' two-wave repeated-measures interaction, the posttest-only comparison, and
#' the optional Stouffer meta-analytic combination.
#'
#' Test I, the Braver & Braver (1988) Stouffer meta-analytic combination, is
#' included for historical replication and teaching. Later work (see
#' Sawilowsky et al., 1994) raised concerns about experiment-wise Type I
#' error when the procedure is used conditionally. Its presence in this function should not be interpreted
#' as a general contemporary recommendation.
#'
#' @param y_post Numeric posttest scores.
#' @param treat Treatment indicator coded 0 = control and 1 = treatment.
#' @param pretested Pretest indicator coded 0 = unpretested and 1 = pretested.
#' @param y_pre Numeric pretest scores. Values should be missing for
#'   participants assigned to the unpretested groups.
#' @param alpha Significance level used to reconstruct the historical decision
#'   pathway. Default is 0.05.
#' @param pretested_test Which historical pretested-group analysis should be
#'   followed in the decision pathway: \code{"ancova"}, \code{"gain"}, or
#'   \code{"repeated"}. All three are still calculated and returned.
#' @param combine_with_stouffer Logical. If \code{TRUE}, include Test I, the
#'   Braver & Braver (1988) Stouffer combination, in the decision pathway when earlier treatment tests are
#'   nonsignificant.
#' @param stouffer_direction Direction of the historical one-tailed treatment
#'   hypothesis used for Test I: \code{"greater"} or \code{"less"}.
#'
#' @param conf_level Confidence level for intervals. Default is 0.95.
#'
#' @section Confidence intervals:
#' Tests A-H are t tests (equivalently, F tests with one numerator degree of
#' freedom) with residual degrees of freedom, and each result carries the
#' matching confidence interval (`conf.low`, `conf.high`). Test I combines
#' p-values and has no interval. The Groups 3-4 standardized mean difference
#' (`g_post`) reports Hedges' g with a noncentral t interval for the
#' population standardized mean difference (Cumming & Finch, 2001; Kelley,
#' 2007).
#'
#' @return An object of class \code{solomon_classic}. The \code{tests}
#'   component contains Tests A-I, while \code{path} records the historical
#'   decision sequence for the observed data.
#'
#' @references
#' Solomon, R. L. (1949). An extension of control group design.
#' *Psychological Bulletin, 46*(2), 137-150.
#'
#' Campbell, D. T., & Stanley, J. C. (1963). *Experimental and
#' quasi-experimental designs for research*. Rand McNally.
#'
#' Huck, S. W., & Sandler, H. M. (1973). A note on the Solomon 4-group
#' design: Appropriate statistical analyses. *The Journal of Experimental
#' Education, 42*(2), 54-55.
#'
#' Braver, M. W., & Braver, S. L. (1988). Statistical treatment of the
#' Solomon four-group design: A meta-analytic approach. *Psychological
#' Bulletin, 104*(1), 150-154.
#'
#' Sawilowsky, S. S., Kelley, D. L., Blair, R. C., & Markman, B. S. (1994).
#' Meta-analysis and the Solomon four-group design. *The Journal of
#' Experimental Education, 62*(4), 361-376.
#'
#' Van Breukelen, G. J. P. (2006). ANCOVA versus change from baseline had
#' more power in randomized studies and more bias in nonrandomized studies.
#' *Journal of Clinical Epidemiology, 59*(9), 920-925.
#'
#' Cumming, G., & Finch, S. (2001). A primer on the understanding, use, and
#' calculation of confidence intervals that are based on central and
#' noncentral distributions. *Educational and Psychological Measurement,
#' 61*(4), 532-574.
#'
#' Kelley, K. (2007). Confidence intervals for standardized effect sizes:
#' Theory, application, and implementation. *Journal of Statistical
#' Software, 20*(8), 1-24.
#'
#' @export
fit_solomon_classic <- function(
    y_post,
    treat,
    pretested,
    y_pre,
    alpha = 0.05,
    pretested_test = c("ancova", "gain", "repeated"),
    combine_with_stouffer = TRUE,
    stouffer_direction = c("greater", "less"),
    conf_level = 0.95
) {

  pretested_test <- match.arg(pretested_test)
  stouffer_direction <- match.arg(stouffer_direction)
  .check_conf_level(conf_level)

  if (length(y_post) != length(treat) ||
      length(y_post) != length(pretested) ||
      length(y_post) != length(y_pre)) {
    stop("y_post, treat, pretested, and y_pre must have equal lengths.")
  }

  treat <- as.integer(treat)
  pretested <- as.integer(pretested)

  if (any(!stats::na.omit(treat) %in% c(0L, 1L))) {
    stop("treat must be coded 0/1.")
  }

  if (any(!stats::na.omit(pretested) %in% c(0L, 1L))) {
    stop("pretested must be coded 0/1.")
  }

  df <- data.frame(
    y_post = y_post,
    treat = treat,
    pretested = pretested,
    y_pre = y_pre
  )

  # ----------------------------------------------------------
  # Tests A-D: four-group posttest model
  #
  # Numeric 0/1 coding lets the historical effects be expressed
  # directly as linear contrasts.
  # ----------------------------------------------------------

  fit_ad <- stats::lm(
    y_post ~ treat * pretested,
    data = df,
    na.action = stats::na.exclude
  )

  # Test A: Treatment x Pretest interaction
  A <- .classic_contrast(
    fit_ad,
    c("treat:pretested" = 1),
    "A"
  )

  # Test B: treatment simple effect among pretested groups
  B <- .classic_contrast(
    fit_ad,
    c(
      "treat" = 1,
      "treat:pretested" = 1
    ),
    "B"
  )

  # Test C: treatment simple effect among unpretested groups
  C <- .classic_contrast(
    fit_ad,
    c("treat" = 1),
    "C"
  )

  # Test D: equal-weighted treatment main effect across
  # pretesting conditions
  D <- .classic_contrast(
    fit_ad,
    c(
      "treat" = 1,
      "treat:pretested" = 0.5
    ),
    "D"
  )

  # Auxiliary pretesting main effect, useful for teaching but
  # not one of the historical A-I treatment tests.
  pretest_main <- .classic_contrast(
    fit_ad,
    c(
      "pretested" = 1,
      "treat:pretested" = 0.5
    ),
    "Pretest main effect"
  )

  # ----------------------------------------------------------
  # Pretested groups: Tests E, F, G
  # ----------------------------------------------------------

  idx_pre <- df$pretested == 1L &
    stats::complete.cases(
      df$y_post,
      df$y_pre,
      df$treat
    )

  pre_df <- df[idx_pre, , drop = FALSE]

  if (length(unique(pre_df$treat)) != 2L) {
    stop(
      "Both treatment conditions must be represented among ",
      "complete pretested cases."
    )
  }

  # Center the pretest so the ANCOVA treatment coefficient is
  # evaluated at the observed mean pretest.
  pre_df$y_pre_c <- pre_df$y_pre - mean(pre_df$y_pre)

  # Test E: ANCOVA
  fit_e <- stats::lm(
    y_post ~ treat + y_pre_c,
    data = pre_df
  )

  E <- .classic_contrast(
    fit_e,
    c("treat" = 1),
    "E"
  )

  # Test F: gain scores
  pre_df$gain <- pre_df$y_post - pre_df$y_pre

  fit_f <- stats::lm(
    gain ~ treat,
    data = pre_df
  )

  F <- .classic_contrast(
    fit_f,
    c("treat" = 1),
    "F"
  )

  # Test G: Treatment x Time interaction in a two-wave
  # repeated-measures model.
  #
  # With exactly two occasions, the Treatment x Time interaction
  # is algebraically equivalent to the between-group comparison of
  # gain scores. We therefore return the same inferential statistic
  # while preserving its historical interpretation.
  G <- F
  G$test <- "G"

  # ----------------------------------------------------------
  # Test H: treatment comparison among unpretested groups
  # ----------------------------------------------------------

  idx_un <- df$pretested == 0L &
    stats::complete.cases(
      df$y_post,
      df$treat
    )

  un_df <- df[idx_un, , drop = FALSE]

  if (length(unique(un_df$treat)) != 2L) {
    stop(
      "Both treatment conditions must be represented among ",
      "complete unpretested cases."
    )
  }

  # The historical independent-samples t test corresponds to the
  # treatment coefficient from this equal-variance linear model.
  fit_h <- stats::lm(
    y_post ~ treat,
    data = un_df
  )

  H <- .classic_contrast(
    fit_h,
    c("treat" = 1),
    "H"
  )

  # ----------------------------------------------------------
  # Test I: Stouffer combinations
  #
  # Calculate all three historical variants. The selected
  # pretested analysis determines the version used in `path`.
  # ----------------------------------------------------------

  I_all <- list(
    E = .classic_stouffer_pair(
      E,
      H,
      stouffer_direction,
      "E + H (ANCOVA + posttest-only)"
    ),
    F = .classic_stouffer_pair(
      F,
      H,
      stouffer_direction,
      "F + H (gain score + posttest-only)"
    ),
    G = .classic_stouffer_pair(
      G,
      H,
      stouffer_direction,
      "G + H (repeated measures + posttest-only)"
    )
  )

  selected_letter <- switch(
    pretested_test,
    ancova = "E",
    gain = "F",
    repeated = "G"
  )

  I_selected <- I_all[[selected_letter]]

  # ----------------------------------------------------------
  # Reconstruct the historical decision pathway
  # ----------------------------------------------------------

  path <- "A"

  if (A$p.value < alpha) {

    path <- c(path, "B", "C")

    b_sig <- B$p.value < alpha
    c_sig <- C$p.value < alpha

    if (b_sig && c_sig) {
      conclusion <- paste(
        "Historical pathway: evidence of pretest sensitization;",
        "the treatment effect is detected in both pretested and",
        "unpretested groups but differs by pretesting condition."
      )
    } else if (b_sig && !c_sig) {
      conclusion <- paste(
        "Historical pathway: evidence of pretest sensitization;",
        "the treatment effect is detected only among pretested groups."
      )
    } else if (!b_sig && c_sig) {
      conclusion <- paste(
        "Historical pathway: evidence of pretest sensitization;",
        "the treatment effect is detected only among unpretested groups."
      )
    } else {
      conclusion <- paste(
        "Historical pathway: the Treatment x Pretest interaction is",
        "significant, but neither simple treatment effect reaches alpha."
      )
    }

  } else {

    path <- c(path, "D")

    if (D$p.value < alpha) {

      conclusion <- paste(
        "Historical pathway: no evidence of pretest sensitization",
        "and the treatment main effect is significant."
      )

    } else {

      selected_results <- list(
        E = E,
        F = F,
        G = G
      )

      selected_result <- selected_results[[selected_letter]]

      path <- c(path, selected_letter)

      if (selected_result$p.value < alpha) {

        conclusion <- paste(
          "Historical pathway: no evidence of pretest sensitization;",
          paste0("Test ", selected_letter),
          "detects a treatment effect in the pretested groups."
        )

      } else {

        path <- c(path, "H")

        if (H$p.value < alpha) {

          conclusion <- paste(
            "Historical pathway: no evidence of pretest sensitization;",
            "Test H detects a treatment effect in the unpretested groups."
          )

        } else if (isTRUE(combine_with_stouffer)) {

          path <- c(path, "I")

          if (I_selected$p.value < alpha) {
            conclusion <- paste(
              "Historical pathway: Test I produces a significant",
              "Stouffer combination. This result is retained for",
              "historical replication and should be interpreted in",
              "light of later Type I error critiques."
            )
          } else {
            conclusion <- paste(
              "Historical pathway: no treatment test in the selected",
              "A-I sequence reaches the specified alpha level."
            )
          }

        } else {

          conclusion <- paste(
            "Historical pathway: Tests A, D,",
            selected_letter,
            "and H are nonsignificant; Test I was not requested."
          )
        }
      }
    }
  }

  # ----------------------------------------------------------
  # Confidence intervals for Tests A-H (t with residual df)
  # ----------------------------------------------------------

  with_ci <- function(res) {
    ci <- .wald_ci(res$estimate, res$std.error, res$df, conf_level)
    res$conf.low <- unname(ci[, "conf.low"])
    res$conf.high <- unname(ci[, "conf.high"])
    res
  }

  A <- with_ci(A)
  B <- with_ci(B)
  C <- with_ci(C)
  D <- with_ci(D)
  E <- with_ci(E)
  F <- with_ci(F)
  G <- with_ci(G)
  H <- with_ci(H)
  pretest_main <- with_ci(pretest_main)

  # ----------------------------------------------------------
  # Hedges g for Groups 3 and 4
  # ----------------------------------------------------------

  x1 <- un_df$y_post[un_df$treat == 1L]
  x0 <- un_df$y_post[un_df$treat == 0L]

  g_post <- hedges_g_ci(
    mean(x1),
    mean(x0),
    stats::sd(x1),
    stats::sd(x0),
    length(x1),
    length(x0),
    conf = conf_level
  )

  # Legacy-friendly objects retained for existing user code.
  aov_legacy <- broom::tidy(
    stats::aov(
      y_post ~ factor(treat) * factor(pretested),
      data = df
    )
  )

  ancova_legacy <- broom::tidy(fit_e)

  t_legacy <- broom::tidy(
    stats::t.test(
      y_post ~ factor(treat),
      data = un_df,
      var.equal = TRUE
    )
  )

  stouffer_legacy <- NULL

  if (isTRUE(combine_with_stouffer)) {
    stouffer_legacy <- list(
      z_meta = I_selected$z,
      p_meta_one_tailed = I_selected$p.value
    )
  }

  tests <- list(
    A = list(
      label = "Pretest x Treatment interaction",
      result = A,
      model = fit_ad
    ),
    B = list(
      label = "Treatment effect among pretested groups",
      result = B,
      model = fit_ad
    ),
    C = list(
      label = "Treatment effect among unpretested groups",
      result = C,
      model = fit_ad
    ),
    D = list(
      label = "Treatment main effect",
      result = D,
      model = fit_ad
    ),
    E = list(
      label = "ANCOVA treatment effect",
      result = E,
      model = fit_e
    ),
    F = list(
      label = "Gain-score treatment effect",
      result = F,
      model = fit_f
    ),
    G = list(
      label = "Repeated-measures Treatment x Time interaction",
      result = G,
      model = fit_f,
      note = "Equivalent to Test F for exactly two measurement occasions."
    ),
    H = list(
      label = "Posttest-only treatment effect",
      result = H,
      model = fit_h
    ),
    I = list(
      label = "Braver & Braver (1988) Stouffer combination",
      result = I_selected,
      all = I_all
    )
  )

  structure(
    list(
      tests = tests,
      path = path,
      path_string = paste(path, collapse = " -> "),
      conclusion = conclusion,
      pretest_main = pretest_main,
      settings = list(
        alpha = alpha,
        pretested_test = pretested_test,
        selected_test = selected_letter,
        combine_with_stouffer = combine_with_stouffer,
        stouffer_direction = stouffer_direction,
        conf_level = conf_level
      ),

      # Legacy fields
      aov = aov_legacy,
      ancova = ancova_legacy,
      t_unpretested = t_legacy,
      stouffer = stouffer_legacy,
      g_post = g_post
    ),
    class = "solomon_classic"
  )
}
