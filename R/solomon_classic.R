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


# Internal helper: Walton Braver & Braver Stouffer combination
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
    procedure = "Walton Braver & Braver (1988)",
    source = source,
    z = combined$z_meta,
    p.value = combined$p_meta_one_tailed,
    p_pretested_one_tailed = p_one[1],
    p_unpretested_one_tailed = p_one[2],
    row.names = NULL
  )
}


# Historical decision pathway from the p-values of Tests A-I (`S` is the
# selected pretested-groups test, E, F, or G).
#
# 1988: Walton Braver and Braver (1988, pp. 151-153). A significant Test A
# leads to Tests B and C and ends the sequence; otherwise testing stops at
# the first significant test among D, the selected test, and H, and Test I
# is reached only when all of them are nonsignificant.
# 1990: Braver and Walton Braver (1990, p. 322). As 1988 through Test D;
# once Tests A and D are nonsignificant, every test through Test I is run
# and Test I is regarded as the most definitive.
# 1995: Walton Braver and Braver (1995, as cited in Sawilowsky, 1996, p. 2)
# removed Test D. The remaining tests follow the 1988 rule, each reached only
# if the ones before it are nonsignificant, as in Sawilowsky's (1996)
# simulation of the revised sequence.
#
# `alpha` is one level for every test, or a named vector of levels for A, B,
# C, D, S, H, and I (see .classic_alpha_levels()).
.classic_path <- function(p, alpha, flow, selected, combine_with_stouffer) {
  if (is.null(names(alpha))) {
    alpha <- stats::setNames(rep(alpha, 7L), c("A", "B", "C", "D", "S", "H", "I"))
  }
  sig <- stats::setNames(p < alpha[names(p)], names(p))

  if (sig[["A"]]) {
    conclusion <- if (sig[["B"]] && sig[["C"]]) {
      paste(
        "Historical pathway: evidence of pretest sensitization;",
        "the treatment effect is detected in both pretested and",
        "unpretested groups but differs by pretesting condition."
      )
    } else if (sig[["B"]]) {
      paste(
        "Historical pathway: evidence of pretest sensitization;",
        "the treatment effect is detected only among pretested groups."
      )
    } else if (sig[["C"]]) {
      paste(
        "Historical pathway: evidence of pretest sensitization;",
        "the treatment effect is detected only among unpretested groups."
      )
    } else {
      paste(
        "Historical pathway: the Treatment x Pretest interaction is",
        "significant, but neither simple treatment effect reaches alpha."
      )
    }
    return(list(path = c("A", "B", "C"), conclusion = conclusion))
  }

  without_d <- flow == "1995"
  first <- if (without_d) "A" else c("A", "D")
  lead <- if (without_d) {
    "Historical pathway (1995 revision, without Test D):"
  } else {
    "Historical pathway:"
  }

  if (!without_d && sig[["D"]]) {
    return(list(
      path = c("A", "D"),
      conclusion = paste(
        "Historical pathway: no evidence of pretest sensitization",
        "and the treatment main effect is significant."
      )
    ))
  }

  if (flow == "1990") {
    return(list(
      path = c("A", "D", selected, "H", "I"),
      conclusion = paste0(
        "Historical pathway (1990 amendment): Tests A and D are ",
        "nonsignificant, so every test through Test I is run and Test I ",
        "is regarded as the most definitive. ",
        if (sig[["I"]]) {
          "Test I produces a significant Stouffer combination"
        } else {
          "Test I does not reach the specified alpha level"
        },
        sprintf(
          " (Test %s p = %.3f, Test H p = %.3f, Test I p = %.3f). ",
          selected, p[["S"]], p[["H"]], p[["I"]]
        ),
        "Interpret this in light of later Type I error critiques."
      )
    ))
  }

  if (sig[["S"]]) {
    return(list(
      path = c(first, selected),
      conclusion = paste(
        lead,
        "no evidence of pretest sensitization;",
        paste0("Test ", selected),
        "detects a treatment effect in the pretested groups."
      )
    ))
  }

  if (sig[["H"]]) {
    return(list(
      path = c(first, selected, "H"),
      conclusion = paste(
        lead,
        "no evidence of pretest sensitization;",
        "Test H detects a treatment effect in the unpretested groups."
      )
    ))
  }

  if (!isTRUE(combine_with_stouffer)) {
    return(list(
      path = c(first, selected, "H"),
      conclusion = paste(
        lead,
        paste0("Tests ", paste(first, collapse = ", "), ","),
        selected,
        "and H are nonsignificant; Test I was not requested."
      )
    ))
  }

  list(
    path = c(first, selected, "H", "I"),
    conclusion = if (sig[["I"]]) {
      paste(
        lead,
        "Test I produces a significant",
        "Stouffer combination. This result is retained for",
        "historical replication and should be interpreted in",
        "light of later Type I error critiques."
      )
    } else {
      paste(
        lead,
        "no treatment test in the selected",
        "A-I sequence reaches the specified alpha level."
      )
    }
  )
}


# Test-wise significance levels for the historical sequence. With an
# allocation, the levels are those Sawilowsky (1996, Table 4) obtained by
# Monte Carlo for Tests A, E, H, and I (the 1995 sequence, nominal alpha =
# .05) under Bradley's (1968, as cited in Sawilowsky, 1996) conservative and
# liberal robustness criteria. Tests B and C, reached only after a
# significant Test A, are not part of that allocation and keep `alpha`.
.classic_alpha_levels <- function(alpha, allocation) {
  levels <- stats::setNames(rep(alpha, 7L), c("A", "B", "C", "D", "S", "H", "I"))
  published <- switch(
    allocation,
    none = NULL,
    method1_conservative = c(A = 0.02, S = 0.02, H = 0.02, I = 0.02),
    method1_liberal = c(A = 0.0275, S = 0.0275, H = 0.0275, I = 0.0275),
    method2_conservative = c(A = 0.05, S = 0.005, H = 0.005, I = 0.005),
    method2_liberal = c(A = 0.05, S = 0.02, H = 0.02, I = 0.02)
  )
  levels[names(published)] <- published
  levels
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
#' Test I, the Walton Braver & Braver (1988) Stouffer meta-analytic combination, is
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
#'   Walton Braver & Braver (1988) Stouffer combination, in the decision pathway when earlier treatment tests are
#'   nonsignificant.
#' @param stouffer_direction Direction of the historical one-tailed treatment
#'   hypothesis used for Test I: \code{"greater"} or \code{"less"}.
#'
#' @param conf_level Confidence level for intervals. Default is 0.95.
#' @param flow Which published version of the decision sequence to follow.
#'   `"1988"` (the default) is the original sequence of Walton Braver and
#'   Braver (1988, pp. 151--153): testing stops at the first significant test
#'   among D, the selected pretested-groups test (E, F, or G), and H, and
#'   Test I is reached only when all are nonsignificant. `"1990"` is the
#'   authors' amendment (Braver & Walton Braver, 1990, p. 322): once Tests A
#'   and D are nonsignificant, every test through Test I is run and Test I is
#'   regarded as the most definitive. `"1995"` is the authors' later
#'   revision, which removed Test D (Walton Braver & Braver, 1995, as cited in
#'   Sawilowsky, 1996, p. 2); the revision itself is unpublished. The
#'   remaining tests follow the 1988 rule, each reached only if the
#'   ones before it are nonsignificant, as in Sawilowsky's (1996) simulation
#'   of the revised sequence. Only the `path` and `conclusion` differ; every
#'   test is always calculated.
#' @param alpha_allocation Test-wise significance levels for the historical
#'   sequence. `"none"` (the default) uses `alpha` for every test. The other
#'   options are the two Bonferroni-type allocations of Sawilowsky (1996,
#'   Table 4), which he calibrated by Monte Carlo so that the experiment-wise
#'   Type I error of the 1995 sequence stays within Bradley's (1968, as cited
#'   in Sawilowsky, 1996) conservative (`"_conservative"`) or liberal
#'   (`"_liberal"`) robustness limit for a nominal alpha of .05:
#'   - `"method1_conservative"` and `"method1_liberal"`: Tests A, E, H, and I
#'     each at .020 or .0275;
#'   - `"method2_conservative"` and `"method2_liberal"`: Test A at .05, as a
#'     preliminary test, and Tests E, H, and I each at .005 or .020.
#'
#'   An allocation requires `flow = "1995"`, `alpha = 0.05`, `pretested_test =
#'   "ancova"`, and `combine_with_stouffer = TRUE`, the conditions for which
#'   the levels were obtained. Tests B and C, reached only after a
#'   significant Test A, are outside the allocation and keep `alpha`; that
#'   choice is solomonR's.
#'
#' @section History and maturation:
#' The `history` component compares the unpretested control posttest (O6)
#' with the pretests of the pretested groups (O1 and O3) by independent-samples
#' t tests, as Mai et al. (2020, p. 8) report. Neither group had received the
#' treatment when these scores were measured, so a difference estimates the
#' combined effect of history and maturation between the two occasions
#' (Campbell & Stanley, 1963/1966, p. 25), provided assignment was random and
#' the measure is comparable at both occasions. Solomon (1949, pp. 146–148)
#' introduced the fourth group for this purpose, attributing its change from
#' the pretest to outside events and the passage of time. It is reported as a
#' historical check, not a test of the treatment.
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
#'   decision sequence for the observed data under the chosen \code{flow}.
#'   The \code{history} component holds the history/maturation comparisons.
#'
#' @references
#' Braver, S. L., & Walton Braver, M. C. (1990). Meta-analysis for Solomon
#' four-group designs reconsidered: A reply to Sawilowsky and Markman.
#' *Perceptual and Motor Skills, 71*(1), 321–322.
#' https://doi.org/10.2466/pms.1990.71.1.321
#'
#' Campbell, D. T., & Stanley, J. C. (1966). *Experimental and
#' quasi-experimental designs for research*. Rand McNally. (Original work
#' published 1963)
#'
#' Cumming, G., & Finch, S. (2001). A primer on the understanding, use, and
#' calculation of confidence intervals that are based on central and noncentral
#' distributions. *Educational and Psychological Measurement, 61*(4), 532–574.
#' https://doi.org/10.1177/00131640121971374
#'
#' Huck, S. W., & Sandler, H. M. (1973). A note on the Solomon 4-group design:
#' Appropriate statistical analyses. *The Journal of Experimental Education,
#' 42*(2), 54–55. https://doi.org/10.1080/00220973.1973.11011460
#'
#' Kelley, K. (2007). Confidence intervals for standardized effect sizes:
#' Theory, application, and implementation. *Journal of Statistical Software,
#' 20*(8), 1–24. https://doi.org/10.18637/jss.v020.i08
#'
#' Mai, N. N., Takahashi, Y., & Oo, M. M. (2020). Testing the effectiveness of
#' transfer interventions using Solomon four-group designs. *Education
#' Sciences, 10*(4), Article 92. https://doi.org/10.3390/educsci10040092
#'
#' Sawilowsky, S. S. (1996, June 23). *Controlling experiment-wise Type I error
#' of meta-analysis in the Solomon four-group design* \[Paper presentation\].
#' First International Conference on Multiple Comparisons, Tel Aviv, Israel.
#' <http://digitalcommons.wayne.edu/coe_tbf/29>
#'
#' Sawilowsky, S. S., Kelley, D. L., Blair, R. C., & Markman, B. S. (1994).
#' Meta-analysis and the Solomon four-group design. *The Journal of Experimental
#' Education, 62*(4), 361–376. https://doi.org/10.1080/00220973.1994.9944140
#'
#' Solomon, R. L. (1949). An extension of control group design. *Psychological
#' Bulletin, 46*(2), 137–150. https://doi.org/10.1037/h0062958
#'
#' Van Breukelen, G. J. P. (2006). ANCOVA versus change from baseline had more
#' power in randomized studies and more bias in nonrandomized studies. *Journal
#' of Clinical Epidemiology, 59*(9), 920–925.
#' https://doi.org/10.1016/j.jclinepi.2006.02.007
#'
#' Walton Braver, M. C., & Braver, S. L. (1988). Statistical treatment of the
#' Solomon four-group design: A meta-analytic approach. *Psychological Bulletin,
#' 104*(1), 150–154. https://doi.org/10.1037/0033-2909.104.1.150
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
    conf_level = 0.95,
    flow = c("1988", "1990", "1995"),
    alpha_allocation = c("none", "method1_conservative", "method1_liberal",
                         "method2_conservative", "method2_liberal")
) {

  pretested_test <- match.arg(pretested_test)
  stouffer_direction <- match.arg(stouffer_direction)
  flow <- match.arg(flow)
  alpha_allocation <- match.arg(alpha_allocation)
  .check_conf_level(conf_level)

  if (alpha_allocation != "none") {
    if (flow != "1995" || !isTRUE(all.equal(alpha, 0.05)) ||
        pretested_test != "ancova" || !isTRUE(combine_with_stouffer)) {
      stop(
        "Sawilowsky's (1996) alpha allocations were obtained for the 1995 ",
        "sequence of Tests A, E, H, and I at a nominal alpha of .05; use ",
        "flow = \"1995\", alpha = 0.05, pretested_test = \"ancova\", and ",
        "combine_with_stouffer = TRUE.",
        call. = FALSE
      )
    }
  }
  alpha_levels <- .classic_alpha_levels(alpha, alpha_allocation)

  if (flow == "1990" && !isTRUE(combine_with_stouffer)) {
    stop(
      "flow = \"1990\" always completes the sequence with Test I ",
      "(Braver & Walton Braver, 1990); use combine_with_stouffer = TRUE.",
      call. = FALSE
    )
  }

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

  pathway <- .classic_path(
    p = c(
      A = A$p.value, B = B$p.value, C = C$p.value, D = D$p.value,
      S = list(E = E, F = F, G = G)[[selected_letter]]$p.value,
      H = H$p.value, I = I_selected$p.value
    ),
    alpha = alpha_levels,
    flow = flow,
    selected = selected_letter,
    combine_with_stouffer = combine_with_stouffer
  )
  path <- pathway$path
  conclusion <- pathway$conclusion

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

  # ----------------------------------------------------------
  # History/maturation check (historical)
  #
  # The unpretested control posttest (O6) is compared with the pretests
  # of the pretested groups (O1, O3) by independent-samples t tests, as
  # Mai et al. (2020, p. 8) report. Neither group received the treatment
  # before these measurements, so a difference estimates the combined
  # effect of history and maturation between the two occasions.
  # ----------------------------------------------------------

  o6 <- stats::na.omit(df$y_post[df$treat == 0L & df$pretested == 0L])
  history_row <- function(pre, label) {
    pre <- stats::na.omit(pre)
    if (length(pre) < 2L || length(o6) < 2L) return(NULL)
    tt <- stats::t.test(o6, pre, var.equal = TRUE, conf.level = conf_level)
    data.frame(
      comparison = label,
      estimate = mean(o6) - mean(pre),
      statistic = unname(tt$statistic),
      df = unname(tt$parameter),
      p.value = tt$p.value,
      conf.low = tt$conf.int[1],
      conf.high = tt$conf.int[2],
      stringsAsFactors = FALSE
    )
  }
  history <- rbind(
    history_row(df$y_pre[df$treat == 1L & df$pretested == 1L],
                "O6 - O1: control posttest vs. treated-group pretest"),
    history_row(df$y_pre[df$treat == 0L & df$pretested == 1L],
                "O6 - O3: control posttest vs. control-group pretest")
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
      label = "Walton Braver & Braver (1988) Stouffer combination",
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
      history = history,
      settings = list(
        alpha = alpha,
        pretested_test = pretested_test,
        selected_test = selected_letter,
        combine_with_stouffer = combine_with_stouffer,
        stouffer_direction = stouffer_direction,
        conf_level = conf_level,
        flow = flow,
        alpha_allocation = alpha_allocation,
        alpha_levels = alpha_levels
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
