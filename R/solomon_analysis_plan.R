# An editable analysis plan for a new Solomon study (issue #81).

.plan_bullets <- function(x) paste0("- ", x)

# Stop when a plan describes a design with several treatments (issue #45):
# group sizes beyond the four cells, or one planning value per treatment.
# plan_solomon() plans four-group designs only, so its results always pass.
.stop_ngroup_plan <- function(plan) {
  if (is.null(plan)) return(invisible(NULL))
  cells <- grep("^n[0-9]+$", names(plan), value = TRUE)
  s <- attr(plan, "settings")
  several <- is.list(s) && (length(s$delta) > 1L || length(s$sens) > 1L)
  if (length(setdiff(cells, c("n1", "n2", "n3", "n4"))) || several) {
    stop(structure(
      class = c("solomonR_ngroup_unsupported", "error", "condition"),
      list(
        message = paste0(
          "analysis_plan_solomon() writes plans for the four-group design; ",
          "plans with several treatments are not yet supported."
        ),
        call = NULL
      )
    ))
  }
  invisible(NULL)
}

#' Write an analysis plan for a Solomon four-group study
#'
#' `r lifecycle::badge("experimental")`
#' Writes an editable Markdown analysis plan, for preregistration or a
#' protocol, from a few choices and, optionally, a [plan_solomon()] result.
#' The plan follows the three sections of van 't Veer and Giner-Sorolla's
#' (2016, Appendix A) preregistration template (hypotheses, methods, and
#' analysis plan) and names the SPIRIT 2013 item each part answers (Chan et
#' al., 2013). The Solomon-specific content comes from the package: the four
#' contrasts, the pretest adjustment in the pretested groups, how to claim
#' that sensitization is negligible, and how missing posttests are handled.
#' Text in square brackets is for the researcher to complete.
#'
#' @details
#' **What the plan covers.**
#' - **Hypotheses** (SPIRIT item 12). Numbered predictions for the
#'   treatment effect and for sensitization. For interactions, the template
#'   asks for "the expected shape" (van 't Veer & Giner-Sorolla, 2016, p.
#'   10), so the sensitization hypothesis states which treatment effect is
#'   expected to be larger. A claim that sensitization is negligible needs
#'   an equivalence test against a smallest effect of interest fixed in
#'   advance (Lakens, 2017).
#' - **Methods** (items 14 and 16a). Randomization to the four groups; the
#'   pretest-posttest interval, since pretest effects can depend on it
#'   (Entwisle, 1961); identical measurement in all groups (French et al.,
#'   2021b, Recommendation 11); and the planned sample, from `plan` when
#'   given, with its planning values.
#' - **Analysis plan** (items 20a-c). The model, the confirmatory contrasts,
#'   the test of sensitization, the handling of missing posttests, with a
#'   tipping-point sensitivity analysis (White et al., 2011), and the
#'   checks run before the model.
#' - **Deviations** (item 25). A section for recording departures from the
#'   plan, which are to be reported with the results (Nosek et al., 2018, p.
#'   2602).
#'
#' **What it does not decide.** The plan states the package's defaults; it
#' does not choose the hypotheses, the sample size, or the smallest effect
#' of interest, and every sentence can be edited. Prediction and postdiction
#' must stay distinguishable (Nosek et al., 2018, p. 2602), so the plan
#' should be registered before the outcomes are seen.
#'
#' **The design.** The plan is for the four-group design: one treatment and a
#' control, each with and without a pretest. Plans for designs with several
#' treatments are not yet supported; see [fit_solomon_glm()] for their
#' analysis.
#'
#' **The round trip.** The returned object records the confirmatory
#' contrasts. Passing it to [report_solomon()] as `design = list(plan = )`
#' makes the report state whether the sensitization analysis was
#' pre-specified, from the plan rather than by guessing.
#'
#' @param plan Optional result of [plan_solomon()]. Its sample sizes and
#'   planning values fill the planned-sample section.
#' @param outcome Description of the primary outcome, measured at posttest.
#' @param treatment Description of the treatment and of the control
#'   condition.
#' @param direction Expected direction of the treatment effect: `"increase"`
#'   (default), `"decrease"`, or `"two-sided"` for no predicted direction.
#' @param sensitization How the Pretest x Treatment contrast is treated:
#'   `"test"` (default), a confirmatory two-sided test of sensitization;
#'   `"equivalence"`, a confirmatory claim that sensitization is negligible,
#'   tested with two one-sided tests against `equivalence_bound`;
#'   `"larger_pretested"` or `"smaller_pretested"`, a confirmatory
#'   prediction of the direction of sensitization; or `"exploratory"`,
#'   estimated and reported but not confirmatory.
#' @param equivalence_bound The smallest sensitization of interest, in
#'   posttest units, for `sensitization = "equivalence"`.
#' @param alpha Significance level of the confirmatory tests. Default 0.05.
#' @param occasions Number of posttest occasions, or their labels. With more
#'   than one, the primary analysis is [fit_solomon_mmrm()], a mixed model for
#'   repeated measures (Mallinckrodt et al., 2008), and the confirmatory
#'   contrasts are those at `primary_occasion`.
#' @param primary_occasion The occasion whose contrasts are confirmatory, when
#'   there are several. Default: the last.
#' @param tipping_groups Groups whose imputed posttests are shifted in the
#'   tipping-point sensitivity analysis, as in [tipping_point_solomon()].
#'   Default `"treatment"`. When sensitization is confirmatory, a second
#'   analysis in the pretested treatment group alone is added.
#' @param title Title of the plan.
#' @param file Optional path. When given, the plan is also written there as
#'   UTF-8 Markdown.
#'
#' @return An object of class `solomon_analysis_plan`: `text` (the Markdown
#'   lines), `settings` (the choices, including `confirmatory`), and `date`.
#'   Printing it shows the Markdown.
#'
#' @references
#' Chan, A.-W., Tetzlaff, J. M., Altman, D. G., Laupacis, A., Gøtzsche, P.
#' C., Krleža-Jerić, K., Hróbjartsson, A., Mann, H., Dickersin, K., Berlin,
#' J. A., Doré, C. J., Parulekar, W. R., Summerskill, W. S. M., Groves, T.,
#' Schulz, K. F., Sox, H. C., Rockhold, F. W., Rennie, D., & Moher, D.
#' (2013). SPIRIT 2013 statement: Defining standard protocol items for
#' clinical trials. *Annals of Internal Medicine, 158*(3), 200–207.
#' https://doi.org/10.7326/0003-4819-158-3-201302050-00583
#'
#' Entwisle, D. R. (1961). Interactive effects of pretesting. *Educational and
#' Psychological Measurement, 21*(3), 607–620.
#' https://doi.org/10.1177/001316446102100307
#'
#' Fitzmaurice, G. M., Laird, N. M., & Ware, J. H. (2011). *Applied
#' longitudinal analysis* (2nd ed.). Wiley. https://doi.org/10.1002/9781119513469
#'
#' French, D. P., Miles, L. M., Elbourne, D., Farmer, A., Gulliford, M.,
#' Locock, L., Sutton, S., McCambridge, J., & MERIT Collaborative Group.
#' (2021b). Reducing bias in trials from reactions to measurement: The MERIT
#' study including developmental work and expert workshop. *Health Technology
#' Assessment, 25*(55), 1–72. https://doi.org/10.3310/hta25550
#'
#' Lakens, D. (2017). Equivalence tests: A practical primer for t tests,
#' correlations, and meta-analyses. *Social Psychological and Personality
#' Science, 8*(4), 355–362. https://doi.org/10.1177/1948550617697177
#'
#' Mallinckrodt, C. H., Lane, P. W., Schnell, D., Peng, Y., & Mancuso, J. P.
#' (2008). Recommendations for the primary analysis of continuous endpoints
#' in longitudinal clinical trials. *Drug Information Journal, 42*(4),
#' 303–319. https://doi.org/10.1177/009286150804200402
#'
#' Nosek, B. A., Ebersole, C. R., DeHaven, A. C., & Mellor, D. T. (2018). The
#' preregistration revolution. *Proceedings of the National Academy of
#' Sciences, 115*(11), 2600–2606. https://doi.org/10.1073/pnas.1708274114
#'
#' van 't Veer, A. E., & Giner-Sorolla, R. (2016). Pre-registration in social
#' psychology—A discussion and suggested template. *Journal of Experimental
#' Social Psychology, 67*, 2–12. https://doi.org/10.1016/j.jesp.2016.03.004
#'
#' White, I. R., Horton, N. J., Carpenter, J., & Pocock, S. J. (2011).
#' Strategy for intention to treat analysis in randomised trials with missing
#' outcome data. *BMJ, 342*, Article d40. https://doi.org/10.1136/bmj.d40
#'
#' @seealso [plan_solomon()] for the sample size, [report_solomon()] for the
#'   results, and the R Markdown template "Solomon four-group study"
#'   (`rmarkdown::draft("study.Rmd", "solomon-study", package = "solomonR")`),
#'   which runs the planned analysis in order.
#'
#' @examples
#' p <- plan_solomon(delta = 0.4, sens = 0, estimand = "ate")
#' analysis_plan_solomon(
#'   plan = p,
#'   outcome = "reading comprehension (0-100)",
#'   treatment = "a six-week tutoring program, against the usual lessons",
#'   sensitization = "equivalence",
#'   equivalence_bound = 3
#' )
#'
#' @export
analysis_plan_solomon <- function(plan = NULL,
                                  outcome = "[the primary outcome, measured at posttest]",
                                  treatment = "[the treatment, and the control condition]",
                                  direction = c("increase", "decrease", "two-sided"),
                                  sensitization = c("test", "equivalence", "larger_pretested",
                                                    "smaller_pretested", "exploratory"),
                                  equivalence_bound = NULL,
                                  alpha = 0.05,
                                  occasions = 1, primary_occasion = NULL,
                                  tipping_groups = "treatment",
                                  title = "Analysis plan for a Solomon four-group study",
                                  file = NULL) {
  direction <- match.arg(direction)
  sensitization <- match.arg(sensitization)
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 0.5) {
    stop("`alpha` must be a single number between 0 and 0.5.", call. = FALSE)
  }
  if (sensitization == "equivalence") {
    if (!is.numeric(equivalence_bound) || length(equivalence_bound) != 1L ||
        !is.finite(equivalence_bound) || equivalence_bound <= 0) {
      stop("`sensitization = \"equivalence\"` needs a positive `equivalence_bound`: the ",
           "smallest sensitization of interest, in posttest units.", call. = FALSE)
    }
  }
  if (!is.null(plan) && !(is.data.frame(plan) && "estimand" %in% names(plan) &&
                          all(c("n1", "n2", "n3", "n4") %in% names(plan)))) {
    stop("`plan` must be a result of plan_solomon().", call. = FALSE)
  }
  .stop_ngroup_plan(plan)
  if (!(is.character(tipping_groups) && length(tipping_groups) == 1L &&
        tipping_groups %in% names(.tipping_groups)) &&
      !(is.numeric(tipping_groups) && length(tipping_groups) &&
        all(tipping_groups %in% 1:4))) {
    stop("`tipping_groups` must be one of ",
         paste(sprintf("\"%s\"", names(.tipping_groups)), collapse = ", "),
         ", or group numbers from 1 to 4.", call. = FALSE)
  }

  if (is.numeric(occasions) && length(occasions) == 1L) {
    if (!is.finite(occasions) || occasions < 1 || occasions != round(occasions)) {
      stop("`occasions` must be a whole number of at least 1, or the occasion labels.",
           call. = FALSE)
    }
    occasion_labels <- as.character(seq_len(occasions))
  } else {
    occasion_labels <- as.character(occasions)
    if (!length(occasion_labels) || anyDuplicated(occasion_labels) || anyNA(occasion_labels)) {
      stop("Occasion labels in `occasions` must be distinct.", call. = FALSE)
    }
  }
  longitudinal <- length(occasion_labels) > 1L
  if (is.null(primary_occasion)) primary_occasion <- occasion_labels[length(occasion_labels)]
  primary_occasion <- as.character(primary_occasion)
  if (length(primary_occasion) != 1L || !primary_occasion %in% occasion_labels) {
    stop("`primary_occasion` must be one of the occasions.", call. = FALSE)
  }
  at <- if (longitudinal) sprintf(" at occasion %s", primary_occasion) else ""

  sens_confirmatory <- sensitization != "exploratory"
  confirmatory <- c("ATE (avg over pretest)", if (sens_confirmatory) "Pretest x Treatment")
  alpha_txt <- sub("^0", "", format(alpha))
  keys <- c("vantveer2016", "chan2013", "solomon1949", "entwisle1961", "french2021b",
            "lin2013", "newman1990", "mackinnon1985", "long2000", "carpenter2023", "white2011",
            "nosek2018")

  # ---- Hypotheses --------------------------------------------------------------
  h1 <- switch(
    direction,
    increase = sprintf("**H1.** Averaged over the pretest conditions, %s will raise %s%s relative to the control condition.", treatment, outcome, at),
    decrease = sprintf("**H1.** Averaged over the pretest conditions, %s will lower %s%s relative to the control condition.", treatment, outcome, at),
    `two-sided` = sprintf("**H1.** Averaged over the pretest conditions, %s will change %s%s relative to the control condition; no direction is predicted.", treatment, outcome, at)
  )
  h2 <- switch(
    sensitization,
    test = "**H2.** The treatment effect may differ between pretested and unpretested participants (pretest sensitization). No direction is predicted; the Pretest x Treatment contrast is tested two-sided.",
    equivalence = sprintf("**H2.** Pretesting will not change the treatment effect by as much as %s points, the smallest sensitization of interest. The Pretest x Treatment contrast is tested for equivalence within -%s and %s.", format(equivalence_bound), format(equivalence_bound), format(equivalence_bound)),
    larger_pretested = "**H2.** The treatment effect will be larger among pretested than among unpretested participants: the Pretest x Treatment contrast will be positive.",
    smaller_pretested = "**H2.** The treatment effect will be smaller among pretested than among unpretested participants: the Pretest x Treatment contrast will be negative.",
    exploratory = "**Exploratory.** The Pretest x Treatment contrast (pretest sensitization) will be estimated and reported with its confidence interval, but no hypothesis about it is tested."
  )
  if (sensitization == "equivalence") keys <- c(keys, "lakens2017", "schuirmann1987")
  if (longitudinal) {
    keys <- setdiff(c(keys, "mallinckrodt2008", "fitzmaurice2011"),
                    c("mackinnon1985", "long2000", "carpenter2023"))
  }

  # ---- Planned sample ----------------------------------------------------------
  sample_lines <- if (is.null(plan)) {
    c("[The number of participants per group, how it was determined, and the sources of the planning values: the expected treatment effect and sensitization, the pretest-posttest correlation, and the posttest standard deviation. State the smallest effect of interest and why (van 't Veer & Giner-Sorolla, 2016, p. 8). `plan_solomon()` computes the group sizes.]")
  } else {
    s <- attr(plan, "settings")
    rows <- vapply(seq_len(nrow(plan)), function(i) {
      r <- plan[i, ]
      sprintf("%s: %d, %d, %d, and %d participants in the pretested treatment, pretested control, unpretested treatment, and unpretested control groups (%d in all), for power %s against a true contrast of %s (%s%s).",
              r$estimand, r$n1, r$n2, r$n3, r$n4, r$total_n, format(round(r$power, 3)),
              format(r$true_effect), if (identical(r$basis, "simulation")) "simulated" else "normal theory",
              if (!is.na(r$mcse)) sprintf(", Monte Carlo SE %s", format(round(r$mcse, 3))) else "")
    }, "")
    c(
      paste0("Group sizes from `plan_solomon()`, at alpha = ", sub("^0", "", format(plan$alpha[1])), ":"),
      "",
      .plan_bullets(rows),
      "",
      if (!is.null(s)) sprintf("Planning values: a treatment effect among unpretested participants of %s, sensitization of %s, a pretest-posttest correlation of %s, and a posttest standard deviation of %s.",
                               format(s$delta), format(s$sens), format(s$rho), format(s$sigma)),
      if (!is.null(s)) "",
      "[The sources of these planning values, and the smallest effect of interest and why (van 't Veer & Giner-Sorolla, 2016, p. 8).]"
    )
  }

  # ---- Analysis ------------------------------------------------------------------
  tipping_label <- if (is.character(tipping_groups)) {
    paste(tipping_groups, "groups")
  } else {
    paste(.solomon_group_labels[sort(tipping_groups)], collapse = ", ")
  }
  sens_analysis <- switch(
    sensitization,
    test = sprintf("H2 is tested with the Pretest x Treatment contrast from the same model, two-sided at alpha = %s.", alpha_txt),
    equivalence = sprintf("H2 is tested with two one-sided tests of the Pretest x Treatment contrast against -%s and %s, each at alpha = %s: `equivalence_solomon(fit, bounds = %s)` (Schuirmann, 1987; Lakens, 2017). Sensitization is called negligible only if both one-sided tests reject.", format(equivalence_bound), format(equivalence_bound), alpha_txt, format(equivalence_bound)),
    larger_pretested = sprintf("H2 is tested with the Pretest x Treatment contrast from the same model, one-sided at alpha = %s in the predicted direction.", alpha_txt),
    smaller_pretested = sprintf("H2 is tested with the Pretest x Treatment contrast from the same model, one-sided at alpha = %s in the predicted direction.", alpha_txt),
    exploratory = "The Pretest x Treatment contrast is reported with its confidence interval as an exploratory result."
  )
  h1_test <- if (direction == "two-sided") {
    sprintf("H1 is tested with the average treatment effect, two-sided at alpha = %s.", alpha_txt)
  } else {
    sprintf("H1 is tested with the average treatment effect at alpha = %s: [two-sided, or one-sided in the predicted direction].", alpha_txt)
  }
  tipping_code <- if (is.character(tipping_groups)) sprintf("\"%s\"", tipping_groups) else
    sprintf("c(%s)", paste(sort(tipping_groups), collapse = ", "))

  text <- c(
    paste("#", title),
    "",
    sprintf("Drafted with solomonR %s on %s. Edit every section before registering it; text in square brackets is for the researcher to complete. The sections follow van 't Veer and Giner-Sorolla's (2016) template, and each names the SPIRIT 2013 item it answers (Chan et al., 2013).",
            as.character(utils::packageVersion("solomonR")), format(Sys.Date())),
    "",
    "## 1. Hypotheses (SPIRIT 12)",
    "",
    sprintf("The primary outcome is %s.", outcome),
    "",
    h1,
    "",
    h2,
    "",
    "## 2. Methods",
    "",
    "### Design and randomization (SPIRIT 16a)",
    "",
    "Participants are randomly assigned to the four groups of a Solomon four-group design (Solomon, 1949), which crosses treatment (treatment or control) with pretesting (pretested or not). The unpretested groups have no pretest by design; that absence is the experimental manipulation, not missing data.",
    "",
    .plan_bullets(c(
      "**Randomization.** [The method of sequence generation, any stratification or blocking, and how allocation is concealed.]",
      "**Pretest-posttest interval.** [The interval.] Pretest effects can depend on it (Entwisle, 1961), so it is fixed in advance.",
      if (longitudinal) sprintf("**Posttest occasions.** %d occasions (%s) [with their times]; the confirmatory contrasts are those at occasion %s.", length(occasion_labels), paste(occasion_labels, collapse = ", "), primary_occasion),
      "**Measurement.** The pretest and posttest are the same instrument, administered in the same way and at the same times in every group (French et al., 2021b, Recommendation 11). [Describe the instrument and its administration.]",
      "**Blinding (SPIRIT 17a).** [Who is blinded to the group assignment, and how.]"
    )),
    "",
    "### Planned sample (SPIRIT 14)",
    "",
    sample_lines,
    "",
    "**Exclusions.** [Criteria for excluding participants or data, fixed in advance.]",
    "",
    "## 3. Analysis plan",
    "",
    "### Primary analysis (SPIRIT 20a)",
    "",
    if (longitudinal) {
      "The four Solomon contrasts at each occasion, and the change in the Pretest x Treatment contrast from the first occasion to the last, are estimated with a mixed model for repeated measures, `fit_solomon_mmrm(y_post, treat, pretested, id, occasion, y_pre)` (Mallinckrodt et al., 2008): occasion, treatment, pretesting, and all their interactions, adjusting for the pretest score among pretested participants at each occasion (Lin, 2013; Newman et al., 1990), with an unstructured covariance estimated separately for pretested and unpretested participants by restricted maximum likelihood, and Kenward-Roger degrees of freedom (Kenward & Roger, 1997, as cited in Fitzmaurice et al., 2011, p. 101). If the unstructured model does not converge, the fallback structures of `fit_solomon_mmrm()` are used, chosen by AIC."
    } else {
      "The four Solomon contrasts are estimated with one linear model for all four groups, `fit_solomon_glm(y_post, treat, pretested, y_pre)`: treatment, pretesting, and their interaction, adjusting for the pretest score among pretested participants (Lin, 2013; Newman et al., 1990), with HC3 standard errors (MacKinnon & White, 1985; Long & Ervin, 2000)."
    },
    "",
    .plan_bullets(c(
      sprintf("**Confirmatory contrasts:** %s%s.", paste(confirmatory, collapse = " and "), at),
      h1_test,
      sens_analysis,
      "**Multiple testing.** [How the confirmatory tests control the error rate, or why no correction is needed.]",
      "**Other contrasts.** The treatment effect within each pretest condition is reported as a secondary result."
    )),
    "",
    "### Checks before the model (SPIRIT 20b)",
    "",
    .plan_bullets(c(
      "`validate_solomon()` confirms the coding of the design and the size of each group.",
      "`check_solomon_missing()` separates pretests absent by design from missing values.",
      "`baseline_solomon()` describes the pretest difference between the two pretested groups. It is reported, not tested for significance, and it does not change the planned model.",
      "[What will be done if an assumption fails, for example markedly unequal variances or a group with too few participants.]"
    )),
    "",
    "### Missing data and the analysis population (SPIRIT 20c)",
    "",
    .plan_bullets(c(
      "**Population.** All randomized participants, analyzed in the groups to which they were assigned.",
      if (longitudinal) {
        "**Main assumption.** Missing posttests, including those of participants who drop out, are missing at random given group, the pretest in the pretested groups, and the posttests observed earlier. The mixed model above is valid under that assumption; per-occasion analyses of the participants still observed are not (Fitzmaurice et al., 2011). Missing pretests among pretested participants are not imputed; [state how they are handled]."
      } else {
        "**Main assumption.** Missing posttests are missing at random given group and, in the pretested groups, the pretest. Under that assumption the model above, fitted to the observed posttests, is the main analysis; multiple imputation under the same assumption agrees with it (Carpenter et al., 2023, p. 256). Missing pretests among pretested participants are not imputed; [state how they are handled]."
      },
      if (longitudinal) {
        "**Sensitivity analysis.** [A sensitivity analysis for departures from missing at random, stated in advance (White et al., 2011).] The package's `tipping_point_solomon()` imputes from the pretest alone, so it does not apply when dropout depends on earlier posttests."
      } else {
        sprintf("**Sensitivity analysis.** `tipping_point_solomon(..., groups = %s)` shifts the imputed posttests of the %s from -1 to 1 standard deviation and reports the offset at which each confirmatory conclusion changes (White et al., 2011; Carpenter et al., 2023).", tipping_code, tipping_label)
      },
      if (sens_confirmatory && !longitudinal) "**Sensitization-specific sensitivity analysis.** `tipping_point_solomon(..., contrast = \"Pretest x Treatment\", groups = 1)` shifts the imputed posttests of the pretested treatment group alone, the departure from missing at random that bears on the sensitization contrast."
    )),
    "",
    "## 4. Deviations from this plan (SPIRIT 25)",
    "",
    "Deviations are recorded here, with the date and the reason, and reported with the results, so that planned and unplanned analyses can be told apart (Nosek et al., 2018, p. 2602). Every confirmatory analysis in this plan is reported, whatever its result.",
    "",
    "| Date | Section | Deviation | Reason |",
    "|---|---|---|---|",
    "| | | | |",
    "",
    "## References",
    ""
  )
  refs <- .apa_sort(unname(.solomon_reference_text[unique(keys)]))
  text <- c(text, unlist(lapply(refs, function(r) c(r, ""))))

  out <- structure(
    list(
      text = text,
      settings = list(
        occasions = occasion_labels, primary_occasion = if (longitudinal) primary_occasion,
        confirmatory = confirmatory, direction = direction, sensitization = sensitization,
        equivalence_bound = equivalence_bound, alpha = alpha, tipping_groups = tipping_groups,
        outcome = outcome, treatment = treatment
      ),
      date = Sys.Date()
    ),
    class = "solomon_analysis_plan"
  )
  if (!is.null(file)) {
    writeLines(enc2utf8(text), file, useBytes = TRUE)
  }
  out
}

#' @export
print.solomon_analysis_plan <- function(x, ...) {
  cat(x$text, sep = "\n")
  invisible(x)
}
