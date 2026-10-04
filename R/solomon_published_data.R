#' Pretest and posttest statistics from El Karkri et al. (2025a)
#'
#' The posttest sample size, mean, and standard deviation of the four Solomon
#' groups in a classroom study of the Cognitive Acceleration through Science
#' Education programme (El Karkri et al., 2025a, Table 8, p. 11), and the
#' pretest mean and standard deviation of the two pretested groups (Table 7,
#' p. 10). Numbers reported in the publication are reused with citation.
#'
#' @details
#' **Design and caveats.** The authors describe the study as
#' quasi-experimental: each Solomon group was one intact class, so class and
#' condition are confounded, and differences between the groups can reflect
#' the classes as well as the treatment and the pretest. [validate_solomon()]
#' flags this design when class membership is supplied. The pretested classes
#' already differed at pretest (9.61 against 7.86), which
#' [baseline_solomon()] reports; the unpretested classes have no pretest.
#' The authors reported a significant Pretest x Treatment interaction.
#'
#' **Known result.** [solomon_from_summary()] reproduces the published
#' two-way ANOVA from these numbers within rounding: interaction F(1, 84) =
#' 11.46 against the published 11.482, treatment 6.78 against 6.794, and
#' pretest 0.18 against 0.186 (pp. 11-12).
#'
#' @format A data frame with 4 rows, one per Solomon group, and 8 variables:
#' \describe{
#'   \item{group}{The Solomon group.}
#'   \item{pretested, treat}{Indicators (1 = yes).}
#'   \item{n, mean, sd}{Posttest sample size, mean, and standard deviation.}
#'   \item{pre_mean, pre_sd}{Pretest mean and standard deviation (pretested
#'     groups only; the pretest sample sizes equal the posttest ones).}
#' }
#'
#' @source El Karkri, M., Quesada, A., & Romero-Ariza, M. (2025a). The dual
#' impact of pretest sensitisation and the cognitive acceleration through
#' science education programme in the Solomon four-group design. *Brain
#' Sciences, 16*(1), Article 64. https://doi.org/10.3390/brainsci16010064
#'
#' @seealso [solomon_from_summary()], [baseline_solomon()], [kvalem1996]
#'
#' @examples
#' with(elkarkri2025a, solomon_from_summary(n, mean, sd))
#' pre <- elkarkri2025a[elkarkri2025a$pretested == 1, ]
#' baseline_solomon(n = pre$n, mean = pre$pre_mean, sd = pre$pre_sd)
"elkarkri2025a"

#' Condom use in the Solomon study of Kvalem et al. (1996)
#'
#' The number of students who used condoms at their most recent intercourse,
#' at the 6-month posttest, among students who had had intercourse before the
#' intervention, in each Solomon group of a school-based sex-education trial
#' (Kvalem et al., 1996, p. 42). Numbers reported in the publication are
#' reused with citation.
#'
#' @details
#' **Design and caveats.** Whole classes were randomized: 30 classes to the
#' intervention, 15 of them pretested, and 94 to control, 47 of them
#' pretested (p. 38). The authors analyzed individuals rather than classes
#' because class sizes varied too much for a class-level analysis (p. 39).
#' These counts are therefore individual-level, and analyses of them ignore
#' the clustering; see [perm_solomon()] for cluster-level inference when
#' class membership is available.
#'
#' **Known result.** Fisher's exact and chi-square comparisons by
#' [fisher_solomon()] reproduce the reported pattern: 70% (51/73) against
#' 51% (76/148) among pretested students, chi-square = 6.85, and 43% (21/49)
#' against 52% (69/133) among unpretested students, chi-square = 1.17
#' (p. 42).
#'
#' @format A data frame with 4 rows, one per Solomon group, and 5 variables:
#' \describe{
#'   \item{group}{The Solomon group.}
#'   \item{pretested, treat}{Indicators (1 = yes).}
#'   \item{events}{Students who used condoms at their most recent intercourse.}
#'   \item{n}{Students in the group with intercourse before the intervention
#'     who answered the 6-month posttest.}
#' }
#'
#' @source Kvalem, I. L., Sundet, J. M., Rivø, K. I., Eilertsen, D. E., &
#' Bakketeig, L. S. (1996). The effect of sex education on adolescents' use
#' of condoms: Applying the Solomon four-group design. *Health Education
#' Quarterly, 23*(1), 34–47. https://doi.org/10.1177/109019819602300103
#'
#' @seealso [fisher_solomon()], [marginal_solomon()], [elkarkri2025a]
#'
#' @examples
#' # One row per student.
#' students <- kvalem1996[rep(seq_len(4), kvalem1996$n), c("pretested", "treat")]
#' students$used <- unlist(Map(function(e, n) rep(1:0, c(e, n - e)),
#'                             kvalem1996$events, kvalem1996$n))
#' with(students, fisher_solomon(used, treat, pretested))
"kvalem1996"

#' Transfer-intervention data from Mai et al. (2020)
#'
#' Individual-level data from a Solomon six-group study of interventions to
#' support the transfer of time-management training: pretesting crossed with
#' relapse prevention, proximal plus distal goal setting, and a control
#' condition (Mai et al., 2020). The data were published as the article's
#' supplementary material and are redistributed under its Creative Commons
#' Attribution 4.0 license (see License).
#'
#' @details
#' **Design.** Final-year students at a management college took a
#' three-hour time-management workshop.
#' - **Pretesting.** Before the workshop, they were "randomly assigned into
#'   either the pretested or unpretested groups based on their roll numbers"
#'   (p. 5). The article does not say how roll numbers were used.
#' - **Interventions.** After it, they were randomly assigned to relapse
#'   prevention (RP), goal setting (GS), or control (p. 6).
#' - **Measures.** Self-reported time-management behavior was measured before
#'   the workshop in the pretested groups and 10 weeks after it in all groups
#'   (p. 7; the note to Table 3 on p. 8 says eight weeks).
#'
#' Every pair of conditions forms a Solomon four-group design. The authors
#' analyzed three such pairs, reusing groups across them. [fit_solomon_glm()]
#' with `control = "Control"` fits one model to all six groups (see the
#' examples and the article "Designs With Several Treatments").
#'
#' **Sample and attrition.** The file has 211 participants, and the article
#' reports 210 (p. 5). 133 answered the posttest, with attrition that
#' differed by group (see the example). The unpretested groups have no
#' pretest by design.
#'
#' **Scale.** `pre_behavior` and `post_behavior` are the authors' scale
#' scores for the 30 self-report items. The article describes the response
#' scale as running from 1 (always) to 5 (never) (p. 7). The value labels in
#' the data file run from 1 (never) to 5 (always). The scores are not the
#' simple means of the 30 items in the file, so the items are not included.
#'
#' **Known results.** The two-by-two ANOVAs of the posttest scores reproduce
#' the authors' Table 4 (p. 8). For relapse prevention against control, the
#' Pretest x Treatment sum of squares is 0.508, the error sum of squares
#' 13.347, and F = 3.461. The history checks of their Table 5 (p. 8) and the
#' ANCOVAs of the first two comparisons in their Table 7 (p. 9) are
#' reproduced as well. The third ANCOVA, relapse prevention against goal
#' setting, gives the published difference of -0.219 with p = .044 against
#' the published .046. The published cell statistics of their Table 6 do not
#' match the data (see issue #53). The article "Worked Example: A Published
#' Solomon Study" carries the relapse-prevention comparison through the
#' package's analyses.
#'
#' @section Acknowledgment:
#' We thank Nu Nu Mai, Yoshi Takahashi, and Mon Mon Oo for making their data
#' publicly available with their article. For the study's design, measures,
#' and findings, read the published article (Mai et al., 2020, cited under
#' Source).
#'
#' @section License:
#' The data are © 2020 by the authors (Nu Nu Mai, Yoshi Takahashi, and Mon Mon
#' Oo) and are published with the
#' article under the Creative Commons Attribution 4.0 International license
#' (<https://creativecommons.org/licenses/by/4.0/>). MDPI's open access
#' policy (<https://www.mdpi.com/openaccess>) states that its articles,
#' "including their data, graphics, and
#' supplementary material", may be reused with attribution. Changes made in
#' solomonR: 6 of the 117 variables in the published file are kept, renamed,
#' and recoded to factors or 0/1 indicators, and an `id` column is added. The
#' original file and the script that builds this data set are in the
#' package's `data-raw` folder on GitHub. See also `inst/COPYRIGHTS`.
#'
#' @format A data frame with 211 rows, one per participant, and 6 variables:
#' \describe{
#'   \item{id}{Row number in the published file.}
#'   \item{gender}{Female or Male.}
#'   \item{pretested}{1 if assigned to be pretested, 0 if not.}
#'   \item{condition}{`RP` (relapse prevention), `GS` (proximal plus distal
#'     goal setting), or `Control`.}
#'   \item{pre_behavior}{Self-reported time-management behavior before the
#'     workshop (pretested groups only).}
#'   \item{post_behavior}{Self-reported time-management behavior after the
#'     workshop (see Details for the timing); missing for participants who
#'     did not answer.}
#' }
#'
#' @source Mai, N. N., Takahashi, Y., & Oo, M. M. (2020). Testing the
#' effectiveness of transfer interventions using Solomon four-group designs.
#' *Education Sciences, 10*(4), Article 92.
#' https://doi.org/10.3390/educsci10040092 (supplementary material, Table S1:
#' dataset.sav)
#'
#' @seealso [baseline_solomon()], [fit_solomon_glm()], [elkarkri2025a],
#'   [kvalem1996]
#'
#' @examples
#' # Posttest completion by group.
#' with(mai2020, table(pretested, condition, posttest = !is.na(post_behavior)))
#'
#' # Relapse prevention against control, as a Solomon four-group design.
#' rp <- subset(mai2020, condition %in% c("RP", "Control"))
#' with(rp, fit_solomon_glm(post_behavior, treat = as.integer(condition == "RP"),
#'                          pretested = pretested, pretest_score = pre_behavior))
#'
#' # One model for all six groups.
#' fit_solomon_glm(post_behavior, condition, pretested, pre_behavior,
#'                 control = "Control", data = mai2020)
"mai2020"

#' Solomon's (1949) spelling experiment
#'
#' Group sizes, means, and standard errors from the first published
#' experiment with the Solomon design: spelling lessons in a fifth-grade and
#' a sixth-grade class, analyzed with the three-group design (Solomon, 1949,
#' Tables II and III, pp. 144–145). Numbers reported in the publication are
#' reused with citation.
#'
#' @details
#' **Design.** In each class, pupils were assigned to three groups that were
#' "roughly equated in spelling ability by means of teachers' judgments"
#' (p. 144), not randomized.
#' - **Experimental group.** Pretested on a list of words, given a standard
#'   spelling lesson, and posttested on the same words.
#' - **Control Group I.** Pretested and posttested, without the lesson.
#' - **Control Group II.** Given the lesson and the posttest, without the
#'   pretest.
#'
#' There was no fourth group; Solomon introduced it later in the article for
#' field studies (p. 147).
#'
#' **Known result.** [fit_solomon_1949()] reproduces the published values
#' from these means: an inferred pretest of 3.0 and 5.7, improvements for
#' Control Group II of 8.2 and 8.7, and interactions I = -2.2 (grade 5) and
#' -3.1 (grade 6). The pretest reduced the effect of the lesson, so the usual
#' two-group design would have underrated it (p. 145).
#'
#' **Caveats.** The groups are small (8 to 10 pupils) and were not randomized.
#' Solomon printed the standard error of Control Group II's improvement as
#' "?", since it depends on the inferred pretest.
#'
#' @format A data frame with 6 rows, one per group and grade, and 11
#'   variables:
#' \describe{
#'   \item{grade}{5 or 6.}
#'   \item{group}{Experimental, Control I, or Control II.}
#'   \item{pretested, treat}{Indicators (1 = yes); `treat` is the spelling
#'     lesson.}
#'   \item{n}{Group size.}
#'   \item{pre_mean, pre_se}{Pretest mean and its standard error (pretested
#'     groups only).}
#'   \item{mean, se}{Posttest mean and its standard error.}
#'   \item{change, change_se}{Mean improvement, as printed, and its standard
#'     error (missing for Control II, printed as "?").}
#' }
#'
#' The standard errors are the values printed with a plus-or-minus sign,
#' which Solomon's notation (sigma_m, Table I) identifies as standard errors
#' of the means.
#'
#' @source Solomon, R. L. (1949). An extension of control group design.
#' *Psychological Bulletin, 46*(2), 137–150. https://doi.org/10.1037/h0062958
#'
#' @seealso [fit_solomon_1949()], [elkarkri2025a], [mai2020]
#'
#' @examples
#' solomon1949
#' g6 <- solomon1949[solomon1949$grade == 6, ]
#' fit_solomon_1949(post_mean = g6$mean, pre_mean = g6$pre_mean[1:2], n = g6$n)
"solomon1949"

#' Posttest statistics from Lana's (1959) attitude experiment
#'
#' The sample size, posttest mean, and standard deviation of the four Solomon
#' groups in one of the first experiments built on the Solomon design: a
#' recorded pro-vivisection talk and a 10-item vivisection questionnaire
#' given to introductory psychology classes (Lana, 1959, Table 2, p. 297).
#' Numbers reported in the publication are reused with citation.
#'
#' @details
#' **Design.** Five classes (156 students) were randomly assigned to five
#' conditions (p. 295). Four form a Solomon design:
#' - **Group I.** Pretest, then the talk 12 days later, then the posttest.
#' - **Group IV.** Pretest, then the posttest 12 days later.
#' - **Group II.** The talk, then the posttest.
#' - **Group III.** The posttest only.
#'
#' A fifth group, with a delayed posttest, is not part of the four-group
#' design and is not included. Because whole classes were assigned, class
#' and condition are confounded, as in [elkarkri2025a].
#'
#' **Known result.** Lana analyzed the posttest means as a two-by-two
#' analysis of variance (Table 3, p. 297). The treatment effect was
#' significant, F = 5.35, and the pretest and the interaction were not.
#' [solomon_from_summary()] reproduces it within rounding (treatment
#' F = 5.36 on 1 and 152 df). Lana's sums of squares are on the scale of the
#' cell means (an unweighted-means analysis). Dividing the package's sums of
#' squares by the harmonic mean of the cell sizes recovers them, and the F
#' tests are the same.
#'
#' **Historical note.** Lana concluded that the pretest did not sensitize
#' participants to the talk (p. 298), the opposite of Solomon's (1949)
#' spelling result. Together they begin the record that Willson and Putnam
#' (1982) later pooled.
#'
#' @format A data frame with 4 rows, one per Solomon group, and 7 variables:
#' \describe{
#'   \item{group}{The Solomon group.}
#'   \item{lana_group}{Lana's group label (I to IV).}
#'   \item{pretested, treat}{Indicators (1 = yes); `treat` is the recorded
#'     talk.}
#'   \item{n, mean, sd}{Posttest sample size, mean, and standard deviation.}
#' }
#'
#' @source Lana, R. E. (1959). Pretest-treatment interaction effects in
#' attitudinal studies. *Psychological Bulletin, 56*(4), 293–300.
#' https://doi.org/10.1037/h0044646
#'
#' @seealso [solomon_from_summary()], [solomon1949], [elkarkri2025a]
#'
#' @examples
#' lana1959
#' with(lana1959, solomon_from_summary(n, mean, sd))
"lana1959"

#' Group statistics from Steyn's (2005) eight-group study
#'
#' The sample size, pretest and posttest means, and standard deviations of the
#' eight groups of a Solomon design with three treatments: a study of how
#' information about one's own ability changes self-efficacy perceptions, with
#' 1,723 police trainees (Steyn, 2005, Table 5.59, pp. 151–152). It is the
#' eight-group study that Steyn (2009) describes, and Steyn and Mynhardt
#' (2008) report it in English. Numbers reported in the publication are
#' reused with citation.
#'
#' @details
#' **Design.** Three treatments and a control, each with and without a
#' pretest (pp. 103–105). Each treatment added a source of information about
#' the participant's ability (p. 102):
#' - **Test.** Completing a 60-item cognitive test.
#' - **Marking.** Completing the test and marking one's own answers.
#' - **Norms.** Completing and marking the test, and receiving the test's
#'   norms.
#'
#' The outcome is the total score on a 45-item questionnaire of self-efficacy
#' perceptions (p. 94). The posttest followed the treatment after a break of
#' 10 minutes (p. 107).
#'
#' **Assignment.** Participants were not randomized individually. Fifty-six
#' existing classes of a police training college were allocated to the eight
#' groups, seven classes to each, in consultation with the college's
#' management, and the thesis states that this allocation was not random
#' (pp. 105–106). The classes had been formed from the order in which
#' trainees reported for training. The thesis is not consistent on this
#' point. Its introduction says that participants were divided at random
#' into eight groups and that the groups were assigned at random to the
#' conditions (p. 9), and it calls the allocation random again on p. 107.
#' So does the English report, which does not mention the classes (Steyn &
#' Mynhardt, 2008, pp. 566–567).
#' The method chapter, followed here, gives the detail: the researcher
#' judged the existing classes to be random groups and made no further
#' random assignment (pp. 105–106), and balanced the groups by sex: each
#' held five men's classes and two women's classes (Steyn, 2005, p. 106).
#' The analyses in the thesis, and those
#' below, treat participants as the units, so they do not allow for the
#' classes. See [baseline_solomon()] for nonrandomized designs and
#' [validate_solomon()] for clustered ones.
#'
#' **Known results.** These statistics reproduce the analyses in the thesis:
#' - **One-way ANOVA of the eight posttest groups** (Table 5.60, p. 152):
#'   F(7, 1715) = 4.545.
#' - **Scheffé's (1953) tests** (Table 5.61, p. 153), within .001. Only the
#'   unpretested Test group differs from the two control groups (p = .002
#'   and p = .011).
#' - **The 2 x 2 ANOVAs of each treatment against the control** (Tables 5.21,
#'   5.34, and 5.47; pp. 128, 135, 142; Steyn & Mynhardt, 2008, Table 2,
#'   p. 569), within rounding. For Test: intervention F = 21.3, pretest
#'   F = 4.5 (p = .033), and interaction F = 2.0 (p = .152). For Marking:
#'   9.3 (p = .002), 0.4 (p = .538), and 0.0 (p = .929). For Norms, the thesis ran this analysis with 213
#'   participants in the unpretested group (Table 5.45, pp. 141–142), one
#'   fewer than in Table 5.59, so its error degrees of freedom are 854 and
#'   those computed from these statistics are 855. The F statistics for
#'   Norms (14.0, 0.9, and 0.1) are the same to one decimal.
#'
#' Steyn and Mynhardt (2008) report only these three analyses, judged at
#' the .01 level (p. 568), so they read the pretest effect for Test
#' (p = .033) as no effect (p. 568); the thesis calls it significant at .05
#' (p. 128). Their means of the treated and untreated groups reproduce as
#' unweighted averages of the two groups' means. Their effect size d
#' divides the difference by what they call the mean standard deviation
#' (p. 568): 13.00, 12.696, and 12.803 (Table 2, p. 569), the "Totaal"
#' rows of the thesis (Tables 5.19, 5.32, and 5.45; pp. 127, 134, 142). For
#' Norms, their means and that standard deviation follow the 213
#' participants of Table 5.45, not the 214 of Table 5.59 and of their own
#' Table 1 (p. 567). Their text gives the Marking pretest p as .583
#' (p. 568); their Table 2 and these statistics give .538.
#'
#' The treatments lowered the scores. The thesis analyzed the design as
#' overlapping four-group designs, one with the treatments pooled and one for
#' each treatment, and then as a one-way analysis of variance of the eight
#' posttests (pp. 103–105). The joint model of [solomon_from_summary()] tests
#' the Pretest x Condition interaction once, F(3, 1715) = 1.00, p = .392.
#'
#' **Language.** The thesis is in Afrikaans. The condition labels are
#' solomonR's.
#'
#' @format A data frame with 8 rows, one per group, and 8 variables:
#' \describe{
#'   \item{group}{The group's label in the thesis: `EG1` to `EG3` are the
#'     pretested treatment groups, `KG1.1` to `KG1.3` the unpretested
#'     treatment groups, `KG2` the pretested control, and `KG3` the
#'     unpretested control.}
#'   \item{condition}{`Norms`, `Marking`, `Test`, or `Control`.}
#'   \item{pretested}{1 if the group was pretested, 0 if not.}
#'   \item{n}{Sample size.}
#'   \item{pre_mean, pre_sd}{Pretest mean and standard deviation (pretested
#'     groups only).}
#'   \item{mean, sd}{Posttest mean and standard deviation.}
#' }
#'
#' @source Steyn, R. (2005). *Self-evaluasie en die vorming van
#' selfdoeltreffendheidspersepsies* \[Self-evaluation and the forming of
#' self-efficacy perceptions\] \[Doctoral thesis, University of South Africa\].
#' Unisa Institutional Repository. https://hdl.handle.net/10500/1745
#'
#' @references
#' Scheffé, H. (1953). A method for judging all contrasts in the analysis of
#' variance. *Biometrika, 40*(1–2), 87–104.
#' https://doi.org/10.1093/biomet/40.1-2.87
#'
#' Steyn, R. (2009). Re-designing the Solomon four-group: Can we improve on
#' this exemplary model? *Design Principles and Practices: An International
#' Journal—Annual Review, 3*(1), 383–394. https://doi.org/10.18848/1833-1874/CGP/v03i01/37588
#'
#' Steyn, R., & Mynhardt, J. (2008). Factors that influence the forming of
#' self-evaluation and self-efficacy perceptions. *South African Journal of
#' Psychology, 38*(3), 563–573. https://doi.org/10.1177/008124630803800310
#'
#' @seealso [solomon_from_summary()], [fit_solomon_steyn()], [mai2020],
#'   [lana1959]
#'
#' @examples
#' steyn2005
#'
#' # One model for all eight groups.
#' with(steyn2005, solomon_from_summary(n, mean, sd, treat = condition,
#'                                      pretested = pretested, control = "Control"))
"steyn2005"

#' Group statistics from Jordaan's (2014) study with three posttest occasions
#'
#' The sample sizes, means, and standard deviations of a Solomon four-group
#' design measured on three posttest occasions: an evaluation of a life
#' skills program for young adult male offenders (Jordaan, 2014, Table 7.3,
#' p. 112). The outcomes are the three subscales of the Coping Strategy
#' Indicator. Numbers reported in the publication are reused with citation.
#'
#' @details
#' **Design.** 120 offenders aged 21 to 25 with long sentences were selected
#' by systematic random sampling in a maximum-security correctional center
#' and assigned at random to the program or the control condition. Half of
#' each condition was assigned at random to be pretested (pp. 86–87). The
#' program ran for six months; the control group followed the center's
#' normal daily activities (p. 87). All groups were tested after the
#' program and again 3 and 6 months later (p. 98).
#'
#' **Attrition.** Transfers removed 17 offenders from the program groups and
#' 7 from the control groups, leaving 96 (p. 109). Those 96 completed every
#' posttest: every analysis in the thesis has 92 error degrees of freedom.
#' Because more were lost from the program groups, the groups analyzed are
#' not guaranteed to be comparable, even though they were randomized.
#'
#' **Measures.** Each subscale of the Coping Strategy Indicator has 11 items
#' scored from 1 (not at all) to 3 (a lot), so scores range from 11 to 33.
#' High scores on problem solving and on seeking social support, and low
#' scores on avoidance, indicate better coping (p. 92).
#'
#' **Known results.** These statistics reproduce, within rounding, the
#' 2 x 2 analyses of variance of the posttests on each occasion (Tables
#' 7.4–7.19, pp. 113–127). The Pretest x Treatment interaction F values are:
#' - **Social support:** 9.678 (p = .002) after the program, 0.266 at 3
#'   months, and 2.306 at 6 months.
#' - **Problem solving:** 0.819, 0.563, and 5.556 (p = .021).
#' - **Avoidance:** 0.373, 0.688, and 0.327.
#'
#' The thesis followed the decision sequence of Walton Braver and Braver
#' (1988) on each occasion separately ([fit_solomon_classic()]). The
#' correlations between occasions are not reported, so analyses of change
#' across occasions cannot be reproduced; see the article "Worked Example:
#' Repeated Posttests".
#'
#' @format A data frame with 42 rows, one per subscale, group, and occasion,
#'   and 8 variables:
#' \describe{
#'   \item{subscale}{`Social support`, `Problem solving`, or `Avoidance`.}
#'   \item{group}{The group's number in the thesis: 1 = program with
#'     pretest, 2 = program without pretest, 3 = control with pretest,
#'     4 = control without pretest.}
#'   \item{treat}{1 for the program groups, 0 for the control groups.}
#'   \item{pretested}{1 if the group was pretested, 0 if not.}
#'   \item{occasion}{`Pretest` (groups 1 and 3 only), `Posttest` (after the
#'     program), `Follow-up 1` (3 months later), or `Follow-up 2` (6 months
#'     later).}
#'   \item{n}{Sample size.}
#'   \item{mean, sd}{Mean and standard deviation.}
#' }
#'
#' @source Jordaan, J. (2014). *The development and evaluation of a life
#' skills programme for young adult prisoners* \[Doctoral thesis, University
#' of the Free State\]. KovsieScholar. https://hdl.handle.net/11660/832
#'
#' @references
#' Walton Braver, M. C., & Braver, S. L. (1988). Statistical treatment of the
#' Solomon four-group design: A meta-analytic approach. *Psychological
#' Bulletin, 104*(1), 150–154. https://doi.org/10.1037/0033-2909.104.1.150
#'
#' @seealso [solomon_from_summary()], [fit_solomon_mmrm()], [steyn2005]
#'
#' @examples
#' # Problem solving, 6 months after the program (Table 7.13, p. 121).
#' ps6 <- subset(jordaan2014, subscale == "Problem solving" & occasion == "Follow-up 2")
#' with(ps6, solomon_from_summary(n, mean, sd, treat = treat, pretested = pretested))
"jordaan2014"

#' The worked example of Walton Braver and Braver (1988)
#'
#' The hypothetical data with which Walton Braver and Braver (1988, Table 3,
#' p. 153) illustrate their sequence of Tests A–I: the posttest means and
#' variances of four groups of 14, and the pretest means and variances and the
#' pretest-posttest correlations of the two pretested groups. Numbers reported
#' in the publication are reused with citation.
#'
#' @details
#' **The table.** Table 3 prints variances, not standard deviations; `sd` and
#' `pre_sd` are their square roots. The two pretest-posttest correlations are
#' printed between the rows of Groups 1 and 2 and of Groups 2 and 3. `r`
#' reads .58 as Group 1's and .62 as Group 2's, the reading that reproduces
#' the published analysis of covariance (Table 5, p. 153); the other reading
#' gives F = 2.96. As the authors note, the means show a treatment effect
#' without pretest sensitization (p. 153).
#'
#' **Known results.** The package reproduces the worked example (p. 153):
#' - **The 2 x 2 analysis of variance** (Table 4): [solomon_from_summary()]
#'   gives the published mean squares (pretest .14, treatment 67.76,
#'   interaction 0, error 20.00 on 52 df). For the treatment (Test D),
#'   F = 3.388 and p = .0714; for the pretest, F = .007 and p = .9336 against
#'   the printed .9337; for the interaction (Test A), F = 0 and p = 1 against
#'   the printed .999.
#' - **Test E** (Table 5): treatment mean square 37.74, error mean square
#'   12.87 on 25 df, F = 2.93, and p = .0992 against the printed .0993.
#' - **Test H:** t(26) = 1.28, p = .2126 against the printed .2127.
#' - **Test I:** the one-tailed p-values of Tests E and H give z = 1.65 and
#'   1.25, and the combined z = 2.05. The printed p = .040 is the two-tailed
#'   p of that z; the unrounded z = 2.047 gives .041.
#' - **The 1988 sequence:** Tests A, D, E, and H are not significant and
#'   Test I is, the path that [fit_solomon_classic()] returns.
#' - **The power remark:** Groups 3 and 4 alone, with 28 participants each,
#'   would give t(54) = 1.807 (two-tailed p = .076), as the authors state
#'   (p. 153).
#'
#' [fit_solomon_classic()] takes individual scores. Scores built to have
#' exactly these means, variances, and correlations give the published
#' values, because Tests A–I depend on the data only through them (see the
#' examples). The article reports no results for Tests B, C, F, and G on
#' these data.
#'
#' **Caveats.** The sequence reaches Test I only when every earlier test is
#' nonsignificant, and the authors read the significance level at each step
#' as conditional on reaching that step (p. 153, note 3). Sawilowsky and
#' Markman (1988, pp. 3–4; Table 2, p. 7) constructed data in which Test H is
#' significant and their computation of Test I is not. Braver and Walton
#' Braver (1990, p. 322) replied and amended the sequence (`flow = "1990"`).
#' In a Monte Carlo study, the 1988 sequence falsely declared an effect about
#' 14% of the time at a nominal 5% per test (Sawilowsky et al., 1994,
#' p. 368). See [fit_solomon_classic()].
#'
#' @format A data frame with 4 rows, one per Solomon group, and 11 variables:
#' \describe{
#'   \item{group}{The Solomon group; the rows are the article's Groups 1 to 4.}
#'   \item{pretested, treat}{Indicators (1 = yes).}
#'   \item{n}{Group size.}
#'   \item{mean, var, sd}{Posttest mean, variance as printed, and its square
#'     root.}
#'   \item{pre_mean, pre_var, pre_sd}{Pretest mean, variance as printed, and
#'     its square root (pretested groups only).}
#'   \item{r}{Pretest-posttest correlation (pretested groups only).}
#' }
#'
#' @source Walton Braver, M. C., & Braver, S. L. (1988). Statistical treatment
#' of the Solomon four-group design: A meta-analytic approach. *Psychological
#' Bulletin, 104*(1), 150–154. https://doi.org/10.1037/0033-2909.104.1.150
#'
#' @references
#' Braver, S. L., & Walton Braver, M. C. (1990). Meta-analysis for Solomon
#' four-group designs reconsidered: A reply to Sawilowsky and Markman.
#' *Perceptual and Motor Skills, 71*(1), 321–322.
#' https://doi.org/10.2466/pms.1990.71.1.321
#'
#' Sawilowsky, S. S., Kelley, D. L., Blair, R. C., & Markman, B. S. (1994).
#' Meta-analysis and the Solomon four-group design. *The Journal of
#' Experimental Education, 62*(4), 361–376.
#' https://doi.org/10.1080/00220973.1994.9944140
#'
#' Sawilowsky, S. S., & Markman, B. S. (1988). *Another look at the power of
#' meta-analysis in the Solomon four-group design* (ED316556). ERIC.
#' https://eric.ed.gov/?id=ED316556
#'
#' @seealso [fit_solomon_classic()], [stouffer_solomon()],
#'   [solomon_from_summary()], [lana1959]
#'
#' @examples
#' waltonbraver1988
#'
#' # Tests A-D from the posttest statistics (Table 4).
#' with(waltonbraver1988, solomon_from_summary(n, mean, sd))
#'
#' # Tests A-I and the 1988 sequence, from scores built to have exactly the
#' # published means, variances, and pretest-posttest correlations.
#' set.seed(1988)
#' scores <- do.call(rbind, lapply(seq_len(4), function(i) {
#'   g <- waltonbraver1988[i, ]
#'   if (g$pretested == 1) {
#'     s <- g$r * g$pre_sd * g$sd
#'     x <- MASS::mvrnorm(g$n, c(g$pre_mean, g$mean),
#'                        matrix(c(g$pre_var, s, s, g$var), 2), empirical = TRUE)
#'   } else {
#'     x <- cbind(NA, g$mean + g$sd * as.vector(scale(stats::rnorm(g$n))))
#'   }
#'   data.frame(treat = g$treat, pretested = g$pretested, y_pre = x[, 1], y_post = x[, 2])
#' }))
#' fit <- fit_solomon_classic(y_post, treat, pretested, y_pre, data = scores)
#' fit$path_string
#' fit$tests$I$result[, c("z", "p.value")]
"waltonbraver1988"
