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
#' analyzed three such pairs, reusing groups across them. A single model
#' for all six groups is planned in issue #45.
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
"mai2020"
