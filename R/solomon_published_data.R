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
