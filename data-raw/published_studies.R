# Published Solomon results stored as small data sets (issue #54).
#
# Numbers reported in the publications, reused with citation. Run from the
# package root:
#
#   Rscript data-raw/published_studies.R

groups <- c("Pretested, treatment", "Pretested, control",
            "Unpretested, treatment", "Unpretested, control")

# El Karkri, Quesada, and Romero-Ariza (2025a), Table 8 (p. 11): posttest
# statistics of the four Solomon groups; Table 7 (p. 10): pretest statistics
# of the two pretested groups.
elkarkri2025a <- data.frame(
  group = factor(groups, levels = groups),
  pretested = c(1L, 1L, 0L, 0L),
  treat = c(1L, 0L, 1L, 0L),
  n = c(9L, 25L, 17L, 37L),
  mean = c(10.94, 7.80, 8.94, 9.35),
  sd = c(2.26, 2.29, 1.98, 2.11),
  pre_mean = c(9.61, 7.86, NA, NA),
  pre_sd = c(2.67, 2.72, NA, NA)
)

# Lana (1959), Table 2 (p. 297): posttest statistics of the four Solomon
# groups in an experiment on attitudes toward vivisection. Lana's Groups I,
# IV, II, and III are, in the package's order, the pretested communication,
# pretested control, unpretested communication, and unpretested control
# groups. Group V (pretest, communication, delayed posttest) is not part of
# the four-group design and is not stored.
lana1959 <- data.frame(
  group = factor(groups, levels = groups),
  lana_group = c("I", "IV", "II", "III"),
  pretested = c(1L, 1L, 0L, 0L),
  treat = c(1L, 0L, 1L, 0L),
  n = c(26L, 32L, 50L, 48L),
  mean = c(42.96, 40.28, 42.90, 40.77),
  sd = c(9.06, 5.48, 5.34, 5.76)
)

# Kvalem, Sundet, Rivo, Eilertsen, and Bakketeig (1996, p. 42): use of condoms
# at the most recent intercourse at the 6-month posttest, among students who
# had had intercourse before the intervention.
kvalem1996 <- data.frame(
  group = factor(groups, levels = groups),
  pretested = c(1L, 1L, 0L, 0L),
  treat = c(1L, 0L, 1L, 0L),
  events = c(51L, 76L, 21L, 69L),
  n = c(73L, 148L, 49L, 133L)
)

# Mai, Takahashi, and Oo (2020): individual-level data published as the
# article's supplementary material ("Table S1: dataset.sav", available at
# https://www.mdpi.com/2227-7102/10/4/92/s1) under the article's Creative
# Commons Attribution 4.0 license. The file is kept in data-raw/. Changes:
# six of its 117 variables are kept, renamed, and recoded as described on the
# help page of mai2020.
mai <- haven::zap_labels(haven::read_sav("data-raw/mai2020-supplementary.sav"))
mai2020 <- data.frame(
  id = seq_len(nrow(mai)),
  gender = factor(mai$Gender, levels = 1:2, labels = c("Female", "Male")),
  pretested = as.integer(mai$TestGroup),
  condition = factor(mai$Conditions, levels = 1:3, labels = c("RP", "GS", "Control")),
  pre_behavior = as.numeric(mai$Prebehavior),
  post_behavior = as.numeric(mai$Postbehavior)
)

# Solomon (1949), Tables II and III (pp. 144-145): the spelling experiment
# with the three-group design, in a fifth-grade and a sixth-grade class. The
# +/- values are printed as standard errors of the means (sigma_m, Table I).
# Control Group II was not pretested; its improvement uses Solomon's inferred
# pretest, and he printed its standard error as "?".
solomon1949 <- data.frame(
  grade = rep(c(5L, 6L), each = 3),
  group = factor(rep(c("Experimental", "Control I", "Control II"), 2),
                 levels = c("Experimental", "Control I", "Control II")),
  pretested = rep(c(1L, 1L, 0L), 2),
  treat = rep(c(1L, 0L, 1L), 2),
  n = c(10L, 10L, 10L, 8L, 9L, 8L),
  pre_mean = c(3.2, 2.8, NA, 5.4, 6.0, NA),
  pre_se = c(0.8, 0.7, NA, 1.4, 1.5, NA),
  mean = c(9.9, 3.5, 11.2, 11.8, 6.8, 14.4),
  se = c(1.6, 0.8, 1.2, 1.0, 1.4, 0.3),
  change = c(6.7, 0.7, 8.2, 6.4, 0.8, 8.7),
  change_se = c(0.9, 0.5, NA, 1.2, 0.6, NA)
)

# Steyn (2005), Table 5.59 (pp. 151-152): the eight groups of a Solomon design
# with three treatments, in a study of self-efficacy perceptions among police
# trainees. The thesis is in Afrikaans; it labels the pretested treatment
# groups EG1 to EG3, the unpretested treatment groups KG1.1 to KG1.3, the
# pretested control KG2, and the unpretested control KG3. The treatments
# (p. 102): EG3 completed a cognitive test, EG2 also marked it, and EG1 also
# received the test's norms. Rows are in the package's group order.
# Reported in English by Steyn and Mynhardt (2008), whose Table 1 (p. 567)
# gives the same group sizes and labels the KG groups CG.
steyn2005 <- data.frame(
  group = c("EG1", "EG2", "EG3", "KG2", "KG1.1", "KG1.2", "KG1.3", "KG3"),
  condition = factor(rep(c("Norms", "Marking", "Test", "Control"), 2),
                     levels = c("Norms", "Marking", "Test", "Control")),
  pretested = rep(c(1L, 0L), each = 4),
  n = c(218L, 214L, 219L, 218L, 214L, 211L, 220L, 209L),
  pre_mean = c(155.422, 155.724, 155.644, 156.991, NA, NA, NA, NA),
  pre_sd = c(12.125, 11.856, 12.784, 11.987, NA, NA, NA, NA),
  mean = c(155.895, 156.196, 156.142, 158.917, 154.827, 155.739, 153.036, 158.306),
  sd = c(13.263, 12.545, 13.396, 13.264, 12.324, 12.830, 12.628, 11.878)
)

# Jordaan (2014), Table 7.3 (p. 112): the three subscales of the Coping
# Strategy Indicator in a Solomon four-group design with three posttest
# occasions: after the six-month program, and 3 and 6 months later (p. 98).
# Groups (Table 7.3 note): 1 = program with pretest, 2 = program without,
# 3 = control with pretest, 4 = control without. Group sizes from Table 7.14
# (p. 122); every analysis in the thesis has 92 error degrees of freedom, so
# all 96 participants have every posttest.
jordaan_one <- function(subscale, pre, post, fu1, fu2) {
  occ <- c("Pretest", "Posttest", "Follow-up 1", "Follow-up 2")
  rows <- rbind(
    data.frame(group = c(1L, 3L), occasion = "Pretest", mean = pre[c(1, 3)], sd = pre[c(2, 4)]),
    data.frame(group = 1:4, occasion = "Posttest", mean = post[c(1, 3, 5, 7)], sd = post[c(2, 4, 6, 8)]),
    data.frame(group = 1:4, occasion = "Follow-up 1", mean = fu1[c(1, 3, 5, 7)], sd = fu1[c(2, 4, 6, 8)]),
    data.frame(group = 1:4, occasion = "Follow-up 2", mean = fu2[c(1, 3, 5, 7)], sd = fu2[c(2, 4, 6, 8)])
  )
  data.frame(subscale = subscale, group = rows$group,
             treat = as.integer(rows$group %in% c(1L, 2L)),
             pretested = as.integer(rows$group %in% c(1L, 3L)),
             occasion = factor(rows$occasion, levels = occ),
             n = c(22L, 21L, 22L, 31L)[rows$group], mean = rows$mean, sd = rows$sd)
}
jordaan2014 <- rbind(
  jordaan_one("Social support", pre = c(27.23, 3.05, 24.09, 4.33),
              post = c(28.86, 3.28, 25.76, 4.89, 24.05, 5.31, 26.74, 4.38),
              fu1 = c(28.59, 2.36, 27.71, 3.51, 26.59, 5.08, 26.55, 4.10),
              fu2 = c(28.00, 3.62, 26.81, 4.08, 26.05, 4.56, 27.39, 3.88)),
  jordaan_one("Problem solving", pre = c(27.00, 3.67, 24.68, 3.67),
              post = c(28.68, 4.10, 29.52, 3.44, 25.82, 6.39, 28.55, 5.47),
              fu1 = c(29.00, 3.69, 27.95, 4.52, 28.00, 5.55, 28.39, 4.57),
              fu2 = c(30.23, 3.49, 28.62, 3.47, 28.64, 4.17, 30.35, 2.65)),
  jordaan_one("Avoidance", pre = c(24.00, 4.34, 23.18, 2.97),
              post = c(21.68, 4.19, 22.10, 3.75, 21.68, 3.53, 23.06, 3.84),
              fu1 = c(23.32, 5.28, 21.24, 3.39, 23.32, 4.39, 22.65, 3.31),
              fu2 = c(21.73, 4.44, 21.90, 3.21, 22.82, 4.31, 22.10, 3.27))
)
jordaan2014$subscale <- factor(jordaan2014$subscale,
                               levels = c("Social support", "Problem solving", "Avoidance"))
rownames(jordaan2014) <- NULL

save(elkarkri2025a, file = "data/elkarkri2025a.rda", compress = "xz")
save(mai2020, file = "data/mai2020.rda", compress = "xz")
save(kvalem1996, file = "data/kvalem1996.rda", compress = "xz")
save(solomon1949, file = "data/solomon1949.rda", compress = "xz")
save(lana1959, file = "data/lana1959.rda", compress = "xz")
save(steyn2005, file = "data/steyn2005.rda", compress = "xz")
save(jordaan2014, file = "data/jordaan2014.rda", compress = "xz")
