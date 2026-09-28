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

save(elkarkri2025a, file = "data/elkarkri2025a.rda", compress = "xz")
save(mai2020, file = "data/mai2020.rda", compress = "xz")
save(kvalem1996, file = "data/kvalem1996.rda", compress = "xz")
