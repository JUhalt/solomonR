# Published Solomon results stored as small data sets (issue #54).
#
# Numbers reported in the publications, reused with citation. Run from the
# package root:
#
#   Rscript data-raw/published_studies.R

groups <- c("Pretested, treatment", "Pretested, control",
            "Unpretested, treatment", "Unpretested, control")

# El Karkri, Quesada, and Romero-Ariza (2025a), Table 8 (p. 11): posttest
# statistics of the four Solomon groups.
elkarkri2025a <- data.frame(
  group = factor(groups, levels = groups),
  pretested = c(1L, 1L, 0L, 0L),
  treat = c(1L, 0L, 1L, 0L),
  n = c(9L, 25L, 17L, 37L),
  mean = c(10.94, 7.80, 8.94, 9.35),
  sd = c(2.26, 2.29, 1.98, 2.11)
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

save(elkarkri2025a, file = "data/elkarkri2025a.rda", compress = "xz")
save(kvalem1996, file = "data/kvalem1996.rda", compress = "xz")
