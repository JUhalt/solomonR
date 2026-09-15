# Generates `solomon_example`, the primary teaching data set for solomonR.
#
# A simulated randomized Solomon four-group study with 30 participants per
# group and scores on a 0-100 scale. The seed and all parameters were fixed and
# posted in issue #21 before the data were generated, and the data were not
# regenerated to obtain particular results.
#
# Data-generating mechanism (true values):
#   - latent baseline ability A ~ N(0, 1) for every participant;
#   - pretest (Groups 1 and 2 only) = 50 + 10 * A;
#   - posttest = 50 + 5 * treat + 2 * pretested + 10 * (0.6 * A + 0.8 * e),
#     with e ~ N(0, 1);
#   - so the treatment effect is 5 points in both pretest conditions (no
#     sensitization), the pretest effect is 2 points, the equal-weighted
#     average treatment effect is 5 points, and the pretest-posttest
#     correlation is 0.6.
# Scores are rounded to whole points and kept within 0-100.
#
# Run from the package root: source("data-raw/solomon_example.R")

seed <- 20260915L
n_cell <- 30L

set.seed(seed)

treat <- rep(c(1L, 0L, 1L, 0L), each = n_cell)
pretested <- rep(c(1L, 1L, 0L, 0L), each = n_cell)

ability <- stats::rnorm(4L * n_cell)
noise <- stats::rnorm(4L * n_cell)

pretest <- 50 + 10 * ability
posttest <- 50 + 5 * treat + 2 * pretested + 10 * (0.6 * ability + 0.8 * noise)

score <- function(x) pmin(pmax(round(x), 0), 100)

solomon_example <- data.frame(
  y_post = score(posttest),
  treat = treat,
  pretested = pretested,
  y_pre = ifelse(pretested == 1L, score(pretest), NA_real_)
)

save(solomon_example, file = "data/solomon_example.rda", compress = "bzip2", version = 2)
