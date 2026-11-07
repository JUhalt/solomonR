# The reference set for perm_solomon() and plot_perm() (issue #134).
#
# perm-reference.rds is written by tools/perm-reference/generate.R, which does
# not load the package. It enumerates the labelings (as compositions of row
# types, with exact multiplicities), so every count is an exact integer, and
# decides each labeling's statistic against the observed one in exact
# rational arithmetic (linear models, saturated cells, cluster scores), or,
# for models fitted by iteration with a covariate in the contrast's cells,
# from fits at the exact maximum-likelihood estimate: from doubles when the
# gap exceeds a margin over the measured error (a margin rule, not a proof),
# otherwise at 50 digits. Equal statistics are found only by exact
# certificates (equal compositions, swaps of the arms, equal sufficient
# statistics, integer relations), never with a tolerance. See
# tools/perm-reference/README.md and the reference-set comment on issue #134.
#
#   Rscript tools/perm-reference/generate.R full <out_dir>
#
# Contents: `cases`, one row per world, observed labeling, contrast,
# statistic, and labeling space, with the expected outcome (`p`,
# `error:<class>`, `excluded:<reason>`, or `no-reference:NB-BND`), the exact
# counts b and M (as text when above 2^53), and the counts under the
# alternative conventions (`_drop`: labelings without a defined statistic
# left out; `_inf`: counted as infinite; `_h1inf`: a leverage-1 labeling whose
# contrast is infinite takes that value), `alt_outcome` (the outcome under a
# pending proposal where it differs, such as D5's), `note`, and, for the
# negative binomial, theta*, glm.nb's theta, and the band over which the
# decisions were certified; and `worlds`, the test tier: for each world its
# data, observed labelings, row types, and per class (composition) the code,
# the statistic as a double, and its relation to the observed class (-1, 0,
# 1; NA where a convention decides), plus `under`, the code a leverage-1 class
# has without the leverage rule, so that a test can classify the labelings
# perm_solomon() draws without fitting.
