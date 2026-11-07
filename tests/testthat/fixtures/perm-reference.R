# The reference set for perm_solomon() and plot_perm() (issue #134).
#
# perm-reference.rds is written by tools/perm-reference/generate.R, which does
# not load the package and computes every reference p-value exactly: by
# enumerating the labelings (as compositions of row types, with exact
# multiplicities), and deciding each labeling's statistic against the
# observed one in exact rational arithmetic (linear models, saturated cells,
# cluster scores), or from fits at the exact maximum-likelihood estimate whose
# order is certified from the measured error or recomputed at 50 digits
# (models fitted by iteration with a covariate in the contrast's cells).
# Equal statistics are found only by exact certificates (equal compositions,
# swaps of the arms, equal sufficient statistics, integer relations), never
# with a tolerance. See tools/perm-reference/README.md and the design posted
# on issue #134.
#
#   Rscript tools/perm-reference/generate.R full <out_dir>
#
# Contents: `cases`, one row per world, observed labeling, contrast,
# statistic, and labeling space, with the expected outcome (`p`, or
# `error:<class>`), the exact counts b and M (as text when above 2^53), and
# the counts under the alternative conventions (`_drop`: labelings without a
# defined statistic left out; `_inf`: counted as infinite); and `worlds`, the
# test tier: for each world its data, observed labelings, row types, and per
# class (composition) the code, the statistic as a double, and its relation
# to the observed class (-1, 0, 1; NA where a convention decides), so that a
# test can classify the labelings perm_solomon() draws without fitting.
