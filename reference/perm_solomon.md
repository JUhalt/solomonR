# Permutation test for a Solomon contrast

Performs a randomization-based test by permuting treatment assignment
within pretest strata. This preserves the Solomon four-group design
while generating the null distribution for a selected treatment
contrast.

## Usage

``` r
perm_solomon(
  object,
  contrast = "ATE (avg over pretest)",
  reps = 5000L,
  seed = NULL,
  return_dist = FALSE
)
```

## Arguments

- object:

  An object returned by
  [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md).

- contrast:

  Character string identifying the contrast to test. One of
  `"ATE (avg over pretest)"`, `"Pretest x Treatment"`,
  `"Treatment | pretested"`, or `"Treatment | unpretested"`.

- reps:

  Number of permutations. Default is 5000.

- seed:

  Optional random-number seed for reproducibility.

- return_dist:

  Logical. If `TRUE`, return the permutation distribution in addition to
  the observed statistic and p-value.

## Value

A list containing the observed studentized statistic (`z_obs`) and
permutation p-value (`p_perm`). If `return_dist = TRUE`, the permutation
distribution (`z_perm`) is also returned.
