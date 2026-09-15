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
  [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md)
  without a clustering variable.

- contrast:

  Character string identifying the contrast to test. One of
  `"ATE (avg over pretest)"`, `"Pretest x Treatment"`,
  `"Treatment | pretested"`, or `"Treatment | unpretested"`.

- reps:

  Number of permutations. Default is 5000.

- seed:

  Optional random-number seed for reproducibility. The global
  random-number state is restored when the function exits.

- return_dist:

  Logical. If `TRUE`, return the permutation distribution in addition to
  the observed statistic and p-value.

## Value

A list containing the observed studentized statistic (`z_obs`) and
permutation p-value (`p_perm`). If `return_dist = TRUE`, the permutation
distribution (`z_perm`) is also returned.

## Details

The test statistic is the HC3-studentized contrast. The permutation
p-value is a valid test of the sharp null hypothesis that treatment has
no effect for any participant; the `+1` correction keeps the Monte Carlo
p-value from being zero (Phipson & Smyth, 2010). Studentizing the
statistic makes permutation tests asymptotically robust when only an
average effect is hypothesized to be zero (DiCiccio & Romano, 2017; Wu &
Ding, 2021); for the Pretest x Treatment contrast that robustness should
be regarded as approximate.

Randomization inference must permute the unit that was randomized.
Because this function permutes individual participants, it refuses fits
that include a clustering variable. For clustered designs, use the CR2
small-sample tests reported by
[`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md).

## References

DiCiccio, C. J., & Romano, J. P. (2017). Robust permutation tests for
correlation and regression coefficients. *Journal of the American
Statistical Association, 112*(519), 1211-1220.

Phipson, B., & Smyth, G. K. (2010). Permutation p-values should never be
zero: Calculating exact p-values when permutations are randomly drawn.
*Statistical Applications in Genetics and Molecular Biology, 9*(1),
Article 39.

Wu, J., & Ding, P. (2021). Randomization tests for weak null hypotheses
in randomized experiments. *Journal of the American Statistical
Association, 116*(536), 1898-1913.
