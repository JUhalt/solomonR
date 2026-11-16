# Plot permutation distribution for a Solomon contrast

**\[stable\]** Draws the permutation distribution of the test statistic
with the observed value marked. The subtitle reports the same p-value as
[`perm_solomon()`](https://juhalt.github.io/solomonR/reference/perm_solomon.md),
and says whether it came from every possible allocation or from sampled
permutations.

## Usage

``` r
plot_perm(perm)
```

## Arguments

- perm:

  An object returned by `perm_solomon(..., return_dist = TRUE)`.

## Value

A ggplot object

## Details

The p-value is the object's `p_perm`. An object without it is recounted
from `z_perm` by the rule of
[`perm_solomon()`](https://juhalt.github.io/solomonR/reference/perm_solomon.md),
with the tie window and the unit's floor that the object records
(`tie_tolerance` and `tie_scale`), over the entries that are not `NA`
(permutations left out): b / M when every allocation was used, and
(b + 1) / (m + 1) otherwise. An object saved by solomonR 0.8.2 or
earlier records neither and is recounted as `plot_perm()` recounted it
then: with a window of a relative 1e-10 and, for the difference
statistic, the largest permuted difference as the unit, with no floor.

Infinite statistics (an emptied or constant arm, or separation; at the
cluster level, constant arms) are counted in the p-value but cannot be
drawn. A caption says how many there are, and whether the observed
statistic is infinite, in which case no line marks it.

## Examples

``` r
fit <- with(solomon_example, fit_solomon_glm(y_post, treat, pretested, y_pre))
perm <- perm_solomon(fit, reps = 199, seed = 1, return_dist = TRUE)
plot_perm(perm)
```
