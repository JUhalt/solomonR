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

## Examples

``` r
fit <- with(solomon_example, fit_solomon_glm(y_post, treat, pretested, y_pre))
perm <- perm_solomon(fit, reps = 199, seed = 1, return_dist = TRUE)
plot_perm(perm)
```
