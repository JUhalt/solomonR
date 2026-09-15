# Plot permutation distribution for a Solomon contrast

Draws the permutation distribution of the studentized statistic with the
observed value marked. The subtitle reports the same Monte Carlo p-value
as
[`perm_solomon()`](https://juhalt.github.io/solomonR/reference/perm_solomon.md).

## Usage

``` r
plot_perm(perm)
```

## Arguments

- perm:

  An object returned by `perm_solomon(..., return_dist = TRUE)`.

## Value

A ggplot object
