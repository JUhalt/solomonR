# Convert p-value to Z (one-tailed) for Stouffer's method

**\[stable\]**

## Usage

``` r
p_to_z(p)
```

## Arguments

- p:

  numeric vector of p-values assumed one-tailed and aligned in the same
  direction

## Value

numeric Z-scores

## References

Stouffer, S. A., Suchman, E. A., DeVinney, L. C., Star, S. A., &
Williams, R. M., Jr. (1949). *The American soldier: Adjustment during
army life* (Vol. 1). Princeton University Press.

## Examples

``` r
p_to_z(c(0.05, 0.025, 0.5))
#> [1] 1.644854 1.959964 0.000000
```
