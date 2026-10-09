# Stouffer's Z combiner (Walton Braver & Braver, 1988, Test I)

**\[stable\]** Combine one-tailed p-values that test the *same
directional* hypothesis into a single Z. This is provided to reproduce
the Walton Braver & Braver (1988) meta-analytic option (Test I) for the
Solomon four-group design. Use cautiously and document assumptions about
homogeneity; see the 1988–1990 exchanges for caveats.

## Usage

``` r
stouffer_solomon(p)
```

## Arguments

- p:

  numeric vector of one-tailed p-values (same direction)

## Value

A list with:

- `z_meta`: the combined z, the sum of the z scores divided by the
  square root of their number.

- `p_meta_two_tailed`: its two-tailed p-value.

- `p_meta_one_tailed`: its one-tailed p-value, in the direction of the
  p-values supplied.

## Details

Walton Braver and Braver (1988, p. 152) convert the one-tailed p-value
of each test to z and refer the combined z to a normal table. Their
worked example reports z = 2.05 with p = .040 (p. 153), the two-tailed
value, so
[`fit_solomon_classic()`](https://juhalt.github.io/solomonR/reference/fit_solomon_classic.md)
judges Test I by `p_meta_two_tailed`. The one-tailed value is also
returned.

## References

Sawilowsky, S. S., Kelley, D. L., Blair, R. C., & Markman, B. S. (1994).
Meta-analysis and the Solomon four-group design. *The Journal of
Experimental Education, 62*(4), 361–376.
https://doi.org/10.1080/00220973.1994.9944140

Stouffer, S. A., Suchman, E. A., DeVinney, L. C., Star, S. A., &
Williams, R. M., Jr. (1949). *The American soldier: Adjustment during
army life* (Vol. 1). Princeton University Press.

Walton Braver, M. C., & Braver, S. L. (1988). Statistical treatment of
the Solomon four-group design: A meta-analytic approach. *Psychological
Bulletin, 104*(1), 150–154. https://doi.org/10.1037/0033-2909.104.1.150

## Examples

``` r
stouffer_solomon(c(0.10, 0.11))
#> $z_meta
#> [1] 1.77348
#> 
#> $p_meta_two_tailed
#> [1] 0.07614918
#> 
#> $p_meta_one_tailed
#> [1] 0.03807459
#> 

# Walton Braver and Braver's (1988, p. 153) example: the ANCOVA (p = .0993)
# and the t test (p = .2127), halved to one-tailed p-values in the
# direction of the effect, give z = 2.05 and p = .040.
stouffer_solomon(c(0.0993, 0.2127) / 2)
#> $z_meta
#> [1] 2.046673
#> 
#> $p_meta_two_tailed
#> [1] 0.04069024
#> 
#> $p_meta_one_tailed
#> [1] 0.02034512
#> 
```
