# Stouffer's Z combiner (Braver & Braver, 1988, Test I)

Combine one-tailed p-values that test the *same directional* hypothesis
into a single Z. This is provided to reproduce the Braver & Braver
(1988) meta-analytic option (Test I) for the Solomon four-group design.
Use cautiously and document assumptions about homogeneity; see the
1988–1990 exchanges for caveats.

## Usage

``` r
stouffer_solomon(p)
```

## Arguments

- p:

  numeric vector of one-tailed p-values (same direction)

## Value

list with z_meta and p_meta (one-tailed)

## References

Stouffer, S. A., Suchman, E. A., DeVinney, L. C., Star, S. A., &
Williams, R. M., Jr. (1949). *The American soldier: Adjustment during
army life* (Vol. 1). Princeton University Press.

Braver, M. W., & Braver, S. L. (1988). Statistical treatment of the
Solomon four-group design: A meta-analytic approach. *Psychological
Bulletin, 104*(1), 150-154.

Sawilowsky, S. S., Kelley, D. L., Blair, R. C., & Markman, B. S. (1994).
Meta-analysis and the Solomon four-group design. *The Journal of
Experimental Education, 62*(4), 361-376.

## Examples

``` r
stouffer_solomon(c(0.10, 0.11))
#> $z_meta
#> [1] 1.77348
#> 
#> $p_meta_one_tailed
#> [1] 0.03807459
#> 
```
