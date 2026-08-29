# Stouffer's Z combiner (a.k.a. "Test I")

Combine one-tailed p-values that test the *same directional* hypothesis
into a single Z. This is provided to reproduce the Braver & Braver
(1988) option for SFGD. Use cautiously and document assumptions about
homogeneity; see the 1988–1990 exchanges for caveats.

## Usage

``` r
stouffer_solomon(p)
```

## Arguments

- p:

  numeric vector of one-tailed p-values (same direction)

## Value

list with z_meta and p_meta (one-tailed)

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
