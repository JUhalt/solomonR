# Condom use in the Solomon study of Kvalem et al. (1996)

The number of students who used condoms at their most recent
intercourse, at the 6-month posttest, among students who had had
intercourse before the intervention, in each Solomon group of a
school-based sex-education trial (Kvalem et al., 1996, p. 42). Numbers
reported in the publication are reused with citation.

## Usage

``` r
kvalem1996
```

## Format

A data frame with 4 rows, one per Solomon group, and 5 variables:

- group:

  The Solomon group.

- pretested, treat:

  Indicators (1 = yes).

- events:

  Students who used condoms at their most recent intercourse.

- n:

  Students in the group with intercourse before the intervention who
  answered the 6-month posttest.

## Source

Kvalem, I. L., Sundet, J. M., Rivø, K. I., Eilertsen, D. E., &
Bakketeig, L. S. (1996). The effect of sex education on adolescents' use
of condoms: Applying the Solomon four-group design. *Health Education
Quarterly, 23*(1), 34–47. https://doi.org/10.1177/109019819602300103

## Details

**Design and caveats.** Whole classes were randomized: 30 classes to the
intervention, 15 of them pretested, and 94 to control, 47 of them
pretested (p. 38). The authors analyzed individuals rather than classes
because class sizes varied too much for a class-level analysis (p. 39).
These counts are therefore individual-level, and analyses of them ignore
the clustering; see
[`perm_solomon()`](https://juhalt.github.io/solomonR/reference/perm_solomon.md)
for cluster-level inference when class membership is available.

**Known result.** Fisher's exact and chi-square comparisons by
[`fisher_solomon()`](https://juhalt.github.io/solomonR/reference/fisher_solomon.md)
reproduce the reported pattern: 70% (51/73) against 51% (76/148) among
pretested students, chi-square = 6.85, and 43% (21/49) against 52%
(69/133) among unpretested students, chi-square = 1.17 (p. 42).

## See also

[`fisher_solomon()`](https://juhalt.github.io/solomonR/reference/fisher_solomon.md),
[`marginal_solomon()`](https://juhalt.github.io/solomonR/reference/marginal_solomon.md),
[elkarkri2025a](https://juhalt.github.io/solomonR/reference/elkarkri2025a.md)

## Examples

``` r
# One row per student.
students <- kvalem1996[rep(seq_len(4), kvalem1996$n), c("pretested", "treat")]
students$used <- unlist(Map(function(e, n) rep(1:0, c(e, n - e)),
                            kvalem1996$events, kvalem1996$n))
with(students, fisher_solomon(used, treat, pretested))
#> Historical categorical Solomon analysis (El Karkri et al., 2025b)
#> 
#>    Condition Events (T) n (T) Events (C) n (C) Chi-square Fisher p
#>    Pretested         51    73         76   148       6.85   0.0095
#>  Unpretested         21    49         69   133       1.17   0.3180
#>     Combined         72   122        145   281       1.88   0.1922
#> 
#> Historical rule (significant among pretested, not among unpretested, alpha = 0.05): sensitization declared
#> Caution: this rule compares significance, not effects, and is not a test of
#> the Pretest x Treatment interaction. See marginal_solomon().
```
