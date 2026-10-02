# Solomon's (1949) spelling experiment

Group sizes, means, and standard errors from the first published
experiment with the Solomon design: spelling lessons in a fifth-grade
and a sixth-grade class, analyzed with the three-group design (Solomon,
1949, Tables II and III, pp. 144–145). Numbers reported in the
publication are reused with citation.

## Usage

``` r
solomon1949
```

## Format

A data frame with 6 rows, one per group and grade, and 11 variables:

- grade:

  5 or 6.

- group:

  Experimental, Control I, or Control II.

- pretested, treat:

  Indicators (1 = yes); `treat` is the spelling lesson.

- n:

  Group size.

- pre_mean, pre_se:

  Pretest mean and its standard error (pretested groups only).

- mean, se:

  Posttest mean and its standard error.

- change, change_se:

  Mean improvement, as printed, and its standard error (missing for
  Control II, printed as "?").

The standard errors are the values printed with a plus-or-minus sign,
which Solomon's notation (sigma_m, Table I) identifies as standard
errors of the means.

## Source

Solomon, R. L. (1949). An extension of control group design.
*Psychological Bulletin, 46*(2), 137–150.
https://doi.org/10.1037/h0062958

## Details

**Design.** In each class, pupils were assigned to three groups that
were "roughly equated in spelling ability by means of teachers'
judgments" (p. 144), not randomized.

- **Experimental group.** Pretested on a list of words, given a standard
  spelling lesson, and posttested on the same words.

- **Control Group I.** Pretested and posttested, without the lesson.

- **Control Group II.** Given the lesson and the posttest, without the
  pretest.

There was no fourth group; Solomon introduced it later in the article
for field studies (p. 147).

**Known result.**
[`fit_solomon_1949()`](https://juhalt.github.io/solomonR/reference/fit_solomon_1949.md)
reproduces the published values from these means: an inferred pretest of
3.0 and 5.7, improvements for Control Group II of 8.2 and 8.7, and
interactions I = -2.2 (grade 5) and -3.1 (grade 6). The pretest reduced
the effect of the lesson, so the usual two-group design would have
underrated it (p. 145).

**Caveats.** The groups are small (8 to 10 pupils) and were not
randomized. Solomon printed the standard error of Control Group II's
improvement as "?", since it depends on the inferred pretest.

## See also

[`fit_solomon_1949()`](https://juhalt.github.io/solomonR/reference/fit_solomon_1949.md),
[elkarkri2025a](https://juhalt.github.io/solomonR/reference/elkarkri2025a.md),
[mai2020](https://juhalt.github.io/solomonR/reference/mai2020.md)

## Examples

``` r
solomon1949
#>   grade        group pretested treat  n pre_mean pre_se mean  se change
#> 1     5 Experimental         1     1 10      3.2    0.8  9.9 1.6    6.7
#> 2     5    Control I         1     0 10      2.8    0.7  3.5 0.8    0.7
#> 3     5   Control II         0     1 10       NA     NA 11.2 1.2    8.2
#> 4     6 Experimental         1     1  8      5.4    1.4 11.8 1.0    6.4
#> 5     6    Control I         1     0  9      6.0    1.5  6.8 1.4    0.8
#> 6     6   Control II         0     1  8       NA     NA 14.4 0.3    8.7
#>   change_se
#> 1       0.9
#> 2       0.5
#> 3        NA
#> 4       1.2
#> 5       0.6
#> 6        NA
g6 <- solomon1949[solomon1949$grade == 6, ]
fit_solomon_1949(post_mean = g6$mean, pre_mean = g6$pre_mean[1:2], n = g6$n)
#> Solomon (1949) improvement-score analysis (historical)
#> ------------------------------------------------------
#> Design: three-group (Solomon, 1949, Table I, p. 142)
#> 
#>  Group        n Pretest Training Pre mean        Post mean Improvement
#>  Experimental 8 yes     yes      5.40            11.80     6.40       
#>  Control I    9 yes     no       6.00            6.80      0.80       
#>  Control II   8 no      yes      5.70 (inferred) 14.40     8.70       
#> 
#> Inferred pretest i = 5.70 (average of the pretested groups' means)
#> Interaction I = d1 - (d2 + d3) = -3.10
#> 
#> Solomon gave no standard error or test for I. Campbell and Stanley
#> (1963/1966, p. 25) judged his gain-score suggestions unacceptable; see
#> fit_solomon_classic() for the tests that followed and fit_solomon_glm()
#> for the recommended analysis.
```
