# Pretest and posttest statistics from El Karkri et al. (2025a)

The posttest sample size, mean, and standard deviation of the four
Solomon groups in a classroom study of the Cognitive Acceleration
through Science Education programme (El Karkri et al., 2025a, Table 8,
p. 11), and the pretest mean and standard deviation of the two pretested
groups (Table 7, p. 10). Numbers reported in the publication are reused
with citation.

## Usage

``` r
elkarkri2025a
```

## Format

A data frame with 4 rows, one per Solomon group, and 8 variables:

- group:

  The Solomon group.

- pretested, treat:

  Indicators (1 = yes).

- n, mean, sd:

  Posttest sample size, mean, and standard deviation.

- pre_mean, pre_sd:

  Pretest mean and standard deviation (pretested groups only; the
  pretest sample sizes equal the posttest ones).

## Source

El Karkri, M., Quesada, A., & Romero-Ariza, M. (2025a). The dual impact
of pretest sensitisation and the cognitive acceleration through science
education programme in the Solomon four-group design. *Brain Sciences,
16*(1), Article 64. https://doi.org/10.3390/brainsci16010064

## Details

**Design and caveats.** The authors describe the study as
quasi-experimental: each Solomon group was one intact class, so class
and condition are confounded, and differences between the groups can
reflect the classes as well as the treatment and the pretest.
[`validate_solomon()`](https://juhalt.github.io/solomonR/reference/validate_solomon.md)
flags this design when class membership is supplied. The pretested
classes already differed at pretest (9.61 against 7.86), which
[`baseline_solomon()`](https://juhalt.github.io/solomonR/reference/baseline_solomon.md)
reports; the unpretested classes have no pretest. The authors reported a
significant Pretest x Treatment interaction.

**Known result.**
[`solomon_from_summary()`](https://juhalt.github.io/solomonR/reference/solomon_from_summary.md)
reproduces the published two-way ANOVA from these numbers within
rounding: interaction F(1, 84) = 11.46 against the published 11.482,
treatment 6.78 against 6.794, and pretest 0.18 against 0.186 (pp.
11-12).

## See also

[`solomon_from_summary()`](https://juhalt.github.io/solomonR/reference/solomon_from_summary.md),
[`baseline_solomon()`](https://juhalt.github.io/solomonR/reference/baseline_solomon.md),
[kvalem1996](https://juhalt.github.io/solomonR/reference/kvalem1996.md)

## Examples

``` r
with(elkarkri2025a, solomon_from_summary(n, mean, sd))
#> Solomon analysis from summary statistics
#> ----------------------------------------
#> Pooled error variance: 4.640 on 84 df (equal variances assumed)
#> 
#> Two-way ANOVA on the posttest (Type III sums of squares)
#>   Treatment            SS =   31.452  df = 1  F = 6.78  p = 0.011
#>   Pretest              SS =    0.855  df = 1  F = 0.18  p = 0.669
#>   Pretest x Treatment  SS =   53.184  df = 1  F = 11.46  p = 0.001
#>   Error                SS =  389.721  df = 84
#> 
#> Contrasts with 95% confidence intervals
#>   Test A: Pretest x Treatment               3.550 [1.465, 5.635], t(84) = 3.39, p = 0.001
#>   Test B: Treatment | pretested             3.140 [1.475, 4.805], t(84) = 3.75, p < .001
#>   Test C: Treatment | unpretested          -0.410 [-1.665, 0.845], t(84) = -0.65, p = 0.518
#>   Test D: ATE (avg over pretest)            1.365 [0.322, 2.408], t(84) = 2.60, p = 0.011
#>   Pretest main effect                       0.225 [-0.818, 1.268], t(84) = 0.43, p = 0.669
pre <- elkarkri2025a[elkarkri2025a$pretested == 1, ]
baseline_solomon(n = pre$n, mean = pre$pre_mean, sd = pre$pre_sd)
#> Baseline comparison of the pretested arms
#> -----------------------------------------
#>   Pretested, treatment   n = 9, M = 9.61, SD = 2.67
#>   Pretested, control     n = 25, M = 7.86, SD = 2.72
#> 
#> Difference: 1.75, 95% CI [-0.39, 3.89], t(32) = 1.66, p = .106
#> Hedges's g: 0.63, 95% CI [-0.14, 1.42] (noncentral t)
#> 
#> The unpretested arms have no pretest, so their baseline cannot be checked
#> or adjusted for with the study's own data.
```
