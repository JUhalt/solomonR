# Solomon analysis from summary statistics

**\[stable\]** Reanalyzes a published Solomon four-group study from the
sample size, mean, and standard deviation of the posttest in each of the
four groups.

## Usage

``` r
solomon_from_summary(n, mean, sd, conf_level = 0.95)
```

## Arguments

- n, mean, sd:

  Posttest sample size, mean, and standard deviation of the four groups,
  in the order above.

- conf_level:

  Confidence level for intervals. Default 0.95.

## Value

An object of class `solomon_summary_fit` with `contrasts` (Tests A-D,
the pretest main effect, and the simple effects: estimate, standard
error, t, degrees of freedom, p-value, confidence interval, F, and Type
III sum of squares), an `anova` table, the `cells`, the error mean
square and degrees of freedom, and the settings.

## Details

The four groups are given in the order used throughout the package:
pretested treated (O2), pretested control (O4), unpretested treated
(O5), and unpretested control (O6). The analysis is the 2 x 2
between-groups model on the posttest with a pooled error variance, so it
reproduces the historical Tests A-D and the simple treatment effects
(Tests B and C) of
[`fit_solomon_classic()`](https://juhalt.github.io/solomonR/reference/fit_solomon_classic.md).
Main effects are contrasts of unweighted cell means, which match the
Type III sums of squares that statistical packages report for unbalanced
cells.

Summary statistics limit the analysis. The pooled error variance assumes
equal variances in the four groups; heteroskedasticity-consistent
standard errors, covariate adjustment, and the analyses of the pretested
groups (Tests E-G) need the individual data. Rounded published
statistics reproduce published tests only to within rounding: for El
Karkri et al. (2025a), the interaction F is 11.46 against the published
11.48.

## References

El Karkri, M., Quesada, A., & Romero-Ariza, M. (2025a). The dual impact
of pretest sensitisation and the cognitive acceleration through science
education programme in the Solomon four-group design. *Brain Sciences,
16*(1), Article 64. https://doi.org/10.3390/brainsci16010064

Walton Braver, M. C., & Braver, S. L. (1988). Statistical treatment of
the Solomon four-group design: A meta-analytic approach. *Psychological
Bulletin, 104*(1), 150–154. https://doi.org/10.1037/0033-2909.104.1.150

## See also

[`solomon_effect_sizes()`](https://juhalt.github.io/solomonR/reference/solomon_effect_sizes.md)
for effect sizes for meta-analysis.

## Examples

``` r
# El Karkri et al. (2025a), Table 8
solomon_from_summary(
  n = c(9, 25, 17, 37),
  mean = c(10.94, 7.80, 8.94, 9.35),
  sd = c(2.26, 2.29, 1.98, 2.11)
)
#> Solomon analysis from summary statistics
#> ----------------------------------------
#> Pooled error variance: 4.640 on 84 df (equal variances assumed)
#> 
#> Two-way ANOVA on the posttest (Type III sums of squares)
#>   Treatment            SS =   31.452  df = 1  F = 6.78  p = 0.011
#>   Pretest              SS =    0.855  df = 1  F = 0.18  p = 0.669
#>   Treatment x Pretest  SS =   53.184  df = 1  F = 11.46  p = 0.001
#>   Error                SS =  389.721  df = 84
#> 
#> Contrasts with 95% confidence intervals
#>   Test A: Pretest x Treatment               3.550 [1.465, 5.635], t(84) = 3.39, p = 0.001
#>   Test B: Treatment | pretested             3.140 [1.475, 4.805], t(84) = 3.75, p < .001
#>   Test C: Treatment | unpretested          -0.410 [-1.665, 0.845], t(84) = -0.65, p = 0.518
#>   Test D: ATE (avg over pretest)            1.365 [0.322, 2.408], t(84) = 2.60, p = 0.011
#>   Pretest main effect                       0.225 [-0.818, 1.268], t(84) = 0.43, p = 0.669
```
