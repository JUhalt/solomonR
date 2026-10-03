# Solomon analysis from summary statistics

**\[stable\]** Reanalyzes a published Solomon four-group study from the
sample size, mean, and standard deviation of the posttest in each of the
four groups. With `treat` and `pretested` it also reanalyzes designs
with several treatments (see "Designs with several treatments").

## Usage

``` r
solomon_from_summary(
  n,
  mean,
  sd,
  conf_level = 0.95,
  treat = NULL,
  pretested = NULL,
  control = NULL
)
```

## Arguments

- n, mean, sd:

  Posttest sample size, mean, and standard deviation of the four groups,
  in the order above. With `treat` and `pretested`, one value for every
  group of the design, in any order.

- conf_level:

  Confidence level for intervals. Default 0.95.

- treat:

  Optional condition of each group, in the order of `n`: a character
  vector or factor, with the control named by `control`, or a 0/1
  treatment indicator. Needed for a design with several treatments. The
  treatments are reported in the order of the factor's levels, or in the
  order they first appear in a character vector.

- pretested:

  1 (or `TRUE`) for each pretested group and 0 for the others, in the
  order of `n`. Needed with `treat`.

- control:

  The control condition, when `treat` is a character vector or factor.

## Value

An object of class `solomon_summary_fit` with `contrasts` (Tests A-D,
the pretest main effect, and the simple effects: estimate, standard
error, t, degrees of freedom, p-value, confidence interval, F, and Type
III sum of squares), an `anova` table, the `cells`, the error mean
square and degrees of freedom, and the settings.

For a design with several treatments, an object of class
`solomon_summary_ngroup` with `contrasts` (for each comparison of a
treatment with the control and each of the four contrasts: estimate,
standard error, t, p-value, Holm-adjusted p-value `p.adjusted`, degrees
of freedom, and confidence interval), the omnibus tests in `anova` (Type
III sum of squares, df, mean square, F, p-value), the `cells` in the
package's group order, the `conditions` (control first), `adjust`
(`"holm"`), the error mean square and degrees of freedom, and the
settings.

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

## Designs with several treatments

A Solomon N-group design has k treatments and a control, each with and
without a pretest: 2(k + 1) groups (Edmonds & Kennedy, 2017; Steyn,
2009). Give `n`, `mean`, and `sd` for every group, in any order, and say
which group each value belongs to with `treat` (its condition) and
`pretested` (1 for a pretested group), naming the control with
`control`.

The analysis is the cell-means model of the posttest with a pooled error
variance on N - 2(k + 1) degrees of freedom. For each treatment against
the control it gives the four Solomon contrasts with t tests and
confidence intervals. They equal those of
`fit_solomon_glm(y_post, treat, pretested, control = , robust = "none")`
on the individual data, without the pretest as a covariate. The p-values
of each contrast are adjusted across the k comparisons by Holm's (1979)
procedure. The confidence intervals are not adjusted.

The omnibus F tests are those of the two-way ANOVA with Type III sums of
squares: Condition (k df), Pretest (1 df), and Pretest x Condition (k
df). Each is a Wald test of contrasts of the unweighted cell means, as
in the four-group analysis. The result has class
`solomon_summary_ngroup`. With two conditions, the result is the
four-group analysis above.

## References

Edmonds, W. A., & Kennedy, T. D. (2017). *An applied guide to research
designs: Quantitative, qualitative, and mixed methods* (2nd ed.). SAGE
Publications. https://doi.org/10.4135/9781071802779

El Karkri, M., Quesada, A., & Romero-Ariza, M. (2025a). The dual impact
of pretest sensitisation and the cognitive acceleration through science
education programme in the Solomon four-group design. *Brain Sciences,
16*(1), Article 64. https://doi.org/10.3390/brainsci16010064

Holm, S. (1979). A simple sequentially rejective multiple test
procedure. *Scandinavian Journal of Statistics, 6*(2), 65–70.
https://www.jstor.org/stable/4615733

Mai, N. N., Takahashi, Y., & Oo, M. M. (2020). Testing the effectiveness
of transfer interventions using Solomon four-group designs. *Education
Sciences, 10*(4), Article 92. https://doi.org/10.3390/educsci10040092

Steyn, R. (2009). Re-designing the Solomon four-group: Can we improve on
this exemplary model? *Design Principles and Practices: An International
Journal—Annual Review, 3*(1), 383–394.
https://doi.org/10.18848/1833-1874/CGP/v03i01/37588

Walton Braver, M. C., & Braver, S. L. (1988). Statistical treatment of
the Solomon four-group design: A meta-analytic approach. *Psychological
Bulletin, 104*(1), 150–154. https://doi.org/10.1037/0033-2909.104.1.150

## See also

[`solomon_effect_sizes()`](https://juhalt.github.io/solomonR/reference/solomon_effect_sizes.md)
for effect sizes for meta-analysis, and
[`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md)
for the analysis of the individual data.

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

# A six-group design: two treatments and a control. Posttest statistics
# computed from the data of Mai et al. (2020); see mai2020.
solomon_from_summary(
  n = c(24, 23, 27, 22, 15, 22),
  mean = c(2.929167, 3.168116, 3.112346, 3.128788, 3.152184, 3.018548),
  sd = c(0.434203, 0.369613, 0.355440, 0.383150, 0.374069, 0.354758),
  treat = c("RP", "GS", "Control", "RP", "GS", "Control"),
  pretested = c(1, 1, 1, 0, 0, 0),
  control = "Control"
)
#> Solomon analysis from summary statistics, N-group design
#> --------------------------------------------------------
#> Conditions: RP, GS; control: Control. With and without a pretest: 6 groups.
#> Pooled error variance: 0.144 on 127 df (equal variances assumed)
#> 
#> Two-way ANOVA on the posttest (Type III sums of squares)
#>   Condition            SS =    0.363  df = 2  F = 1.26  p = 0.288
#>   Pretest              SS =    0.029  df = 1  F = 0.20  p = 0.655
#>   Pretest x Condition  SS =    0.535  df = 2  F = 1.86  p = 0.161
#>   Error                SS =   18.311  df = 127
#> 
#> Contrasts with 95% confidence intervals
#> Comparison     Contrast                       Est (SE)      t   df      p  p adj.           95% CI
#> RP vs Control  ATE (avg over pretest)   -0.036 (0.078)  -0.47  127  0.642   0.642  [-0.191, 0.118]
#> GS vs Control  ATE (avg over pretest)    0.095 (0.083)   1.14  127  0.258   0.516  [-0.070, 0.260]
#> RP vs Control  Pretest x Treatment      -0.293 (0.156)  -1.88  127  0.063   0.126  [-0.603, 0.016]
#> GS vs Control  Pretest x Treatment      -0.078 (0.167)  -0.47  127  0.641   0.641  [-0.408, 0.252]
#> RP vs Control  Treatment | pretested    -0.183 (0.107)  -1.72  127  0.088   0.176  [-0.394, 0.028]
#> GS vs Control  Treatment | pretested     0.056 (0.108)   0.52  127  0.606   0.606  [-0.157, 0.269]
#> RP vs Control  Treatment | unpretested   0.110 (0.114)   0.96  127  0.337   0.590  [-0.116, 0.337]
#> GS vs Control  Treatment | unpretested   0.134 (0.127)   1.05  127  0.295   0.590  [-0.118, 0.385]
#> 
#> p adj.: adjusted by Holm's (1979) procedure within each contrast, across the 2 comparisons.
#> Confidence intervals are not adjusted.
```
