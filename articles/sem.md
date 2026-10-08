# Structural Equation Models for Solomon Designs

When a Solomon study’s outcome is measured by several items, a
structural equation model can compare latent means instead of scale
scores. This article shows the package’s two SEM functions:

- **The observed-variable model** of the four groups, and its ANCOVA
  version for the pretested groups.
- **The latent-variable model**, and the measurement invariance its
  comparisons require.

The SEM functions use lavaan (Rosseel, 2012), and their tests and
intervals are lavaan’s large-sample Wald results. Latent-variable
analysis of a Solomon design has a published precedent, Dukes et
al. (1995), which the last section compares with the package’s approach.
The four-group models of the SEM functions, and the observed-variable
model, are solomonR extensions built on established SEM methods.

## Observed outcomes

[`fit_solomon_sem()`](https://juhalt.github.io/solomonR/reference/fit_solomon_sem.md)
estimates the posttest mean of each Solomon group in a four-group
mean-structure model and reports the four Solomon contrasts:

``` r

sem <- with(solomon_example, fit_solomon_sem(y_post, treat, pretested))
sem
#> Solomon SEM (multi-group mean structure)
#> Fit: saturated four-group mean structure (df = 0)
#> Global CFI/RMSEA/SRMR are not diagnostic for this model.
#> 
#> Key contrasts  Est (SE)            z      p           95% CI
#> ATE            2.683 (1.726)    1.55  0.120  [-0.699, 6.066]
#> Sens           -1.900 (3.451)  -0.55  0.582  [-8.665, 4.865]
#> Pre_Eff        1.733 (2.696)    0.64  0.520  [-3.551, 7.018]
#> Unpre_Eff      3.633 (2.155)    1.69  0.092  [-0.590, 7.857]
```

- **No global fit test.** The four-group mean-structure model is
  saturated, so its global fit indices are not diagnostic, and the
  function does not report them.
- **Unadjusted contrasts.** The contrasts are unadjusted posttest mean
  differences, which is the SEM counterpart of the two-by-two analysis.

With `ancova = TRUE`, the model is restricted to the pretested groups
and adjusts for the pretest. That is the analysis of covariance Campbell
and Stanley (1963/1966, p. 25) and Huck and Sandler (1973) describe for
those groups:

``` r

with(solomon_example, fit_solomon_sem(y_post, treat, pretested, y_pre, ancova = TRUE))
#> Solomon SEM (ANCOVA in pretested groups)
#> Fit: CFI=0.929, RMSEA=0.247, SRMR=0.125; df=1
#> 
#> Key contrasts  Est (SE)          z      p           95% CI
#> Pre_Eff        1.694 (2.142)  0.79  0.429  [-2.504, 5.893]
```

The unpretested groups have no pretest by design, so the adjusted model
uses only the pretested groups. It does not treat the absent pretests as
missing data.

## Latent outcomes

[`fit_solomon_sem_latent()`](https://juhalt.github.io/solomonR/reference/fit_solomon_sem_latent.md)
defines a latent posttest from several indicators and compares its means
across the four groups.

- **Identification.** The latent mean of the unpretested control group
  is fixed at 0; the contrasts are differences between latent means, so
  they do not depend on that choice. With the default `std_lv = TRUE`,
  the latent variance is fixed at 1 in the pretested treated group, so
  the contrasts are in standard deviations of the latent posttest there.
- **Optional ANCOVA.** With `ancova = TRUE` and pretest indicators, a
  latent ANCOVA in the pretested groups is added. With `std_lv = TRUE`,
  its adjusted effect is in residual standard deviations of the latent
  posttest given the latent pretest, a different unit from that of the
  four-group contrasts.

The example data below are simulated for this article. Four indicators
measure one construct in every group, and the treatment raises the
latent mean by 0.4 standard deviations:

``` r

set.seed(2049)
n <- 80
g <- rep(1:4, each = n)
treat <- c(1, 0, 1, 0)[g]
pretested <- c(1, 1, 0, 0)[g]
f <- rnorm(4 * n, 0.4 * treat)
items <- data.frame(
  y1 = 0.8 * f + rnorm(4 * n, 0, 0.6),
  y2 = 0.9 * f + rnorm(4 * n, 0, 0.6),
  y3 = 0.7 * f + rnorm(4 * n, 0, 0.6),
  y4 = 0.8 * f + rnorm(4 * n, 0, 0.6)
)
latent <- fit_solomon_sem_latent(items, c("y1", "y2", "y3", "y4"), treat, pretested)
latent
#> Solomon SEM (latent)
#> POST measurement invariance: scalar
#> Latent mean reference: U0 (unpretested control) fixed at 0
#> 
#> POST model fit (4 groups): CFI=1.000, RMSEA=0.000, SRMR=0.056; df=26
#> Invariance check: scalar invariance supported by both criteria
#> 
#> Key contrasts (latent POST)  Est (SE)          z      p           95% CI
#> ATE                          0.399 (0.130)  3.08  0.002   [0.145, 0.653]
#> Sens                         0.060 (0.252)  0.24  0.812  [-0.434, 0.554]
#> Pre_Eff                      0.429 (0.163)  2.63  0.008   [0.110, 0.748]
#> Unpre_Eff                    0.369 (0.197)  1.87  0.061  [-0.017, 0.755]
```

## Measurement invariance

Latent means can be compared only if the indicators measure the
construct in the same way in every group. The loadings and the
intercepts must be equal across groups, which is scalar invariance
(Meredith, 1993).

- **Why this matters for Solomon designs.** Taking the pretest may
  change how participants answer some items at the posttest, for example
  by making them familiar with the items.
- **Where the damage falls.** A shift shared by both pretested groups is
  absorbed into their latent means and looks like an effect of
  pretesting. It cancels from the sensitization contrast, which compares
  treatment effects within each pretest condition. A shift in one group
  only, such as the treated pretested participants reading an item
  differently, biases the sensitization contrast itself. The package’s
  simulation study confirmed both patterns, with the one-group shift
  placed in the unpretested control group (see “Latent Contrasts:
  Validating the Invariance Check”).

[`invariance_solomon()`](https://juhalt.github.io/solomonR/reference/invariance_solomon.md)
fits the configural, metric, and scalar models across the four groups,
in the sequence Vandenberg and Lance (2000, pp. 56–57) recommend:

``` r

inv <- invariance_solomon(items, c("y1", "y2", "y3", "y4"), treat, pretested)
inv
#> Measurement invariance across the four Solomon groups
#> Indicators: y1, y2, y3, y4; estimator: MLR; group sizes (P1, P0, U1, U0): 80, 80, 80, 80
#> 
#>       model chisq df cfi rmsea  srmr
#>  configural  1.93  8   1     0 0.007
#>      metric  8.51 17   1     0 0.041
#>      scalar 23.33 26   1     0 0.056
#> 
#>             comparison chisq_diff df_diff p.value delta_cfi delta_rmsea
#>  metric vs. configural       6.71       9   0.667         0           0
#>      scalar vs. metric      14.25       9   0.114         0           0
#>  delta_srmr noninvariant_chisq noninvariant_chen
#>       0.034              FALSE             FALSE
#>       0.015              FALSE             FALSE
#> 
#> Chi-square difference test (scaled; Satorra & Bentler, 2001) at alpha = 0.05 (Vandenberg & Lance, 2000, p. 46): scalar supported.
#> Change in fit (Chen, 2007, pp. 501-502; total N > 300 and equal group sizes): scalar supported.
```

Each step reports two published criteria:

- **The chi-square difference test.** Scaled for the robust MLR
  estimator (Satorra & Bentler, 2001). Vandenberg and Lance (2000,
  p. 46) recommend it as the primary criterion.
- **The change in fit.** The drop in CFI and the rises in RMSEA and
  SRMR, against Chen’s (2007, pp. 501–502) cutoffs. Cheung and Rensvold
  (2002, pp. 234–235) favor such changes over the chi-square difference
  because the latter depends on sample size.

Here both criteria support scalar invariance, so the latent contrasts
above can be interpreted.

### When an intercept is not invariant

Suppose the pretest makes item `y4` read higher in the two pretested
groups at the same latent level:

``` r

shifted <- items
shifted$y4[pretested == 1] <- shifted$y4[pretested == 1] + 0.6
invariance_solomon(shifted, c("y1", "y2", "y3", "y4"), treat, pretested)
#> Measurement invariance across the four Solomon groups
#> Indicators: y1, y2, y3, y4; estimator: MLR; group sizes (P1, P0, U1, U0): 80, 80, 80, 80
#> 
#>       model chisq df   cfi rmsea  srmr
#>  configural  1.93  8 1.000  0.00 0.007
#>      metric  8.51 17 1.000  0.00 0.041
#>      scalar 86.32 26 0.905  0.17 0.100
#> 
#>             comparison chisq_diff df_diff  p.value delta_cfi delta_rmsea
#>  metric vs. configural       6.71       9 6.67e-01     0.000        0.00
#>      scalar vs. metric      78.16       9 3.74e-13    -0.095        0.17
#>  delta_srmr noninvariant_chisq noninvariant_chen
#>       0.034              FALSE             FALSE
#>       0.059               TRUE              TRUE
#> 
#> Chi-square difference test (scaled; Satorra & Bentler, 2001) at alpha = 0.05 (Vandenberg & Lance, 2000, p. 46): metric supported.
#> Change in fit (Chen, 2007, pp. 501-502; total N > 300 and equal group sizes): metric supported.
```

The scalar step fails. If theory identifies the noninvariant parameter,
freeing it gives a partial-invariance model, in which the latent means
are compared through the invariant indicators (Byrne et al., 1989).
Byrne et al. caution against relaxing constraints that do not make sense
substantively (p. 465). Vandenberg and Lance (2000, p. 38) set three
conditions:

- the freed parameters involve only a minority of the indicators;
- the choice rests on as strong a theoretical basis as possible;
- cross-validation evidence supports it.

The package enforces the first condition and never chooses the
parameters itself:

``` r

fit_partial <- fit_solomon_sem_latent(shifted, c("y1", "y2", "y3", "y4"), treat, pretested,
                                      partial_post = "y4 ~ 1")
fit_partial$effects_post[, c("contrast", "estimate", "conf.low", "conf.high")]
#>    contrast estimate conf.low conf.high
#> 1       ATE    0.422    0.160     0.684
#> 2      Sens   -0.039   -0.540     0.463
#> 3   Pre_Eff    0.403    0.072     0.734
#> 4 Unpre_Eff    0.441    0.050     0.833
```

## What the package does when invariance fails

[`fit_solomon_sem_latent()`](https://juhalt.github.io/solomonR/reference/fit_solomon_sem_latent.md)
runs
[`invariance_solomon()`](https://juhalt.github.io/solomonR/reference/invariance_solomon.md)
on the POST indicators before it fits the latent means. When a criterion
does not support scalar invariance (or partial scalar invariance, with
`partial_post`), it warns. It does not refuse the contrasts:

``` r

fit_shifted <- withCallingHandlers(
  fit_solomon_sem_latent(shifted, c("y1", "y2", "y3", "y4"), treat, pretested),
  solomonR_invariance_warning = function(w) {
    message("Warning: ", conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)
#> Warning: Latent mean contrasts assume scalar invariance of the POST indicators, which the invariance check did not support: the chi-square difference test supports metric invariance only; the change in fit (Chen, 2007) supports metric invariance only. In Solomon-sized groups these criteria can also reject invariance that holds; see 'Invariance check' in ?fit_solomon_sem_latent and the fit's `invariance` element.
fit_shifted$invariance_status
#> [1] "scalar invariance not supported: the chi-square difference test supports metric invariance only; the change in fit (Chen, 2007) supports metric invariance only"
fit_shifted$effects_post[, c("contrast", "estimate", "conf.low", "conf.high")]
#>    contrast estimate conf.low conf.high
#> 1       ATE    0.410    0.151     0.669
#> 2      Sens    0.033   -0.463     0.530
#> 3   Pre_Eff    0.426    0.100     0.753
#> 4 Unpre_Eff    0.393    0.004     0.782
```

The package’s simulation study settled that choice. Its decision rules
were posted on issue \#55 before any run (see “Latent Contrasts:
Validating the Invariance Check”). A criterion could govern an automatic
refusal only if it falsely rejected invariance at a rate of .060 or less
in Solomon-sized groups, and none did. The warning is therefore a prompt
to examine the invariance results, not a verdict.

The shift in this example is common to both pretested groups, so the
sensitization contrast is close to the one from the partial model above;
the shift is absorbed into the effect of pretesting. A shift in one
group only would bias the sensitization contrast, and there the partial
model matters. `check_invariance = FALSE` skips the check when
invariance was established elsewhere.

## Reporting latent models

[`report_solomon()`](https://juhalt.github.io/solomonR/reference/report_solomon.md)
drafts the methods and results of a latent fit. It gives what the
reporting standards for structural equation models ask for (Appelbaum et
al., 2018):

- **The estimation.** The software and its version, the estimator, and
  how missing values were handled.
- **The model.** The constraints across groups, any freed parameters,
  and how the latent scale was identified.
- **The fit.** The chi-square test, CFI, RMSEA, and SRMR.

Its statement about invariance follows the fit’s own check, so for the
partial model above it names the freed intercept:

``` r

report_solomon(fit_partial)$method
#> [1] "Latent posttest means were compared across the four Solomon groups in a multiple-group structural equation model in which four posttest indicators (y1, y2, y3, and y4) measured one latent posttest, fitted with lavaan (Version 0.7-2; Rosseel, 2012) by maximum likelihood with Huber-White robust standard errors and a scaled test statistic asymptotically equal to the Yuan-Bentler statistic (MLR; Yuan & Bentler, 2000), using full-information maximum likelihood for missing indicator values. The loadings and intercepts of the indicators were constrained to be equal across the four groups, the scalar measurement invariance that latent mean comparisons require (Meredith, 1993; Vandenberg & Lance, 2000), except the intercept of y4, which was estimated separately in each group (partial invariance; Byrne et al., 1989). For identification, the latent posttest mean of the unpretested control group was fixed at 0, and the latent variance at 1 in the pretested treated group; the contrasts are therefore in standard deviations of the latent posttest in that group. Measurement invariance was tested by comparing configural, metric, and scalar models in sequence, with the same parameters freed (Vandenberg & Lance, 2000); the scaled chi-square difference test (Satorra & Bentler, 2001) and the change in fit (Chen, 2007) both supported partial scalar invariance. The contrasts between latent means were tested with Wald z tests."
```

For
[`invariance_solomon()`](https://juhalt.github.io/solomonR/reference/invariance_solomon.md)
results, it reports the fit of each model, the comparisons between them,
and the decision under each criterion, the minimal information Putnick
and Bornstein (2016) propose for invariance tests:

``` r

report_solomon(inv)$results
#> [1] "The configural model gave a scaled χ²(8) = 1.98, p = .982 (unscaled χ²(8) = 1.93), CFI = 1.000, RMSEA = .000, and SRMR = .007; the metric model gave a scaled χ²(17) = 8.71, p = .949 (unscaled χ²(17) = 8.51), CFI = 1.000, RMSEA = .000, and SRMR = .041; and the scalar model gave a scaled χ²(26) = 23.35, p = .613 (unscaled χ²(26) = 23.33), CFI = 1.000, RMSEA = .000, and SRMR = .056."
#> [2] "Constraining the loadings to be equal (metric against configural) gave a scaled Δχ²(9) = 6.71, p = .667, ΔCFI = .000, ΔRMSEA = .000, and ΔSRMR = .034; constraining the intercepts as well (scalar against metric) gave a scaled Δχ²(9) = 14.25, p = .114, ΔCFI = .000, ΔRMSEA = .000, and ΔSRMR = .015."                                                                                      
#> [3] "The scaled chi-square difference test and the change in fit both supported scalar invariance."
```

## Dukes et al. (1995) and solomonR

Dukes et al. (1995) evaluated a school drug-prevention program with a
Solomon design and four latent outcomes measured by self-report items
(pp. 412, 414). They compared latent means in multisample structural
equation models.

- **Three two-group models.** They fitted three separate two-group
  models rather than one model of the four groups (p. 420): the
  pretested control’s pretest against the unpretested control’s
  posttest, for maturation, and the program against the control within
  the pretested and within the unpretested groups.
- **Loadings tested, intercepts imposed.** In each comparison, equal
  factor loadings were tested by chi-square difference tests. The latent
  means were then compared with the indicators’ means held equal, a
  constraint that was imposed, not tested (pp. 420–426).
- **No test of the interaction.** Pretest effects were judged by
  comparing which program effects were significant in the two program
  comparisons (p. 426), the comparison Gelman and Stern (2006) caution
  against. The abstract calls the pretest reactive for one outcome,
  resistance to peer pressure (p. 409). The results are more consistent
  with a pretest main effect, seen in one of that outcome’s two
  indicators, than with sensitization. The paper reports no test of the
  interaction, and solomonR’s approximate re-analyses find none, either
  in the difference between the two latent program effects (Table 4,
  p. 425) or in two-by-two analyses of the two indicators’ classroom
  posttest means (Table 2, pp. 415–416).
- **Classroom means.** They analyzed 440 classroom means (pp. 413, 417).
  Such means are more stable and precise (p. 412), and the program was
  delivered to intact classes, although the school, the unit assigned to
  the program, was the appropriate unit (p. 432, note 3).

solomonR shares the core of their method: a multiple-group model with a
mean structure (pp. 414, 420), equal loadings tested before latent means
are compared, and latent means compared with the indicators’ means held
equal (pp. 420–426). For the pretested groups, Dukes et al. also held
the path from each pretest factor to its posttest factor equal across
the two groups (p. 423). That is the common slope of the latent analysis
of covariance in `fit_solomon_sem_latent(ancova = TRUE)`; they tested
the equality, and the package imposes it.

solomonR fits the four groups in one model. Its latent mean structure is
identified by fixing the unpretested control’s mean at 0
([\#16](https://github.com/JUhalt/solomonR/issues/16)), so the
sensitization contrast, the average treatment effect, and both simple
effects have their own Wald tests on one scale.
[`invariance_solomon()`](https://juhalt.github.io/solomonR/reference/invariance_solomon.md)
tests configural, metric, and scalar invariance across all four groups
([\#55](https://github.com/JUhalt/solomonR/issues/55)), including the
posttests of pretested against unpretested groups, which Dukes et
al. did not compare. The SEM functions make no adjustment for clusters.
When clusters such as schools were assigned to conditions, analyze at
the level of assignment, or analyze observed scores with
`fit_solomon_glm(robust = "CR2", cluster = )`.

Four elements of their analyses are planned for the SEM functions
([\#117](https://github.com/JUhalt/solomonR/issues/117)):

- a latent maturation contrast, comparing the pretested control’s
  pretest with the unpretested control’s posttest (pp. 420–422);
- a pretest main-effect contrast, the effect they describe (p. 426);
- the latent pretest difference between the pretested groups, as a
  baseline check in the latent analysis of covariance (p. 423; Table 4,
  p. 425);
- a standardized latent contrast, using the pooled latent standard
  deviation (p. 426).

All works cited in solomonR are listed, with notes on how the package
uses them, on the
[References](https://juhalt.github.io/solomonR/articles/references.md)
page.

## References

Appelbaum, M., Cooper, H., Kline, R. B., Mayo-Wilson, E., Nezu, A. M., &
Rao, S. M. (2018). Journal article reporting standards for quantitative
research in psychology: The APA Publications and Communications Board
task force report. *American Psychologist, 73*(1), 3–25.
<https://doi.org/10.1037/amp0000191>

Byrne, B. M., Shavelson, R. J., & Muthén, B. (1989). Testing for the
equivalence of factor covariance and mean structures: The issue of
partial measurement invariance. *Psychological Bulletin, 105*(3),
456–466. <https://doi.org/10.1037/0033-2909.105.3.456>

Campbell, D. T., & Stanley, J. C. (1966). *Experimental and
quasi-experimental designs for research*. Rand McNally. (Original work
published 1963)

Chen, F. F. (2007). Sensitivity of goodness of fit indexes to lack of
measurement invariance. *Structural Equation Modeling: A
Multidisciplinary Journal, 14*(3), 464–504.
<https://doi.org/10.1080/10705510701301834>

Cheung, G. W., & Rensvold, R. B. (2002). Evaluating goodness-of-fit
indexes for testing measurement invariance. *Structural Equation
Modeling: A Multidisciplinary Journal, 9*(2), 233–255.
<https://doi.org/10.1207/S15328007SEM0902_5>

Dukes, R. L., Ullman, J. B., & Stein, J. A. (1995). An evaluation of
D.A.R.E. (Drug Abuse Resistance Education), using a Solomon four-group
design with latent variables. *Evaluation Review, 19*(4), 409–435.
<https://doi.org/10.1177/0193841X9501900404>

Gelman, A., & Stern, H. (2006). The difference between “significant” and
“not significant” is not itself statistically significant. *The American
Statistician, 60*(4), 328–331.
<https://doi.org/10.1198/000313006X152649>

Huck, S. W., & Sandler, H. M. (1973). A note on the Solomon 4-group
design: Appropriate statistical analyses. *The Journal of Experimental
Education, 42*(2), 54–55.
<https://doi.org/10.1080/00220973.1973.11011460>

Meredith, W. (1993). Measurement invariance, factor analysis and
factorial invariance. *Psychometrika, 58*(4), 525–543.
<https://doi.org/10.1007/BF02294825>

Putnick, D. L., & Bornstein, M. H. (2016). Measurement invariance
conventions and reporting: The state of the art and future directions
for psychological research. *Developmental Review, 41*, 71–90.
<https://doi.org/10.1016/j.dr.2016.06.004>

Rosseel, Y. (2012). lavaan: An R package for structural equation
modeling. *Journal of Statistical Software, 48*(2), 1–36.
<https://doi.org/10.18637/jss.v048.i02>

Satorra, A., & Bentler, P. M. (2001). A scaled difference chi-square
test statistic for moment structure analysis. *Psychometrika, 66*(4),
507–514. <https://doi.org/10.1007/BF02296192>

Vandenberg, R. J., & Lance, C. E. (2000). A review and synthesis of the
measurement invariance literature: Suggestions, practices, and
recommendations for organizational research. *Organizational Research
Methods, 3*(1), 4–70. <https://doi.org/10.1177/109442810031002>
