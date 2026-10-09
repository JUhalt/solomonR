# Mixed model for repeated measures in a longitudinal Solomon design

**\[experimental\]** Analyzes a Solomon four-group design with several
posttest occasions by a mixed model for repeated measures (MMRM), and
estimates the four Solomon contrasts and the pretest effects at each
occasion, and the change in pretest sensitization from the first
occasion to the last. The model is likelihood-based, so it is valid when
posttests are missing at random, for example when participants drop out
depending on their earlier scores (Fitzmaurice et al., 2011, pp. 497,
505). Per-occasion analyses of the participants still observed, and
unweighted generalized estimating equations, require the stronger
assumption that the posttests are missing completely at random (Liang &
Zeger, 1986; Fitzmaurice et al., 2011, pp. 358, 498).

## Usage

``` r
fit_solomon_mmrm(
  y_post,
  treat,
  pretested,
  id,
  occasion,
  y_pre = NULL,
  df = c("kenward-roger", "satterthwaite"),
  conf_level = 0.95,
  data = NULL
)
```

## Arguments

- y_post:

  Numeric posttest scores, one per participant and occasion (long
  format), with `NA` for missing posttests.

- treat:

  Treatment indicator coded 0/1 (or logical), the same in every row of a
  participant. Designs with several treatments are not supported; see
  [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md).

- pretested:

  Pretest indicator coded 0/1 (or logical), the same in every row of a
  participant.

- id:

  Participant identifier.

- occasion:

  Posttest occasion. A factor keeps its level order; otherwise the
  sorted unique values are used.

- y_pre:

  Optional pretest score, the same in every row of a participant and
  missing by design for unpretested participants.

- df:

  Degrees of freedom for tests and intervals: `"kenward-roger"`
  (default) or `"satterthwaite"`.

- conf_level:

  Confidence level. Default 0.95.

- data:

  Optional data frame. When supplied, the other data arguments are
  looked up in it first, as bare column names (`y_post = post`) or as
  strings (`y_post = "post"`).

## Value

An object of class `solomon_mmrm` with `effects` (the four Solomon
contrasts and the three pretest effects at each occasion, then the
change in sensitization), `model` (the mmrm fit), `covariance` (the
structure used), `observed` (observed posttests by group and occasion),
`pretest_mean` (the center of the pretest; `NA` without `y_pre`), and
the settings used.

## Model

The specification follows the example of Mallinckrodt et al. (2008, p.
312), adapted to the Solomon design:

- **Fixed effects.** Occasion x Treatment x Pretested, all categorical.

- **Pretest adjustment.** The pretest enters as in
  [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md):
  in the pretested groups, the score minus the mean pretest of the
  pretested participants in the model (each counted once; returned as
  `pretest_mean`), and 0 in the unpretested groups, with a separate
  slope at each occasion ("a full interaction of the covariate with
  time", p. 312). Pretests absent by design are never imputed.

- **Pretest effects.** At each occasion, the pretest effects among
  controls, among treated participants, and their average compare
  pretested and unpretested participants at that mean pretest, as
  described in the section "The pretest effect" of
  [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md).
  Pretest effects may fade with time (Entwisle, 1961, p. 610), and these
  estimates show whether they do. Their standard errors include the
  sampling variance of the mean pretest, s^2 / n for n pretested
  participants, times the square of that occasion's pretest slope, and
  their degrees of freedom combine those of the contrast with n - 1 by
  the Welch-Satterthwaite formula (Satterthwaite, 1946; Welch, 1947).
  The pretest effects were not part of the validation study below.

- **Covariance.** Unstructured within participant, estimated by
  restricted maximum likelihood (Laird & Ware, 1982). With a pretest, it
  is estimated separately for pretested and unpretested participants:
  the pretest adjustment makes the residual covariance conditional on
  the pretest in the pretested groups only, so a single shared matrix
  would be misspecified.

- **Fallback.** If the unstructured model does not converge, the
  heterogeneous Toeplitz, heterogeneous AR(1), and heterogeneous
  compound symmetry structures are fitted, and the converged one with
  the smallest AIC is used (Mallinckrodt et al., 2008, p. 312). The
  output names the structure used.

- **Inference.** t tests and intervals with Kenward-Roger degrees of
  freedom and adjusted covariance (Kenward & Roger, 1997, as cited in
  Fitzmaurice et al., 2011, p. 101), as Mallinckrodt et al.
  (2008, p. 312) specify, or with Satterthwaite (1946) degrees of
  freedom. The Kenward-Roger adjustment depends on how the covariance is
  parameterized; it is computed with the matrix's own elements as the
  parameters, the form that reproduces SAS PROC MIXED (mmrm's
  "Kenward-Roger-Linear"; Sabanes Bove et al., 2026). Fitzmaurice et al.
  (2011, p. 101) note that the small-sample properties of both
  approximations in longitudinal models "have not been extensively
  studied"; the simulation study posted on issue \#57 compares them in
  Solomon designs.

The model is fitted with the mmrm package (Sabanes Bove et al., 2026),
which must be installed.

## Validation

A simulation study under a protocol posted on issue \#57 before any run
(12 scenarios, 2,000 replications each; see the article "Longitudinal
Designs: Validating the Repeated-Measures Analysis") generated monotone
dropout that depended on the last posttest. With 30 to 120 participants
per group:

- **Coverage and Type I error.** With Kenward-Roger degrees of freedom,
  95% intervals covered from 0.937 to 0.962 of the time, and Type I
  error where the true value was zero ranged from 0.041 to 0.0615.
  Satterthwaite degrees of freedom performed alike.

- **Bias.** The contrasts were unbiased; the largest bias was 0.029 SD.

- **Normal reference distribution.** It gave coverage as low as 0.935 at
  30 per group.

- **Shared covariance.** A covariance shared by all four groups
  misstated the standard errors by up to 12% for the unpretested
  contrasts and 17% for the pretested ones.

- **Complete-case analyses.** Per-occasion analyses of the participants
  still observed were biased by up to 0.09 SD at the last occasion.

## Lifecycle

Experimental. The study's pre-specified rule for choosing the default
degrees of freedom, every tolerance met in 90% of cells and in every
cell with 60 or more per group, was not met by either approximation (145
of 156 cells each; 97 and 96 of 104), although Monte Carlo error alone
makes the second part unlikely to be met even by a correctly calibrated
method. Kenward-Roger remains the default, as Mallinckrodt et al. (2008)
specify.

## References

Entwisle, D. R. (1961). Interactive effects of pretesting. *Educational
and Psychological Measurement, 21*(3), 607–620.
https://doi.org/10.1177/001316446102100307

Fitzmaurice, G. M., Laird, N. M., & Ware, J. H. (2011). *Applied
longitudinal analysis* (2nd ed.). Wiley.
https://doi.org/10.1002/9781119513469

Laird, N. M., & Ware, J. H. (1982). Random-effects models for
longitudinal data. *Biometrics, 38*(4), 963–974.
https://doi.org/10.2307/2529876

Liang, K.-Y., & Zeger, S. L. (1986). Longitudinal data analysis using
generalized linear models. *Biometrika, 73*(1), 13–22.
https://doi.org/10.1093/biomet/73.1.13

Mallinckrodt, C. H., Lane, P. W., Schnell, D., Peng, Y., & Mancuso, J.
P. (2008). Recommendations for the primary analysis of continuous
endpoints in longitudinal clinical trials. *Drug Information Journal,
42*(4), 303–319. https://doi.org/10.1177/009286150804200402

Sabanes Bove, D., Li, L., Dedic, J., Kelkhoff, D., Kunzmann, K., Lang,
B. M., Stock, C., Wang, Y., James, D., Sidi, J., Leibovitz, D., Sjoberg,
D. D., Krieger, N. I., Panagos, A., & Jones, J. (2026). *mmrm: Mixed
models for repeated measures* (Version 0.3.18) \[R package\].
https://doi.org/10.32614/CRAN.package.mmrm

Satterthwaite, F. E. (1946). An approximate distribution of estimates of
variance components. *Biometrics Bulletin, 2*(6), 110–114.
https://doi.org/10.2307/3002019

Welch, B. L. (1947). The generalization of "Student's" problem when
several different population variances are involved. *Biometrika,
34*(1–2), 28–35. https://doi.org/10.1093/biomet/34.1-2.28

## See also

[`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md)
for a single posttest occasion;
[jordaan2014](https://juhalt.github.io/solomonR/reference/jordaan2014.md)
and the article "Worked Example: Repeated Posttests" for a published
study with three posttest occasions.

## Examples

``` r
# \donttest{
if (requireNamespace("mmrm", quietly = TRUE)) {
  # Three posttest occasions for the example participants, with dropout.
  set.seed(57)
  d <- solomon_example
  d$id <- seq_len(nrow(d))
  long <- do.call(rbind, lapply(1:3, function(t) {
    x <- d
    x$occasion <- t
    x$y_post <- d$y_post - 2 * (t - 1) + rnorm(nrow(d), 0, 4)
    x
  }))
  long$y_post[long$occasion == 3 & runif(nrow(long)) < 0.3] <- NA
  fit_solomon_mmrm(y_post, treat, pretested, id, occasion, y_pre, data = long)
}
#> Solomon MMRM: mixed model for repeated measures (Mallinckrodt et al., 2008)
#> Occasions: 1, 2, 3
#> Covariance: unstructured, separately for pretested and unpretested participants (REML)
#> Degrees of freedom: Kenward-Roger
#> Pretest centered at the pretested participants' mean: 49.600
#> Pretest effects: standard errors include the sampling variance of that mean.
#> 
#> Observed posttests by group and occasion:
#>                        1  2  3
#> pretested treatment   30 30 22
#> pretested control     30 30 26
#> unpretested treatment 30 30 23
#> unpretested control   30 30 20
#> 
#>  Occasion Contrast                      Estimate SE    df    t     p    
#>  1        ATE (avg over pretest)         3.027   1.586 114.6  1.91 0.059
#>  1        Pretest x Treatment           -3.301   3.172 114.6 -1.04 0.300
#>  1        Treatment | pretested          1.377   2.298  57.0  0.60 0.552
#>  1        Treatment | unpretested        4.677   2.187  58.0  2.14 0.037
#>  1        Pretest effect | control       4.639   2.444 151.1  1.90 0.060
#>  1        Pretest effect | treated       1.338   2.444 151.1  0.55 0.585
#>  1        Pretest main effect            2.988   1.859 170.1  1.61 0.110
#>  2        ATE (avg over pretest)         3.107   1.730 115.0  1.80 0.075
#>  2        Pretest x Treatment           -3.606   3.460 115.0 -1.04 0.300
#>  2        Treatment | pretested          1.305   2.447  57.0  0.53 0.596
#>  2        Treatment | unpretested        4.910   2.447  58.0  2.01 0.049
#>  2        Pretest effect | control       3.446   2.600 141.9  1.33 0.187
#>  2        Pretest effect | treated      -0.160   2.600 141.9 -0.06 0.951
#>  2        Pretest main effect            1.643   1.940 161.0  0.85 0.398
#>  3        ATE (avg over pretest)         3.918   1.770 109.1  2.21 0.029
#>  3        Pretest x Treatment           -0.023   3.540 109.1 -0.01 0.995
#>  3        Treatment | pretested          3.906   2.559  56.7  1.53 0.132
#>  3        Treatment | unpretested        3.929   2.446  52.4  1.61 0.114
#>  3        Pretest effect | control       2.897   2.602 127.3  1.11 0.268
#>  3        Pretest effect | treated       2.874   2.633 129.2  1.09 0.277
#>  3        Pretest main effect            2.885   1.928 144.4  1.50 0.137
#>  3 vs 1   Change in Pretest x Treatment  3.277   2.353  89.9  1.39 0.167
#>  95% CI          
#>  [ -0.115, 6.169]
#>  [ -9.584, 2.983]
#>  [ -3.225, 5.978]
#>  [  0.299, 9.055]
#>  [ -0.189, 9.467]
#>  [ -3.490, 6.166]
#>  [ -0.681, 6.658]
#>  [ -0.320, 6.535]
#>  [-10.460, 3.249]
#>  [ -3.596, 6.205]
#>  [  0.013, 9.808]
#>  [ -1.693, 8.585]
#>  [ -5.299, 4.979]
#>  [ -2.189, 5.475]
#>  [  0.409, 7.426]
#>  [ -7.039, 6.993]
#>  [ -1.219, 9.031]
#>  [ -0.979, 8.837]
#>  [ -2.251, 8.045]
#>  [ -2.336, 8.084]
#>  [ -0.926, 6.697]
#>  [ -1.398, 7.953]
# }
```
