# Measurement invariance across the four Solomon groups

**\[experimental\]** Tests whether a set of indicators measures the same
construct in the same way in the four Solomon groups, the condition
latent mean contrasts require (Meredith, 1993). It fits the configural,
metric (equal loadings), and scalar (equal loadings and intercepts)
models in the sequence Vandenberg and Lance (2000, pp. 56–57) recommend
and compares each with the one before it.

## Usage

``` r
invariance_solomon(
  data,
  items,
  treat,
  pretested,
  estimator = "MLR",
  partial = NULL,
  alpha = 0.05
)
```

## Arguments

- data:

  A data frame with the indicators.

- items:

  Names of at least three indicators.

- treat, pretested:

  Treatment and pretest indicators (0/1).

- estimator:

  lavaan estimator, default `"MLR"`.

- partial:

  Optional character vector of freed parameters (see Partial
  invariance).

- alpha:

  Significance level of the chi-square difference test. Default 0.05.

## Value

An object of class `solomon_invariance`: the three `fits`, the fit
indexes in `models`, the step `tests` with the decision under each
criterion, the most constrained level `supported` under each
(`"configural"`, `"metric"`, `"scalar"`, or `"partial scalar"`), and the
`cutoffs` applied.

## Two published criteria

Each step compares the more constrained model with the one before it in
two ways, and the output gives the decision under each.

- **Chi-square difference test** (`noninvariant_chisq`): noninvariance
  when the difference test rejects at `alpha`. For the robust estimators
  MLR, MLM, and MLMV the difference is scaled (Satorra & Bentler, 2001).
  Vandenberg and Lance (2000, p. 46) recommend this test as the primary
  criterion, with the change in CFI as a supplement.

- **Change in fit** (`noninvariant_chen`): noninvariance when CFI drops
  by at least .005 and, in addition, RMSEA rises by at least .010 or
  SRMR by at least .025 (loadings) or .005 (intercepts), the cutoffs
  Chen (2007, pp. 501–502) gives for a total N of 300 or less with
  unequal group sizes. For a total N above 300 with equal group sizes
  they are .010, .015, and .030 or .010. Chen does not cover the two
  mixed cases, for which solomonR uses the small-sample values. Reading
  "a change in CFI, supplemented by a change in RMSEA or SRMR" as
  requiring both is solomonR's; Chen chose CFI as the main criterion (p.
  502). Cheung and Rensvold (2002, pp. 234–235) favor such changes over
  the chi-square difference, which depends on sample size.

The criteria can disagree, and both rest on simulations with two groups
and maximum likelihood estimation of multivariate normal data (Cheung &
Rensvold, 2002, p. 251; Chen, 2007). Chen (2007, p. 502) notes that
RMSEA and SRMR tend to over-reject invariant models when samples are
small, as Solomon groups often are. A simulation study under a protocol
posted on issue \#55 found that neither criterion, nor the two together,
kept false rejections of invariance at or below .060 in Solomon-sized
groups (see the article "Latent Contrasts: Validating the Invariance
Check"). This function therefore reports both and decides nothing for
the user, and
[`fit_solomon_sem_latent()`](https://juhalt.github.io/solomonR/reference/fit_solomon_sem_latent.md)
runs it and warns, rather than refuses, when a criterion flags
noninvariance. The fit indexes are those of the maximum likelihood fit,
whose estimates MLR shares.

## Partial invariance

When scalar invariance fails, latent means can still be compared if the
noninvariant parameters are freed and enough indicators stay invariant
(Byrne et al., 1989, p. 458). `partial` names the freed parameters in
lavaan syntax (for example, `"item3 ~ 1"` for an intercept). They must
be chosen on substantive grounds, not by searching the data (Byrne et
al., 1989, p. 465), and must involve only a minority of the indicators
(Vandenberg & Lance, 2000, p. 38); at least two indicators must stay
fully invariant. The function never chooses them.

## Lifecycle

Experimental. The issue \#55 study found no criterion that holds its
false-rejection rate in Solomon-sized groups (see Two published
criteria), so the criteria reported may change as better small-sample
criteria are published.

## References

Byrne, B. M., Shavelson, R. J., & Muthén, B. (1989). Testing for the
equivalence of factor covariance and mean structures: The issue of
partial measurement invariance. *Psychological Bulletin, 105*(3),
456–466. https://doi.org/10.1037/0033-2909.105.3.456

Chen, F. F. (2007). Sensitivity of goodness of fit indexes to lack of
measurement invariance. *Structural Equation Modeling: A
Multidisciplinary Journal, 14*(3), 464–504.
https://doi.org/10.1080/10705510701301834

Cheung, G. W., & Rensvold, R. B. (2002). Evaluating goodness-of-fit
indexes for testing measurement invariance. *Structural Equation
Modeling: A Multidisciplinary Journal, 9*(2), 233–255.
https://doi.org/10.1207/S15328007SEM0902_5

Meredith, W. (1993). Measurement invariance, factor analysis and
factorial invariance. *Psychometrika, 58*(4), 525–543.
https://doi.org/10.1007/BF02294825

Satorra, A., & Bentler, P. M. (2001). A scaled difference chi-square
test statistic for moment structure analysis. *Psychometrika, 66*(4),
507–514. https://doi.org/10.1007/BF02296192

Vandenberg, R. J., & Lance, C. E. (2000). A review and synthesis of the
measurement invariance literature: Suggestions, practices, and
recommendations for organizational research. *Organizational Research
Methods, 3*(1), 4–70. https://doi.org/10.1177/109442810031002

## See also

[`fit_solomon_sem_latent()`](https://juhalt.github.io/solomonR/reference/fit_solomon_sem_latent.md)

## Examples

``` r
# \donttest{
if (requireNamespace("lavaan", quietly = TRUE)) {
  set.seed(55)
  g <- rep(1:4, each = 80)
  treat <- c(1, 0, 1, 0)[g]
  pretested <- c(1, 1, 0, 0)[g]
  f <- stats::rnorm(320, 0.4 * treat)
  items <- data.frame(y1 = f + stats::rnorm(320, 0, 0.6),
                      y2 = 0.9 * f + stats::rnorm(320, 0, 0.6),
                      y3 = 0.8 * f + stats::rnorm(320, 0, 0.6),
                      y4 = 0.7 * f + stats::rnorm(320, 0, 0.6))
  invariance_solomon(items, c("y1", "y2", "y3", "y4"), treat, pretested)
}
#> Measurement invariance across the four Solomon groups
#> Indicators: y1, y2, y3, y4; estimator: MLR; group sizes (P1, P0, U1, U0): 80, 80, 80, 80
#> 
#>       model chisq df   cfi rmsea  srmr
#>  configural  6.31  8 1.000 0.000 0.012
#>      metric 20.60 17 0.995 0.051 0.063
#>      scalar 32.80 26 0.991 0.057 0.069
#> 
#>             comparison chisq_diff df_diff p.value delta_cfi delta_rmsea
#>  metric vs. configural      15.17       9  0.0863    -0.005       0.051
#>      scalar vs. metric      11.67       9  0.2320    -0.004       0.006
#>  delta_srmr noninvariant_chisq noninvariant_chen
#>       0.051              FALSE             FALSE
#>       0.007              FALSE             FALSE
#> 
#> Chi-square difference test (scaled; Satorra & Bentler, 2001) at alpha = 0.05 (Vandenberg & Lance, 2000, p. 46): scalar supported.
#> Change in fit (Chen, 2007, pp. 501-502; total N > 300 and equal group sizes): scalar supported.
# }
```
