# Equivalence test for a Solomon contrast

Tests whether a Solomon contrast, by default the Pretest x Treatment
(sensitization) contrast, is small enough to be considered negligible,
using the two one-sided tests (TOST) procedure (Schuirmann, 1987;
Lakens, 2017). A nonsignificant sensitization test is not evidence that
sensitization is absent; an equivalence test against a prespecified
smallest effect size of interest (SESOI) can provide that evidence.

## Usage

``` r
equivalence_solomon(
  object,
  bounds,
  contrast = "Pretest x Treatment",
  alpha = 0.05
)
```

## Arguments

- object:

  A fit from
  [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md)
  or
  [`fit_solomon_ml()`](https://juhalt.github.io/solomonR/reference/fit_solomon_ml.md).

- bounds:

  Equivalence bounds on the raw posttest scale: one positive number
  `delta`, giving `c(-delta, delta)`, or `c(lower, upper)` with
  `lower < 0 < upper`. There is no default; bounds must be chosen in
  advance.

- contrast:

  Solomon contrast to test. Default is `"Pretest x Treatment"`.

- alpha:

  Significance level for each one-sided test. Default is 0.05.

## Value

An object of class `solomon_equivalence` containing the estimate,
standard error, degrees of freedom, both one-sided tests (`t_lower`,
`p_lower`, `t_upper`, `p_upper`), the equivalence p-value
(`p_equivalence`), the test against zero (`statistic`, `p_zero`), both
confidence intervals, the logical results `equivalent`, `different`, and
`exceeds_bounds`, the `outcome`, and a plain-language `interpretation`.

## Choosing equivalence bounds

The bounds define the smallest effect size of interest and should be set
before the data are examined, for example in a preregistration (Lakens,
2017; Lakens et al., 2018). Justify them substantively, such as the
smallest change in posttest scores that would alter a conclusion, or
from prior research.

Bounds are on the raw posttest scale. To use a standardized SESOI (for
example, d = 0.2), multiply it by a standard deviation fixed in advance,
such as one reported in prior studies. Do not use the standard deviation
of the current data: the same standardized bound then implies different
raw bounds in different samples (Lakens, 2017).

## Inference

Each one-sided test uses the estimate, standard error, and reference
distribution of the fitted model: t with the model's degrees of freedom
for
[`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md)
(Satterthwaite degrees of freedom with CR2); for
[`fit_solomon_ml()`](https://juhalt.github.io/solomonR/reference/fit_solomon_ml.md),
the normal distribution with its default Wald inference or t with
Welch-Satterthwaite degrees of freedom with
`inference = "satterthwaite"`. The equivalence p-value is the larger of
the two one-sided p-values, and the matching interval has confidence
level 1 - 2 `alpha` (90\\ The conventional two-sided test against zero
is reported alongside, with its 1 - `alpha` interval.

## Outcomes

Combining the equivalence test with the test against zero gives four
outcomes (Lakens, 2017):

- `"equivalent"`: statistically equivalent and not different from zero;

- `"trivial"`: different from zero but statistically equivalent, so
  smaller than the smallest effect size of interest;

- `"different"`: different from zero and not statistically equivalent;

- `"inconclusive"`: neither different from zero nor statistically
  equivalent.

`exceeds_bounds` is `TRUE` when the 1 - 2 `alpha` interval lies entirely
beyond one bound, which rejects effects no larger than the smallest
effect size of interest in that direction (a minimum-effect test; Murphy
& Myors, 1999).

## References

Lakens, D. (2017). Equivalence tests: A practical primer for t tests,
correlations, and meta-analyses. *Social Psychological and Personality
Science, 8*(4), 355-362.

Lakens, D., Scheel, A. M., & Isager, P. M. (2018). Equivalence testing
for psychological research: A tutorial. *Advances in Methods and
Practices in Psychological Science, 1*(2), 259-269.

Murphy, K. R., & Myors, B. (1999). Testing the hypothesis that
treatments have negligible effects: Minimum-effect tests in the general
linear model. *Journal of Applied Psychology, 84*(2), 234-248.

Schuirmann, D. J. (1987). A comparison of the two one-sided tests
procedure and the power approach for assessing the equivalence of
average bioavailability. *Journal of Pharmacokinetics and
Biopharmaceutics, 15*(6), 657-680.

## Examples

``` r
data(solomon_example)
fit <- with(solomon_example, fit_solomon_glm(y_post, treat, pretested, y_pre))

# Illustrative bounds only: half a standard deviation (5 points) on this
# test. In a real study, justify the smallest effect size of interest and
# fix the bounds before examining the data.
equivalence_solomon(fit, bounds = 5)
#> Solomon equivalence test (TOST)
#> Contrast: Pretest x Treatment
#> Equivalence bounds (raw scale): [-5.000, 5.000]; alpha = 0.05
#> Inference: HC3 heteroskedasticity-consistent; t tests (df = 115)
#> 
#> Estimate = -1.940 (SE = 3.168)
#> 90% CI [-7.193, 3.313] (equivalence); 95% CI [-8.214, 4.335] (test against zero)
#> 
#> Lower bound test:   t(115) = 0.97, p = 0.168
#> Upper bound test:   t(115) = -2.19, p = 0.015
#> Equivalence (TOST): p = 0.168
#> Test against zero:  t(115) = -0.61, p = 0.541
#> 
#> Conclusion: Inconclusive: the contrast is neither different from zero nor
#>   statistically equivalent.
#> Equivalence bounds must be justified and fixed before the data are examined;
#> see ?equivalence_solomon.
```
