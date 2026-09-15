# Full-information ML analysis for a Solomon Four-Group Design

Fits the maximum-likelihood regression model described by van
Engelenburg (1999). Pretest information is incorporated for the
pretested groups while structurally missing pretests in the unpretested
groups are handled through a separate residual variance.

## Usage

``` r
fit_solomon_ml(
  y_post,
  treat,
  pretested,
  y_pre,
  weights = c("equal"),
  control = list(),
  conf_level = 0.95,
  inference = c("wald", "satterthwaite")
)
```

## Arguments

- y_post:

  Numeric posttest scores.

- treat:

  Treatment indicator coded 0 = control and 1 = treatment.

- pretested:

  Pretest indicator coded 0 = unpretested and 1 = pretested.

- y_pre:

  Numeric pretest scores. These should be missing by design for
  participants assigned to the unpretested groups.

- weights:

  Character. How to define the average treatment effect across pretest
  conditions. Currently `"equal"` gives equal weight to the pretested
  and unpretested treatment effects.

- control:

  Optional list passed to
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html).

- conf_level:

  Confidence level for intervals. Default is 0.95.

- inference:

  How standard errors, tests, and intervals are computed: `"wald"`
  (default) for van Engelenburg's (1999) large-sample Wald inference, or
  `"satterthwaite"` for the small-sample option. The point estimates are
  the same. See the Inference options section.

## Value

An object of class `solomon_ml`.

## Details

The model estimates the treatment effect, pretest effect, Treatment x
Pretest interaction, pretest-posttest slope, and separate residual
standard deviations for pretested and unpretested participants.

## Inference options

The point estimates are maximum-likelihood estimates, which coincide
with separate regressions in the pretested and unpretested groups.

- `inference = "wald"` (default) follows van Engelenburg (1999):
  standard errors come from the observed information matrix, and tests
  and intervals use a normal reference distribution, the usual
  large-sample basis for maximum-likelihood inference.

- `inference = "satterthwaite"` is a small-sample option. Standard
  errors use unbiased residual variances within each pretest condition.
  Contrasts within one condition use t tests with that condition's
  residual degrees of freedom, and contrasts that combine the conditions
  (the ATE and Pretest x Treatment) use Welch-Satterthwaite degrees of
  freedom (Satterthwaite, 1946; Welch, 1947).

In the package's simulation validation (issues \#10 and \#22; 84
scenarios with 2,000 replications each, reported in the article
"Validating fit_solomon_ml()" on the package website), both options
recovered the Solomon contrasts without bias. Wald intervals were too
narrow in small samples: mean coverage of nominal 95% intervals was
0.893 with 6 participants per cell, 0.920 with 10, 0.936 with 20, 0.941
with 30, and 0.948 with 100, and the Pretest x Treatment test rejected a
true null hypothesis in 9.9% of samples with 6 per cell and 5.9% with
30. The small-sample option had mean coverage of 0.949 to 0.950 and Type
I error of 0.050 to 0.053 at every cell size studied, from 6 to 100 per
cell.

When the smallest cell has fewer than 40 participants and `inference` is
not supplied, `fit_solomon_ml()` issues a warning of class
`solomonR_small_sample_warning` that suggests the small-sample option.
Supplying `inference = "wald"` explicitly keeps the default without the
warning. The threshold follows a rule set before the validation results
were examined: 40 is the smallest cell size at which Wald inference had
mean coverage of at least 0.940 and Type I error of at most 0.060, with
equal and unequal residual variances, at that size and every larger size
studied.

## References

Satterthwaite, F. E. (1946). An approximate distribution of estimates of
variance components. *Biometrics Bulletin, 2*(6), 110-114.

van Engelenburg, G. (1999). Statistical analysis for the Solomon
four-group design. University of Twente Research Report 99-06.

Welch, B. L. (1947). The generalization of "Student's" problem when
several different population variances are involved. *Biometrika,
34*(1/2), 28-35.
