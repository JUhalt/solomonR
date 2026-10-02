# Count Outcomes: Validating Poisson Fits and marginal_solomon()

This article reports the simulation validation of the Solomon contrasts
for count outcomes. The aims, data-generating mechanisms, estimands,
methods, performance measures, tolerances, and a decision rule about a
negative-binomial option were posted on [issue
\#44](https://github.com/JUhalt/solomonR/issues/44) before any run,
following the ADEMP structure of Morris et al. (2019). The results come
from package commit c9e8da8 and 2000 replications per scenario.

## Design

**Data.** Each participant has a latent baseline X ~ N(0, 1), observed
only in the pretested groups, an exposure E ~ Uniform(0.5, 1.5), and a
count with

log E\[Y\] = log E + b0 + bX X + bT T + bP P + bTP T P.

Counts are Poisson or overdispersed: negative binomial with Var(Y) = μ +
α μ², generated as Poisson counts with a gamma multiplier. The 96
scenarios cross participants per cell (20, 50, 100), the control rate
(0.5, 3), overdispersion (α = 0 or 1), the association of the pretest
with the outcome (bX = 0 or 0.5), and four effect patterns: no effects,
a treatment effect only (bT = log 1.5), a treatment and a pretest effect
(bP = log 1.25), and sensitization only (bTP = log 1.5).

**Estimands.** The model’s contrasts are the conditional log rate
ratios. The marginal rates per unit of exposure have a closed form, from
which come the Solomon contrasts as rate differences and as marginal log
rate ratios. Log-link rate ratios are collapsible (Daniel et al., 2021),
so the model’s treatment and sensitization contrasts estimate the
marginal ones.

**Methods.** `fit_solomon_glm(family = poisson(), exposure = )` with the
default HC3 covariance and with model-based covariance, and
[`marginal_solomon()`](https://juhalt.github.io/solomonR/reference/marginal_solomon.md)
with delta-method intervals.

**Tolerances.** Bias within 2 Monte Carlo standard errors (MCSE);
coverage of 95% intervals from 0.940 to 0.960; model standard error
within 10% of the empirical standard error; and, where the true value is
zero, a rejection rate from 0.040 to 0.060.

## Failures

| Per cell | Control rate | Failed fits (max over scenarios) |
|---------:|-------------:|---------------------------------:|
|       20 |          0.5 |                                6 |
|       50 |          0.5 |                                0 |
|      100 |          0.5 |                                0 |
|       20 |          3.0 |                                0 |
|       50 |          3.0 |                                0 |
|      100 |          3.0 |                                0 |

A fit fails when it does not converge or when a Solomon cell has no
counts at all. 0 of 96 scenarios had more than 5% failed fits.

## Results by method

| Overdispersion | Method | Scale | Meeting every tolerance | Coverage within tolerance |
|:---|:---|:---|:---|:---|
| none (Poisson) | Poisson fit, HC3 | Log rate ratio (model contrast) | 148 of 192 | 175 of 192 |
| none (Poisson) | marginal_solomon(), delta method | Log rate ratio (marginal) | 156 of 192 | 177 of 192 |
| none (Poisson) | marginal_solomon(), delta method | Rate difference | 178 of 192 | 186 of 192 |
| none (Poisson) | Poisson fit, model-based | Log rate ratio (model contrast) | 104 of 192 | 123 of 192 |
| alpha = 1 | Poisson fit, HC3 | Log rate ratio (model contrast) | 129 of 192 | 145 of 192 |
| alpha = 1 | marginal_solomon(), delta method | Log rate ratio (marginal) | 143 of 192 | 153 of 192 |
| alpha = 1 | marginal_solomon(), delta method | Rate difference | 162 of 192 | 170 of 192 |
| alpha = 1 | Poisson fit, model-based | Log rate ratio (model contrast) | 0 of 192 | 0 of 192 |

| Overdispersion | Method | Scale | Per cell | Coverage (range) |
|:---|:---|:---|---:|:---|
| none | Poisson fit, HC3 | Log rate ratio (model contrast) | 20 | 0.941 to 0.966 |
| none | Poisson fit, HC3 | Log rate ratio (model contrast) | 50 | 0.938 to 0.967 |
| none | Poisson fit, HC3 | Log rate ratio (model contrast) | 100 | 0.940 to 0.962 |
| none | marginal_solomon(), delta method | Log rate ratio (marginal) | 20 | 0.941 to 0.972 |
| none | marginal_solomon(), delta method | Log rate ratio (marginal) | 50 | 0.938 to 0.967 |
| none | marginal_solomon(), delta method | Log rate ratio (marginal) | 100 | 0.940 to 0.962 |
| none | marginal_solomon(), delta method | Rate difference | 20 | 0.941 to 0.971 |
| none | marginal_solomon(), delta method | Rate difference | 50 | 0.937 to 0.962 |
| none | marginal_solomon(), delta method | Rate difference | 100 | 0.939 to 0.961 |
| none | Poisson fit, model-based | Log rate ratio (model contrast) | 20 | 0.820 to 0.964 |
| none | Poisson fit, model-based | Log rate ratio (model contrast) | 50 | 0.816 to 0.964 |
| none | Poisson fit, model-based | Log rate ratio (model contrast) | 100 | 0.804 to 0.962 |
| alpha = 1 | Poisson fit, HC3 | Log rate ratio (model contrast) | 20 | 0.913 to 0.954 |
| alpha = 1 | Poisson fit, HC3 | Log rate ratio (model contrast) | 50 | 0.928 to 0.956 |
| alpha = 1 | Poisson fit, HC3 | Log rate ratio (model contrast) | 100 | 0.939 to 0.961 |
| alpha = 1 | marginal_solomon(), delta method | Log rate ratio (marginal) | 20 | 0.913 to 0.954 |
| alpha = 1 | marginal_solomon(), delta method | Log rate ratio (marginal) | 50 | 0.930 to 0.959 |
| alpha = 1 | marginal_solomon(), delta method | Log rate ratio (marginal) | 100 | 0.939 to 0.961 |
| alpha = 1 | marginal_solomon(), delta method | Rate difference | 20 | 0.933 to 0.968 |
| alpha = 1 | marginal_solomon(), delta method | Rate difference | 50 | 0.932 to 0.966 |
| alpha = 1 | marginal_solomon(), delta method | Rate difference | 100 | 0.941 to 0.963 |
| alpha = 1 | Poisson fit, model-based | Log rate ratio (model contrast) | 20 | 0.526 to 0.916 |
| alpha = 1 | Poisson fit, model-based | Log rate ratio (model contrast) | 50 | 0.513 to 0.899 |
| alpha = 1 | Poisson fit, model-based | Log rate ratio (model contrast) | 100 | 0.513 to 0.902 |

## Overdispersion and the decision rule

The protocol fixed a rule before the run: if HC3 inference from the
Poisson fit met the coverage and Type I tolerances in at least 90% of
the supported overdispersed scenario-contrasts with 50 or 100
participants per cell, no negative-binomial option would be added. The
share was 82% (105 of 128), so a negative-binomial option is proposed,
with its own protocol.

Where HC3 fell short, it was with strong overdispersion and small cells.
With a control rate of 3 (mean Pearson dispersion of about 4 to 6), its
coverage ranged from 0.913 to 0.949 with 20 participants per cell and
0.928 to 0.954 with 50, and reached 0.939 to 0.961 with 100. With a
control rate of 0.5, overdispersion was milder and HC3 met the
tolerances in most contrasts at every cell size.

| Control rate | Per cell | HC3 contrasts within the coverage and Type I tolerances |
|---:|---:|:---|
| 0.5 | 20 | 31 of 32 |
| 0.5 | 50 | 27 of 32 |
| 0.5 | 100 | 30 of 32 |
| 3.0 | 20 | 9 of 32 |
| 3.0 | 50 | 20 of 32 |
| 3.0 | 100 | 28 of 32 |

Overdispersed counts (alpha = 1). {.table}

Model-based standard errors understate uncertainty whenever the counts
are overdispersed: their coverage ranged from 0.513 to 0.916 with α = 1.
This is the failure Cameron and Trivedi (2013, p. 72) describe. They
also failed with Poisson counts when the pretest predicts the outcome
(coverage from 0.804 to 0.959), which the protocol did not anticipate.
The unpretested participants have no pretest in the model, so their
counts mix over the unobserved baseline and are overdispersed even when
each person’s count is Poisson. With no pretest effect, model-based
coverage ranged from 0.944 to 0.964. This is the reason the package’s
default is the robust covariance.

## The negative-binomial option

Because robust Poisson inference fell short of the rule above, a
negative-binomial (NB2) option was added under its own protocol, posted
on [issue \#62](https://github.com/JUhalt/solomonR/issues/62) before
implementation and before any run.
`fit_solomon_glm(family = "negative_binomial")` fits NB2 by maximum
likelihood (Venables & Ripley, 2002) with HC3 standard errors by
default, as Cameron and Trivedi (2013, pp. 84–85) advise.

The study reused the datasets of the count study: same generating code,
same seeds. The Poisson results it recomputed reproduce the results
above exactly, so the comparison between the two models is paired.

| Overdispersion | Per cell | Model | Covariance | Contrasts within the coverage and Type I tolerances |
|:---|---:|:---|:---|:---|
| none (Poisson) | 20 | Negative binomial | HC3 | 58 of 64 |
| none (Poisson) | 20 | Negative binomial | model-based | 39 of 64 |
| none (Poisson) | 20 | Poisson | HC3 | 57 of 64 |
| none (Poisson) | 50 | Negative binomial | HC3 | 59 of 64 |
| none (Poisson) | 50 | Negative binomial | model-based | 41 of 64 |
| none (Poisson) | 50 | Poisson | HC3 | 57 of 64 |
| none (Poisson) | 100 | Negative binomial | HC3 | 59 of 64 |
| none (Poisson) | 100 | Negative binomial | model-based | 45 of 64 |
| none (Poisson) | 100 | Poisson | HC3 | 61 of 64 |
| alpha = 1 | 20 | Negative binomial | HC3 | 38 of 64 |
| alpha = 1 | 20 | Negative binomial | model-based | 18 of 64 |
| alpha = 1 | 20 | Poisson | HC3 | 40 of 64 |
| alpha = 1 | 50 | Negative binomial | HC3 | 49 of 64 |
| alpha = 1 | 50 | Negative binomial | model-based | 34 of 64 |
| alpha = 1 | 50 | Poisson | HC3 | 47 of 64 |
| alpha = 1 | 100 | Negative binomial | HC3 | 60 of 64 |
| alpha = 1 | 100 | Negative binomial | model-based | 45 of 64 |
| alpha = 1 | 100 | Poisson | HC3 | 58 of 64 |

Log rate-ratio contrasts, supported scenarios. {.table}

The protocol fixed five rules before the run.

1.  **Overdispersed counts.** With α = 1 and 50 or 100 participants per
    cell, the NB2 fit with HC3 met the coverage and Type I tolerances in
    85% (109 of 128) of the contrasts, and robust Poisson in 82% (105 of
    128). The 90% threshold was not reached, so robust Poisson remains
    the recommendation and the negative-binomial option is an
    alternative. With strong overdispersion (a control rate of 3), the
    NB2 estimates were more precise: their empirical standard errors
    averaged 3.5% smaller than Poisson’s, as Cameron and Trivedi (2013,
    p. 88) found for NB2 data. The small-sample shortfall of the robust
    standard errors remained.
2.  **Model-based standard errors.** They fell outside the coverage
    tolerance in 162 of 384 contrasts (coverage 0.876 to 0.988). As with
    Poisson fits, they fail in the unpretested groups whenever the
    pretest predicts the outcome, even with Poisson counts, because
    those groups’ counts mix over the unobserved baseline. The default
    stays HC3.
3.  **Poisson counts.** With α = 0, the NB2 fit with HC3 fell outside a
    tolerance that robust Poisson met in 2 of 192 contrasts, so the help
    page recommends Poisson when there is no evidence of overdispersion.
    With Poisson counts, θ did not converge in up to 1464 of 2000 fits
    per scenario, which the classed warning
    `solomonR_theta_boundary_warning` reports.
4.  **Designs neither method supports.** With 20 participants per cell
    and strong overdispersion, neither model’s HC3 intervals met the
    coverage tolerance consistently (coverage from 0.913 to 0.949).
5.  **Guidance.** Where both methods fell short (α = 1, control rate 3),
    the Poisson fit’s mean Pearson dispersion ranged from 3.8 to 6.8;
    with a control rate of 0.5 it ranged from 1.4 to 2.0. No dispersion
    threshold was fixed in advance, so none is recommended.

## Summary

- With the default HC3 covariance, Poisson fits give valid inference for
  the Solomon contrasts on the log rate-ratio scale when the counts are
  Poisson or mildly overdispersed.
- With strong overdispersion, HC3 intervals are somewhat too narrow with
  20 to 50 participants per cell and close to nominal with 100, for
  Poisson and negative-binomial fits alike.
- The negative-binomial option gives slightly more precise estimates
  under strong overdispersion but did not reach the threshold for being
  recommended; robust Poisson remains the recommendation.
- [`marginal_solomon()`](https://juhalt.github.io/solomonR/reference/marginal_solomon.md)
  adds rate differences, with delta-method intervals, for both models.
- Model-based standard errors should not be used for counts, from either
  model: they fail under overdispersion, and in the unpretested groups
  whenever the pretest predicts the outcome.

## Reproducibility

The scripts, `count-validation/count-simulation.R` and
`count-validation/nb-simulation.R`, and their results
(`performance.csv`, `run-information.csv`, `nb-performance.csv`, and
`nb-run-information.csv`) are in the package repository. The count study
ran from 2026-09-26 13:36:08 UTC to 2026-09-26 13:59:23 UTC and the
negative-binomial study from 2026-09-27 05:19:25 UTC to 2026-09-27
06:21:12 UTC (package commit e8e374b), each on 11 workers (R version
4.6.1 (2026-06-24 ucrt)).

All works cited in solomonR are listed, with notes on how the package
uses them, on the
[References](https://juhalt.github.io/solomonR/articles/references.md)
page.

## References

Cameron, A. C., & Trivedi, P. K. (2013). *Regression analysis of count
data* (2nd ed.). Cambridge University Press.
<https://doi.org/10.1017/CBO9781139013567>

Daniel, R., Zhang, J., & Farewell, D. (2021). Making apples from
oranges: Comparing noncollapsible effect estimators and their standard
errors after adjustment for different covariate sets. *Biometrical
Journal, 63*(3), 528–557. <https://doi.org/10.1002/bimj.201900297>

Morris, T. P., White, I. R., & Crowther, M. J. (2019). Using simulation
studies to evaluate statistical methods. *Statistics in Medicine,
38*(11), 2074–2102. <https://doi.org/10.1002/sim.8086>

Venables, W. N., & Ripley, B. D. (2002). *Modern applied statistics with
S* (4th ed.). Springer. <https://doi.org/10.1007/978-0-387-21706-2>
