# Latent Contrasts: Validating the Invariance Check

[`fit_solomon_sem_latent()`](https://juhalt.github.io/solomonR/reference/fit_solomon_sem_latent.md)
compares latent means across the four Solomon groups. That comparison
assumes scalar measurement invariance: equal loadings and intercepts in
every group (Meredith, 1993). This study asked which published
criterion, if any, should decide whether the package refuses the latent
contrasts when invariance fails.

The aims, data-generating mechanisms, estimands, methods, performance
measures, and decision rules were posted on [issue
\#55](https://github.com/JUhalt/solomonR/issues/55) before any run,
following the ADEMP structure of Morris et al. (2019). The results come
from package commit 27bf616, with 1000 replications per scenario.

## Design

**Data.** One latent factor with k = 3, 4, or 6 normally distributed
indicators:

- loadings of 0.8, 0.9, 0.7, 0.8, 0.9, and 0.7 (the first k);
- residual SD 0.6;
- latent variance 1 in each group;
- latent means of 0.4, 0, 0.4, and 0 in the pretested treatment (P1),
  pretested control (P0), unpretested treatment (U1), and unpretested
  control (U0) groups, so there is a treatment effect and no
  sensitization.

There were 30, 60, or 120 participants per group.

**Noninvariance patterns,** all on the last indicator:

| Pattern                             | Step that should fail |
|:------------------------------------|:----------------------|
| 1: none                             | none                  |
| 2: intercept +0.3, pretested groups | scalar                |
| 3: intercept +0.6, pretested groups | scalar                |
| 4: intercept +0.6, U0 only          | scalar                |
| 5: loading -0.3, pretested groups   | metric                |

Patterns 2 and 3 are the response shift a pretest could cause in both
pretested groups. The full factorial has 45 scenarios.

**Criteria.**
[`invariance_solomon()`](https://juhalt.github.io/solomonR/reference/invariance_solomon.md)
with MLR estimation, under three decision criteria:

- **C1:** the scaled chi-square difference test at alpha = .05 (Satorra
  & Bentler, 2001), which Vandenberg and Lance (2000, p. 46) recommend
  as the primary criterion;
- **C2:** the change in CFI, RMSEA, and SRMR, judged against Chen’s
  (2007, pp. 501–502) cutoffs;
- **C3:** both C1 and C2 indicate noninvariance.

**Decision rules** (fixed before the run):

1.  **Eligibility.** A criterion is eligible to govern refusals if it
    falsely rejects metric or scalar invariance at a rate of .060 or
    less in at least 8 of the 9 pattern-1 scenarios, and in every
    pattern-1 scenario with 60 or more per group.
2.  **Choice.** The eligible criterion that detects noninvariance most
    often becomes the refusal rule.
3.  **No eligible criterion.** If none is eligible, the fit runs the
    check, reports it, and warns, without refusing.
4.  **Small samples.** Group sizes at which every criterion detects
    noninvariance less than half the time are named as too small for the
    check to protect the contrasts.

**Amendment.** The simulation script recorded the rejection rate of each
step separately, not the rate of rejecting either step. The rate of
rejecting either step is at least the larger of the two, so eligibility
below uses the larger rate. A criterion that fails with the larger rate
fails with the combined rate too. Here every criterion failed on the
larger rate alone, so the missing combined rate does not change any
decision.

## Failures

24 of 45,000 replications failed, in 6 scenarios, all with 30 or 60 per
group. The most in one scenario was 17 (1.7%), so every scenario is
within the supported range. Most failures were a model that did not
converge, usually the three-indicator configural model, which is just
identified within each group. At least two came from the scaled
difference test, which lavaan could not compute for those data. The
package now reports that criterion as undetermined instead of stopping.
One such replication was reproduced exactly and is kept as a test.

## False rejection when invariance holds

The table gives, for each pattern-1 scenario, the larger of the two
steps’ false-rejection rates. Values above .060 are outside the
tolerance.

| Indicators | n per group | C3: both | C2: change in fit (Chen, 2007) | C1: scaled chi-square difference |
|---:|---:|---:|---:|---:|
| 3 | 30 | 0.130 | 0.400 | 0.130 |
| 4 | 30 | 0.135 | 0.314 | 0.139 |
| 6 | 30 | 0.135 | 0.341 | 0.144 |
| 3 | 60 | 0.097 | 0.301 | 0.097 |
| 4 | 60 | 0.080 | 0.234 | 0.084 |
| 6 | 60 | 0.090 | 0.178 | 0.100 |
| 3 | 120 | 0.063 | 0.075 | 0.068 |
| 4 | 120 | 0.030 | 0.030 | 0.080 |
| 6 | 120 | 0.006 | 0.006 | 0.072 |

No criterion was eligible.

- **The chi-square difference test** falsely rejected at rates of 0.068
  to 0.144.
- **Chen’s cutoffs** ranged from 0.006 to 0.400. They rejected far too
  often with 30 or 60 per group (up to 0.400), as Chen (2007, p. 502)
  anticipated for RMSEA and SRMR in small samples. With 120 per group
  they stayed within the tolerance except with three indicators (0.075).
- **Requiring both** ranged from 0.006 to 0.135.

## Detection of noninvariance

Rates at which each criterion rejected the step that should fail
(patterns 2–5):

| Pattern | Indicators | n per group | C3: both | C2: change in fit (Chen, 2007) | C1: scaled chi-square difference |
|:---|---:|---:|---:|---:|---:|
| 2: intercept +0.3, pretested groups | 3 | 30 | 0.425 | 0.722 | 0.425 |
| 2: intercept +0.3, pretested groups | 4 | 30 | 0.351 | 0.639 | 0.360 |
| 2: intercept +0.3, pretested groups | 6 | 30 | 0.350 | 0.656 | 0.358 |
| 2: intercept +0.3, pretested groups | 3 | 60 | 0.741 | 0.894 | 0.741 |
| 2: intercept +0.3, pretested groups | 4 | 60 | 0.640 | 0.781 | 0.660 |
| 2: intercept +0.3, pretested groups | 6 | 60 | 0.582 | 0.741 | 0.617 |
| 2: intercept +0.3, pretested groups | 3 | 120 | 0.944 | 0.949 | 0.962 |
| 2: intercept +0.3, pretested groups | 4 | 120 | 0.836 | 0.836 | 0.949 |
| 2: intercept +0.3, pretested groups | 6 | 120 | 0.643 | 0.643 | 0.938 |
| 3: intercept +0.6, pretested groups | 3 | 30 | 0.968 | 0.994 | 0.968 |
| 3: intercept +0.6, pretested groups | 4 | 30 | 0.932 | 0.987 | 0.934 |
| 3: intercept +0.6, pretested groups | 6 | 30 | 0.937 | 0.981 | 0.943 |
| 3: intercept +0.6, pretested groups | 3 | 60 | 1.000 | 1.000 | 1.000 |
| 3: intercept +0.6, pretested groups | 4 | 60 | 1.000 | 1.000 | 1.000 |
| 3: intercept +0.6, pretested groups | 6 | 60 | 1.000 | 1.000 | 1.000 |
| 3: intercept +0.6, pretested groups | 3 | 120 | 1.000 | 1.000 | 1.000 |
| 3: intercept +0.6, pretested groups | 4 | 120 | 1.000 | 1.000 | 1.000 |
| 3: intercept +0.6, pretested groups | 6 | 120 | 1.000 | 1.000 | 1.000 |
| 4: intercept +0.6, U0 only | 3 | 30 | 0.891 | 0.973 | 0.892 |
| 4: intercept +0.6, U0 only | 4 | 30 | 0.801 | 0.928 | 0.809 |
| 4: intercept +0.6, U0 only | 6 | 30 | 0.778 | 0.934 | 0.794 |
| 4: intercept +0.6, U0 only | 3 | 60 | 0.995 | 1.000 | 0.995 |
| 4: intercept +0.6, U0 only | 4 | 60 | 0.991 | 0.997 | 0.992 |
| 4: intercept +0.6, U0 only | 6 | 60 | 0.981 | 0.985 | 0.991 |
| 4: intercept +0.6, U0 only | 3 | 120 | 1.000 | 1.000 | 1.000 |
| 4: intercept +0.6, U0 only | 4 | 120 | 1.000 | 1.000 | 1.000 |
| 4: intercept +0.6, U0 only | 6 | 120 | 0.999 | 0.999 | 1.000 |
| 5: loading -0.3, pretested groups | 3 | 30 | 0.388 | 0.721 | 0.388 |
| 5: loading -0.3, pretested groups | 4 | 30 | 0.399 | 0.664 | 0.403 |
| 5: loading -0.3, pretested groups | 6 | 30 | 0.388 | 0.657 | 0.397 |
| 5: loading -0.3, pretested groups | 3 | 60 | 0.623 | 0.874 | 0.623 |
| 5: loading -0.3, pretested groups | 4 | 60 | 0.582 | 0.783 | 0.586 |
| 5: loading -0.3, pretested groups | 6 | 60 | 0.563 | 0.706 | 0.598 |
| 5: loading -0.3, pretested groups | 3 | 120 | 0.922 | 0.945 | 0.923 |
| 5: loading -0.3, pretested groups | 4 | 120 | 0.844 | 0.845 | 0.904 |
| 5: loading -0.3, pretested groups | 6 | 120 | 0.645 | 0.645 | 0.915 |

Mean detection across patterns 2–5 was 0.807 for C1, 0.874 for C2, and
0.782 for C3. In every scenario at least one criterion detected
noninvariance at least half the time (rule 4).

The small shifts were hard to detect with 30 per group. For the
intercept shift of 0.3 (pattern 2) and the loading shift (pattern 5),
the chi-square test detected noninvariance in 0.36 to 0.42 of
replications. Chen’s cutoffs detected it in 0.64 to 0.72, at the cost of
the false rejections above.

## Bias in the sensitization contrast

The sensitization contrast is 0 in every pattern. The table gives its
bias and the coverage of its 95% interval under the scalar model, which
assumes full invariance, and under the partial model, which frees the
last indicator’s intercept.

| Pattern | Indicators | n per group | Model | Bias | MCSE | Coverage |
|:---|---:|---:|:---|---:|---:|---:|
| 1: none | 3 | 30 | partial | -0.013 | 0.014 | 0.952 |
| 1: none | 3 | 30 | scalar | -0.010 | 0.014 | 0.952 |
| 1: none | 4 | 30 | partial | -0.018 | 0.013 | 0.953 |
| 1: none | 4 | 30 | scalar | -0.017 | 0.013 | 0.953 |
| 1: none | 6 | 30 | partial | -0.011 | 0.013 | 0.951 |
| 1: none | 6 | 30 | scalar | -0.011 | 0.013 | 0.944 |
| 1: none | 3 | 60 | partial | 0.003 | 0.010 | 0.941 |
| 1: none | 3 | 60 | scalar | 0.002 | 0.010 | 0.943 |
| 1: none | 4 | 60 | partial | -0.005 | 0.009 | 0.949 |
| 1: none | 4 | 60 | scalar | -0.005 | 0.009 | 0.949 |
| 1: none | 6 | 60 | partial | -0.005 | 0.009 | 0.962 |
| 1: none | 6 | 60 | scalar | -0.004 | 0.009 | 0.963 |
| 1: none | 3 | 120 | partial | -0.005 | 0.007 | 0.947 |
| 1: none | 3 | 120 | scalar | -0.004 | 0.007 | 0.943 |
| 1: none | 4 | 120 | partial | -0.007 | 0.007 | 0.942 |
| 1: none | 4 | 120 | scalar | -0.008 | 0.006 | 0.946 |
| 1: none | 6 | 120 | partial | -0.010 | 0.006 | 0.948 |
| 1: none | 6 | 120 | scalar | -0.011 | 0.006 | 0.945 |
| 2: intercept +0.3, pretested groups | 3 | 30 | partial | 0.012 | 0.014 | 0.953 |
| 2: intercept +0.3, pretested groups | 3 | 30 | scalar | 0.011 | 0.014 | 0.949 |
| 2: intercept +0.3, pretested groups | 4 | 30 | partial | -0.008 | 0.014 | 0.936 |
| 2: intercept +0.3, pretested groups | 4 | 30 | scalar | -0.010 | 0.014 | 0.941 |
| 2: intercept +0.3, pretested groups | 6 | 30 | partial | -0.008 | 0.013 | 0.962 |
| 2: intercept +0.3, pretested groups | 6 | 30 | scalar | -0.008 | 0.013 | 0.963 |
| 2: intercept +0.3, pretested groups | 3 | 60 | partial | -0.005 | 0.009 | 0.956 |
| 2: intercept +0.3, pretested groups | 3 | 60 | scalar | -0.006 | 0.009 | 0.954 |
| 2: intercept +0.3, pretested groups | 4 | 60 | partial | -0.003 | 0.009 | 0.954 |
| 2: intercept +0.3, pretested groups | 4 | 60 | scalar | -0.002 | 0.009 | 0.957 |
| 2: intercept +0.3, pretested groups | 6 | 60 | partial | -0.010 | 0.009 | 0.954 |
| 2: intercept +0.3, pretested groups | 6 | 60 | scalar | -0.012 | 0.009 | 0.952 |
| 2: intercept +0.3, pretested groups | 3 | 120 | partial | -0.009 | 0.006 | 0.952 |
| 2: intercept +0.3, pretested groups | 3 | 120 | scalar | -0.008 | 0.006 | 0.952 |
| 2: intercept +0.3, pretested groups | 4 | 120 | partial | -0.002 | 0.006 | 0.956 |
| 2: intercept +0.3, pretested groups | 4 | 120 | scalar | -0.002 | 0.006 | 0.958 |
| 2: intercept +0.3, pretested groups | 6 | 120 | partial | -0.006 | 0.006 | 0.951 |
| 2: intercept +0.3, pretested groups | 6 | 120 | scalar | -0.005 | 0.006 | 0.949 |
| 3: intercept +0.6, pretested groups | 3 | 30 | partial | -0.013 | 0.014 | 0.944 |
| 3: intercept +0.6, pretested groups | 3 | 30 | scalar | -0.017 | 0.014 | 0.941 |
| 3: intercept +0.6, pretested groups | 4 | 30 | partial | -0.008 | 0.014 | 0.955 |
| 3: intercept +0.6, pretested groups | 4 | 30 | scalar | -0.006 | 0.014 | 0.952 |
| 3: intercept +0.6, pretested groups | 6 | 30 | partial | 0.015 | 0.013 | 0.953 |
| 3: intercept +0.6, pretested groups | 6 | 30 | scalar | 0.017 | 0.013 | 0.957 |
| 3: intercept +0.6, pretested groups | 3 | 60 | partial | -0.010 | 0.009 | 0.958 |
| 3: intercept +0.6, pretested groups | 3 | 60 | scalar | -0.016 | 0.009 | 0.956 |
| 3: intercept +0.6, pretested groups | 4 | 60 | partial | 0.022 | 0.009 | 0.965 |
| 3: intercept +0.6, pretested groups | 4 | 60 | scalar | 0.022 | 0.009 | 0.966 |
| 3: intercept +0.6, pretested groups | 6 | 60 | partial | -0.009 | 0.009 | 0.948 |
| 3: intercept +0.6, pretested groups | 6 | 60 | scalar | -0.010 | 0.009 | 0.947 |
| 3: intercept +0.6, pretested groups | 3 | 120 | partial | -0.005 | 0.006 | 0.969 |
| 3: intercept +0.6, pretested groups | 3 | 120 | scalar | -0.006 | 0.006 | 0.968 |
| 3: intercept +0.6, pretested groups | 4 | 120 | partial | 0.001 | 0.007 | 0.936 |
| 3: intercept +0.6, pretested groups | 4 | 120 | scalar | 0.002 | 0.007 | 0.942 |
| 3: intercept +0.6, pretested groups | 6 | 120 | partial | -0.004 | 0.006 | 0.948 |
| 3: intercept +0.6, pretested groups | 6 | 120 | scalar | -0.004 | 0.006 | 0.943 |
| 4: intercept +0.6, U0 only | 3 | 30 | partial | -0.010 | 0.015 | 0.947 |
| 4: intercept +0.6, U0 only | 3 | 30 | scalar | 0.156 | 0.015 | 0.936 |
| 4: intercept +0.6, U0 only | 4 | 30 | partial | -0.023 | 0.014 | 0.955 |
| 4: intercept +0.6, U0 only | 4 | 30 | scalar | 0.120 | 0.014 | 0.947 |
| 4: intercept +0.6, U0 only | 6 | 30 | partial | 0.004 | 0.013 | 0.950 |
| 4: intercept +0.6, U0 only | 6 | 30 | scalar | 0.080 | 0.013 | 0.949 |
| 4: intercept +0.6, U0 only | 3 | 60 | partial | 0.001 | 0.009 | 0.953 |
| 4: intercept +0.6, U0 only | 3 | 60 | scalar | 0.162 | 0.009 | 0.927 |
| 4: intercept +0.6, U0 only | 4 | 60 | partial | 0.003 | 0.009 | 0.949 |
| 4: intercept +0.6, U0 only | 4 | 60 | scalar | 0.142 | 0.009 | 0.921 |
| 4: intercept +0.6, U0 only | 6 | 60 | partial | 0.005 | 0.009 | 0.946 |
| 4: intercept +0.6, U0 only | 6 | 60 | scalar | 0.081 | 0.009 | 0.938 |
| 4: intercept +0.6, U0 only | 3 | 120 | partial | -0.007 | 0.007 | 0.950 |
| 4: intercept +0.6, U0 only | 3 | 120 | scalar | 0.153 | 0.006 | 0.902 |
| 4: intercept +0.6, U0 only | 4 | 120 | partial | -0.008 | 0.006 | 0.949 |
| 4: intercept +0.6, U0 only | 4 | 120 | scalar | 0.129 | 0.006 | 0.907 |
| 4: intercept +0.6, U0 only | 6 | 120 | partial | -0.003 | 0.006 | 0.961 |
| 4: intercept +0.6, U0 only | 6 | 120 | scalar | 0.072 | 0.006 | 0.942 |

Two patterns of results stand out:

- **A shift common to both pretested groups** (patterns 2 and 3) left
  the sensitization contrast unbiased under the scalar model. The bias
  ranged from -0.017 to 0.022. The largest was 2.5 Monte Carlo standard
  errors from zero, within what chance alone produces across 18 unbiased
  scenarios (95th percentile of the largest deviation: 3.0). The
  contrast is a difference of differences within each pretest condition,
  so a shift shared by P1 and P0 cancels.
- **A shift in one group** (pattern 4, U0 only) biased the contrast
  under the scalar model, by 0.072 to 0.162 latent standard deviations,
  with coverage from 0.902 to 0.949. Freeing the shifted intercept
  removed the bias: under the partial model the largest bias was 1.7
  Monte Carlo standard errors from zero.

## Decision

Under rule 3,
[`fit_solomon_sem_latent()`](https://juhalt.github.io/solomonR/reference/fit_solomon_sem_latent.md)
does not refuse latent contrasts. It runs
[`invariance_solomon()`](https://juhalt.github.io/solomonR/reference/invariance_solomon.md)
on the POST indicators, stores the result in the fit’s `invariance`
element, prints a one-line summary, and warns (class
`solomonR_invariance_warning`) when a criterion does not support scalar
(or partial scalar) invariance. Researchers who established invariance
elsewhere can set `check_invariance = FALSE`.

The warning is a prompt to look at the invariance results, not a
verdict. In groups of 30 to 60, the criteria often flag invariance that
holds. Freeing parameters for a partial-invariance model must rest on
substantive grounds, not on the data at hand (Byrne et al., 1989,
p. 465; Vandenberg & Lance, 2000, p. 38).

## Limitations

- **The data-generating mechanisms.** Normal data, one factor, equal
  group sizes, and noninvariance in one indicator only.
- **The latent contrasts.** Bias and coverage were studied for the
  sensitization contrast only, under patterns 1–4.
- **The criteria.** Other criteria, such as modification indices or
  Bayesian approaches, were not studied.
- **Improper solutions.** With three indicators and 30 per group, lavaan
  often estimated a negative residual variance (a Heywood case). When
  the loading-pattern scenario was rerun to reproduce its failures,
  lavaan gave that warning 663 times across the three models of 1,000
  replications. The study counted such fits like any other, so its rates
  describe the criteria as they behave on these data.

The script is
`vignettes/articles/invariance-validation/invariance-simulation.R` in
the package repository.

All works cited in solomonR are listed, with notes on how the package
uses them, on the
[References](https://juhalt.github.io/solomonR/articles/references.md)
page.

## References

Byrne, B. M., Shavelson, R. J., & Muthén, B. (1989). Testing for the
equivalence of factor covariance and mean structures: The issue of
partial measurement invariance. *Psychological Bulletin, 105*(3),
456–466. <https://doi.org/10.1037/0033-2909.105.3.456>

Chen, F. F. (2007). Sensitivity of goodness of fit indexes to lack of
measurement invariance. *Structural Equation Modeling: A
Multidisciplinary Journal, 14*(3), 464–504.
<https://doi.org/10.1080/10705510701301834>

Meredith, W. (1993). Measurement invariance, factor analysis and
factorial invariance. *Psychometrika, 58*(4), 525–543.
<https://doi.org/10.1007/BF02294825>

Morris, T. P., White, I. R., & Crowther, M. J. (2019). Using simulation
studies to evaluate statistical methods. *Statistics in Medicine,
38*(11), 2074–2102. <https://doi.org/10.1002/sim.8086>

Satorra, A., & Bentler, P. M. (2001). A scaled difference chi-square
test statistic for moment structure analysis. *Psychometrika, 66*(4),
507–514. <https://doi.org/10.1007/BF02296192>

Vandenberg, R. J., & Lance, C. E. (2000). A review and synthesis of the
measurement invariance literature: Suggestions, practices, and
recommendations for organizational research. *Organizational Research
Methods, 3*(1), 4–70. <https://doi.org/10.1177/109442810031002>
