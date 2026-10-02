# Historical Tests: Replicating the Published Error Rates

The historical analysis of the Solomon design is a sequence of tests,
Tests A to I, proposed by Walton Braver and Braver (1988). Two Monte
Carlo studies estimated how often the sequence declares an effect when
there is none (Sawilowsky et al., 1994; Sawilowsky, 1996). This article
reports solomonR’s replication of those published rates with
[`fit_solomon_classic()`](https://juhalt.github.io/solomonR/reference/fit_solomon_classic.md).
The protocol was posted on [issue
\#51](https://github.com/JUhalt/solomonR/issues/51) before any run,
following the ADEMP structure of Morris et al. (2019). It covers the
1988 sequence, the 1995 revision without Test D, and Sawilowsky’s (1996)
alpha allocations. The results come from package commit 0d970c4 and
20,000 replications per condition.

## Design

**Data.** All four groups have the same distribution, so there is no
treatment effect and no pretest effect. The published studies are
replicated in three sets of conditions:

- **Set 1.** Eight distributions with 30 per group and independent
  pretest and posttest: Sawilowsky et al. (1994, Table 2) and Sawilowsky
  (1996, Tables 1 and 4).
- **Set 2.** Normal data with pretest–posttest correlations from .05 to
  .95 (Sawilowsky, 1996, Table 3).
- **Set 3.** Normal data with 3, 10, or 20 per group (Sawilowsky, 1996,
  Table 3).

**Methods.** Each dataset was analyzed once with
[`fit_solomon_classic()`](https://juhalt.github.io/solomonR/reference/fit_solomon_classic.md).
The paths of every flow and allocation were derived from its p-values.
For the first 200 datasets of each condition, the fitted function was
also run with each flow and allocation, and its path had to match; all
42,000 checks matched. Test I, the Stouffer combination, was judged by
two criteria:

- the one-tailed p of the combined z, which the package used when the
  study ran;
- the two-tailed p that Walton Braver and Braver report in their worked
  example (1988, p. 153: z = 2.05, p = .040), which the package uses now
  (see the decisions below).

**Agreement.** A replicated rate agrees with a published one if they
differ by at most two standard errors of the difference. That standard
error allows for the Monte Carlo error of both studies; the published
studies used 10,000 replications. A published table is reproduced if at
least 90% of its values agree and none differs by more than four
standard errors.

## Tests A to H agree with the published rates

For normal data, the rate at which each classical test is reached and
rejects agrees with the published values:

| Source | Flow | Test | Published | Replicated | Difference (SE) |
|:---|:---|:--:|---:|---:|---:|
| Sawilowsky (1996) Table 1 | 1988 | A | 0.0514 | 0.0491 | -0.9 |
| Sawilowsky (1996) Table 1 | 1988 | D | 0.0463 | 0.0500 | 1.4 |
| Sawilowsky (1996) Table 1 | 1988 | E | 0.0205 | 0.0180 | -1.5 |
| Sawilowsky (1996) Table 1 | 1988 | H | 0.0184 | 0.0181 | -0.2 |
| Sawilowsky (1996) Table 1 | 1995 | A | 0.0514 | 0.0491 | -0.9 |
| Sawilowsky (1996) Table 1 | 1995 | E | 0.0337 | 0.0346 | 0.4 |
| Sawilowsky (1996) Table 1 | 1995 | H | 0.0326 | 0.0347 | 0.9 |
| Sawilowsky et al. (1994) Table 2 | 1988 | A | 0.0516 | 0.0491 | -0.9 |
| Sawilowsky et al. (1994) Table 2 | 1988 | D | 0.0463 | 0.0500 | 1.4 |
| Sawilowsky et al. (1994) Table 2 | 1988 | E | 0.0215 | 0.0180 | -2.0 |
| Sawilowsky et al. (1994) Table 2 | 1988 | H | 0.0184 | 0.0181 | -0.2 |

Rate at which each test is reached and rejects, normal data, 30 per
group. {.table}

Across the eight distributions, Tests A to H agreed in 42 of 56
comparisons. Of the 14 disagreements, 10 were in the uniform and gamma
conditions. The published Test A rates for uniform (0.0308) and gamma
(0.1047) data could not be reproduced. The replicated rates for those
distributions, 0.0510 and 0.0472, are close to the nominal .05 and to
the rates for the other distributions, and no explanation for the
published values was found. Because the later tests are reached only
when Test A is nonsignificant, their published rates in those two
conditions differ as well.

## Test I does not agree, and why

The published conditional rates of Test I are far lower than the
replication’s under either criterion:

| Source | Flow | Allocation | Published | One-tailed | Two-tailed |
|:---|:---|:---|---:|---:|---:|
| Sawilowsky (1996) Table 1 | 1988 | none | 0.0002 | 0.0159 | 0.0019 |
| Sawilowsky (1996) Table 4 | 1995 | method1_conservative | 0.0042 | 0.0123 | 0.0089 |
| Sawilowsky (1996) Table 4 | 1995 | method1_liberal | 0.0045 | 0.0159 | 0.0107 |
| Sawilowsky (1996) Table 4 | 1995 | method2_conservative | 0.0020 | 0.0040 | 0.0032 |
| Sawilowsky (1996) Table 4 | 1995 | method2_liberal | 0.0041 | 0.0123 | 0.0089 |
| Sawilowsky (1996) Table 1 | 1995 | none | 0.0071 | 0.0239 | 0.0164 |
| Sawilowsky et al. (1994) Table 2 | 1988 | none | 0.0002 | 0.0159 | 0.0019 |

Rate at which Test I is reached and rejects, normal data, 30 per group.
{.table}

Test I agreed with the published value in 0 of 21 comparisons under the
one-tailed criterion and 0 of 21 under the two-tailed criterion. Under
the protocol’s decision rule 2, the published rates support neither
criterion, so the rule does not decide it (see the decisions below).

**Investigation (post hoc).** The disagreement was investigated after
the run. The investigation was not part of the protocol. A separate
implementation of the sequence, written without solomonR, reproduced the
package’s Test I rates. It also reproduced the published ones, but only
under a different reading of the Stouffer step:

- **Walton Braver and Braver (1988, p. 152).** Their Stouffer step
  combines the z values of the *one-tailed* p-values of Tests E and H,
  in the direction of the hypothesized effect.
- **The reading that fits.** Each test’s *two-sided* p-value is
  converted to z as if it were one-tailed, so the sign of the effect is
  ignored, and the sum is compared with the one-tailed critical value.

| Flow and allocation | Measure | Published | Directional, one-tailed | Directional, two-tailed | Two-sided p read as one-tailed |
|:---|:---|---:|---:|---:|---:|
| 1988 none | Test I | 0.0002 | 0.0164 | 0.0021 | 0.0001 |
| 1995 method1_conservative | Test I | 0.0042 | 0.0119 | 0.0092 | 0.0039 |
| 1995 method1_liberal | Test I | 0.0045 | 0.0158 | 0.0120 | 0.0048 |
| 1995 method2_conservative | Test I | 0.0020 | 0.0034 | 0.0030 | 0.0014 |
| 1995 method2_liberal | Test I | 0.0041 | 0.0119 | 0.0092 | 0.0038 |
| 1995 none | Test I | 0.0071 | 0.0255 | 0.0173 | 0.0063 |
| 1988 none | Any rejection | 0.1368 | 0.1497 | 0.1354 | 0.1335 |
| 1995 method1_conservative | Any rejection | 0.0536 | 0.0608 | 0.0581 | 0.0527 |
| 1995 method1_liberal | Any rejection | 0.0722 | 0.0822 | 0.0784 | 0.0712 |
| 1995 method2_conservative | Any rejection | 0.0569 | 0.0571 | 0.0567 | 0.0552 |
| 1995 method2_liberal | Any rejection | 0.0763 | 0.0838 | 0.0810 | 0.0756 |
| 1995 none | Any rejection | 0.1248 | 0.1409 | 0.1328 | 0.1218 |

Independent implementation, normal data, 30 per group, 1e+05
replications. {.table}

That reading reproduces every published normal-data rate of Test I, and
the experiment-wise rates, within Monte Carlo error. This suggests that
the published studies computed Test I differently from its definition.
The suggestion is an inference from the numbers: the published reports
do not describe their computation in that detail.

## The error rate of the sequence as defined

For the procedure as Walton Braver and Braver define it, the error rates
are higher than those published:

| Sequence | Published | Replicated, one-tailed Test I | Replicated, two-tailed Test I |
|----|----|----|----|
| 1988 (with Test D) | 0.1368 | 0.1512 | 0.1371 |
| 1995 (without Test D) | 0.1248 | 0.1422 | 0.1348 |

Both sequences falsely declare an effect nearly three times as often as
the nominal 5%. The conclusion of Sawilowsky et al. (1994) stands: the
sequence is not a safe analysis.

Across pretest–posttest correlations and group sizes (Sets 2 and 3), the
1995 sequence’s error rate ranged from 0.1302 to 0.1459 with the
two-tailed criterion. The published values from Sawilowsky’s (1996)
Table 3 agreed in 0 of 22 comparisons. Those values are totals that
include Test I.

## Sawilowsky’s alpha allocations

Sawilowsky (1996, Table 4) chose test-wise levels so that the
experiment-wise error rate of the 1995 sequence stays within the
conservative (.055) or liberal (.075) robustness limit of Bradley (1968,
as cited in Sawilowsky, 1996):

| Allocation | Limit | Published | Replicated, one-tailed | Replicated, two-tailed |
|:---|---:|---:|---:|---:|
| method1_conservative | 0.055 | 0.0536 | 0.0618 | 0.0585 |
| method1_liberal | 0.075 | 0.0722 | 0.0840 | 0.0788 |
| method2_conservative | 0.055 | 0.0569 | 0.0568 | 0.0559 |
| method2_liberal | 0.075 | 0.0763 | 0.0843 | 0.0810 |

Experiment-wise error rate under each allocation, normal data, 30 per
group. {.table}

The allocations reproduce the published rates only under the reading of
Test I described above. With Test I as defined, 4 of the 4 allocations
exceed their limit under the two-tailed criterion. As the protocol’s
decision rule 3 requires,
[`fit_solomon_classic()`](https://juhalt.github.io/solomonR/reference/fit_solomon_classic.md)
implements Sawilowsky’s published levels, and its help page reports
these replicated rates next to his.

## The 1990 amendment

The 1990 flow has no published error rates. With Test I regarded as
definitive, a significant Test A, Test D, or Test I declared an effect
in 12.2% of datasets (one-tailed Test I). Any significant test on the
path did so in 15.1%.

## Decisions under the protocol

1.  **Flows.** The published tables were not reproduced as a whole: 60
    of 136 published values agreed under the one-tailed criterion and 72
    of 136 under the two-tailed criterion. By rule 1, the flows are not
    described as validated against the published tables. The
    investigation locates the disagreement:
    - Tests A to H and the flow logic agree for normal data, and the
      fitted function’s paths matched in every check.
    - The published Test I rates follow a computation that differs from
      Walton Braver and Braver’s definition.
    - The published Test A rates for uniform and gamma data, and the
      rates of the tests that follow them, remain unexplained.
2.  **Test I criterion.** Neither criterion reproduces the published
    Test I rates, so rule 2 did not decide between them, and the
    question went to the maintainer. The maintainer’s decision was to
    follow the literature. Walton Braver and Braver’s worked example
    reports the two-tailed p of the combined z (1988, p. 153), so
    [`fit_solomon_classic()`](https://juhalt.github.io/solomonR/reference/fit_solomon_classic.md)
    and
    [`power_solomon()`](https://juhalt.github.io/solomonR/reference/power_solomon.md)
    now judge Test I by it. The one-tailed value is still returned. The
    replicated error rates above are given under both criteria.
3.  **Allocations.** The option implements the published levels. The
    help page reports the replicated rates and that 4 of the 4 exceed
    their limits with Test I as defined.

## Reproducibility

The files are in the package repository:

- the pre-registered script, `classic-validation/classic-simulation.R`;
- the published targets, `published.csv`, with table and page;
- the results, `performance.csv` and `run-information.csv`;
- the post hoc investigation, `test-i-investigation.R` and its output.

The run took from 2026-09-27 17:43:25 UTC to 2026-09-27 19:09:28 UTC on
11 workers (R version 4.6.1 (2026-06-24 ucrt)).

All works cited in solomonR are listed, with notes on how the package
uses them, on the
[References](https://juhalt.github.io/solomonR/articles/references.md)
page.

## References

Morris, T. P., White, I. R., & Crowther, M. J. (2019). Using simulation
studies to evaluate statistical methods. *Statistics in Medicine,
38*(11), 2074–2102. <https://doi.org/10.1002/sim.8086>

Sawilowsky, S. S. (1996, June 23). *Controlling experiment-wise Type I
error of meta-analysis in the Solomon four-group design* \[Paper
presentation\]. First International Conference on Multiple Comparisons,
Tel Aviv, Israel. <https://digitalcommons.wayne.edu/coe_tbf/29/>

Sawilowsky, S. S., Kelley, D. L., Blair, R. C., & Markman, B. S. (1994).
Meta-analysis and the Solomon four-group design. *The Journal of
Experimental Education, 62*(4), 361–376.
<https://doi.org/10.1080/00220973.1994.9944140>

Walton Braver, M. C., & Braver, S. L. (1988). Statistical treatment of
the Solomon four-group design: A meta-analytic approach. *Psychological
Bulletin, 104*(1), 150–154.
<https://doi.org/10.1037/0033-2909.104.1.150>
