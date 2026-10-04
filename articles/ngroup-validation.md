# Designs With Several Treatments: Validating the Joint Model

This article reports the simulation validation of
[`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md)
for Solomon designs with several treatments. The aims, data-generating
mechanisms, estimands, methods, performance measures, and decision rules
were posted on [issue
\#45](https://github.com/JUhalt/solomonR/issues/45) before the study was
run, following the ADEMP structure of Morris et al. (2019). One
amendment was posted before the run. It changed the post hoc test of
Steyn’s sequence to Scheffé’s (1953) test, which Steyn (2005) used, and
kept the Holm-adjusted version as a second variant. The study has 5,000
replications in each of 112 scenarios.

## Design

**Data.** Normal outcomes, generated as
[`simulate_solomon()`](https://juhalt.github.io/solomonR/reference/simulate_solomon.md)
does, with a posttest standard deviation of 10 in every group when the
spread is equal. A design has k treatments and a control, each with and
without a pretest.

**Scenarios.** They cross:

- **treatments:** two (six groups) or three (eight groups);
- **participants per group:** 10, 25, or 50, and, for two treatments,
  the unequal posttest counts of Mai et al. (2020), 15 to 27 per group;
- **pretest-posttest correlation:** 0.3 or 0.7;
- **residual spread:** the part of the posttest that the pretest does
  not predict has the same standard deviation in every group, or one 1.5
  times larger in the first treatment’s groups. Those groups’ posttest
  standard deviation is then 14.6 at a correlation of 0.3 and 12.8 at
  0.7;
- **effects:**
  - S0, none;
  - S1, a pretest effect of 5 points and nothing else;
  - S2, a treatment effect of 5 points for every treatment, with no
    sensitization;
  - S3, a treatment effect and sensitization of 5 points for the first
    treatment only.

**Methods.**

- **M1:** `fit_solomon_glm(control = )` with the pretest adjustment and
  HC3 standard errors, the package’s default.
- **M2:** the same without the pretest adjustment.
- **M3:** M2 with conventional standard errors, the joint two-way
  analysis of variance.
- **M4:** overlapping four-group analyses, a separate two-by-two
  analysis of variance for each pair of conditions, each interaction
  judged at .05. The analyses of each treatment against the control
  follow McCarthy and Tucker (2002); all pairs follow Mai et al. (2020).
- **M5:** the tests of Steyn’s (2009) sequence, computed as
  [`fit_solomon_steyn()`](https://juhalt.github.io/solomonR/reference/fit_solomon_steyn.md)
  computes them, with Scheffé’s (1953) post hoc tests. **M5h** uses
  Holm-adjusted pairwise t tests. As the protocol defined them, the
  outcomes of steps E4 and E5 are counted in every replication, whether
  or not the function’s decision path reaches those steps.

**Decision rules,** for the HC3 fits (M1 and M2). Each allows for Monte
Carlo error by a Bonferroni adjustment across the rates it checks:

1.  No omnibus Type I error and no Holm family’s familywise error rate
    is significantly above .05. The families are the comparisons of each
    treatment with the control and the comparisons of all pairs.
2.  No contrast’s coverage is significantly below .95, for each
    treatment against the control.
3.  No contrast’s bias is significantly different from zero, for each
    treatment against the control.

If the three rules hold,
[`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md)
keeps its stable lifecycle stage for designs with several treatments. If
any fails, that analysis is labeled experimental and the scenarios that
failed are named.

## Agreement with the package

For speed the study computed every method directly. It checked those
computations against the package’s functions on the first 100
replications of every scenario for M1 to M3 and the first 20 for the
sequence. The protocol required the values to agree within 1e-8 and the
outcomes counted for Steyn’s sequence to be identical. The largest
difference in any estimate, standard error, p-value, or adjusted p-value
was 2.6e-12. The outcomes counted for Steyn’s sequence (E1 significant,
a treatment declared effective, any E4 test significant, any E5 pair
significant), rebuilt from the package’s p-values by the protocol’s
rules, were identical in 4480 of 4480 checks. The agreement criterion is
**met**.

## Rule 1: error control

| Rate | Per group | Rates | Mean | Highest | Above .05 (unadjusted test) | Fails rule 1 |
|:---|:---|---:|---:|---:|---:|---:|
| Omnibus test | 10 | 144 | 0.0495 | 0.0688 | 43 | 16 |
| Omnibus test | 25 | 144 | 0.0503 | 0.0614 | 19 | 0 |
| Omnibus test | 50 | 144 | 0.0496 | 0.0584 | 8 | 0 |
| Omnibus test | 15 to 27 | 72 | 0.0494 | 0.0588 | 9 | 0 |
| Holm family, against the control | 10 | 208 | 0.0387 | 0.0570 | 3 | 0 |
| Holm family, against the control | 25 | 208 | 0.0428 | 0.0530 | 0 | 0 |
| Holm family, against the control | 50 | 208 | 0.0443 | 0.0556 | 1 | 0 |
| Holm family, against the control | 15 to 27 | 104 | 0.0426 | 0.0522 | 0 | 0 |
| Holm family, all pairs | 10 | 256 | 0.0329 | 0.0586 | 1 | 0 |
| Holm family, all pairs | 25 | 256 | 0.0365 | 0.0490 | 0 | 0 |
| Holm family, all pairs | 50 | 256 | 0.0385 | 0.0552 | 0 | 0 |
| Holm family, all pairs | 15 to 27 | 128 | 0.0356 | 0.0520 | 0 | 0 |

Rejection rates of the HC3 fits (M1 and M2) where the null hypothesis
holds. An unadjusted one-sided test at .05 flags about 5% of valid rates
by chance. {.table}

Rule 1 checked 2,128 rates. 16 failed it, so the rule is **not met**.
With an unadjusted one-sided test at .05, 84 rates lay above .05, where
about 106 would be expected by chance for a valid method.

| Method | Treatments | Per group | Correlation | Spread | Effects | Rate | Test | Estimate |
|:---|---:|:---|---:|:---|:---|:---|:---|:---|
| M2 | 3 | 10 | 0.3 | unequal | S0 | Omnibus test | Condition \| pretested | 0.0688 |
| M2 | 3 | 10 | 0.7 | unequal | S1 | Omnibus test | Condition \| pretested | 0.0688 |
| M1 | 3 | 10 | 0.7 | equal | S1 | Omnibus test | Condition \| unpretested | 0.0662 |
| M2 | 3 | 10 | 0.7 | equal | S1 | Omnibus test | Condition \| unpretested | 0.0662 |
| M2 | 3 | 10 | 0.3 | equal | S1 | Omnibus test | Condition \| pretested | 0.0658 |
| M2 | 3 | 10 | 0.7 | unequal | S0 | Omnibus test | Condition \| pretested | 0.0652 |
| M1 | 3 | 10 | 0.3 | equal | S0 | Omnibus test | Condition \| unpretested | 0.0648 |
| M2 | 3 | 10 | 0.3 | equal | S0 | Omnibus test | Condition \| unpretested | 0.0648 |
| M1 | 3 | 10 | 0.3 | unequal | S1 | Omnibus test | Condition \| unpretested | 0.0648 |
| M2 | 3 | 10 | 0.3 | unequal | S1 | Omnibus test | Condition \| unpretested | 0.0648 |
| M2 | 3 | 10 | 0.7 | equal | S0 | Omnibus test | Condition \| pretested | 0.0638 |
| M2 | 3 | 10 | 0.3 | unequal | S1 | Omnibus test | Condition \| pretested | 0.0636 |
| M1 | 3 | 10 | 0.3 | unequal | S0 | Omnibus test | Condition \| unpretested | 0.0634 |
| M2 | 3 | 10 | 0.3 | unequal | S0 | Omnibus test | Condition \| unpretested | 0.0634 |
| M1 | 3 | 10 | 0.7 | unequal | S1 | Omnibus test | Condition \| unpretested | 0.0632 |
| M2 | 3 | 10 | 0.7 | unequal | S1 | Omnibus test | Condition \| unpretested | 0.0632 |

Rates that fail rule 1. {.table style="width:100%;"}

Every failing rate is an omnibus test of Condition \| pretested or
Condition \| unpretested, with 3 treatments and 10 participants per
group; the failing rates run from 0.063 to 0.069. No Holm family failed.

The two omnibus tests that carry the design’s main questions, Pretest x
Condition and Condition averaged over pretest, had rejection rates of
0.036 to 0.055, and 0 of their 280 rates failed the rule.

| Omnibus test | 10 | 25 | 50 | 15 to 27 |
|:---|:---|:---|:---|:---|
| Condition (avg over pretest) | 0.0420 (0.0478) | 0.0467 (0.0502) | 0.0495 (0.0544) | 0.0472 (0.0528) |
| Condition \| pretested | 0.0575 (0.0688) | 0.0537 (0.0614) | 0.0516 (0.0584) | 0.0511 (0.0588) |
| Condition \| unpretested | 0.0594 (0.0662) | 0.0539 (0.0590) | 0.0510 (0.0582) | 0.0540 (0.0584) |
| Pretest x Condition | 0.0425 (0.0526) | 0.0481 (0.0548) | 0.0473 (0.0534) | 0.0469 (0.0514) |

Type I error of each omnibus test by participants per group: mean
(highest) over the scenarios and the two HC3 fits. {.table
style="width:100%;"}

| Factor | Level | Checked | Flagged at an unadjusted .05 | Expected by chance | Fails rule 1 |
|:---|:---|---:|---:|---:|---:|
| Treatments | 2 | 1216 | 29 | 60.8 | 0 |
| Treatments | 3 | 912 | 55 | 45.6 | 16 |
| Per group | 10 | 608 | 47 | 30.4 | 16 |
| Per group | 25 | 608 | 19 | 30.4 | 0 |
| Per group | 50 | 608 | 9 | 30.4 | 0 |
| Per group | 15 to 27 | 304 | 9 | 15.2 | 0 |
| Contrast | Condition (avg over pretest) | 112 | 0 | 5.6 | 0 |
| Contrast | Pretest x Condition | 168 | 0 | 8.4 | 0 |
| Contrast | Condition \| pretested | 112 | 35 | 5.6 | 6 |
| Contrast | Condition \| unpretested | 112 | 44 | 5.6 | 10 |
| Contrast | ATE (avg over pretest) | 392 | 0 | 19.6 | 0 |
| Contrast | Pretest x Treatment | 448 | 0 | 22.4 | 0 |
| Contrast | Treatment \| pretested | 392 | 3 | 19.6 | 0 |
| Contrast | Treatment \| unpretested | 392 | 2 | 19.6 | 0 |
| Spread | equal | 1064 | 36 | 53.2 | 6 |
| Spread | unequal | 1064 | 48 | 53.2 | 10 |

Rule 1 by number of treatments, participants per group, test or
contrast, and spread. {.table}

The Holm families were conservative: their familywise error rates were
0.012 to 0.059, with a mean of 0.039. With 10 participants per group,
509 of 608 rates were below .05 (mean 0.0388), and 432 were
significantly below it at an unadjusted .05, where about 30 would be
expected by chance. Rates below .05 are conservatism, which the protocol
does not count as a failure.

## Rule 2: coverage

| Method | Per group | Intervals |   Mean | Lowest | Highest | Fails rule 2 |
|:-------|:----------|----------:|-------:|-------:|--------:|-------------:|
| M1     | 10        |       320 | 0.9557 | 0.9454 |  0.9664 |            0 |
| M1     | 25        |       320 | 0.9519 | 0.9432 |  0.9620 |            0 |
| M1     | 50        |       320 | 0.9512 | 0.9390 |  0.9580 |            0 |
| M1     | 15 to 27  |       128 | 0.9529 | 0.9398 |  0.9614 |            0 |
| M2     | 10        |       320 | 0.9543 | 0.9432 |  0.9668 |            0 |
| M2     | 25        |       320 | 0.9515 | 0.9432 |  0.9606 |            0 |
| M2     | 50        |       320 | 0.9507 | 0.9390 |  0.9578 |            0 |
| M2     | 15 to 27  |       128 | 0.9524 | 0.9398 |  0.9622 |            0 |

Coverage of nominal 95% intervals for every contrast of each treatment
against the control. {.table}

Rule 2 checked 2,176 intervals. 0 failed it, so the rule is **met**.
With an unadjusted one-sided test at .05, 43 intervals lay below .95,
where about 109 would be expected by chance for a valid method.

| Factor | Level | Checked | Flagged at an unadjusted .05 | Expected by chance | Fails rule 2 |
|:---|:---|---:|---:|---:|---:|
| Per group | 10 | 640 | 3 | 32.0 | 0 |
| Per group | 25 | 640 | 15 | 32.0 | 0 |
| Per group | 50 | 640 | 20 | 32.0 | 0 |
| Per group | 15 to 27 | 256 | 5 | 12.8 | 0 |
| Contrast | ATE (avg over pretest) | 544 | 2 | 27.2 | 0 |
| Contrast | Pretest x Treatment | 544 | 3 | 27.2 | 0 |
| Contrast | Treatment \| pretested | 544 | 16 | 27.2 | 0 |
| Contrast | Treatment \| unpretested | 544 | 22 | 27.2 | 0 |
| Spread | equal | 1088 | 23 | 54.4 | 0 |
| Spread | unequal | 1088 | 20 | 54.4 | 0 |

Rule 2 by participants per group, contrast, and spread. {.table}

The 1,664 intervals for the treatment-against-treatment comparisons,
which no rule covers, had coverage of 0.9400 to 0.9678; 0 were
significantly below .95 at .05/1664.

## Rule 3: bias

Rule 3 checked 2,176 contrasts. 0 failed it, so the rule is **met**. The
largest bias was 0.241 points, on a scale whose control-group standard
deviation is 10, and 117 contrasts differed from zero at an unadjusted
.05, where about 109 would be expected by chance.

| Factor | Level | Checked | Flagged at an unadjusted .05 | Expected by chance | Fails rule 3 |
|:---|:---|---:|---:|---:|---:|
| Per group | 10 | 640 | 42 | 32.0 | 0 |
| Per group | 25 | 640 | 38 | 32.0 | 0 |
| Per group | 50 | 640 | 32 | 32.0 | 0 |
| Per group | 15 to 27 | 256 | 5 | 12.8 | 0 |
| Contrast | ATE (avg over pretest) | 544 | 39 | 27.2 | 0 |
| Contrast | Pretest x Treatment | 544 | 24 | 27.2 | 0 |
| Contrast | Treatment \| pretested | 544 | 22 | 27.2 | 0 |
| Contrast | Treatment \| unpretested | 544 | 32 | 27.2 | 0 |
| Spread | equal | 1088 | 56 | 54.4 | 0 |
| Spread | unequal | 1088 | 61 | 54.4 | 0 |

Rule 3 by participants per group, contrast, and spread. {.table}

For the treatment-against-treatment comparisons, which no rule covers,
the largest bias was 0.225 points, and 0 of 1,664 differed from zero at
.05/1664.

## The conventional analysis

M3, the joint analysis of variance with conventional standard errors,
assumes equal variances in the groups. No rule applies to it.

| Per group | Equal spread    | Unequal spread  |
|:----------|:----------------|:----------------|
| 10        | 0.0437 (0.0530) | 0.0420 (0.0602) |
| 25        | 0.0467 (0.0580) | 0.0425 (0.0594) |
| 50        | 0.0465 (0.0560) | 0.0433 (0.0604) |
| 15 to 27  | 0.0460 (0.0588) | 0.0421 (0.0542) |

M3: omnibus Type I error and Holm familywise error, mean (highest) over
the scenarios. {.table}

## Overlapping four-group analyses

Where no treatment was sensitized (S0 to S2), the chance that at least
one of the separate interaction tests was significant at .05 was:

| Treatments | Analyses | At least one significant interaction: mean (range) |
|---:|:---|:---|
| 2 | Each treatment against the control | 0.092 (0.084 to 0.104) |
| 2 | Every pair of conditions | 0.121 (0.111 to 0.131) |
| 3 | Each treatment against the control | 0.128 (0.119 to 0.138) |
| 3 | Every pair of conditions | 0.203 (0.192 to 0.217) |

The joint model’s single Pretest x Condition test rejected in 0.036 to
0.055 of replications in the same scenarios (M2, which like M4 uses the
posttests alone).

## The tests of Steyn’s sequence

| Outcome | Scenarios | Scheffé (M5) | Holm (M5h) |
|:---|:---|:---|:---|
| The first test (E1) is significant, no effects at all | S0 | 0.052 (0.042 to 0.061) | 0.052 (0.042 to 0.061) |
| The first test (E1) is significant, pretest effect only | S1 | 0.576 (0.193 to 0.974) | 0.576 (0.193 to 0.974) |
| A treatment with no effect is declared effective | S3 | 0.000 (0.000 to 0.000) | 0.001 (0.000 to 0.002) |
| Any E5 post hoc pair significant when no E4 test is (equal treatments) | S2 | 0.040 (0.029 to 0.050) | 0.042 (0.031 to 0.050) |
| Any E4 test significant, no effects at all | S0 | 0.117 (0.090 to 0.155) | 0.117 (0.090 to 0.155) |
| Any E4 test significant, pretest effect only | S1 | 0.638 (0.282 to 0.970) | 0.638 (0.282 to 0.970) |

Tests of Steyn’s sequence: the probability of each outcome as the
protocol defined it, mean (range) over the scenarios. {.table}

E1 compares all the posttest groups, so it responds to a pretest effect
as well as to the treatments. With no effects at all (S0) its null
hypothesis holds, and the rate is its Type I error. With a pretest
effect only (S1) the pretested groups’ means are 5 points above the
unpretested groups’, so the null hypothesis is false: the rate is the
test’s power to detect the pretest effect, and the sequence continues
although no treatment has an effect.

E4 compares a treatment’s pretested and unpretested groups, and Steyn
cautions about internal validity when they differ. With a pretest effect
and no sensitization (S1), at least one E4 test was significant in 0.28
to 0.97 of replications, depending on the scenario. Averaged over the
scenarios with equal group sizes, with two and three treatments, the
rate was 0.37 with 10 per group, 0.69 with 25 per group, 0.91 with 50
per group (rising with the group size). At the unequal sizes of Mai et
al. (2020), 15 to 27 per group with two treatments, it was 0.45 to 0.57.
That comparison picks up the pretest’s own effect, which is not
sensitization. In the same scenarios the joint model’s Pretest x
Condition test (M1 and M2) rejected in 0.038 to 0.054 of replications,
and 0 of its 56 rates failed rule 1.

These rates count the E4 tests in every replication, as the protocol set
out.
[`fit_solomon_steyn()`](https://juhalt.github.io/solomonR/reference/fit_solomon_steyn.md)
acts on E4 only when E1 is significant and a treatment group differs
from both control groups. With a pretest effect only, no group is
expected to: the pretested groups match the pretested control, and the
unpretested groups match the unpretested control. Its decision path then
usually stops at step E2 and its conclusion does not raise the caution.
The E4 tests are still printed, marked as off the decision path. The
agreement check ran the function on 20 replications of each S1 scenario,
with each post hoc test: E4 was on its decision path in 1 of those 1,120
fits. This count was not planned in the protocol and is descriptive. In
the same way, the S2 rate counts significant E5 post hoc pairs whenever
no E4 test is significant, whereas the function’s conclusion reads the
E5 one-way analysis of variance, and only when E5 is on its path.

## Verdict

Rule 1 is not met, rule 2 is met, and rule 3 is met. As the protocol set
out, the analysis of designs with several treatments is therefore
labeled experimental. The tables above list every rate, interval, and
contrast that failed.

The failures are confined to the omnibus tests named under rule 1. The
Holm-adjusted comparisons, which the package reports for every contrast,
and the Pretest x Condition test did not fail. With few participants per
group and three treatments, judge whether the conditions differ among
pretested, or among unpretested, participants by the adjusted
comparisons rather than by those two omnibus tests.

## Limits of this study

- Outcomes were normal. The four-group studies of binary and count
  outcomes and of clustered designs have not been repeated for several
  treatments.
- The comparisons were each treatment against the control and all pairs.
  Planned comparisons given as weights use the same estimates and
  standard errors, and were not simulated separately.
- Before the study, runs with other seeds were made while the simulation
  code was written. Short runs timed the code and checked it against the
  package. One rate in the timing run looked high, so a further run of
  2,000 replications in each of three scenarios with three treatments
  (S0 and S1; 10 or 25 per group) looked at the omnibus tests. With 10
  per group, the HC3 omnibus test of Condition \| unpretested rejected
  in 6.1% and 6.4% of replications. Nothing in the methods, scenarios,
  or rules was changed in response, and rule 1 judged that test as
  posted. These runs are described in the amendment on issue \#45 and
  are not part of these results.

## Run information

| item                      | value                                    |
|:--------------------------|:-----------------------------------------|
| simulation_commit         | 7905c41e19828755c44f91a3d6d1f3f2cfffd165 |
| agreement_commit          | 459d17dfbbd91053caca52e4a940a4bae34e56ff |
| started                   | 2026-10-02 16:25:21 UTC                  |
| finished                  | 2026-10-02 19:55:40 UTC                  |
| replications_per_scenario | 5000                                     |
| scenarios                 | 112                                      |
| workers                   | 9                                        |
| R_version                 | R version 4.6.1 (2026-06-24 ucrt)        |
| platform                  | x86_64-w64-mingw32                       |
| seed                      | 45045                                    |
| rng                       | L’Ecuyer-CMRG                            |

The replications were run from commit 7905c41. They do not call the
package. The agreement check compared them with the package at commit
459d17d, after the summaries were rebuilt from the saved replications.
The run was stopped once, by a time limit of the session that started
it, after 303 of its 1,120 tasks, and was resumed from the saved tasks.
A resumed run gives the same results as an uninterrupted one, because
each task has its own random-number substream. The start time of this
run was taken from its log file. The replications finished at 17:47 UTC;
the finish time above is that of the rebuilt summaries.

The script, `ngroup-simulation.R`, and the files `performance.csv`,
`agreement.csv`, and `run-information.csv` are in the package repository
under `vignettes/articles/ngroup-validation/`.

## References

Holm, S. (1979). A simple sequentially rejective multiple test
procedure. *Scandinavian Journal of Statistics, 6*(2), 65–70.
<https://www.jstor.org/stable/4615733>

Mai, N. N., Takahashi, Y., & Oo, M. M. (2020). Testing the effectiveness
of transfer interventions using Solomon four-group designs. *Education
Sciences, 10*(4), Article 92. <https://doi.org/10.3390/educsci10040092>

McCarthy, A. M., & Tucker, M. L. (2002). Encouraging community service
through service learning. *Journal of Management Education, 26*(6),
629–647. <https://doi.org/10.1177/1052562902238322>

Morris, T. P., White, I. R., & Crowther, M. J. (2019). Using simulation
studies to evaluate statistical methods. *Statistics in Medicine,
38*(11), 2074–2102. <https://doi.org/10.1002/sim.8086>

Scheffé, H. (1953). A method for judging all contrasts in the analysis
of variance. *Biometrika, 40*(1–2), 87–104.
<https://doi.org/10.1093/biomet/40.1-2.87>

Steyn, R. (2005). *Self-evaluasie en die vorming van
selfdoeltreffendheidspersepsies* \[Self-evaluation and the forming of
self-efficacy perceptions\] \[Doctoral thesis, University of South
Africa\]. Unisa Institutional Repository.
<https://hdl.handle.net/10500/1745>

Steyn, R. (2009). Re-designing the Solomon four-group: Can we improve on
this exemplary model? *Design Principles and Practices: An International
Journal—Annual Review, 3*(1), 383–394.
<https://doi.org/10.18848/1833-1874/CGP/v03i01/37588>
