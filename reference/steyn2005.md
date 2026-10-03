# Group statistics from Steyn's (2005) eight-group study

The sample size, pretest and posttest means, and standard deviations of
the eight groups of a Solomon design with three treatments: a study of
how information about one's own ability changes self-efficacy
perceptions, with 1,723 police trainees (Steyn, 2005, Table 5.59, pp.
151–152). It is the eight-group study that Steyn (2009) describes.
Numbers reported in the publication are reused with citation.

## Usage

``` r
steyn2005
```

## Format

A data frame with 8 rows, one per group, and 8 variables:

- group:

  The group's label in the thesis: `EG1` to `EG3` are the pretested
  treatment groups, `KG1.1` to `KG1.3` the unpretested treatment groups,
  `KG2` the pretested control, and `KG3` the unpretested control.

- condition:

  `Norms`, `Marking`, `Test`, or `Control`.

- pretested:

  1 if the group was pretested, 0 if not.

- n:

  Sample size.

- pre_mean, pre_sd:

  Pretest mean and standard deviation (pretested groups only).

- mean, sd:

  Posttest mean and standard deviation.

## Source

Steyn, R. (2005). *Self-evaluasie en die vorming van
selfdoeltreffendheidspersepsies* \[Self-evaluation and the forming of
self-efficacy perceptions\] \[Doctoral thesis, University of South
Africa\]. Unisa Institutional Repository.
https://hdl.handle.net/10500/1745

## Details

**Design.** Three treatments and a control, each with and without a
pretest (pp. 103–105). Each treatment added a source of information
about the participant's ability (p. 102):

- **Test.** Completing a 60-item cognitive test.

- **Marking.** Completing the test and marking one's own answers.

- **Norms.** Completing and marking the test, and receiving the test's
  norms.

The outcome is the total score on a 45-item questionnaire of
self-efficacy perceptions (p. 94). The posttest followed the treatment
after a break of 10 minutes (p. 107).

**Assignment.** Participants were not randomized individually. Fifty-six
existing classes of a police training college were allocated to the
eight groups, seven classes to each, in consultation with the college's
management, and the thesis states that this allocation was not random
(pp. 105–106). The classes had been formed from the order in which
trainees reported for training. The thesis is not consistent on this
point. Its introduction says that participants were divided at random
into eight groups and that the groups were assigned at random to the
conditions (p. 9), and it calls the allocation random again on p. 107.
The method chapter, followed here, gives the detail: the researcher
judged the existing classes to be random groups and made no further
random assignment (pp. 105–106). The analyses in the thesis, and those
below, treat participants as the units, so they do not allow for the
classes. See
[`baseline_solomon()`](https://juhalt.github.io/solomonR/reference/baseline_solomon.md)
for nonrandomized designs and
[`validate_solomon()`](https://juhalt.github.io/solomonR/reference/validate_solomon.md)
for clustered ones.

**Known results.** These statistics reproduce the analyses in the
thesis:

- **One-way ANOVA of the eight posttest groups** (Table 5.60, p. 152):
  F(7, 1715) = 4.545.

- **Scheffé tests** (Table 5.61, p. 153), within .001. Only the
  unpretested Test group differs from the two control groups (p = .002
  and p = .011).

- **The 2 x 2 ANOVAs of each treatment against the control** (Tables
  5.21, 5.34, and 5.47; pp. 128, 135, 142), within rounding. For Test:
  intervention F = 21.3, pretest F = 4.5 (p = .033), and interaction F =
  2.0 (p = .152). For Norms, the thesis ran this analysis with 213
  participants in the unpretested group (Table 5.45, pp. 141–142), one
  fewer than in Table 5.59, so its error degrees of freedom are 854 and
  those computed from these statistics are 855. The F statistics for
  Norms (14.0, 0.9, and 0.1) are the same to one decimal.

The treatments lowered the scores. The thesis analyzed the design as
overlapping four-group designs, one with the treatments pooled and one
for each treatment, and then as a one-way analysis of variance of the
eight posttests (pp. 103–105). The joint model of
[`solomon_from_summary()`](https://juhalt.github.io/solomonR/reference/solomon_from_summary.md)
tests the Pretest x Condition interaction once, F(3, 1715) = 1.00, p =
.392.

**Language.** The thesis is in Afrikaans. The condition labels are
solomonR's.

## References

Steyn, R. (2009). Re-designing the Solomon four-group: Can we improve on
this exemplary model? *Design Principles and Practices: An International
Journal, 3*(1), 383–394.
https://doi.org/10.18848/1833-1874/CGP/v03i01/37588

## See also

[`solomon_from_summary()`](https://juhalt.github.io/solomonR/reference/solomon_from_summary.md),
[`fit_solomon_steyn()`](https://juhalt.github.io/solomonR/reference/fit_solomon_steyn.md),
[mai2020](https://juhalt.github.io/solomonR/reference/mai2020.md),
[lana1959](https://juhalt.github.io/solomonR/reference/lana1959.md)

## Examples

``` r
steyn2005
#>   group condition pretested   n pre_mean pre_sd    mean     sd
#> 1   EG1     Norms         1 218  155.422 12.125 155.895 13.263
#> 2   EG2   Marking         1 214  155.724 11.856 156.196 12.545
#> 3   EG3      Test         1 219  155.644 12.784 156.142 13.396
#> 4   KG2   Control         1 218  156.991 11.987 158.917 13.264
#> 5 KG1.1     Norms         0 214       NA     NA 154.827 12.324
#> 6 KG1.2   Marking         0 211       NA     NA 155.739 12.830
#> 7 KG1.3      Test         0 220       NA     NA 153.036 12.628
#> 8   KG3   Control         0 209       NA     NA 158.306 11.878

# One model for all eight groups.
with(steyn2005, solomon_from_summary(n, mean, sd, treat = condition,
                                     pretested = pretested, control = "Control"))
#> Solomon analysis from summary statistics, N-group design
#> --------------------------------------------------------
#> Conditions: Norms, Marking, Test; control: Control. With and without a pretest: 8 groups.
#> Pooled error variance: 163.370 on 1715 df (equal variances assumed)
#> 
#> Two-way ANOVA on the posttest (Type III sums of squares)
#>   Condition            SS = 3937.177  df = 3  F = 8.03  p < .001
#>   Pretest              SS =  739.550  df = 1  F = 4.53  p = 0.034
#>   Pretest x Condition  SS =  489.696  df = 3  F = 1.00  p = 0.392
#>   Error                SS = 280179.215  df = 1715
#> 
#> Contrasts with 95% confidence intervals
#> Comparison          Contrast                       Est (SE)      t    df      p  p adj.            95% CI
#> Norms vs Control    ATE (avg over pretest)   -3.251 (0.872)  -3.73  1715  <.001   <.001  [-4.961, -1.540]
#> Marking vs Control  ATE (avg over pretest)   -2.644 (0.876)  -3.02  1715  0.003   0.003  [-4.362, -0.926]
#> Test vs Control     ATE (avg over pretest)   -4.023 (0.869)  -4.63  1715  <.001   <.001  [-5.727, -2.318]
#> Norms vs Control    Pretest x Treatment       0.457 (1.745)   0.26  1715  0.793   1.000   [-2.965, 3.879]
#> Marking vs Control  Pretest x Treatment      -0.154 (1.752)  -0.09  1715  0.930   1.000   [-3.590, 3.282]
#> Test vs Control     Pretest x Treatment       2.495 (1.738)   1.44  1715  0.151   0.454   [-0.913, 5.903]
#> Norms vs Control    Treatment | pretested    -3.022 (1.224)  -2.47  1715  0.014   0.041  [-5.423, -0.621]
#> Marking vs Control  Treatment | pretested    -2.721 (1.230)  -2.21  1715  0.027   0.047  [-5.133, -0.309]
#> Test vs Control     Treatment | pretested    -2.775 (1.223)  -2.27  1715  0.023   0.047  [-5.173, -0.377]
#> Norms vs Control    Treatment | unpretested  -3.479 (1.243)  -2.80  1715  0.005   0.010  [-5.917, -1.041]
#> Marking vs Control  Treatment | unpretested  -2.567 (1.247)  -2.06  1715  0.040   0.040  [-5.014, -0.120]
#> Test vs Control     Treatment | unpretested  -5.270 (1.235)  -4.27  1715  <.001   <.001  [-7.692, -2.848]
#> 
#> p adj.: adjusted by Holm's (1979) procedure within each contrast, across the 3 comparisons.
#> Confidence intervals are not adjusted.
```
