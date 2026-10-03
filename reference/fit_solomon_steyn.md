# Steyn's (2009) analysis of the extended Solomon design

**\[experimental\]** Carries out the sequence of tests Steyn (2009)
proposed for the Solomon four-group design and for its extension to
several interventions, the extended Solomon design: k interventions and
a control, each with and without a pretest, giving 2(k + 1) groups. The
sequence checks the threats to internal validity the design can detect
(nonequivalent groups, history and maturation, testing, the
pretest-intervention interaction, instrumentation, regression to the
mean, and attrition) before it compares the interventions. It is a
published proposal, kept for replication and teaching. The package's
recommended analysis is
[`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md)
with `control =`, which estimates the Solomon contrasts of every
intervention in one model.

## Usage

``` r
fit_solomon_steyn(
  y_post,
  treat,
  pretested,
  y_pre = NULL,
  control = NULL,
  alpha = 0.05,
  include_unpretested_control = FALSE,
  posthoc = c("scheffe", "holm"),
  data = NULL
)
```

## Arguments

- y_post:

  Numeric posttest scores, with `NA` for participants who dropped out.

- treat:

  The condition: a 0/1 (or logical) indicator for one intervention (1 =
  intervention), or a factor or character vector of conditions, with the
  control named by `control`.

- pretested:

  Pretest indicator coded 0/1 (or logical).

- y_pre:

  Optional numeric pretest scores, `NA` for unpretested participants.
  Without them the steps that use the pretest are skipped.

- control:

  The control condition, when `treat` is a factor or character vector.

- alpha:

  Significance level of every test in the sequence. Default 0.05.

- include_unpretested_control:

  Logical. If `TRUE`, the equivalence step adds a one-way ANOVA with the
  unpretested control's posttests (`Of`) as a further group. Default
  `FALSE`.

- posthoc:

  The post hoc tests of pairs of groups: `"scheffe"` (the default),
  Scheffé tests, as in Steyn (2005); or `"holm"`, pairwise t tests with
  Holm's (1979) adjustment. Both use a pooled SD. See the
  Operationalization section.

- data:

  Optional data frame. When supplied, the other data arguments are
  looked up in it first, as bare column names (`y_post = post`) or as
  strings (`y_post = "post"`).

## Value

An object of class `solomon_steyn`, a list with:

- `equivalence`, `testing`, `reliability`, `regression`: tables of tests
  for steps 1, 3, 5, and 6. `testing` adds `comparison`, `term`
  (`"Pretest"`, `"Intervention"`, `"Interaction"`), and `sum_sq`;
  `regression` adds `var_pre`, `var_post`, and `direction`.

- `history`: `tests`, and `pattern`, the classification of step 2.

- `classic`: `summary` (one row per intervention: Test A, the path, and
  the conclusion) and `fits`, the
  [`fit_solomon_classic()`](https://juhalt.github.io/solomonR/reference/fit_solomon_classic.md)
  results.

- `attrition`: `counts` (randomized, observed, and missing posttests in
  each group), `dropouts` (the two-by-two table of dropouts), and
  `tests`.

- `effects`: `tests` (E1 to E5), `groups` (for several interventions,
  whether each intervention group differs from `Od` and from `Of`), and
  `highest` (the intervention with the highest combined posttest mean in
  E5, or `NA`; see the Operationalization section).

- `path`: the steps of step 8 that Steyn's sequence acts on.

- `conclusions`: one plain-language summary for each step.

- `notes`: skipped steps and data notes.

- `n`: counts of participants, posttests, and pretests in each group.

- `conditions`: the control and the interventions.

- `settings`: `alpha`, `include_unpretested_control`, `posthoc`, `k`,
  and `pretest` (whether pretest scores were supplied).

The tables of tests share the columns `step`, `test`, `groups` (Steyn's
labels with the condition names), `n`, `estimate` (the first group's
mean minus the second's for a t test; the mean change for a paired t
test; r for a correlation; the variance ratio; or the difference in
dropout proportions), `statistic`, `reference` (`"t"`, `"F"`, `"chisq"`,
or `"z"`), `df1`, `df2`, `p.value`, `p.adjusted` (the post hoc p-value
of a pair of groups: Scheffé's, or Holm-adjusted, as `posthoc` sets;
`NA` for the other tests), and `interpretation` (Steyn's reading, where
he gives one). For a post hoc test, `statistic` and `p.value` are the
pairwise t with a pooled SD and its unadjusted p-value; the decisions of
the sequence read `p.adjusted`.

## Steyn's sequence

Steyn labels the scores by group: `Oa1`, `Ob1` are the pretest and
posttest of the pretested group of intervention 1 (EG1); `Oc`, `Od` the
pretest and posttest of the pretested control (CG1); `Oe1` the posttest
of the unpretested group of intervention 1 (CG2.1); and `Of` the
posttest of the unpretested control (CG3). With one intervention the
number is dropped (`Oa`, `Ob`, `Oe`; EG, CG2). The steps below follow a
pre-publication draft of the article, dated March 2, 2009; they will be
checked against the published version. In Steyn's order:

1.  **Equivalence after randomization.** The pretests of the pretested
    groups are compared: a t test for one intervention, a one-way ANOVA
    otherwise. If they differ, Steyn advises post hoc tests to locate
    the difference, an explanation, and reconsidering whether to
    continue. With `include_unpretested_control = TRUE`, a further
    one-way ANOVA adds `Of`, which Steyn suggests when the study is
    short or history and maturation are believed negligible.

2.  **History and maturation.** A paired t test of `Oc` and `Od`; a t
    test of all the pretests (the Time 1 value for the unpretested
    control) against `Of`; and a one-way ANOVA of `Oc`, `Od`, and `Of`.
    If `Oc` differs from `Od` and `Of`, Steyn reads history or
    maturation; if `Od` differs from `Oc` and `Of`, the pretest (a
    testing effect). Steyn (2009) presents this reading of the three
    sets of scores (`Oc`, `Od`, and `Of`) as new.

3.  **Testing effect.** The pretest main effect in the two-way
    between-groups ANOVA of the four posttest groups of each
    intervention and the control.

4.  **Pretest-intervention interaction.** The decision sequence of
    Walton Braver and Braver (1988),
    [`fit_solomon_classic()`](https://juhalt.github.io/solomonR/reference/fit_solomon_classic.md),
    for each intervention against the control.

5.  **Test-retest reliability and instrumentation.** The correlation of
    `Oc` and `Od`; the paired t test of `Oc` and `Od`; and a t test of
    `Oc` and `Of`.

6.  **Regression to the mean.** A chi-square test for the variance of
    `Od` against that of `Oc`. Steyn reads an increase in variance as a
    sign of regression to the mean, a threat when groups were selected
    for extreme scores.

7.  **Attrition.** A z test of the dropout proportions of the
    intervention groups against the non-intervention groups (CG1 and
    CG3), and Steyn's chi-square on the two-by-two table of dropouts,
    intervention or not by pretested or not (df = 1). That chi-square
    asks whether intervention and pretesting are associated among the
    dropouts: whether the dropouts of the intervention groups were
    pretested more, or less, often than the dropouts of the
    non-intervention groups. It counts dropouts only, so it does not
    compare dropout rates, and it depends on how many participants each
    group started with. To compare dropout rates by group, see
    `attrition$counts` and
    [`check_solomon_missing()`](https://juhalt.github.io/solomonR/reference/check_solomon_missing.md),
    which counts the missing posttests in each group.

8.  **Effects of the interventions.** For one intervention: E1, a
    one-way ANOVA of the four posttest groups, and E2, a t test of `Ob`
    plus `Oe` against `Od` plus `Of`. For several:

    - E1, a one-way ANOVA of all 2(k + 1) posttest groups;

    - E2, post hoc tests of every pair of groups, summarized by whether
      each intervention group differs from `Od` and from `Of`;

    - E3, a one-way ANOVA of the 2k intervention groups;

    - E4, for each intervention, a t test of `Ob` against `Oe`;

    - E5, when no E4 test is significant, a one-way ANOVA of the k
      interventions with their two groups combined, with post hoc tests,
      which Steyn uses to find which intervention had the greatest
      effect. When an E4 test is significant, Steyn cautions that
      internal validity is in question and the groups are not combined.

Every test is computed when the data allow it; `path` records the steps
Steyn's sequence acts on at `alpha`. Steps 1, 2, 4, 5, and 6 use the
pretest and are skipped, with a note, when `y_pre` is `NULL`.

## Operationalization

Steyn's description leaves these choices open. solomonR made them as
follows.

- **t tests.** Independent-samples t tests are Student's pooled-variance
  tests, matching the equal-variance one-way ANOVAs.

- **Post hoc tests.** Steyn (2009) names no post hoc test. In the study
  the article describes, Steyn (2005, Table 5.61, p. 153) used Scheffé
  tests after the one-way ANOVA of the eight posttest groups, so they
  are the default (`posthoc = "scheffe"`). `posthoc = "holm"` gives
  pairwise t tests with a pooled SD and Holm's (1979) adjustment, the
  default of
  [`stats::pairwise.t.test()`](https://rdrr.io/r/stats/pairwise.t.test.html),
  which is less conservative for pairwise comparisons. Both use the
  pooled error variance of the groups in the ANOVA. For a pair of groups
  among G groups with N scores in all, the Scheffé statistic is the
  squared pairwise t divided by G - 1, on G - 1 and N - G degrees of
  freedom. In step 1 the post hoc tests follow a significant ANOVA of
  more than two groups; in step 2, and in step 8 with several
  interventions (E2 and E5), they are always computed.

- **Two-way ANOVA.** Type III sums of squares, with effect coding, one
  analysis for each intervention against the control.

- **History classification.** Two sets "differ" when their post hoc
  p-value (Scheffé's, or Holm-adjusted) is below `alpha`. A
  nonsignificant ANOVA gives "no evidence of history, maturation, or a
  testing effect"; `Oc` differing from `Od` and from `Of`, which do not
  differ, gives "history or maturation"; `Od` differing from `Oc` and
  from `Of`, which do not differ, gives "the pretest (a testing
  effect)"; any other pattern is reported as one Steyn's rule does not
  cover.

- **Scores used.** The one-way ANOVA of `Oc`, `Od`, and `Of` and the t
  test of `Oc` and `Of` use every available score in each set. The
  paired t test, the correlation, and the variance test use the
  pretested controls with both scores.

- **Variance test.** The statistic is (n - 1) times the variance of `Od`
  divided by the variance of `Oc`, on n - 1 degrees of freedom, with a
  two-sided p-value (twice the smaller tail).

- **Attrition.** A participant with a missing posttest is a dropout.
  "Intervention" means any intervention. The z test pools the two
  proportions; neither test has a continuity correction. The chi-square
  is computed on the dropout counts, as Steyn describes; a note is
  recorded when an expected count is below 5.

- **"Differs from both control groups"** (E2). An intervention group
  differs from both when its post hoc p-values against `Od` and against
  `Of` are both below `alpha`. The sequence continues to E3 when at
  least one intervention group does.

- **"Greatest effect"** (E5). Steyn does not say how the intervention
  with the greatest effect is identified. solomonR reports the
  intervention with the highest combined posttest mean (`highest`) and
  the post hoc tests of every pair. The highest mean is the greatest
  effect only when the interventions raise the scores; when they lower
  them, as in Steyn (2005), the greatest effect is the lowest mean. Read
  `highest` with the direction of the outcome and the post hoc tests.

- **Missing scores.** Participants with a missing posttest stay in the
  data for the attrition step and are left out of the posttest analyses.
  The pretest analyses use the pretested participants with a pretest.

- **Tests that cannot be computed.** A test of groups that are empty,
  too small, or without variation is left out, with a note in `notes`.

- **Step 4** uses `fit_solomon_classic(flow = "1988")` with its defaults
  (ANCOVA as the pretested-groups test, and Test I).

## Cautions

- Steyn notes that multiple ANOVAs capitalize on chance. The sequence
  runs many tests, each at `alpha`, with no control of the error rate
  across them.

- A nonsignificant test is not evidence that groups are equivalent, or
  that a threat is absent.

- The one-way ANOVA of `Oc`, `Od`, and `Of` treats the paired scores of
  `Oc` and `Od`, which come from the same participants, as independent.

- The chi-square test for the variance treats the variance of `Oc` as a
  known value and ignores the pairing of `Oc` and `Od`.

- A difference between `Ob` and `Oe` (E4) mixes a pretest main effect
  with pretest sensitization. Step 4 tests the interaction for each
  intervention separately; the joint model
  ([`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md))
  tests it for all the interventions in one model, with adjusted
  p-values for the comparisons.

- The package's recommended analysis is `fit_solomon_glm(control = )`.

## Lifecycle

Experimental. This function follows a pre-publication draft of Steyn's
(2009) article, dated March 2, 2009, and will be checked against the
published version.

## References

Holm, S. (1979). A simple sequentially rejective multiple test
procedure. *Scandinavian Journal of Statistics, 6*(2), 65–70.
https://www.jstor.org/stable/4615733

Steyn, R. (2005). *Self-evaluasie en die vorming van
selfdoeltreffendheidspersepsies* \[Self-evaluation and the forming of
self-efficacy perceptions\] \[Doctoral thesis, University of South
Africa\]. Unisa Institutional Repository.
https://hdl.handle.net/10500/1745

Steyn, R. (2009). Re-designing the Solomon four-group: Can we improve on
this exemplary model? *Design Principles and Practices: An International
Journal, 3*(1), 383–394.
https://doi.org/10.18848/1833-1874/CGP/v03i01/37588

Walton Braver, M. C., & Braver, S. L. (1988). Statistical treatment of
the Solomon four-group design: A meta-analytic approach. *Psychological
Bulletin, 104*(1), 150–154. https://doi.org/10.1037/0033-2909.104.1.150

## See also

[`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md)
for the recommended analysis,
[`fit_solomon_classic()`](https://juhalt.github.io/solomonR/reference/fit_solomon_classic.md),
[`check_solomon_missing()`](https://juhalt.github.io/solomonR/reference/check_solomon_missing.md),
[`report_solomon()`](https://juhalt.github.io/solomonR/reference/report_solomon.md)
for APA 7 text of the results, and
[steyn2005](https://juhalt.github.io/solomonR/reference/steyn2005.md)
for the summary statistics of Steyn's (2005) eight-group study.

## Examples

``` r
# Six groups: two interventions and a control.
fit <- fit_solomon_steyn(post_behavior, condition, pretested, pre_behavior,
                         control = "Control", data = mai2020)
fit
#> Steyn's (2009) analysis of the extended Solomon design (a published proposal; the package's recommended analysis is fit_solomon_glm())
#> Follows a pre-publication draft of the article.
#> Interventions: RP, GS; control: Control. With and without a pretest: 6 groups. alpha = 0.05.
#> 
#> Groups
#> Group  Condition  Pretested  Pretest  Posttest   n  Posttests
#> EG1    RP         yes        Oa1      Ob1       35         24
#> EG2    GS         yes        Oa2      Ob2       33         23
#> CG1    Control    yes        Oc       Od        50         27
#> CG2.1  RP         no         -        Oe1       31         22
#> CG2.2  GS         no         -        Oe2       27         15
#> CG3    Control    no         -        Of        35         22
#> 
#> 1. Equivalence after randomization
#> Groups        Test           Estimate         Statistic      p
#> Oa1, Oa2, Oc  One-way ANOVA            F(2, 115) = 0.78  0.459
#>   No significant difference between the pretested groups at pretest.
#> 
#> 2. History and maturation
#> Groups                Test                                   Estimate        Statistic      p  p adj.
#> Od vs Oc, paired      Paired t test                            -0.010    t(26) = -0.13  0.901        
#> Of vs Oa1 + Oa2 + Oc  t test (pooled variance)                 -0.142   t(138) = -1.75  0.083        
#> Oc, Od, Of            One-way ANOVA (scores as independent)            F(2, 96) = 0.78  0.462        
#> Oc vs Od              Pairwise (pooled SD, Scheffe)             0.016     t(96) = 0.19  0.851   0.982
#> Oc vs Of              Pairwise (pooled SD, Scheffe)             0.109     t(96) = 1.23  0.223   0.474
#> Od vs Of              Pairwise (pooled SD, Scheffe)             0.094     t(96) = 0.94  0.351   0.646
#>   Oc and Od do not differ (paired t test). The pretests and Of do not differ.
#>   The one-way ANOVA of Oc, Od, and Of finds no evidence of history,
#>   maturation, or a testing effect.
#> 
#> 3. Testing effect: two-way between-groups ANOVA (Type III)
#> Comparison     Term             SS        Statistic      p
#> RP vs Control  Pretest       0.066  F(1, 91) = 0.45  0.504
#> RP vs Control  Intervention  0.031  F(1, 91) = 0.21  0.645
#> RP vs Control  Interaction   0.508  F(1, 91) = 3.46  0.066
#> GS vs Control  Pretest       0.063  F(1, 83) = 0.48  0.492
#> GS vs Control  Intervention  0.186  F(1, 83) = 1.42  0.237
#> GS vs Control  Interaction   0.031  F(1, 83) = 0.24  0.626
#>   RP vs Control: no pretest main effect. GS vs Control: no pretest main
#>   effect.
#> 
#> 4. Pretest-intervention interaction: Tests A-I (Walton Braver & Braver, 1988)
#> Comparison              Test A      p                   Path
#> RP vs Control  F(1, 91) = 3.46  0.066  A -> D -> E -> H -> I
#> GS vs Control  F(1, 83) = 0.24  0.626  A -> D -> E -> H -> I
#>   RP vs Control (path A -> D -> E -> H -> I): Historical pathway: no
#>   treatment test in the selected A-I sequence reaches the specified alpha
#>   level. GS vs Control (path A -> D -> E -> H -> I): Historical pathway: no
#>   treatment test in the selected A-I sequence reaches the specified alpha
#>   level.
#> 
#> 5. Test-retest reliability and instrumentation
#> Step             Groups            Test                      Estimate      Statistic      p
#> Reliability      Oc with Od        Pearson r (test-retest)      0.247   t(25) = 1.28  0.213
#> Instrumentation  Od vs Oc, paired  Paired t test               -0.010  t(26) = -0.13  0.901
#> Instrumentation  Of vs Oc          t test (pooled variance)    -0.109  t(70) = -1.24  0.221
#>   Test-retest r = 0.25 (n = 27). Oc and Od do not differ. Oc and Of do not
#>   differ.
#> 
#> 6. Regression to the mean: chi-square test for the variance
#> Groups                     Var(Oc)  Var(Od)  Ratio         Statistic      p
#> Od vs Oc, paired controls    0.095    0.126   1.33  chi2(26) = 34.59  0.242
#>   The variance changed from 0.09496 (Oc) to 0.1263 (Od), a ratio of 1.33; the
#>   change is not significant.
#> 
#> 7. Attrition
#> Group  Condition  Randomized  Observed  Missing   Rate
#> EG1    RP                 35        24       11  31.4%
#> EG2    GS                 33        23       10  30.3%
#> CG1    Control            50        27       23  46.0%
#> CG2.1  RP                 31        22        9  29.0%
#> CG2.2  GS                 27        15       12  44.4%
#> CG3    Control            35        22       13  37.1%
#> 
#> Groups                                  Test                         Estimate       Statistic      p
#> EG1 + EG2 + CG2.1 + CG2.2 vs CG1 + CG3  Two-proportion z test          -0.090       z = -1.33  0.183
#> Dropouts: intervention by pretested     Chi-square (2 x 2 dropouts)            chi2(1) = 1.52  0.218
#>   78 of 211 participants (37.0%) have no posttest. Dropout was 33.3% with an
#>   intervention and 42.4% without; the difference is not significant. Steyn's
#>   chi-square: intervention and pretesting are not significantly associated
#>   among the dropouts.
#> 
#> 8. Effects of the interventions
#> Step  Groups                      Test                            Estimate         Statistic      p  p adj.
#> E1    Ob1, Ob2, Od, Oe1, Oe2, Of  One-way ANOVA                             F(5, 127) = 1.34  0.251        
#> E3    Ob1, Ob2, Oe1, Oe2          One-way ANOVA                              F(3, 80) = 1.83  0.148        
#> E4    Ob1 vs Oe1                  t test (pooled variance)          -0.200     t(44) = -1.65  0.107        
#> E4    Ob2 vs Oe2                  t test (pooled variance)           0.016      t(36) = 0.13  0.898        
#> E5    Ob1 + Oe1, Ob2 + Oe2        One-way ANOVA, groups combined             F(1, 82) = 2.50  0.118        
#> E5    Ob1 + Oe1 vs Ob2 + Oe2      Pairwise (pooled SD, Scheffe)     -0.137     t(82) = -1.58  0.118   0.118
#> Steps not on the decision path are shown for completeness.
#> 
#> E2: intervention groups against Od and Of (Scheffé tests)
#> Group  Condition   Mean  p vs Od  p vs Of  Differs from both
#> Ob1    RP         2.929    0.707    0.986                 no
#> Ob2    GS         3.168    0.998    0.882                 no
#> Oe1    RP         3.129    1.000    0.968                 no
#> Oe2    GS         3.152    1.000    0.953                 no
#> 
#> Decision path: E1
#>   E1: the posttest groups do not differ, so Steyn's sequence finds no
#>   evidence that the interventions had an effect.

# Pairwise t tests with Holm's adjustment as the post hoc tests.
fit_solomon_steyn(post_behavior, condition, pretested, pre_behavior,
                  control = "Control", posthoc = "holm", data = mai2020)
#> Steyn's (2009) analysis of the extended Solomon design (a published proposal; the package's recommended analysis is fit_solomon_glm())
#> Follows a pre-publication draft of the article.
#> Interventions: RP, GS; control: Control. With and without a pretest: 6 groups. alpha = 0.05.
#> 
#> Groups
#> Group  Condition  Pretested  Pretest  Posttest   n  Posttests
#> EG1    RP         yes        Oa1      Ob1       35         24
#> EG2    GS         yes        Oa2      Ob2       33         23
#> CG1    Control    yes        Oc       Od        50         27
#> CG2.1  RP         no         -        Oe1       31         22
#> CG2.2  GS         no         -        Oe2       27         15
#> CG3    Control    no         -        Of        35         22
#> 
#> 1. Equivalence after randomization
#> Groups        Test           Estimate         Statistic      p
#> Oa1, Oa2, Oc  One-way ANOVA            F(2, 115) = 0.78  0.459
#>   No significant difference between the pretested groups at pretest.
#> 
#> 2. History and maturation
#> Groups                Test                                   Estimate        Statistic      p  p adj.
#> Od vs Oc, paired      Paired t test                            -0.010    t(26) = -0.13  0.901        
#> Of vs Oa1 + Oa2 + Oc  t test (pooled variance)                 -0.142   t(138) = -1.75  0.083        
#> Oc, Od, Of            One-way ANOVA (scores as independent)            F(2, 96) = 0.78  0.462        
#> Oc vs Od              Pairwise t (pooled SD, Holm)              0.016     t(96) = 0.19  0.851   0.851
#> Oc vs Of              Pairwise t (pooled SD, Holm)              0.109     t(96) = 1.23  0.223   0.668
#> Od vs Of              Pairwise t (pooled SD, Holm)              0.094     t(96) = 0.94  0.351   0.703
#>   Oc and Od do not differ (paired t test). The pretests and Of do not differ.
#>   The one-way ANOVA of Oc, Od, and Of finds no evidence of history,
#>   maturation, or a testing effect.
#> 
#> 3. Testing effect: two-way between-groups ANOVA (Type III)
#> Comparison     Term             SS        Statistic      p
#> RP vs Control  Pretest       0.066  F(1, 91) = 0.45  0.504
#> RP vs Control  Intervention  0.031  F(1, 91) = 0.21  0.645
#> RP vs Control  Interaction   0.508  F(1, 91) = 3.46  0.066
#> GS vs Control  Pretest       0.063  F(1, 83) = 0.48  0.492
#> GS vs Control  Intervention  0.186  F(1, 83) = 1.42  0.237
#> GS vs Control  Interaction   0.031  F(1, 83) = 0.24  0.626
#>   RP vs Control: no pretest main effect. GS vs Control: no pretest main
#>   effect.
#> 
#> 4. Pretest-intervention interaction: Tests A-I (Walton Braver & Braver, 1988)
#> Comparison              Test A      p                   Path
#> RP vs Control  F(1, 91) = 3.46  0.066  A -> D -> E -> H -> I
#> GS vs Control  F(1, 83) = 0.24  0.626  A -> D -> E -> H -> I
#>   RP vs Control (path A -> D -> E -> H -> I): Historical pathway: no
#>   treatment test in the selected A-I sequence reaches the specified alpha
#>   level. GS vs Control (path A -> D -> E -> H -> I): Historical pathway: no
#>   treatment test in the selected A-I sequence reaches the specified alpha
#>   level.
#> 
#> 5. Test-retest reliability and instrumentation
#> Step             Groups            Test                      Estimate      Statistic      p
#> Reliability      Oc with Od        Pearson r (test-retest)      0.247   t(25) = 1.28  0.213
#> Instrumentation  Od vs Oc, paired  Paired t test               -0.010  t(26) = -0.13  0.901
#> Instrumentation  Of vs Oc          t test (pooled variance)    -0.109  t(70) = -1.24  0.221
#>   Test-retest r = 0.25 (n = 27). Oc and Od do not differ. Oc and Of do not
#>   differ.
#> 
#> 6. Regression to the mean: chi-square test for the variance
#> Groups                     Var(Oc)  Var(Od)  Ratio         Statistic      p
#> Od vs Oc, paired controls    0.095    0.126   1.33  chi2(26) = 34.59  0.242
#>   The variance changed from 0.09496 (Oc) to 0.1263 (Od), a ratio of 1.33; the
#>   change is not significant.
#> 
#> 7. Attrition
#> Group  Condition  Randomized  Observed  Missing   Rate
#> EG1    RP                 35        24       11  31.4%
#> EG2    GS                 33        23       10  30.3%
#> CG1    Control            50        27       23  46.0%
#> CG2.1  RP                 31        22        9  29.0%
#> CG2.2  GS                 27        15       12  44.4%
#> CG3    Control            35        22       13  37.1%
#> 
#> Groups                                  Test                         Estimate       Statistic      p
#> EG1 + EG2 + CG2.1 + CG2.2 vs CG1 + CG3  Two-proportion z test          -0.090       z = -1.33  0.183
#> Dropouts: intervention by pretested     Chi-square (2 x 2 dropouts)            chi2(1) = 1.52  0.218
#>   78 of 211 participants (37.0%) have no posttest. Dropout was 33.3% with an
#>   intervention and 42.4% without; the difference is not significant. Steyn's
#>   chi-square: intervention and pretesting are not significantly associated
#>   among the dropouts.
#> 
#> 8. Effects of the interventions
#> Step  Groups                      Test                            Estimate         Statistic      p  p adj.
#> E1    Ob1, Ob2, Od, Oe1, Oe2, Of  One-way ANOVA                             F(5, 127) = 1.34  0.251        
#> E3    Ob1, Ob2, Oe1, Oe2          One-way ANOVA                              F(3, 80) = 1.83  0.148        
#> E4    Ob1 vs Oe1                  t test (pooled variance)          -0.200     t(44) = -1.65  0.107        
#> E4    Ob2 vs Oe2                  t test (pooled variance)           0.016      t(36) = 0.13  0.898        
#> E5    Ob1 + Oe1, Ob2 + Oe2        One-way ANOVA, groups combined             F(1, 82) = 2.50  0.118        
#> E5    Ob1 + Oe1 vs Ob2 + Oe2      Pairwise t (pooled SD, Holm)      -0.137     t(82) = -1.58  0.118   0.118
#> Steps not on the decision path are shown for completeness.
#> 
#> E2: intervention groups against Od and Of (Holm-adjusted pairwise t tests)
#> Group  Condition   Mean  p vs Od  p vs Of  Differs from both
#> Ob1    RP         2.929    1.000    1.000                 no
#> Ob2    GS         3.168    1.000    1.000                 no
#> Oe1    RP         3.129    1.000    1.000                 no
#> Oe2    GS         3.152    1.000    1.000                 no
#> 
#> Decision path: E1
#>   E1: the posttest groups do not differ, so Steyn's sequence finds no
#>   evidence that the interventions had an effect.

# The four-group design: one intervention.
fit_solomon_steyn(y_post, treat, pretested, y_pre, data = solomon_example)
#> Steyn's (2009) analysis of the extended Solomon design (a published proposal; the package's recommended analysis is fit_solomon_glm())
#> Follows a pre-publication draft of the article.
#> Intervention: Treatment; control: Control. With and without a pretest: 4 groups. alpha = 0.05.
#> 
#> Groups
#> Group  Condition  Pretested  Pretest  Posttest   n  Posttests
#> EG     Treatment  yes        Oa       Ob        30         30
#> CG1    Control    yes        Oc       Od        30         30
#> CG2    Treatment  no         -        Oe        30         30
#> CG3    Control    no         -        Of        30         30
#> 
#> 1. Equivalence after randomization
#> Groups    Test                      Estimate     Statistic      p
#> Oa vs Oc  t test (pooled variance)     0.067  t(58) = 0.02  0.981
#>   No significant difference between the pretested groups at pretest.
#> 
#> 2. History and maturation
#> Groups            Test                                   Estimate        Statistic      p  p adj.
#> Od vs Oc, paired  Paired t test                             4.933     t(29) = 2.72  0.011        
#> Of vs Oa + Oc     t test (pooled variance)                  1.500     t(88) = 0.65  0.516        
#> Oc, Od, Of        One-way ANOVA (scores as independent)            F(2, 87) = 2.14  0.124        
#> Oc vs Od          Pairwise (pooled SD, Scheffe)            -4.933    t(87) = -2.02  0.046   0.136
#> Oc vs Of          Pairwise (pooled SD, Scheffe)            -1.533    t(87) = -0.63  0.531   0.821
#> Od vs Of          Pairwise (pooled SD, Scheffe)             3.400     t(87) = 1.39  0.167   0.383
#>   Oc and Od differ (paired t test). The pretests and Of do not differ. The
#>   one-way ANOVA of Oc, Od, and Of finds no evidence of history, maturation,
#>   or a testing effect.
#> 
#> 3. Testing effect: two-way between-groups ANOVA (Type III)
#> Comparison            Term               SS         Statistic      p
#> Treatment vs Control  Pretest       180.075  F(1, 116) = 1.95  0.165
#> Treatment vs Control  Intervention  216.008  F(1, 116) = 2.34  0.129
#> Treatment vs Control  Interaction    27.075  F(1, 116) = 0.29  0.589
#>   Treatment vs Control: no pretest main effect.
#> 
#> 4. Pretest-intervention interaction: Tests A-I (Walton Braver & Braver, 1988)
#> Comparison                      Test A      p                   Path
#> Treatment vs Control  F(1, 116) = 0.29  0.589  A -> D -> E -> H -> I
#>   Treatment vs Control (path A -> D -> E -> H -> I): Historical pathway: no
#>   treatment test in the selected A-I sequence reaches the specified alpha
#>   level.
#> 
#> 5. Test-retest reliability and instrumentation
#> Step             Groups            Test                      Estimate     Statistic      p
#> Reliability      Oc with Od        Pearson r (test-retest)      0.456  t(28) = 2.71  0.011
#> Instrumentation  Od vs Oc, paired  Paired t test                4.933  t(29) = 2.72  0.011
#> Instrumentation  Of vs Oc          t test (pooled variance)     1.533  t(58) = 0.61  0.543
#>   Test-retest r = 0.46 (n = 30). Oc and Od differ, a possible instrumentation
#>   effect. Oc and Of do not differ.
#> 
#> 6. Regression to the mean: chi-square test for the variance
#> Groups                     Var(Oc)  Var(Od)  Ratio         Statistic      p
#> Od vs Oc, paired controls  100.530   79.569   0.79  chi2(29) = 22.95  0.443
#>   The variance changed from 100.5 (Oc) to 79.57 (Od), a ratio of 0.79; the
#>   change is not significant.
#> 
#> 7. Attrition
#> Group  Condition  Randomized  Observed  Missing  Rate
#> EG     Treatment          30        30        0  0.0%
#> CG1    Control            30        30        0  0.0%
#> CG2    Treatment          30        30        0  0.0%
#> CG3    Control            30        30        0  0.0%
#>   No dropouts: every posttest was observed.
#> 
#> 8. Effect of the intervention
#> Step  Groups              Test                      Estimate         Statistic      p
#> E1    Ob, Od, Oe, Of      One-way ANOVA                       F(3, 116) = 1.53  0.211
#> E2    Ob + Oe vs Od + Of  t test (pooled variance)     2.683     t(118) = 1.53  0.129
#> 
#> Decision path: E1 -> E2
#>   E1: the four posttest groups do not differ. E2: the intervention groups
#>   (Ob, Oe) do not differ from the non-intervention groups (Od, Of).
#> 
#> Notes
#> - Every posttest was observed, so the attrition tests were skipped.
```
