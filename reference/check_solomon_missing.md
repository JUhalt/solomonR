# Distinguish structural and incidental missingness in a Solomon design

**\[stable\]** Classifies missing values in Solomon design data and
explains the supported response to each kind. Pretest scores are
*structurally absent* for participants assigned to the unpretested
groups: withholding the pretest is the experimental manipulation
(Solomon, 1949), so those values must never be imputed. Other missing
values are *incidental* and are handled according to the missing-data
literature (Rubin, 1976; Little & Rubin, 2019).

## Usage

``` r
check_solomon_missing(
  y_post,
  treat,
  pretested,
  y_pre = NULL,
  control = NULL,
  data = NULL,
  covariates = NULL
)
```

## Arguments

- y_post:

  Numeric posttest scores.

- treat:

  Treatment indicator coded 0/1 (or logical); or a factor or character
  vector of conditions, with the control named by `control`.

- pretested:

  Pretest indicator coded 0/1 (or logical).

- y_pre:

  Optional numeric pretest scores, missing by design for unpretested
  participants.

- control:

  The control condition when `treat` is a factor or character vector
  with more than two conditions, a Solomon N-group design; see
  [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md).
  With two conditions the result is the same as with a 0/1 `treat`.

- data:

  Optional data frame. When supplied, the other data arguments are
  looked up in it first, as bare column names (`y_post = post`) or as
  strings (`y_post = "post"`).

- covariates:

  Optional covariates recorded in all groups, such as gender: a data
  frame, or, with `data`, the names of its columns. They enter the fit
  behind the gap and the delta line, and loss is shown by each of them
  (see "Missing posttests and the pretest effects").

## Value

An object of class `solomon_missing`, a list with:

- `by_cell`: the counts in each Solomon cell, in the columns `group`,
  `cell`, `pretested`, `treat`, `n`, `posttest_missing`,
  `pretest_structural`, `pretest_incidental`, and `pretest_unexpected`.

- `counts`: the totals by category.

- `pattern`: `"none"`, `"structural"`, `"incidental"`, or `"mixed"`.

- `guidance`: for each category present, its count, the interpretation,
  the supported response, and the sources (`category`, `n`,
  `interpretation`, `response`, and `sources`).

- `pretest_supplied`: whether pretest scores were supplied.

The checks of missing posttests for the pretest effects (see "Missing
posttests and the pretest effects"), the elements that need a pretest
being `NULL` without a numeric `y_pre`:

- `attrition`: `groups`, each group's `condition`, `pretested`, `n`,
  `posttest_missing`, and `share`; `overall`, the share missing;
  `by_pretest`, the shares of the pretested and the unpretested;
  `pretest_difference`, the first minus the second; `by_condition`, the
  share missing in each condition; `condition_difference`, the largest
  difference between conditions; `p.value`, Fisher's exact test of loss
  by pretest condition; and `differ`, whether it is below .05.

- `pretest_scores`: for each pretested condition, `n` participants with
  a pretest and their mean pretest `mean_all`; `n_with` and `mean_with`
  for those with a posttest; `n_without` and `mean_without` for those
  without; `difference` (without minus with); `shift` (with minus all);
  and `shift_se`, its standard error.

- `selection_test`: the Wald test of equal shifts, its `statistic`
  (chi-square), `df`, and `p.value`; and `selection_differ`, whether the
  p-value is below .05.

- `gap`: `slope`, `m1`, the mean pretest of the pretested participants
  `analyzed` and `not_analyzed`, the `gap` between the centers, and
  `predicted`, what the pretested dropouts' pretests predict for them.

- `unpretested`: for each unpretested group, `condition`, `n`,
  `posttest_missing`, `m0`, and `delta`, the offset the default assumes.

- `covariates`: one table for each covariate: by pretest condition and
  level, `n`, `posttest_missing`, `share`, `analyzed_share`, `p.value`,
  and `test` (`"exact"` or `"simulated"`); or, for a numeric covariate
  with more than two values, its mean among those with (`mean_with`) and
  without (`mean_without`) a posttest. `NULL` without `covariates`.

- `missing_pretest`: for each pretested condition, the participants with
  a posttest and no pretest (`n_without_pretest`,
  `mean_without_pretest`, their mean posttest) and those with both
  scores (`n_with_pretest`, `mean_with_pretest`).

- `raw_pretest_effects`: the pretest effects of a fit without the
  pretest, in the columns of an effects table, when the selection
  differs or a pretested participant lacks a pretest; otherwise `NULL`.

- `unverifiable`: whether any unpretested posttest is missing.

For a design with several treatments, `by_cell` has one row for each of
the 2(k + 1) groups, its `treat` column holds the condition, and
`conditions` names the control and the treatments.

## Details

Categories:

- **Structural pretest absence**: unpretested participants have no
  pretest. Unlike planned missing-data designs, in which unmeasured
  values exist and can be imputed (Graham et al., 2006), an imputed
  pretest here would describe a measurement that never occurred.

- **Incidental pretest missingness**: pretested participants without a
  pretest score. solomonR analyses currently use complete cases and
  warn. Because treatment is randomized, deterministic mean imputation
  of the pretest (without using treatment or outcome), with a
  missingness indicator when pretests may not be missing completely at
  random, retains these participants without biasing the treatment
  effect (White & Thompson, 2005; Groenwold et al., 2012). This concerns
  the treatment effect only: every fit with the pretest as a covariate
  drops these participants from the pretest effects, at either pretest
  center, and when their loss relates to baseline those effects are
  biased while a fit without the pretest is not (issue \#142).

- **Incidental posttest missingness**: missing outcomes. For the
  treatment contrasts, complete-case analysis is unbiased when
  missingness is unrelated to the outcome given the variables in the
  model (Little & Rubin, 2019); attrition that differs across the groups
  should be reported. The pretest effects also depend on the mean
  pretest at which the pretested groups are compared; see "Missing
  posttests and the pretest effects" below.

- **Unexpected pretest scores**: pretest values recorded for unpretested
  participants, which usually indicate a coding or assignment error.
  They are ignored by solomonR analyses.

- **Unassigned participants**: missing treatment or pretest assignment.

Pretest categories are `NA` when `y_pre` is not supplied.

Designs with several treatments: give `treat` as a factor or character
vector of conditions and name the control with `control`. The counts are
then given for each of the 2(k + 1) groups of a design with k treatments
(Steyn, 2009): each treatment and the control, with and without a
pretest.

Not supported: this function does not test the missingness mechanism,
perform imputation, or provide sensitivity analyses for outcome
missingness that depends on unobserved values.

## Missing posttests and the pretest effects

With a pretest covariate, the pretest effects compare the pretested
groups, at a mean pretest, with the unpretested groups, and when
posttests are missing that mean can be taken over the pretested
participants analyzed (the default of
[`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md))
or over every observed pretest (`pretest_center = "observed"`); see "The
pretest effect when posttests are missing" in
[`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md).
The data cannot tell which is right. The result gives the checks of the
decision record on issue \#135 (section 5), and the print method
summarizes them:

1.  **Posttests missing by group**: the overall rate, and the
    differential attrition between the pretest conditions and between
    the conditions (What Works Clearinghouse, 2022, p. 40). Equal rates
    do not show that the same kinds of participants were lost (Campbell
    & Stanley, 1963/1966, pp. 15–16). When Fisher's exact test of loss
    by pretest condition gives p \< .05, the print says that the loss
    rates differ, which suggests that taking the pretest changed
    attrition, so that neither center is assured. That rule is the
    package's: no source fixes one for the pretest conditions.

2.  **Pretest scores of those lost**: for each pretested group, the mean
    pretest of the participants with and without a posttest (Campbell,
    1957, pp. 301–302; Campbell & Stanley, 1963/1966, pp. 15–16).

3.  **Selection by condition**: each pretested condition's selection
    shift, the mean pretest of those with a posttest minus that of every
    participant with a pretest, with its standard error from the shift's
    influence function (Stefanski & Boos, 2002), and a Wald chi-square
    test of equal shifts across the conditions. A condition that lost no
    posttest has a shift of exactly 0. When the test gives p \< .05 (the
    package's rule; no source fixes one), the conditions may have lost
    different kinds of participants. Neither center then gives both
    conditional pretest effects: in the decision record's case S3
    (losses by baseline in the control arms and at random in the treated
    arms, alike in both pretest conditions) the default misstated both,
    although not their average when the pretested groups kept equal
    numbers, and the observed center misstated one of them and their
    average. Pretest x Treatment is biased at every center (issue
    \#143). The pretest effects of a fit without the pretest, the raw
    comparison of posttest means in each condition, are then given. The
    test does not depend on the units of the pretest. When every group
    loses by the same rule it rejects at about its nominal rate: in a
    check on issue \#135 (1,000 data sets for each setting; part 8 of
    `tools/pretest-center-check.R` in the package's source), in 5.7% of
    data sets with 200 participants per group and 4.4% with 50 when loss
    depended on baseline, and in 5.8% and 4.5% when it was completely at
    random.

4.  **The gap between the centers**: the pretest effects at the analyzed
    center minus those at the observed center, b m1 (completers' minus
    dropouts' mean pretest), with b the pretest slope and m1 the share
    of pretested participants with a pretest who are not analyzed, from
    a fit on the identity link with `covariates`. This is the adjustment
    of What Works Clearinghouse (2022, Equation H.2, p. 207).

5.  **The unpretested groups**: for each, its share without a posttest,
    m0, and the offset the default assumes for its dropouts, the delta
    line of
    [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md).
    A group that lost no posttests is named: if the data hold everyone
    randomized to it, the default's pretest effect in that condition is
    right only if the pretested groups' losses, taken together, were
    unrelated to baseline, because the default's center pools the
    pretested groups.

6.  **Covariates recorded in all groups** (`covariates`): loss by
    covariate within each pretest condition, and the covariate's
    distribution among the participants with a posttest (Lee, 2005,
    Remark 4, p. 16; Groenwold et al., 2014). For a covariate with
    levels, loss by level is tested with Fisher's exact test when there
    are two levels, and with a p-value simulated from 10,000 tables with
    the same margins, under a fixed seed, when there are more, because
    the exact test of a larger table can exceed its workspace at
    ordinary sample sizes; that choice is the package's, and no source
    fixes one. Dukes et al. (1995, p. 413) compared the course grades of
    the students lost in this way. A covariate that predicts the
    posttest belongs in the model; that it weakens both centers'
    assumptions is the record's inference, not a checked result.

7.  **Pretested participants with a posttest and no pretest**: their
    number and mean posttest beside those of the pretested participants
    with both scores, and the raw comparison, which keeps them (issue
    \#142).

8.  **What cannot be checked**: the unpretested groups' losses cannot be
    checked against baseline (Campbell, 1957, p. 305). When only those
    groups lost posttests, the two centers coincide, but neither is
    protected; only the delta line reaches that case. If an unpretested
    group's dropouts would have scored Delta above its completers, the
    pretest effect in its condition is m0 Delta lower, m0 being that
    group's share without a posttest (see "The pretest effect when
    posttests are missing" in
    [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md)).

The checks on the pretest (items 2 to 5 and 7) need numeric pretest
scores and are left out otherwise.
[`validate_solomon()`](https://juhalt.github.io/solomonR/reference/validate_solomon.md)
leaves out all of them.

## References

Campbell, D. T. (1957). Factors relevant to the validity of experiments
in social settings. *Psychological Bulletin, 54*(4), 297–312.
https://doi.org/10.1037/h0040950

Campbell, D. T., & Stanley, J. C. (1966). *Experimental and
quasi-experimental designs for research*. Rand McNally. (Original work
published 1963)

Dukes, R. L., Ullman, J. B., & Stein, J. A. (1995). An evaluation of
D.A.R.E. (Drug Abuse Resistance Education), using a Solomon four-group
design with latent variables. *Evaluation Review, 19*(4), 409–435.
https://doi.org/10.1177/0193841X9501900404

Graham, J. W., Taylor, B. J., Olchowski, A. E., & Cumsille, P. E.
(2006). Planned missing data designs in psychological research.
*Psychological Methods, 11*(4), 323–343.
https://doi.org/10.1037/1082-989X.11.4.323

Groenwold, R. H. H., Moons, K. G. M., & Vandenbroucke, J. P. (2014).
Randomized trials with missing outcome data: How to analyze and what to
report. *Canadian Medical Association Journal, 186*(15), 1153–1157.
https://doi.org/10.1503/cmaj.131353

Groenwold, R. H. H., White, I. R., Donders, A. R. T., Carpenter, J. R.,
Altman, D. G., & Moons, K. G. M. (2012). Missing covariate data in
clinical research: When and when not to use the missing-indicator method
for analysis. *Canadian Medical Association Journal, 184*(11),
1265–1269. https://doi.org/10.1503/cmaj.110977

Lee, D. S. (2005). *Training, wages, and sample selection: Estimating
sharp bounds on treatment effects* (NBER Working Paper No. 11721).
National Bureau of Economic Research. https://doi.org/10.3386/w11721

Little, R. J. A., & Rubin, D. B. (2019). *Statistical analysis with
missing data* (3rd ed.). Wiley. https://doi.org/10.1002/9781119482260

Mai, N. N., Takahashi, Y., & Oo, M. M. (2020). Testing the effectiveness
of transfer interventions using Solomon four-group designs. *Education
Sciences, 10*(4), Article 92. https://doi.org/10.3390/educsci10040092

Rubin, D. B. (1976). Inference and missing data. *Biometrika, 63*(3),
581–592. https://doi.org/10.1093/biomet/63.3.581

Solomon, R. L. (1949). An extension of control group design.
*Psychological Bulletin, 46*(2), 137–150.
https://doi.org/10.1037/h0062958

Stefanski, L. A., & Boos, D. D. (2002). The calculus of M-estimation.
*The American Statistician, 56*(1), 29–38.
https://doi.org/10.1198/000313002753631330

Steyn, R. (2009). Re-designing the Solomon four-group: Can we improve on
this exemplary model? *Design Principles and Practices: An International
Journal—Annual Review, 3*(1), 383–394.
https://doi.org/10.18848/1833-1874/CGP/v03i01/37588

What Works Clearinghouse. (2022). *What Works Clearinghouse procedures
and standards handbook, version 5.0* (WWC 2022008REV). U.S. Department
of Education, Institute of Education Sciences, National Center for
Education Evaluation and Regional Assistance.
https://ies.ed.gov/ncee/wwc/Handbooks

White, I. R., & Thompson, S. G. (2005). Adjusting for partially missing
baseline measurements in randomized trials. *Statistics in Medicine,
24*(7), 993–1007. https://doi.org/10.1002/sim.1981

## See also

[`validate_solomon()`](https://juhalt.github.io/solomonR/reference/validate_solomon.md)

## Examples

``` r
data(solomon_example)

# Structural absence only: pretests are absent by design in Groups 3 and 4.
with(solomon_example, check_solomon_missing(y_post, treat, pretested, y_pre))
#> Solomon missingness check
#> Pattern: structural pretest absence only (expected in a Solomon design)
#> 
#>  Group                   Cell  n Post missing Pre absent (design) Pre missing
#>      1   Pretested, treatment 30            0                   0           0
#>      2     Pretested, control 30            0                   0           0
#>      3 Unpretested, treatment 30            0                  30           0
#>      4   Unpretested, control 30            0                  30           0
#>  Pre unexpected
#>               0
#>               0
#>               0
#>               0
#> 
#> Structural pretest absence (n = 60)
#>   Participants assigned to the unpretested groups were never pretested; the
#>   absence of a pretest is the experimental manipulation.
#>   Response: Do not impute. Use analyses that respect the design, such as
#>     fit_solomon_glm(), fit_solomon_ml(), fit_solomon_classic(), or SEM.
#>     Unlike planned missing-data designs, where unmeasured values exist and
#>     can be imputed, an imputed pretest here would describe a measurement that
#>     never occurred.
#>   Sources: Solomon (1949); Graham et al. (2006)

# Mixed: add a lost pretest and a missing posttest.
d <- solomon_example
d$y_pre[which(d$pretested == 1)[1]] <- NA
d$y_post[which(d$pretested == 0)[1]] <- NA
with(d, check_solomon_missing(y_post, treat, pretested, y_pre))
#> Solomon missingness check
#> Pattern: structural pretest absence and incidental missingness
#> 
#>  Group                   Cell  n Post missing Pre absent (design) Pre missing
#>      1   Pretested, treatment 30            0                   0           1
#>      2     Pretested, control 30            0                   0           0
#>      3 Unpretested, treatment 30            1                  30           0
#>      4   Unpretested, control 30            0                  30           0
#>  Pre unexpected
#>               0
#>               0
#>               0
#>               0
#> 
#> Structural pretest absence (n = 60)
#>   Participants assigned to the unpretested groups were never pretested; the
#>   absence of a pretest is the experimental manipulation.
#>   Response: Do not impute. Use analyses that respect the design, such as
#>     fit_solomon_glm(), fit_solomon_ml(), fit_solomon_classic(), or SEM.
#>     Unlike planned missing-data designs, where unmeasured values exist and
#>     can be imputed, an imputed pretest here would describe a measurement that
#>     never occurred.
#>   Sources: Solomon (1949); Graham et al. (2006)
#> 
#> Incidental pretest missingness (n = 1)
#>   Participants assigned to be pretested have no pretest score. Determine
#>   whether the pretest was administered but the score was lost, or never
#>   administered, which is a departure from the assigned pretest condition.
#>   Response: solomonR analyses currently use complete cases and warn. Because
#>     treatment is randomized, deterministic mean imputation of the pretest
#>     (without using treatment or outcome), with a missingness indicator when
#>     pretests may not be missing completely at random, retains these
#>     participants without biasing the treatment effect. That holds for the
#>     treatment effect only: the pretest effects of a fit with the pretest drop
#>     these participants at either pretest center, and are biased when their
#>     loss relates to baseline (issue #142). If the pretest was never
#>     administered, report the departure and analyze participants as assigned.
#>   Sources: White & Thompson (2005); Groenwold et al. (2012)
#> 
#> Incidental posttest missingness (n = 1)
#>   Posttest scores are missing, so these participants are excluded from
#>   complete-case analyses.
#>   Response: For the treatment contrasts, complete-case analysis is unbiased
#>     when missingness is unrelated to the outcome given the variables in the
#>     model. The pretest effects also depend on the mean pretest at which they
#>     are evaluated; see `pretest_center` in fit_solomon_glm() and the
#>     diagnostics below. Report missingness by cell, because attrition that
#>     differs across the four groups can undermine the randomized comparisons,
#>     and consider sensitivity analyses if missingness may depend on the
#>     unobserved outcome.
#>   Sources: Rubin (1976); Little & Rubin (2019)
#> 
#> Missing posttests and the pretest effects (see ?fit_solomon_glm)
#>   Posttests missing: 0.8% overall; pretested 0.0%, unpretested 1.7% (a
#>   difference of -1.7 points; Fisher's exact test, p = 1.000); by condition,
#>   control 0.0% and treatment 1.7% (the largest difference 1.7 points).
#>   Equal loss rates do not show that the same kinds of participants were lost.
#>   Pretested participants with a posttest and no pretest: 1. Their mean posttest
#>   was 60.00, against 55.29 for the pretested participants with both scores.
#>   Every fit with the pretest drops them, at either center (issue #142); the raw
#>   comparison below keeps them.
#>   Raw comparison (a fit without the pretest, HC3; 95% CI):
#>     Pretest effect | control   3.40  [-1.36,  8.16]
#>     Pretest effect | treated   0.99  [-4.19,  6.17]
#>     Pretest main effect        2.20  [-1.32,  5.71]
#>   The unpretested groups' losses cannot be checked against baseline (Campbell,
#>   1957, p. 305); with only those groups losing posttests, the two centers
#>   coincide, but neither is protected. Only the delta line reaches that case: if
#>   an unpretested group's dropouts would have scored Delta above its completers,
#>   the pretest effect in its condition is m0 x Delta lower, m0 being that
#>   group's share without a posttest (see ?fit_solomon_glm).

# A six-group design: two treatments and a control (Mai et al., 2020),
# with loss by gender in each pretest condition.
check_solomon_missing(post_behavior, condition, pretested, pre_behavior,
                      control = "Control", data = mai2020, covariates = "gender")
#> Solomon missingness check
#> Solomon N-group design: two treatments (RP, GS) and a control (Control), six groups
#> Pattern: structural pretest absence and incidental missingness
#> 
#>  Group                 Cell  n Post missing Pre absent (design) Pre missing
#>      1        Pretested, RP 35           11                   0           0
#>      2        Pretested, GS 33           10                   0           0
#>      3   Pretested, Control 50           23                   0           0
#>      4      Unpretested, RP 31            9                  31           0
#>      5      Unpretested, GS 27           12                  27           0
#>      6 Unpretested, Control 35           13                  35           0
#>  Pre unexpected
#>               0
#>               0
#>               0
#>               0
#>               0
#>               0
#> 
#> Structural pretest absence (n = 93)
#>   Participants assigned to the unpretested groups were never pretested; the
#>   absence of a pretest is the experimental manipulation.
#>   Response: Do not impute. Use analyses that respect the design, such as
#>     fit_solomon_glm(), fit_solomon_ml(), fit_solomon_classic(), or SEM.
#>     Unlike planned missing-data designs, where unmeasured values exist and
#>     can be imputed, an imputed pretest here would describe a measurement that
#>     never occurred.
#>   Sources: Solomon (1949); Graham et al. (2006)
#> 
#> Incidental posttest missingness (n = 78)
#>   Posttest scores are missing, so these participants are excluded from
#>   complete-case analyses.
#>   Response: For the treatment contrasts, complete-case analysis is unbiased
#>     when missingness is unrelated to the outcome given the variables in the
#>     model. The pretest effects also depend on the mean pretest at which they
#>     are evaluated; see `pretest_center` in fit_solomon_glm() and the
#>     diagnostics below. Report missingness by cell, because attrition that
#>     differs across the six groups can undermine the randomized comparisons,
#>     and consider sensitivity analyses if missingness may depend on the
#>     unobserved outcome.
#>   Sources: Rubin (1976); Little & Rubin (2019)
#> 
#> Missing posttests and the pretest effects (see ?fit_solomon_glm)
#>   Posttests missing: 37.0% overall; pretested 37.3%, unpretested 36.6% (a
#>   difference of 0.7 points; Fisher's exact test, p = 1.000); by condition,
#>   Control 42.4%, RP 30.3%, and GS 36.7% (the largest difference 12.0 points).
#>   Equal loss rates do not show that the same kinds of participants were lost.
#>   Mean pretest of the pretested participants without and with a posttest:
#>   Control 3.13 (23) and 3.12 (27); RP 3.22 (11) and 3.11 (24); GS 3.39 (10) and
#>   3.15 (23).
#>   Selection shift in each pretested condition, the mean pretest with a posttest
#>   minus every pretest's: Control -0.01 (SE 0.05), RP -0.03 (SE 0.04), and GS
#>   -0.07 (SE 0.04); Wald test of equal shifts, chi-square(2) = 1.10, p = .577.
#>   Gap between the centers: b m1 (completers - dropouts) = 0.45 x 0.37 x (3.13 -
#>   3.21) = -0.014, the pretest effects at the analyzed center minus those at the
#>   observed center (identity link, with the covariates).
#>   Unpretested groups, the share without a posttest and the offset the default
#>   assumes for their dropouts: Control 37.1% and 0.04; RP 29.0% and 0.05; GS
#>   44.4% and 0.03; the pretested dropouts' pretests predict 0.04 for them.
#>   Recorded in all groups, gender: posttests missing by level, pretested 33%
#>   Female and 57% Male (Fisher's exact test, p = .048); unpretested 32% Female
#>   and 64% Male (Fisher's exact test, p = .033). Among those with a posttest:
#>   pretested 88% Female and 12% Male; unpretested 92% Female and 8% Male.
#>   The unpretested groups' losses cannot be checked against baseline (Campbell,
#>   1957, p. 305).
```
