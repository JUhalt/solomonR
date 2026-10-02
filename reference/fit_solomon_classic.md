# Historical Solomon Four-Group Analysis

**\[stable\]** Implements the historical Test A-I framework associated
with analysis of the Solomon four-group design. All tests are calculated
when possible, while the historical decision pathway is stored
separately.

## Usage

``` r
fit_solomon_classic(
  y_post,
  treat,
  pretested,
  y_pre,
  alpha = 0.05,
  pretested_test = c("ancova", "gain", "repeated"),
  stouffer = TRUE,
  stouffer_direction = c("greater", "less"),
  conf_level = 0.95,
  flow = c("1988", "1990", "1995"),
  alpha_allocation = c("none", "method1_conservative", "method1_liberal",
    "method2_conservative", "method2_liberal"),
  data = NULL,
  combine_with_stouffer = deprecated()
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

  Numeric pretest scores. Values should be missing for participants
  assigned to the unpretested groups.

- alpha:

  Significance level used to reconstruct the historical decision
  pathway. Default is 0.05.

- pretested_test:

  Which historical pretested-group analysis should be followed in the
  decision pathway: `"ancova"`, `"gain"`, or `"repeated"`. All three are
  still calculated and returned.

- stouffer:

  Logical. If `TRUE`, include Test I, the Walton Braver & Braver (1988)
  Stouffer combination, in the decision pathway when earlier treatment
  tests are nonsignificant.

- stouffer_direction:

  Direction of the historical one-tailed treatment hypothesis used for
  Test I: `"greater"` or `"less"`.

- conf_level:

  Confidence level for intervals. Default is 0.95.

- flow:

  Which published version of the decision sequence to follow. `"1988"`
  (the default) is the original sequence of Walton Braver and Braver
  (1988, pp. 151–153): testing stops at the first significant test among
  D, the selected pretested-groups test (E, F, or G), and H, and Test I
  is reached only when all are nonsignificant. `"1990"` is the authors'
  amendment (Braver & Walton Braver, 1990, p. 322): once Tests A and D
  are nonsignificant, every test through Test I is run and Test I is
  regarded as the most definitive. `"1995"` is the authors' later
  revision, which removed Test D (Walton Braver & Braver, 1995, as cited
  in Sawilowsky, 1996, p. 2); the revision itself is unpublished. The
  remaining tests follow the 1988 rule, each reached only if the ones
  before it are nonsignificant, as in Sawilowsky's (1996) simulation of
  the revised sequence. Only the `path` and `conclusion` differ; every
  test is always calculated.

- alpha_allocation:

  Test-wise significance levels for the historical sequence. `"none"`
  (the default) uses `alpha` for every test. The other options are the
  two Bonferroni-type allocations of Sawilowsky (1996, Table 4), which
  he calibrated by Monte Carlo so that the experiment-wise Type I error
  of the 1995 sequence stays within Bradley's (1968, as cited in
  Sawilowsky, 1996) conservative (`"_conservative"`) or liberal
  (`"_liberal"`) robustness limit for a nominal alpha of .05:

  - `"method1_conservative"` and `"method1_liberal"`: Tests A, E, H, and
    I each at .020 or .0275;

  - `"method2_conservative"` and `"method2_liberal"`: Test A at .05, as
    a preliminary test, and Tests E, H, and I each at .005 or .020.

  An allocation requires `flow = "1995"`, `alpha = 0.05`,
  `pretested_test = "ancova"`, and `stouffer = TRUE`, the conditions for
  which the levels were obtained. Tests B and C, reached only after a
  significant Test A, are outside the allocation and keep `alpha`; that
  choice is solomonR's.

- data:

  Optional data frame. When supplied, the other data arguments are
  looked up in it first, as bare column names (`y_post = post`) or as
  strings (`y_post = "post"`).

- combine_with_stouffer:

  **\[deprecated\]** Use `stouffer`.

## Value

An object of class `solomon_classic`. The `tests` component contains
Tests A-I, while `path` records the historical decision sequence for the
observed data under the chosen `flow`. The `history` component holds the
history/maturation comparisons.

## Details

The historical sequence includes the four-group posttest factorial
model, simple treatment effects, ANCOVA, gain-score analysis, the
equivalent two-wave repeated-measures interaction, the posttest-only
comparison, and the optional Stouffer meta-analytic combination.

Test I, the Walton Braver & Braver (1988) Stouffer meta-analytic
combination, is included for historical replication and teaching. It
converts the one-tailed p-values of the pretested-groups test and Test
H, in the direction of `stouffer_direction`, to z and combines them (p.
152). The combined z is judged by its two-tailed p-value, as in their
worked example (p. 153: z = 2.05, p = .040); the one-tailed value is
kept as `p_one_tailed`. Later work (see Sawilowsky et al., 1994) raised
concerns about experiment-wise Type I error when the procedure is used
conditionally. Its presence in this function should not be interpreted
as a general contemporary recommendation.

## Error rates and the published replication

The package's replication of the published Type I error rates (issue
\#51; normal data, 30 per group, 20,000 datasets per condition) found:

- **Tests A to H** reached and rejected at the rates Sawilowsky et al.
  (1994, Table 2) and Sawilowsky (1996, Table 1) report.

- **Test I** rejected far more often than published. With its two-tailed
  criterion, the experiment-wise error rates were .137 for the 1988
  sequence and .135 for the 1995 sequence, against the published .137
  and .125. A post hoc investigation reproduced the published Test I
  rates only when the two-sided p-values of Tests E and H were converted
  to z as if one-tailed, a reading that differs from Walton Braver and
  Braver's (1988, p. 152) definition.

- **The alpha allocations**, with Test I as defined here, gave
  experiment-wise error rates of .059, .079, .056, and .081 for
  `"method1_conservative"`, `"method1_liberal"`,
  `"method2_conservative"`, and `"method2_liberal"`, against
  Sawilowsky's .054, .072, .057, and .076. All four exceed their limits
  of .055 and .075.

Until this version, Test I was judged by the one-tailed p of the
combined z. It now follows Walton Braver and Braver's worked example (p.
153), which reports the two-tailed p.

See the article "Historical Tests: Replicating the Published Error
Rates".

## History and maturation

The `history` component compares the unpretested control posttest (O6)
with the pretests of the pretested groups (O1 and O3) by
independent-samples t tests, as Mai et al. (2020, p. 8) report. Neither
group had received the treatment when these scores were measured, so a
difference estimates the combined effect of history and maturation
between the two occasions (Campbell & Stanley, 1963/1966, p. 25),
provided assignment was random and the measure is comparable at both
occasions. Solomon (1949, pp. 146–148) introduced the fourth group for
this purpose, attributing its change from the pretest to outside events
and the passage of time, and Campbell (1957, p. 303) first proposed the
t test, together with the 2 x 2 analysis of variance of the posttests
that Tests A–D carry out. It is reported as a historical check, not a
test of the treatment.

## Confidence intervals

Tests A-H are t tests (equivalently, F tests with one numerator degree
of freedom) with residual degrees of freedom, and each result carries
the matching confidence interval (`conf.low`, `conf.high`). Test I
combines p-values and has no interval. The Groups 3-4 standardized mean
difference (`g_post`) reports Hedges' g with a noncentral t interval for
the population standardized mean difference (Cumming & Finch, 2001;
Kelley, 2007).

## References

Braver, S. L., & Walton Braver, M. C. (1990). Meta-analysis for Solomon
four-group designs reconsidered: A reply to Sawilowsky and Markman.
*Perceptual and Motor Skills, 71*(1), 321–322.
https://doi.org/10.2466/pms.1990.71.1.321

Campbell, D. T. (1957). Factors relevant to the validity of experiments
in social settings. *Psychological Bulletin, 54*(4), 297–312.
https://doi.org/10.1037/h0040950

Campbell, D. T., & Stanley, J. C. (1966). *Experimental and
quasi-experimental designs for research*. Rand McNally. (Original work
published 1963)

Cumming, G., & Finch, S. (2001). A primer on the understanding, use, and
calculation of confidence intervals that are based on central and
noncentral distributions. *Educational and Psychological Measurement,
61*(4), 532–574. https://doi.org/10.1177/00131640121971374

Huck, S. W., & Sandler, H. M. (1973). A note on the Solomon 4-group
design: Appropriate statistical analyses. *The Journal of Experimental
Education, 42*(2), 54–55. https://doi.org/10.1080/00220973.1973.11011460

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24. https://doi.org/10.18637/jss.v020.i08

Mai, N. N., Takahashi, Y., & Oo, M. M. (2020). Testing the effectiveness
of transfer interventions using Solomon four-group designs. *Education
Sciences, 10*(4), Article 92. https://doi.org/10.3390/educsci10040092

Sawilowsky, S. S. (1996, June 23). *Controlling experiment-wise Type I
error of meta-analysis in the Solomon four-group design* \[Paper
presentation\]. First International Conference on Multiple Comparisons,
Tel Aviv, Israel. <https://digitalcommons.wayne.edu/coe_tbf/29/>

Sawilowsky, S. S., Kelley, D. L., Blair, R. C., & Markman, B. S. (1994).
Meta-analysis and the Solomon four-group design. *The Journal of
Experimental Education, 62*(4), 361–376.
https://doi.org/10.1080/00220973.1994.9944140

Solomon, R. L. (1949). An extension of control group design.
*Psychological Bulletin, 46*(2), 137–150.
https://doi.org/10.1037/h0062958

Van Breukelen, G. J. P. (2006). ANCOVA versus change from baseline had
more power in randomized studies and more bias in nonrandomized studies.
*Journal of Clinical Epidemiology, 59*(9), 920–925.
https://doi.org/10.1016/j.jclinepi.2006.02.007

Walton Braver, M. C., & Braver, S. L. (1988). Statistical treatment of
the Solomon four-group design: A meta-analytic approach. *Psychological
Bulletin, 104*(1), 150–154. https://doi.org/10.1037/0033-2909.104.1.150

## Examples

``` r
classic <- with(solomon_example, fit_solomon_classic(y_post, treat, pretested, y_pre))
classic
#> Classic Solomon analysis (historical teaching workflow)
#> -------------------------------------------------------
#> Selected pretested-group method: Test E (ancova)
#> Historical decision path (1988 flow): A -> D -> E -> H -> I
#> 
#> Historical Tests A-I
#> --------------------
#> All tests are shown below. Tests marked [PATH] were reached by
#> the historical decision sequence for these data.
#> 
#> [PATH] Test A: Pretest x Treatment interaction               F(1, 116) = 0.29, p = 0.589
#>        Test B: Treatment effect among pretested groups       F(1, 116) = 0.49, p = 0.486
#>        Test C: Treatment effect among unpretested groups     F(1, 116) = 2.14, p = 0.146
#> [PATH] Test D: Treatment main effect                         F(1, 116) = 2.34, p = 0.129
#> [PATH] Test E: ANCOVA treatment effect                       F(1, 57) = 0.59, p = 0.444
#>        Test F: Gain-score treatment effect                   F(1, 58) = 0.46, p = 0.499
#>        Test G: Repeated-measures Treatment x Time interaction F(1, 58) = 0.46, p = 0.499
#> [PATH] Test H: Posttest-only treatment effect                t(58) = 1.66, p = 0.103
#> [PATH] Test I: Walton Braver & Braver (1988) Stouffer combination Z = 1.69, p = 0.090 [E + H (ANCOVA + posttest-only)]
#> 
#> Historical interpretation
#> -------------------------
#> Historical pathway: no treatment test in the selected A-I sequence reaches the specified alpha level. 
#> 
#> Groups 3-4 effect size: Hedges g = 0.423, 95% CI [-0.086, 0.938] (noncentral t)
#> 
#> History/maturation check (historical; Mai et al., 2020)
#>   O6 - O1: control posttest vs. treated-group pretest: difference = 1.467, t(58) = 0.54, p = 0.590
#>   O6 - O3: control posttest vs. control-group pretest: difference = 1.533, t(58) = 0.61, p = 0.543
#> 
#> Caution: Test I, the Walton Braver & Braver (1988) Stouffer combination, is
#> reproduced for historical teaching and replication. Later simulation
#> work (see Sawilowsky et al., 1994) raised concerns about Type I error
#> for the conditional meta-analytic sequence; it is not the default
#> modern inferential recommendation in solomonR.

# The 1995 revision of the sequence, without Test D (Sawilowsky, 1996).
with(solomon_example, fit_solomon_classic(y_post, treat, pretested, y_pre, flow = "1995"))
#> Classic Solomon analysis (historical teaching workflow)
#> -------------------------------------------------------
#> Selected pretested-group method: Test E (ancova)
#> Historical decision path (1995 flow): A -> E -> H -> I
#> 
#> Historical Tests A-I
#> --------------------
#> All tests are shown below. Tests marked [PATH] were reached by
#> the historical decision sequence for these data.
#> 
#> [PATH] Test A: Pretest x Treatment interaction               F(1, 116) = 0.29, p = 0.589
#>        Test B: Treatment effect among pretested groups       F(1, 116) = 0.49, p = 0.486
#>        Test C: Treatment effect among unpretested groups     F(1, 116) = 2.14, p = 0.146
#>        Test D: Treatment main effect                         F(1, 116) = 2.34, p = 0.129
#> [PATH] Test E: ANCOVA treatment effect                       F(1, 57) = 0.59, p = 0.444
#>        Test F: Gain-score treatment effect                   F(1, 58) = 0.46, p = 0.499
#>        Test G: Repeated-measures Treatment x Time interaction F(1, 58) = 0.46, p = 0.499
#> [PATH] Test H: Posttest-only treatment effect                t(58) = 1.66, p = 0.103
#> [PATH] Test I: Walton Braver & Braver (1988) Stouffer combination Z = 1.69, p = 0.090 [E + H (ANCOVA + posttest-only)]
#> 
#> Historical interpretation
#> -------------------------
#> Historical pathway (1995 revision, without Test D): no treatment test in the selected A-I sequence reaches the specified alpha level. 
#> 
#> Groups 3-4 effect size: Hedges g = 0.423, 95% CI [-0.086, 0.938] (noncentral t)
#> 
#> History/maturation check (historical; Mai et al., 2020)
#>   O6 - O1: control posttest vs. treated-group pretest: difference = 1.467, t(58) = 0.54, p = 0.590
#>   O6 - O3: control posttest vs. control-group pretest: difference = 1.533, t(58) = 0.61, p = 0.543
#> 
#> Caution: Test I, the Walton Braver & Braver (1988) Stouffer combination, is
#> reproduced for historical teaching and replication. Later simulation
#> work (see Sawilowsky et al., 1994) raised concerns about Type I error
#> for the conditional meta-analytic sequence; it is not the default
#> modern inferential recommendation in solomonR.
```
