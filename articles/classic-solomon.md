# The Classic Solomon Four-Group Analysis

## Why preserve the classic analysis?

The Solomon four-group design has a long methodological history. Earlier
approaches analyzed different pieces of the design sequentially, using
tests that came to be described as Tests A through I.

`solomonR` reproduces this historical workflow for two reasons:

1.  it remains useful for teaching how the Solomon design developed; and
2.  it allows researchers to reproduce analyses based on historical
    recommendations.

The classic workflow should not be interpreted as the default modern
analysis merely because it is historically important. `solomonR`
therefore distinguishes the historical procedure from the package’s
modern model-based and randomization-based approaches.

The implementation draws especially on the historical discussions of
Huck and Sandler (1973), Walton Braver and Braver (1988), and later
simulation work examining the operating characteristics of conditional
Solomon procedures.

## The four groups

A Solomon four-group design contains:

| Group | Pretest | Treatment | Posttest |
|:-----:|:-------:|:---------:|:--------:|
|   1   |   Yes   |    Yes    |   Yes    |
|   2   |   Yes   |    No     |   Yes    |
|   3   |   No    |    Yes    |   Yes    |
|   4   |   No    |    No     |   Yes    |

Groups 1 and 2 permit a pretest-adjusted treatment comparison.

Groups 3 and 4 provide a treatment comparison that cannot itself have
been influenced by administration of the pretest.

Together, the four groups permit investigation of both the treatment
effect and possible **pretest sensitization**.

## Example data

`solomonR` includes an example dataset with all four Solomon groups.

``` r

library(solomonR)

data(solomon_example)

demo_preview <- head(
  solomon_example,
  6
)

demo_preview$y_post <- round(
  demo_preview$y_post,
  2
)

demo_preview$y_pre <- round(
  demo_preview$y_pre,
  2
)

knitr::kable(
  demo_preview,
  align = "rrrr"
)
```

| y_post | treat | pretested | y_pre |
|-------:|------:|----------:|------:|
|     60 |     1 |         1 |    55 |
|     71 |     1 |         1 |    42 |
|     58 |     1 |         1 |    58 |
|     54 |     1 |         1 |    34 |
|     57 |     1 |         1 |    42 |
|     45 |     1 |         1 |    38 |

The pretest variable is structurally absent for participants who were
not pretested. Those missing values are therefore part of the design,
not ordinary missing observations.

## Running the historical analysis

The complete historical analysis is available through
[`fit_solomon_classic()`](https://juhalt.github.io/solomonR/reference/fit_solomon_classic.md):

``` r

classic <- with(
  solomon_example,
  fit_solomon_classic(
    y_post,
    treat,
    pretested,
    y_pre
  )
)

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
```

The printed output reports all of the historical tests but marks the
tests actually reached by the historical decision sequence with
`[PATH]`.

This distinction is important.

A test can be **calculated** for teaching or inspection without having
been **reached** under the historical conditional procedure.

[`plot_classic_flow()`](https://juhalt.github.io/solomonR/reference/plot_classic_flow.md)
draws the conditional sequence as a decision tree and highlights the
route these data took, with the p-value of each test reached:

``` r

plot_classic_flow(classic)
```

![Decision tree of the historical Test A to I sequence, with the path
taken by the example data
highlighted.](classic-solomon_files/figure-html/classic-flow-1.png)

The caption repeats the caution discussed below: the sequence ending in
Test I is kept for teaching and replication, not as a recommended
analysis.

## Tests A-D: the four-group posttest model

The first part of the classic workflow is based on the four-group
posttest model.

### Test A: Pretest x Treatment interaction

Test A asks whether the treatment effect differs between pretested and
unpretested participants.

In modern language, this is the pretest-by-treatment interaction and is
the primary statistical representation of pretest sensitization.

### Tests B and C: simple treatment effects

Test B estimates the treatment effect among pretested participants.

Test C estimates the treatment effect among unpretested participants.

### Test D: treatment main effect

Test D evaluates the treatment effect averaged equally across the
pretested and unpretested conditions.

In `solomonR`, this is represented by the contrast

``` math
\beta_{\mathrm{treat}} +
\frac{1}{2}\beta_{\mathrm{treat}\times\mathrm{pretest}}.
```

This equal-weighted contrast is important because the treatment
coefficient by itself represents the treatment effect only in the
reference pretest condition.

The individual historical results are available from the fitted object.
For example:

``` r

classic$tests$A$result
classic$tests$D$result
```

`classic$effects` holds Tests A-H and the pretest main effect in one
table, with the contrast labels and the columns of
[`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md).

### The analysis of variance table

`classic$anova` is the two-by-two analysis of variance of the posttest
that Tests A-D come from, with Type III sums of squares:

``` r

classic$anova
#>                source      sumsq  df    meansq         F   p.value
#> 1           Treatment   216.0083   1 216.00833 2.3370911 0.1290476
#> 2             Pretest   180.0750   1 180.07500 1.9483123 0.1654347
#> 3 Pretest x Treatment    27.0750   1  27.07500 0.2929366 0.5893828
#> 4               Error 10721.4333 116  92.42615        NA        NA
```

Each effect has one degree of freedom, and its F is the square of a t
test above: the treatment row is Test D, the interaction row is Test A,
and the pretest row is the pretest main effect. With unequal group sizes
these are not the sequential sums of squares that
[`aov()`](https://rdrr.io/r/stats/aov.html) prints, whose treatment row
ignores pretesting.

## Tests E-G: analyses within the pretested groups

Historical approaches proposed several ways to make use of the baseline
measurement available in Groups 1 and 2.

### Test E: ANCOVA

Test E compares treatment and control among pretested participants while
adjusting posttest scores for the pretest.

This is the default pretested-group method in
[`fit_solomon_classic()`](https://juhalt.github.io/solomonR/reference/fit_solomon_classic.md).

``` r

classic$tests$E$result
```

### Test F: gain-score analysis

Test F compares change scores between the two pretested groups.

``` r

classic$tests$F$result
```

The change being compared is visible in the pretested groups’ pretest
and posttest means.
[`plot_solomon_change()`](https://juhalt.github.io/solomonR/reference/plot_solomon_change.md)
shows them beside the unpretested groups, which are observed at posttest
only:

``` r

with(solomon_example, plot_solomon_change(y_post, treat, pretested, y_pre))
```

![Pretest and posttest means of the two pretested groups, beside the
posttest means of the two unpretested
groups.](classic-solomon_files/figure-html/change-plot-1.png)

### Test G: repeated-measures formulation

Test G represents the corresponding two-wave treatment-by-time
comparison.

For two measurement occasions, the inferential contrast is algebraically
connected to the gain-score comparison. `solomonR` retains Test G as a
separate historical label because the distinction is useful
pedagogically.

``` r

classic$tests$G$result
```

The pretested-group analysis can also be selected explicitly when
reproducing a particular historical strategy.

``` r

fit_solomon_classic(
  y_post,
  treat,
  pretested,
  y_pre,
  pretested_test = "gain"
)
```

## Test H: posttest-only comparison

Test H compares Groups 3 and 4, which did not receive the pretest.

``` r

classic$tests$H$result
```

Because these groups were never pretested, this comparison provides a
direct treatment contrast uncontaminated by administration of a baseline
measure.

`solomonR` also reports Hedges’ $`g`$ for this posttest-only comparison
when available.

``` r

cat(
  sprintf(
    "Hedges' g = %.2f, 95%% CI [%.2f, %.2f]",
    classic$g_post["g"],
    classic$g_post["lower"],
    classic$g_post["upper"]
  )
)
#> Hedges' g = 0.42, 95% CI [-0.09, 0.94]
```

## Test I: Walton Braver & Braver (1988) Stouffer combination

Walton Braver and Braver (1988) proposed combining evidence from the
pretested and unpretested treatment comparisons using Stouffer’s method.

When requested, `solomonR` reproduces this historical analysis using
**directional one-tailed p-values aligned with the same treatment
direction**.

``` r

cat(
  sprintf(
    "Stouffer Z = %.2f, p = %s",
    classic$tests$I$result$z,
    if (
      classic$tests$I$result$p.value < .001
    ) {
      "< .001"
    } else {
      sub(
        "^0",
        "",
        sprintf(
          "%.3f",
          classic$tests$I$result$p.value
        )
      )
    }
  )
)
#> Stouffer Z = 1.69, p = .090
```

This detail matters. A directional one-tailed p-value is not obtained
correctly by mechanically dividing every two-sided p-value by two
without considering the sign of the effect. The package’s replication of
the published error rates found that the published simulations appear to
have combined the tests without regard to sign, converting each
two-sided p-value to z as if it were one-tailed. That reading reproduces
the published Test I rates, which are far lower than the procedure as
defined produces (see [Historical Tests: Replicating the Published Error
Rates](https://juhalt.github.io/solomonR/articles/classic-validation.html)).
`solomonR` refers the combined z to its two-tailed p-value, as Walton
Braver and Braver report it in their worked example (1988, p. 153).

## Versions of the decision sequence

The sequence changed after it was published, and
[`fit_solomon_classic()`](https://juhalt.github.io/solomonR/reference/fit_solomon_classic.md)
reproduces each version through its `flow` argument.

- **1988, the default.** In the original sequence, testing stops at the
  first significant test among D, the selected pretested-groups test,
  and H; Test I is reached only when all of them are nonsignificant
  (Walton Braver & Braver, 1988, pp. 151–153).
- **1990.** Replying to Sawilowsky and Markman (1990a), the authors
  agreed that using the meta-analysis selectively introduces selection
  bias. Once Tests A and D are nonsignificant, the analyst “should
  complete all the tests through Test I” (Braver & Walton Braver, 1990,
  p. 322), and Test I is regarded as the most definitive test.
- **1995.** A later revision removed Test D. It is unpublished and known
  from Sawilowsky (1996, p. 2), so it is cited as Walton Braver and
  Braver (1995, as cited in Sawilowsky, 1996). The remaining tests
  follow the 1988 rule, as in Sawilowsky’s simulation of the revised
  sequence.

``` r

classic_1990 <- with(
  solomon_example,
  fit_solomon_classic(y_post, treat, pretested, y_pre, flow = "1990")
)

classic_1990$path_string
#> [1] "A -> D -> E -> H -> I"
classic_1990$conclusion
#> [1] "Historical pathway (1990 amendment): Tests A and D are nonsignificant, so every test through Test I is run and Test I is regarded as the most definitive. Test I does not reach the specified alpha level (Test E p = 0.444, Test H p = 0.103, Test I p = 0.090). Interpret this in light of later Type I error critiques."
```

The tests themselves are identical in every version; only the path and
the conclusion differ. For the 1995 version, `alpha_allocation` applies
the test-wise significance levels that Sawilowsky (1996, Table 4)
proposed to control the experiment-wise error rate:

``` r

classic_1995 <- with(
  solomon_example,
  fit_solomon_classic(y_post, treat, pretested, y_pre, flow = "1995",
                      alpha_allocation = "method2_conservative")
)

classic_1995$path_string
#> [1] "A -> E -> H -> I"
```

## Published worked examples

Walton Braver and Braver (1988) illustrated the sequence with
hypothetical data: four groups of 14, with the means and variances and
the pretest-posttest correlations printed in their Table 3 (p. 153). The
data set `waltonbraver1988` holds those numbers. Table 3 prints
variances, so `sd` and `pre_sd` are their square roots.

``` r

waltonbraver1988[, c("group", "n", "mean", "var", "pre_mean", "pre_var", "r")]
#>                    group  n mean  var pre_mean pre_var    r
#> 1   Pretested, treatment 14 12.4 22.0     10.5    19.3 0.58
#> 2     Pretested, control 14 10.2 16.5     10.7    17.2 0.62
#> 3 Unpretested, treatment 14 12.5 19.0       NA      NA   NA
#> 4   Unpretested, control 14 10.3 22.5       NA      NA   NA
```

Tests A to D need only the posttest statistics.
[`solomon_from_summary()`](https://juhalt.github.io/solomonR/reference/solomon_from_summary.md)
gives the published two-by-two analysis of variance (Table 4, p. 153):
no interaction at all (Test A), and a treatment effect (Test D) that is
not quite significant, F(1, 52) = 3.388, p = .0714.

``` r

with(waltonbraver1988, solomon_from_summary(n, mean, sd))
#> Solomon analysis from summary statistics
#> ----------------------------------------
#> Pooled error variance: 20.000 on 52 df (equal variances assumed)
#> 
#> Two-way ANOVA on the posttest (Type III sums of squares)
#>   Treatment            SS =   67.760  df = 1  F = 3.39  p = 0.071
#>   Pretest              SS =    0.140  df = 1  F = 0.01  p = 0.934
#>   Pretest x Treatment  SS =    0.000  df = 1  F = 0.00  p = 1.000
#>   Error                SS = 1040.000  df = 52
#> 
#> Contrasts with 95% confidence intervals
#>   Test A: Pretest x Treatment               0.000 [-4.797, 4.797], t(52) = 0.00, p = 1.000
#>   Test B: Treatment | pretested             2.200 [-1.192, 5.592], t(52) = 1.30, p = 0.199
#>   Test C: Treatment | unpretested           2.200 [-1.192, 5.592], t(52) = 1.30, p = 0.199
#>   Test D: ATE (avg over pretest)            2.200 [-0.198, 4.598], t(52) = 1.84, p = 0.071
#>   Pretest main effect                      -0.100 [-2.498, 2.298], t(52) = -0.08, p = 0.934
```

[`fit_solomon_classic()`](https://juhalt.github.io/solomonR/reference/fit_solomon_classic.md)
takes individual scores. Scores built to have exactly the published
means, variances, and correlations give the same tests as the original
scores would, because Tests A to I depend on the data only through those
statistics:

``` r

set.seed(1988)
wbb_scores <- do.call(rbind, lapply(seq_len(4), function(i) {
  g <- waltonbraver1988[i, ]
  if (g$pretested == 1) {
    s <- g$r * g$pre_sd * g$sd
    x <- MASS::mvrnorm(g$n, c(g$pre_mean, g$mean),
                       matrix(c(g$pre_var, s, s, g$var), 2), empirical = TRUE)
  } else {
    x <- cbind(NA, g$mean + g$sd * as.vector(scale(stats::rnorm(g$n))))
  }
  data.frame(treat = g$treat, pretested = g$pretested, y_pre = x[, 1], y_post = x[, 2])
}))
wbb <- fit_solomon_classic(y_post, treat, pretested, y_pre, data = wbb_scores)
wbb$path_string
#> [1] "A -> D -> E -> H -> I"
wbb$tests$E$result[, c("test", "F", "df", "p.value")]
#>   test        F df    p.value
#> 1    E 2.931688 25 0.09923573
wbb$tests$H$result[, c("test", "statistic", "df", "p.value")]
#>   test statistic df   p.value
#> 1    H  1.277799 26 0.2126117
wbb$tests$I$result[, c("z", "p.value")]
#>          z    p.value
#> 1 2.047064 0.04065177
```

The path is the published one: Tests A, D, E, and H are not significant,
and Test I is. The published values are F(1, 25) = 2.93, p = .0993, for
the analysis of covariance (Test E; Table 5); t(26) = 1.28, p = .2127,
for Test H; and z = 2.05, p = .040, for Test I (p. 153). The published p
of Test I is the two-tailed p of the rounded z.
[`?waltonbraver1988`](https://juhalt.github.io/solomonR/reference/waltonbraver1988.md)
lists each published value beside its reproduction.

Sawilowsky and Markman (1988) answered with fabricated scores (Table 2,
p. 7) for which Test H is significant, t(26) = 2.07, p = .048 (Table 5,
p. 10), and Test I, as they computed it, is not; they judged the
combined z one-tailed (pp. 3–4). The manuscript does not say how the two
z values were obtained. Their values, 1.98 and .08, correspond to
halving each two-sided p-value without regard to the sign of Test E’s
effect, which is negative in their data; this is solomonR’s
reconstruction from the numbers. That halving differs from the reading
that reproduces the published simulation rates, described under Test I
above. The package’s tests reproduce their Tests E and H from the
scores. Its own Test I is directional, as Walton Braver and Braver
(1988, p. 152) define it, so it gives z = 1.34, not their 1.46. The two
worked examples also judge Test I differently, two-tailed in Walton
Braver and Braver’s and one-tailed in Sawilowsky and Markman’s; the
package follows the former.

## History and maturation

The four-group design also allows a check that does not involve the
treatment. The unpretested control group’s posttest (O6) was measured at
the same time as the other posttests, but its participants received
neither the pretest nor the treatment. Comparing it with the pretests of
the pretested groups (O1 and O3), which were measured before any
treatment, estimates the combined effect of history and maturation
between the two occasions. Mai et al. (2020, p. 8) report these
comparisons as independent-samples t tests, which
[`fit_solomon_classic()`](https://juhalt.github.io/solomonR/reference/fit_solomon_classic.md)
returns as `history`:

``` r

classic$history
#>                                            comparison estimate statistic df
#> 1 O6 - O1: control posttest vs. treated-group pretest 1.466667 0.5415318 58
#> 2 O6 - O3: control posttest vs. control-group pretest 1.533333 0.6121767 58
#>     p.value  conf.low conf.high
#> 1 0.5902153 -3.954718  6.888051
#> 2 0.5428135 -3.480416  6.547082
```

The check assumes random assignment and a measure that is comparable at
the two occasions.

## Why Test I is treated cautiously

The Stouffer procedure is included because it is part of the history of
the Solomon design and may be required when reproducing historical
analyses.

It is not the default modern recommendation in `solomonR`.

Later simulation work showed that conditional analysis sequences can
have undesirable experiment-wise Type I error properties. In particular,
choosing later tests based on the statistical significance of earlier
tests changes the operating characteristics of the overall procedure.
Sawilowsky et al. (1994) found that the 1988 sequence falsely declared
an effect about 14% of the time at a nominal 5% per test.

The package’s own replication found similar rates for the sequence as
Walton Braver and Braver define it. With Test I judged two-tailed, as in
their worked example, the rate was 13.5% to 13.7% for the 1988 and 1995
versions, and Sawilowsky’s (1996) alpha allocations stayed above their
targets.

For contemporary applied work, researchers should therefore consider a
single prespecified model or a randomization-based analysis rather than
automatically following a significance-driven testing tree.

See the unified GLM vignette for the primary modern observed-variable
workflow:

``` r

vignette("glm-solomon", package = "solomonR")
```

## Historical analysis versus modern analysis

The two approaches answer related questions but serve different
purposes.

| Goal | Suggested approach |
|----|----|
| Teach the historical Solomon procedure | [`fit_solomon_classic()`](https://juhalt.github.io/solomonR/reference/fit_solomon_classic.md) |
| Reproduce a historical Tests A-I analysis | [`fit_solomon_classic()`](https://juhalt.github.io/solomonR/reference/fit_solomon_classic.md) |
| Estimate Solomon effects in one model | [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md) |
| Conduct randomization-based inference | [`perm_solomon()`](https://juhalt.github.io/solomonR/reference/perm_solomon.md) |
| Use full-information likelihood | [`fit_solomon_ml()`](https://juhalt.github.io/solomonR/reference/fit_solomon_ml.md) |
| Model latent outcomes | [`fit_solomon_sem_latent()`](https://juhalt.github.io/solomonR/reference/fit_solomon_sem_latent.md) |

The historical procedure is therefore preserved rather than erased,
while modern methods are available alongside it. [A History of the
Solomon Design and Its
Analysis](https://juhalt.github.io/solomonR/articles/history.html) tells
how each version arose and how the debate over it unfolded.

All works cited in solomonR are listed, with notes on how the package
uses them, on the
[References](https://juhalt.github.io/solomonR/articles/references.html)
page.

## References

Braver, S. L., & Walton Braver, M. C. (1990). Meta-analysis for Solomon
four-group designs reconsidered: A reply to Sawilowsky and Markman.
*Perceptual and Motor Skills, 71*(1), 321–322.
<https://doi.org/10.2466/pms.1990.71.1.321>

Huck, S. W., & Sandler, H. M. (1973). A note on the Solomon 4-group
design: Appropriate statistical analyses. *The Journal of Experimental
Education, 42*(2), 54–55.
<https://doi.org/10.1080/00220973.1973.11011460>

Mai, N. N., Takahashi, Y., & Oo, M. M. (2020). Testing the effectiveness
of transfer interventions using Solomon four-group designs. *Education
Sciences, 10*(4), Article 92. <https://doi.org/10.3390/educsci10040092>

Sawilowsky, S. S. (1996, June 23). *Controlling experiment-wise Type I
error of meta-analysis in the Solomon four-group design* \[Paper
presentation\]. First International Conference on Multiple Comparisons,
Tel Aviv, Israel. <https://digitalcommons.wayne.edu/coe_tbf/29/>

Sawilowsky, S. S., Kelley, D. L., Blair, R. C., & Markman, B. S. (1994).
Meta-analysis and the Solomon four-group design. *The Journal of
Experimental Education, 62*(4), 361–376.
<https://doi.org/10.1080/00220973.1994.9944140>

Sawilowsky, S. S., & Markman, B. S. (1988). *Another look at the power
of meta-analysis in the Solomon four-group design* (ED316556). ERIC.
<https://eric.ed.gov/?id=ED316556>

Sawilowsky, S. S., & Markman, B. S. (1990a). Another look at the power
of meta-analysis in the Solomon four-group design. *Perceptual and Motor
Skills, 71*(1), 177–178. <https://doi.org/10.2466/pms.1990.71.1.177>

Walton Braver, M. C., & Braver, S. L. (1988). Statistical treatment of
the Solomon four-group design: A meta-analytic approach. *Psychological
Bulletin, 104*(1), 150–154.
<https://doi.org/10.1037/0033-2909.104.1.150>
