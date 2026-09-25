# Solomon Methods: History, Recommendations, and Extensions

## How to read this guide

This guide describes every analysis in `solomonR`: what it estimates,
what it assumes, where it falls short, and where it comes from. Each
method carries one of four labels.

| Label | Meaning |
|----|----|
| **Historical procedure** | Proposed in the Solomon literature and kept for teaching and replication. |
| **Contemporary recommendation** | A general methodological recommendation from current literature, applied here to the Solomon design. |
| **Published Solomon proposal** | A Solomon-specific method proposed in the literature but not yet widely established. |
| **solomonR extension** | An implementation choice or combination introduced by this package. |

A label describes where a method comes from, not whether it is correct.
Historical procedures can still give valid estimates, and package
extensions are tested against the methods they build on.

## The design and its questions

Solomon (1949) added two posttest-only groups to the pretest-posttest
control-group design so that the effect of taking a pretest could be
separated from the effect of the treatment. Campbell and Stanley (1963)
later classified it among the true experimental designs.

| Group | Pretest | Treatment | Posttest |
|:-----:|:-------:|:---------:|:--------:|
|   1   |   Yes   |    Yes    |   Yes    |
|   2   |   Yes   |    No     |   Yes    |
|   3   |   No    |    Yes    |   Yes    |
|   4   |   No    |    No     |   Yes    |

Every analysis in `solomonR` estimates some of four contrasts in mean
posttest scores, always as treatment minus control:

| Contrast | Question |
|----|----|
| Treatment \| pretested | Does the treatment work among participants who took the pretest (Group 1 vs. 2)? |
| Treatment \| unpretested | Does the treatment work among participants who did not (Group 3 vs. 4)? |
| Pretest x Treatment | Does the treatment effect depend on pretesting? This is pretest sensitization. |
| ATE (avg over pretest) | What is the treatment effect averaged equally over both pretest conditions? |

Naming the target quantity before choosing an estimator makes it clear
when different analyses answer the same question and when they do not
(Lundberg et al., 2021).

## Historical procedures

The historical analyses are available through
[`fit_solomon_classic()`](https://juhalt.github.io/solomonR/reference/fit_solomon_classic.md)
and are described in more detail in
[`vignette("classic-solomon")`](https://juhalt.github.io/solomonR/articles/classic-solomon.md).

### Four-group analysis of posttest scores (Tests A-D)

**Label:** historical procedure.

- **Estimates:** all four contrasts from a two-by-two model of posttest
  scores (Test A: Pretest x Treatment; B and C: the simple treatment
  effects; D: the equal-weighted average).
- **Assumptions:** independent observations, normally distributed
  errors, and equal variances across the four groups.
- **Limitations:** the pretest scores in Groups 1 and 2 are ignored,
  which gives up precision. Campbell and Stanley (1963) suggested
  analyzing the posttest scores with a two-by-two analysis of variance
  and, when neither the pretesting effect nor the interaction was
  significant, reanalyzing the pretested groups with ANCOVA for greater
  power. Choosing later tests according to earlier significance tests
  makes the overall procedure conditional.
- **Sources:** Campbell and Stanley (1963); Huck and Sandler (1973).

### Analyses within the pretested groups (Tests E-G)

**Label:** historical procedure.

- **Estimates:** Treatment \| pretested, using the pretest information.
  Test E is analysis of covariance (ANCOVA), Test F compares gain
  scores, and Test G is the two-occasion repeated-measures interaction,
  which is algebraically equivalent to Test F.
- **Assumptions:** for ANCOVA, a common pretest-posttest slope in both
  groups; a gain score implicitly assumes that slope equals one.
- **Limitations:** only Groups 1 and 2 are used. In randomized studies,
  ANCOVA has more power than change scores; in nonrandomized
  comparisons, ANCOVA can be more biased (Van Breukelen, 2006). Because
  Solomon groups are randomized, ANCOVA is generally the stronger
  choice.
- **Sources:** Huck and Sandler (1973); Van Breukelen (2006).

### Posttest-only comparison (Test H)

**Label:** historical procedure.

- **Estimates:** Treatment \| unpretested from Groups 3 and 4.
- **Assumptions:** independent observations, normal errors, equal
  variances in the two groups.
- **Limitations:** half of the sample is ignored.
- **Sources:** Huck and Sandler (1973).

### Meta-analytic combination (Test I)

**Label:** historical procedure.

- **Estimates:** nothing directly. Test I combines one-tailed p-values
  from the pretested and unpretested treatment comparisons with
  Stouffer’s method (Stouffer et al., 1949), as proposed by Braver and
  Braver (1988).
- **Assumptions:** both tests address the same directional hypothesis,
  and the combination is specified in advance.
- **Limitations:** in a Monte Carlo study, the conditional sequence that
  ends in Test I had an experiment-wise Type I error rate nearly three
  times the nominal alpha, and traditional procedures used without
  conditioning were more powerful (Sawilowsky et al., 1994). The
  proposal was debated in published exchanges between the two research
  groups. `solomonR` reproduces Test I for replication but does not
  recommend it as a default analysis.
- **Sources:** Braver and Braver (1988); Sawilowsky et al. (1994).

## Contemporary recommendations

### One prespecified analysis with explicit contrasts

**Label:** contemporary recommendation.

Rather than following a sequence of significance-driven decisions,
contemporary practice specifies the analysis and its target quantities
in advance (Lundberg et al., 2021). In randomized experiments, adjusting
for a baseline covariate such as the pretest is justified for estimating
treatment effects, and heteroskedasticity-robust standard errors protect
that inference (Lin, 2013).

### Heteroskedasticity-consistent and small-sample inference

**Label:** contemporary recommendation.

- **HC3 standard errors** are recommended as a routine choice,
  especially below about 250 observations (MacKinnon & White, 1985; Long
  & Ervin, 2000; Hayes & Cai, 2007).
  [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md)
  uses HC3 by default.
- **t reference distributions** with residual degrees of freedom are
  used for Gaussian models. Normal-reference robust intervals can
  under-cover in small samples (Imbens & Kolesár, 2016), and t with n -
  p degrees of freedom is the choice evaluated in recent guidance for
  behavioral research (Rajh-Weber et al., 2025).
- **Clustered data** use the CR2 estimator with Satterthwaite degrees of
  freedom (Bell & McCaffrey, 2002; Pustejovsky & Tipton, 2018).

### Randomization inference

**Label:** contemporary recommendation.

- **Tests:**
  [`perm_solomon()`](https://juhalt.github.io/solomonR/reference/perm_solomon.md)
  permutes treatment labels within pretest conditions, which reproduces
  how treatment was assigned. Its p-value is valid for the sharp null
  hypothesis of no effect for any participant, and the +1 correction
  keeps it from being zero (Phipson & Smyth, 2010). Studentized
  statistics make such tests robust when only an average effect is
  hypothesized to be zero (DiCiccio & Romano, 2017; Wu & Ding, 2021).
- **Limitations:** the permuted unit must be the randomized unit, so
  cluster-randomized studies are not supported yet.

### Equivalence testing for sensitization

**Label:** contemporary recommendation.

A nonsignificant Pretest x Treatment test does not show that
sensitization is absent.
[`equivalence_solomon()`](https://juhalt.github.io/solomonR/reference/equivalence_solomon.md)
applies the two one-sided tests procedure (Schuirmann, 1987) against
bounds set in advance at the smallest effect size of interest (Lakens,
2017; Lakens et al., 2018).

### Effect sizes with confidence intervals

**Label:** contemporary recommendation.

Standardized mean differences use noncentral t intervals (Cumming &
Finch, 2001; Kelley, 2007), and partial R-squared uses noncentral F
intervals when conventional covariance is used (Steiger, 2004). Hedges’
g corrects the small-sample bias of the standardized mean difference
(Hedges, 1981).

### Missing data

**Label:** contemporary recommendation, applied to the Solomon design.

- **Structural absence:** Groups 3 and 4 have no pretest by design. That
  absence is the manipulation, so it is never imputed. This differs from
  planned missing-data designs, in which unmeasured values exist (Graham
  et al., 2006).
- **Incidental missing pretests:** because treatment is randomized,
  simple approaches such as mean imputation with a missingness indicator
  retain participants without biasing the treatment effect (White &
  Thompson, 2005; Groenwold et al., 2012).
- **Missing outcomes:** the classification of missing-data mechanisms
  follows Rubin (1976) and Little and Rubin (2019).
- **Tools:**
  [`check_solomon_missing()`](https://juhalt.github.io/solomonR/reference/check_solomon_missing.md)
  and
  [`validate_solomon()`](https://juhalt.github.io/solomonR/reference/validate_solomon.md)
  report each category.

### Assumption diagnostics as description

**Label:** contemporary recommendation.

Choosing an analysis because a preliminary test of an assumption was or
was not significant can distort Type I error rates (Zimmerman, 2004).
[`check_solomon_assumptions()`](https://juhalt.github.io/solomonR/reference/check_solomon_assumptions.md)
therefore reports Brown-Forsythe tests (Brown & Forsythe, 1974),
normality tests, and a slope-homogeneity test as descriptions rather
than gates.

## Published Solomon proposals

### Full-information maximum likelihood

**Label:** published Solomon proposal.

- **Estimates:** all four contrasts, using the pretest in Groups 1 and 2
  and separate residual variances for pretested and unpretested
  participants (van Engelenburg, 1999).
- **Assumptions:** normal errors; a common pretest-posttest slope in the
  pretested groups; equal variances within each pretest condition.
- **Limitations:** its point estimates coincide with separate
  regressions in the two pretest conditions. Its default Wald intervals
  use a normal reference distribution and were too narrow with small
  groups in the package’s simulation validation. A small-sample option,
  `inference = "satterthwaite"`, uses unbiased residual variances and
  Welch-Satterthwaite degrees of freedom (Satterthwaite, 1946; Welch,
  1947); the function warns when groups are small and no option has been
  chosen.
- **Function:**
  [`fit_solomon_ml()`](https://juhalt.github.io/solomonR/reference/fit_solomon_ml.md).

## solomonR extensions

### Unified GLM with a design-coded pretest

**Label:** solomonR extension, built on contemporary recommendations.

- **Estimates:** all four contrasts from one model,
  `y ~ treat * pretested + pre_obs`. `pre_obs` equals the pretest in
  Groups 1 and 2 and 0 in Groups 3 and 4, so no group is dropped.
- **Assumptions:** a common pretest slope in the pretested groups;
  independent observations (or clusters with CR2).
- **Limitations:** because adjusting for the pretest reduces residual
  variance only in the pretested groups, the model is heteroskedastic
  whenever the pretest predicts the posttest. HC3 or CR2 covariance
  addresses this.
- **Checks:** its estimate of Treatment \| pretested equals the ANCOVA
  estimate (Test E), and its estimate of Treatment \| unpretested equals
  Test C. The package tests these identities.
- **Function:**
  [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md).

### Structural equation models

**Label:** solomonR extension, built on established SEM methods.

- **Estimates:** the four contrasts from a multi-group mean structure
  ([`fit_solomon_sem()`](https://juhalt.github.io/solomonR/reference/fit_solomon_sem.md)),
  or on a latent outcome measured by several items
  ([`fit_solomon_sem_latent()`](https://juhalt.github.io/solomonR/reference/fit_solomon_sem_latent.md)).
- **Assumptions:** latent mean comparisons require scalar measurement
  invariance across groups (Meredith, 1993; Vandenberg & Lance, 2000).
  Models are estimated with lavaan (Rosseel, 2012).
- **Limitations:** the four-group observed mean model is saturated, so
  its fit indices are not diagnostic. Latent contrasts are on the
  latent-variable scale and are not directly comparable with
  observed-score contrasts.

### Method comparison and design checks

**Label:** solomonR extension.

[`compare_solomon_methods()`](https://juhalt.github.io/solomonR/reference/compare_solomon_methods.md)
places the estimates from several analyses side by side and explains
which analyses are not comparable.
[`validate_solomon()`](https://juhalt.github.io/solomonR/reference/validate_solomon.md)
and
[`check_solomon_missing()`](https://juhalt.github.io/solomonR/reference/check_solomon_missing.md)
check the design before analysis.

### Power simulation

**Label:** solomonR extension, validated by simulation.

[`power_solomon()`](https://juhalt.github.io/solomonR/reference/power_solomon.md)
reports the rejection rate of each Solomon test with its Monte Carlo
standard error, naming the estimand and true effect behind every row. A
pre-specified study of 126 scenarios found that the 2x2 ANOVA
interaction matched its analytic benchmark in every scenario, that the
historical Test I held its nominal size under the complete null, and
that no fit failed. Rejection rates for the GLM tests are conservative
with small cells, following the HC3 standard errors used by default, so
treat them as a lower bound with 20 or fewer participants per cell. The
article “Validating power_solomon()” reports the study in full.
[`plan_solomon()`](https://juhalt.github.io/solomonR/reference/plan_solomon.md)
inverts the same calculations to find the smallest design that reaches a
target power, either analytically or by simulating the package’s own GLM
test.

## Summary

| Method | Label | Contrasts | Pretest used | Key assumption |
|----|----|----|----|----|
| Tests A-D | Historical | All four | No | Equal variances across four groups |
| Tests E-G | Historical | Treatment \| pretested | Yes | Common slope (E); slope of one (F, G) |
| Test H | Historical | Treatment \| unpretested | No | Equal variances in Groups 3-4 |
| Test I | Historical | None (combined p-value) | Yes | Same directional hypothesis |
| Unified GLM | solomonR extension | All four | Yes | Common slope; HC3 or CR2 inference |
| Randomization test | Contemporary | One contrast (test only) | Optional | Permute the randomized unit |
| Equivalence test | Contemporary | One contrast | As fitted | Bounds fixed in advance |
| Maximum likelihood | Published proposal | All four | Yes | Normal errors; separate variances by pretest condition |
| SEM (observed) | solomonR extension | All four | Optional | Multi-group mean structure |
| SEM (latent) | solomonR extension | All four (latent scale) | Optional | Scalar invariance |

## References

Bell, R. M., & McCaffrey, D. F. (2002). Bias reduction in standard
errors for linear regression with multi-stage samples. *Survey
Methodology, 28*, 169-181.

Braver, M. W., & Braver, S. L. (1988). Statistical treatment of the
Solomon four-group design: A meta-analytic approach. *Psychological
Bulletin, 104*, 150-154.

Brown, M. B., & Forsythe, A. B. (1974). Robust tests for the equality of
variances. *Journal of the American Statistical Association, 69*,
364-367.

Campbell, D. T., & Stanley, J. C. (1963). *Experimental and
quasi-experimental designs for research*. Rand McNally.

Cumming, G., & Finch, S. (2001). A primer on the understanding, use, and
calculation of confidence intervals that are based on central and
noncentral distributions. *Educational and Psychological Measurement,
61*, 532-574.

DiCiccio, C. J., & Romano, J. P. (2017). Robust permutation tests for
correlation and regression coefficients. *Journal of the American
Statistical Association, 112*, 1211-1220.

Graham, J. W., Taylor, B. J., Olchowski, A. E., & Cumsille, P. E.
(2006). Planned missing data designs in psychological research.
*Psychological Methods, 11*, 323-343.

Groenwold, R. H. H., White, I. R., Donders, A. R. T., Carpenter, J. R.,
Altman, D. G., & Moons, K. G. M. (2012). Missing covariate data in
clinical research: When and when not to use the missing-indicator method
for analysis. *Canadian Medical Association Journal, 184*, 1265-1269.

Hayes, A. F., & Cai, L. (2007). Using heteroskedasticity-consistent
standard error estimators in OLS regression: An introduction and
software implementation. *Behavior Research Methods, 39*, 709-722.

Hedges, L. V. (1981). Distribution theory for Glass’s estimator of
effect size and related estimators. *Journal of Educational Statistics,
6*, 107-128.

Huck, S. W., & Sandler, H. M. (1973). A note on the Solomon 4-group
design: Appropriate statistical analyses. *The Journal of Experimental
Education, 42*, 54-55.

Imbens, G. W., & Kolesár, M. (2016). Robust standard errors in small
samples: Some practical advice. *The Review of Economics and Statistics,
98*, 701-712.

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1-24.

Lakens, D. (2017). Equivalence tests: A practical primer for t tests,
correlations, and meta-analyses. *Social Psychological and Personality
Science, 8*, 355-362.

Lakens, D., Scheel, A. M., & Isager, P. M. (2018). Equivalence testing
for psychological research: A tutorial. *Advances in Methods and
Practices in Psychological Science, 1*, 259-269.

Lin, W. (2013). Agnostic notes on regression adjustments to experimental
data: Reexamining Freedman’s critique. *The Annals of Applied
Statistics, 7*, 295-318.

Little, R. J. A., & Rubin, D. B. (2019). *Statistical analysis with
missing data* (3rd ed.). Wiley.

Long, J. S., & Ervin, L. H. (2000). Using heteroscedasticity consistent
standard errors in the linear regression model. *The American
Statistician, 54*, 217-224.

Lundberg, I., Johnson, R., & Stewart, B. M. (2021). What is your
estimand? Defining the target quantity connects statistical evidence to
theory. *American Sociological Review, 86*, 532-565.

MacKinnon, J. G., & White, H. (1985). Some heteroskedasticity-consistent
covariance matrix estimators with improved finite sample properties.
*Journal of Econometrics, 29*, 305-325.

Meredith, W. (1993). Measurement invariance, factor analysis and
factorial invariance. *Psychometrika, 58*, 525-543.

Phipson, B., & Smyth, G. K. (2010). Permutation p-values should never be
zero: Calculating exact p-values when permutations are randomly drawn.
*Statistical Applications in Genetics and Molecular Biology, 9*(1),
Article 39.

Pustejovsky, J. E., & Tipton, E. (2018). Small-sample methods for
cluster-robust variance estimation and hypothesis testing in fixed
effects models. *Journal of Business & Economic Statistics, 36*,
672-683.

Rajh-Weber, H., Huber, S. E., & Arendasy, M. (2025). A practice-oriented
guide to statistical inference in linear modeling for non-normal or
heteroskedastic error distributions. *Behavior Research Methods, 57*,
Article 338.

Rosseel, Y. (2012). lavaan: An R package for structural equation
modeling. *Journal of Statistical Software, 48*(2), 1-36.

Rubin, D. B. (1976). Inference and missing data. *Biometrika, 63*,
581-592.

Satterthwaite, F. E. (1946). An approximate distribution of estimates of
variance components. *Biometrics Bulletin, 2*, 110-114.

Sawilowsky, S. S., Kelley, D. L., Blair, R. C., & Markman, B. S. (1994).
Meta-analysis and the Solomon four-group design. *The Journal of
Experimental Education, 62*, 361-376.

Schuirmann, D. J. (1987). A comparison of the two one-sided tests
procedure and the power approach for assessing the equivalence of
average bioavailability. *Journal of Pharmacokinetics and
Biopharmaceutics, 15*, 657-680.

Solomon, R. L. (1949). An extension of control group design.
*Psychological Bulletin, 46*, 137-150.

Steiger, J. H. (2004). Beyond the F test: Effect size confidence
intervals and tests of close fit in the analysis of variance and
contrast analysis. *Psychological Methods, 9*, 164-182.

Stouffer, S. A., Suchman, E. A., DeVinney, L. C., Star, S. A., &
Williams, R. M., Jr. (1949). *The American soldier: Adjustment during
army life* (Vol. 1). Princeton University Press.

Van Breukelen, G. J. P. (2006). ANCOVA versus change from baseline had
more power in randomized studies and more bias in nonrandomized studies.
*Journal of Clinical Epidemiology, 59*, 920-925.

van Engelenburg, G. (1999). *Statistical analysis for the Solomon
four-group design* (Research Report 99-06). University of Twente.

Vandenberg, R. J., & Lance, C. E. (2000). A review and synthesis of the
measurement invariance literature: Suggestions, practices, and
recommendations for organizational research. *Organizational Research
Methods, 3*, 4-70.

Welch, B. L. (1947). The generalization of “Student’s” problem when
several different population variances are involved. *Biometrika, 34*,
28-35.

White, I. R., & Thompson, S. G. (2005). Adjusting for partially missing
baseline measurements in randomized trials. *Statistics in Medicine,
24*, 993-1007.

Wu, J., & Ding, P. (2021). Randomization tests for weak null hypotheses
in randomized experiments. *Journal of the American Statistical
Association, 116*, 1898-1913.

Zimmerman, D. W. (2004). A note on preliminary tests of equality of
variances. *British Journal of Mathematical and Statistical Psychology,
57*, 173-181.
