# Simulate a Solomon four-group study with known effects

**\[stable\]** Generates individual data for a randomized Solomon
four-group study with chosen treatment, pretest, and sensitization
effects, and attaches the true value of every Solomon estimand. It is a
solomonR teaching tool: instructors can show that each analysis recovers
the effects that were built in, and researchers can check an analysis
plan on data whose answer is known. With several treatment effects it
simulates a Solomon N-group design (see "Designs with several
treatments").

## Usage

``` r
simulate_solomon(
  n = 30,
  delta = 0,
  sens = 0,
  pretest_effect = 0,
  rho = 0.5,
  sigma = 1,
  mean = 0,
  digits = NULL,
  limits = NULL,
  seed = NULL
)
```

## Arguments

- n:

  Participants per group: one number for all groups, or one number per
  group. For the four-group design, four numbers in the order pretested
  treatment, pretested control, unpretested treatment, unpretested
  control. For k treatments, 2(k + 1) numbers in the order given in
  "Designs with several treatments". Default 30 per group.

- delta:

  Treatment effect among unpretested participants. For a design with
  several treatments, one effect per treatment, named by treatment.
  Default 0.

- sens:

  Pretest x Treatment interaction (sensitization): the extra treatment
  effect among pretested participants. For a design with several
  treatments, one number for every treatment or one per treatment.
  Default 0.

- pretest_effect:

  Effect of taking the pretest on the posttest, among controls. Default
  0.

- rho:

  Pretest-posttest correlation within the pretested groups. Default 0.5.

- sigma:

  Standard deviation of the scores within each group. Default 1.

- mean:

  Mean pretest and control-group posttest score. Default 0.

- digits:

  Optional number of decimal places to round scores to.

- limits:

  Optional lower and upper limits for scores, applied after rounding.

- seed:

  Optional integer seed. The global random number state is left
  unchanged.

## Value

A data frame with `y_post`, `treat`, `pretested`, and `y_pre` (missing
by design in the unpretested groups), in the form the analysis functions
take. `treat` is 0/1 for the four-group design and a factor of
conditions for a design with several treatments. Its `"truth"` attribute
is a data frame of the true estimands, and its `"settings"` attribute
holds the arguments, with the group sizes named by group.

## Details

**Data-generating model.** It is the model of
[`power_solomon()`](https://juhalt.github.io/solomonR/reference/power_solomon.md),
whose rebuild was validated under a protocol posted on issue \#18, with
a pretesting main effect and a location added. Every participant has a
latent baseline A from a standard normal distribution. Participants in
the pretested groups have a pretest score of `mean + sigma * A`. Every
participant's posttest score is

`mean + delta * treat + pretest_effect * pretested + sens * treat * pretested + sigma * (rho * A + sqrt(1 - rho^2) * e)`,

where e is independent standard normal error. So `rho` is the
pretest-posttest correlation within each pretested group, and `sigma` is
the standard deviation of both scores within each group. The four groups
are generated in the package's order (pretested treatment, pretested
control, unpretested treatment, unpretested control), and the baselines
are drawn before the errors.

**Relation to the package's other simulations.** With `mean = 0`,
`sigma = 1`, and `pretest_effect = 0`, the data are the first data set
that
[`power_solomon()`](https://juhalt.github.io/solomonR/reference/power_solomon.md)
simulates with the same seed. With `n = 30`, `delta = 5`,
`pretest_effect = 2`, `rho = 0.6`, `sigma = 10`, `mean = 50`,
`seed = 20260915`, `digits = 0`, and `limits = c(0, 100)`, they are the
bundled
[solomon_example](https://juhalt.github.io/solomonR/reference/solomon_example.md)
data.

**True values.** The `"truth"` attribute gives each estimand:

- the treatment effect among unpretested participants, `delta`;

- among pretested participants, `delta + sens`;

- the Pretest x Treatment interaction, `sens`;

- the equal-weighted average treatment effect, `delta + sens / 2`;

- the pretest effect among controls, `pretest_effect`.

Rounding and limits (`digits`, `limits`) make the scores look like test
scores, but they shift the true values slightly, and more so when many
scores reach a limit.

**Choosing values.** Willson and Putnam (1982) found an average pretest
effect of about a fifth of a standard deviation in randomized studies;
the article "Planning a Solomon Study" discusses planning values.

## Designs with several treatments

A Solomon N-group design crosses k treatments and a control with
pretesting, giving 2(k + 1) groups: six for two treatments and eight for
three (Edmonds & Kennedy, 2017; Steyn, 2009). Give `delta` one effect
per treatment, named by treatment, for example
`delta = c(A = 0.5, B = 0.2)`. Without names the treatments are called
`"T1"`, `"T2"`, and so on; the control is `"Control"`. `sens` is one
number for every treatment or one per treatment, in the order of
`delta`; a named `sens` is matched to the treatments by name.

The model is the one above, with the effects `delta[j]` and `sens[j]` of
treatment j: each treatment differs from the control only by its own
effects. `n` is one size for every group, or 2(k + 1) sizes in this
order: the pretested treatments (in the order of `delta`), the pretested
control, the unpretested treatments, and the unpretested control. Sizes
named by group, as in the `"settings"` attribute (such as
`"Pretested, A"`), or named `n1`, `n2`, and so on for the groups in that
order, may be given in any order. Other names are refused. The groups
are generated in the order above, and the baselines are drawn before the
errors.

The `treat` column is then a factor whose levels are the treatments and
`"Control"`, so the data go directly to
`fit_solomon_glm(y_post, treat, pretested, y_pre, control = "Control")`.
The `"truth"` attribute has one row for each treatment-control
comparison and contrast, with the columns `comparison` (such as
`"A vs Control"`), `contrast`, and `true_value`, in the order of the
`effects` table of that fit; a last row gives the pretest effect among
controls. With one treatment effect the data and the `"truth"` attribute
are those described above.

## References

Edmonds, W. A., & Kennedy, T. D. (2017). *An applied guide to research
designs: Quantitative, qualitative, and mixed methods* (2nd ed.). SAGE
Publications. https://doi.org/10.4135/9781071802779

Steyn, R. (2009). Re-designing the Solomon four-group: Can we improve on
this exemplary model? *Design Principles and Practices: An International
Journal—Annual Review, 3*(1), 383–394.
https://doi.org/10.18848/1833-1874/CGP/v03i01/37588

Willson, V. L., & Putnam, R. R. (1982). A meta-analysis of pretest
sensitization effects in experimental design. *American Educational
Research Journal, 19*(2), 249–258.
https://doi.org/10.3102/00028312019002249

## See also

[`power_solomon()`](https://juhalt.github.io/solomonR/reference/power_solomon.md)
and
[`plan_solomon()`](https://juhalt.github.io/solomonR/reference/plan_solomon.md)
for planning,
[`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md)
for the recommended analysis, and the article "Teaching with solomonR".

## Examples

``` r
d <- simulate_solomon(n = 40, delta = 0.5, sens = 0.4, rho = 0.6, seed = 1)
attr(d, "truth")
#>                   estimand true_value
#> 1   ATE (avg over pretest)        0.7
#> 2      Pretest x Treatment        0.4
#> 3    Treatment | pretested        0.9
#> 4  Treatment | unpretested        0.5
#> 5 Pretest effect | control        0.0
with(d, fit_solomon_glm(y_post, treat, pretested, y_pre))
#> Solomon GLM (unified model)
#> Formula: y ~ treat * pretested + pre_obs
#> Covariance: HC3 heteroskedasticity-consistent; t tests (df = 155)
#> 
#> Term             Est (SE)            t   df      p           95% CI
#> (Intercept)      -0.045 (0.136)  -0.33  155  0.739  [-0.314, 0.223]
#> treat            0.656 (0.204)    3.22  155  0.002   [0.254, 1.058]
#> pretested        -0.011 (0.200)  -0.05  155  0.958  [-0.406, 0.384]
#> pre_obs          0.756 (0.125)    6.05  155  <.001   [0.509, 1.002]
#> treat:pretested  0.413 (0.280)    1.48  155  0.142  [-0.140, 0.966]
#> 
#> Key contrasts            Est (SE)          t   df      p           95% CI  Wald R2
#> ATE (avg over pretest)   0.862 (0.140)  6.16  155  <.001   [0.586, 1.139]    0.197
#> Pretest x Treatment      0.413 (0.280)  1.48  155  0.142  [-0.140, 0.966]    0.014
#> Treatment | pretested    1.069 (0.192)  5.56  155  <.001   [0.689, 1.449]    0.166
#> Treatment | unpretested  0.656 (0.204)  3.22  155  0.002   [0.254, 1.058]    0.063
#> 
#> Wald R2: partial R-squared for conventional Gaussian OLS;
#> a Wald-based descriptive approximation when robust covariance is used.

# The bundled teaching data.
ex <- simulate_solomon(n = 30, delta = 5, pretest_effect = 2, rho = 0.6,
                       sigma = 10, mean = 50, digits = 0, limits = c(0, 100),
                       seed = 20260915)
attr(ex, "truth") <- attr(ex, "settings") <- NULL
identical(ex, solomon_example)
#> [1] TRUE

# A six-group design: two treatments and a control.
d6 <- simulate_solomon(n = 50, delta = c(A = 0.5, B = 0.2), sens = c(0.3, 0),
                       rho = 0.6, seed = 1)
attr(d6, "truth")
#>     comparison                 contrast true_value
#> 1 A vs Control   ATE (avg over pretest)       0.65
#> 2 B vs Control   ATE (avg over pretest)       0.20
#> 3 A vs Control      Pretest x Treatment       0.30
#> 4 B vs Control      Pretest x Treatment       0.00
#> 5 A vs Control    Treatment | pretested       0.80
#> 6 B vs Control    Treatment | pretested       0.20
#> 7 A vs Control  Treatment | unpretested       0.50
#> 8 B vs Control  Treatment | unpretested       0.20
#> 9      Control Pretest effect | control       0.00
fit_solomon_glm(y_post, treat, pretested, y_pre, control = "Control", data = d6)
#> Solomon GLM (unified model), N-group design
#> Conditions: A, B; control: Control. With and without a pretest: 6 groups.
#> Formula: y ~ (treat_A + treat_B) * pretested + pre_obs
#> Covariance: HC3 heteroskedasticity-consistent; t tests (df = 293)
#> 
#> Omnibus tests
#> Test                                  Statistic      p
#> Condition (avg over pretest)  F(2, 293) = 15.00  <.001
#> Pretest x Condition            F(2, 293) = 1.15  0.318
#> Condition | pretested         F(2, 293) = 17.50  <.001
#> Condition | unpretested        F(2, 293) = 3.07  0.048
#> 
#> Contrasts
#> Comparison    Contrast                      Est (SE)     t   df      p  p adj.           95% CI
#> A vs Control  ATE (avg over pretest)   0.767 (0.144)  5.32  293  <.001   <.001   [0.483, 1.052]
#> B vs Control  ATE (avg over pretest)   0.197 (0.130)  1.52  293  0.130   0.130  [-0.059, 0.453]
#> A vs Control  Pretest x Treatment      0.438 (0.289)  1.52  293  0.130   0.261  [-0.130, 1.006]
#> B vs Control  Pretest x Treatment      0.194 (0.260)  0.74  293  0.457   0.457  [-0.318, 0.705]
#> A vs Control  Treatment | pretested    0.986 (0.174)  5.66  293  <.001   <.001   [0.644, 1.329]
#> B vs Control  Treatment | pretested    0.294 (0.170)  1.73  293  0.084   0.084  [-0.040, 0.628]
#> A vs Control  Treatment | unpretested  0.549 (0.230)  2.38  293  0.018   0.036   [0.096, 1.002]
#> B vs Control  Treatment | unpretested  0.100 (0.197)  0.51  293  0.611   0.611  [-0.287, 0.488]
#> 
#> p adj.: adjusted by Holm's (1979) procedure within each contrast, across the 2 comparisons.
#> Confidence intervals are not adjusted.
#> 
#> Experimental: in the package's simulation study (issue #45), the omnibus
#> tests of Condition | pretested and Condition | unpretested rejected in up
#> to 6.9% of replications at the .05 level with three treatments and 10
#> participants per group. No other test, and no family of adjusted
#> comparisons, failed the study's rule. See ?fit_solomon_glm.
```
