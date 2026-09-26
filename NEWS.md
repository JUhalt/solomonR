# solomonR (development version)

## Binary outcomes (#43)

* New `marginal_solomon()` estimates the Solomon contrasts for binary
  outcomes as risk differences, risk ratios, or odds ratios from marginal
  (standardized) risks, following Daniel et al. (2021) and Localio et al.
  (2007). Intervals come from a bootstrap that resamples within the four
  Solomon cells, or from the delta method.
* New `fisher_solomon()` reproduces the historical categorical analysis of
  El Karkri et al. (2025b), with a caution that its rule compares
  significance, not effects.
* `fit_solomon_glm()` now warns (`solomonR_noncollapsible_warning`) when a
  noncollapsible link such as the logit is combined with `pretest_score`:
  its Pretest x Treatment contrast then compares a conditional with a
  marginal effect and is nonzero without sensitization.
* New article "Binary Outcomes: Validating marginal_solomon()" reports the
  simulation study registered on #43 (48 scenarios, 2,000 replications
  each), and the "Validation Evidence" tables include it. Risk differences
  are validated across the supported range; risk ratios and odds ratios are
  conservative with 20 to 50 participants per cell.

## Clustered designs (#46)

* `validate_solomon()` gains a `cluster` argument. It reports the number of
  clusters and cluster sizes in each cell, and whether treatment and
  pretesting were assigned to whole clusters or within them. A cell made up
  of a single cluster is an error, because cluster and condition are then
  completely confounded, as when each Solomon condition is one intact class.
* `fit_solomon_glm(robust = "CR2")` refuses such designs with a classed
  error (`solomonR_confounded_clusters`) instead of reporting cluster-robust
  standard errors that cannot be estimated.
* When whole clusters are randomized, a cell with two or three clusters is a
  warning, following the rule of thumb that four clusters per arm is an
  absolute minimum (Hayes & Moulton, 2017, p. 128).
* CR2 fits give a classed warning (`solomonR_small_df_warning`) when a
  Solomon contrast has Satterthwaite degrees of freedom below 4, where Tipton
  (2015) advises that p-values not be trusted.

## Sensitization figure for maximum-likelihood fits (#47)

* `plot_sensitization()` now accepts `fit_solomon_ml()` fits. Intervals for
  the adjusted cell means use the fit's own inference: the normal reference
  under the default Wald inference, or Welch-Satterthwaite t under
  `inference = "satterthwaite"`. `fit_solomon_ml()` now stores what those
  intervals need; its estimates, standard errors, and intervals are
  unchanged.
* The figure's caption now puts the adjustment and the interval method on
  separate lines, so it is no longer cut off at common figure widths.

## Attribution (#42)

* Every reference is now in APA Style (7th ed.) with its DOI, verified
  against Crossref, and a new "References and the Solomon Literature" article
  is the package's canonical reference list. It adds the Solomon literature
  the package draws on or plans to, each with a note on its contribution.
* Corrected the attribution of the 1988 meta-analytic procedure (Test I) to
  Walton Braver and Braver (1988), as the authors cite their own article.
  Earlier documentation and output wrote "Braver & Braver (1988)" and
  "Braver, M. W.". Printed labels change accordingly, including the `test`
  column of `power_solomon()` and the committed validation tables; no
  numeric result changes.
* The roadmap and README now link each planned item to its issue, with a new
  v0.7.0 milestone for longitudinal and quasi-experimental designs.

## Validation evidence

* The "Validation Evidence" article and its shared tables now include a
  re-simulation check of `plan_solomon()` (#24). Fifty-six analytically
  planned designs were re-simulated with the package's GLM test, with no
  failed fits. From 30 participants per cell, 37 of 39 plans reached the 80%
  target within 0.02, and the two exceptions were above it. A new article,
  "Checking plan_solomon()", reports every plan.

# solomonR 0.4.0

## Validation evidence in one format (#11)

* New article "Validation Evidence" gathers the package's simulation studies
  in one place, with each study's protocol, scenarios, replications, methods,
  failed fits, and results against its pre-registered tolerances.
* Two shared tables, rebuilt by a committed script from each study's results,
  can be downloaded from the repository. `studies.csv` has one row per study
  and `benchmarks.csv` has one row per scenario, method, estimand, and
  performance measure, with Monte Carlo standard errors. Failed fits are
  counted separately, and measures a method cannot have are listed with a
  reason rather than omitted.

## Power curves (#30)

* New `plot_power_solomon()` draws power against the size of the smallest
  cell for each Solomon estimand, with the target power marked and the
  design `plan_solomon()` returns for that target. Curves use the validated
  normal-theory power by default; `method = "simulation"` uses the package's
  GLM test through `power_solomon()` and shows two-MCSE bands. One of
  `delta`, `sens`, or `rho` may vary across curves.

## Design, change, and historical-path figures (#25, #28, #29)

* New `plot_solomon_design()` draws the four-group design in Campbell and
  Stanley's (1963) notation. Without data it gives the teaching schematic;
  with data or a GLM fit it labels each group with its size and posttest
  mean and flags empty or sparse groups, using the rule in
  `validate_solomon()`.
* New `plot_solomon_change()` shows pretest-to-posttest change for the
  pretested groups, with t intervals, beside posttest means for the
  unpretested groups, which are labeled as unpretested by design.
  Pretested participants with incidentally missing pretests are excluded
  from the trajectories and counted in the caption, never imputed.
* New `plot_classic_flow()` draws the historical Tests A-I sequence as a
  decision tree. Given a `fit_solomon_classic()` result, it highlights the
  path taken with each test's p-value. Every version carries the caution
  that the conditional sequence inflates Type I error (Sawilowsky et al.,
  1994).

## Pretest sensitization figure (#26)

* New `plot_sensitization()` draws the Pretest x Treatment interaction as
  model-adjusted group means with confidence intervals. Pretested groups are
  evaluated at the mean pretest among pretested participants and unpretested
  groups without a pretest, so the difference of differences among the
  plotted means equals the fitted sensitization contrast exactly. Observed
  means are overlaid for comparison.
* Its intervals use the fit's own covariance matrix and reference
  distribution, the subtitle reports the fitted contrast, and optional bounds
  add the outcome of `equivalence_solomon()`.
* `plot_solomon()` and `plot_solomon_gg()` now use t intervals with n - 1
  degrees of freedom for cell means instead of the normal quantile, which was
  too narrow with small cells.

## Forest plot of the Solomon contrasts (#27)

* New `plot_solomon_effects()` draws the four Solomon contrasts with their
  confidence intervals from a fit by `fit_solomon_glm()`, `fit_solomon_ml()`,
  `fit_solomon_sem()`, or `fit_solomon_sem_latent()`. Estimates and intervals
  are taken unchanged from the fitted object, and the caption states the
  confidence level, the inference used, and the reference distribution.
* Contrasts a model does not estimate are named in the caption rather than
  drawn as zero, and optional equivalence bounds are shaded on the
  sensitization row, matching `equivalence_solomon()`.

## Sample-size planning (#24)

* New `plan_solomon()` finds the smallest Solomon design, at a chosen
  allocation across the four cells, whose power for each estimand reaches a
  target. Analytic planning uses the normal-theory benchmarks validated in
  #18 and returns the exact minimum. `method = "simulation"` plans for the
  package's own GLM test through `power_solomon()`, which matters with small
  cells, where HC3 standard errors are conservative.
* The help page shows that the sensitization contrast always has four times
  the sampling variance of the average treatment effect, so detecting
  sensitization as large as the average effect needs about four times as many
  participants.

## Power simulation rebuilt and validated (#18)

* `power_solomon()` is rebuilt. It reports the rejection rate of each Solomon
  test with its Monte Carlo standard error, names the estimand, the test, and
  the true effect behind every row, and gains `alpha` and `seed` arguments.
  The global random number state is restored after use.
* The simulator draws a latent baseline that only pretested participants
  observe, so structural pretest absence is part of the design, and `sigma`
  applies to every cell.
* The experimental warning is removed. A pre-specified simulation study of
  126 scenarios and 315,000 replications, with the protocol and its amendment
  posted on #18 before any results were examined, found exact agreement with
  the analytic benchmark for the 2x2 ANOVA interaction, nominal size for
  Test I under the complete null, exact scale invariance, and no fit
  failures.
* Rejection rates from the unified GLM are conservative with small cells,
  following the HC3 standard errors used by default: Type I error averaged
  0.041 with 10 participants per cell and 0.049 with 100. The help page and
  the new article "Validating power_solomon()" report the findings.
* Test I rows report `NA` rather than zero power when `stouffer = FALSE`.

# solomonR 0.3.0

## Inference corrections

* `fit_solomon_glm(robust = "CR2")` now uses Satterthwaite degrees of freedom
  for coefficient and contrast tests (Pustejovsky & Tipton, 2018), matching
  clubSandwich. The previous normal-reference CR2 tests over-rejected with few
  clusters (#14). Coefficient and contrast tables gain a `df` column (`Inf` for
  normal-reference tests).
* `fit_solomon_glm(robust = "CR2")` no longer fails when outcomes are missing;
  the clustering variable is aligned with the rows used by the model (#14).
* `perm_solomon()` refuses fits that include a clustering variable, because
  permuting individuals is not a valid randomization test when treatment was
  assigned to clusters (#15). Cluster-level randomization inference is planned
  (#19).
* `fit_solomon_sem_latent()` fixes a reference latent mean (U0; P0 in the
  latent ANCOVA) so the latent mean structure is identified. Solomon contrasts
  are unchanged; model degrees of freedom and fit indices are now correct
  (#16).

## Confidence intervals and reference distributions (#8)

* `fit_solomon_glm()` now defaults to `robust = "HC3"`, following Long & Ervin
  (2000) and Hayes & Cai (2007).
* Tests and intervals share one reference distribution: t with residual
  degrees of freedom for Gaussian GLMs (conventional and HC3), the normal
  distribution for binomial and Poisson models, and Satterthwaite t for CR2.
* Effect summaries from `fit_solomon_glm()`, `fit_solomon_classic()` (Tests
  A-H), `fit_solomon_ml()`, `fit_solomon_sem()`, and `fit_solomon_sem_latent()`
  gain `conf.low` and `conf.high` columns and a `conf_level` argument. ML and
  SEM intervals are large-sample Wald intervals.
* The Wald partial R-squared has a noncentral F confidence interval for
  conventional Gaussian fits (Steiger, 2004); none is reported under robust
  covariance.
* The Groups 3-4 Hedges' g interval now uses the noncentral t method (Cumming
  & Finch, 2001; Kelley, 2007) instead of a normal approximation.
* Printed output shows each interval with its confidence level.

## Teaching data and getting-started guide (#17, #21)

* New `solomon_example` data set: a simulated Solomon four-group study with 30
  participants per group, generated by a documented script with a seed and
  parameters fixed in advance, so users can compare estimates with known true
  effects. It is now the example in the README, help pages, and vignettes.
* `solomon_demo` is kept as a second example. Its documentation now describes
  the negative pretest-posttest correlation and chance pretest imbalance that
  make unadjusted, ANCOVA, and gain-score estimates disagree.
* New article `vignette("getting-started")` follows `solomon_example` from design
  checks through the historical workflow, the recommended analysis, an
  equivalence test for sensitization, and a method comparison to a reporting
  example, with an annotated reading list.

## Small-sample inference for maximum likelihood (#22)

* `fit_solomon_ml()` gains `inference = c("wald", "satterthwaite")`. The
  default keeps van Engelenburg's (1999) Wald inference. The small-sample
  option keeps the maximum-likelihood point estimates but uses unbiased
  residual variances within each pretest condition, t tests within a
  condition, and Welch-Satterthwaite degrees of freedom for contrasts that
  combine conditions (Satterthwaite, 1946; Welch, 1947). Coefficient and
  effect tables gain a `df` column.
* When the smallest cell has fewer than 40 participants and `inference` is
  not supplied, `fit_solomon_ml()` warns with class
  `solomonR_small_sample_warning`, and printed output notes the small cells.
  Supplying `inference = "wald"` keeps the default without the warning.
* `compare_solomon_methods()` reports both inference options for maximum
  likelihood.

## Maximum-likelihood validation (#10, #22)

* A pre-specified simulation study (Morris, White & Crowther, 2019), extended
  for #22 to 84 scenarios with 2,000 replications each, validated estimand
  recovery by `fit_solomon_ml()`. Its default Wald intervals were too narrow
  in small samples (mean coverage 0.893 with 6 participants per cell and
  0.936 with 20; sensitization Type I error 0.099 with 6 per cell). The
  small-sample option was calibrated at every cell size studied (mean
  coverage 0.949 to 0.950; Type I error 0.050 to 0.053). The unified GLM with
  HC3 was conservative at 10 or fewer per cell and close to nominal from 20.
* The warning threshold of 40 participants per cell comes from a rule posted
  on #22 before the extended results were examined.
* Results, script, and scenario definitions are in the article "Validating
  fit_solomon_ml()", and the help pages report the findings.

## Method guide (#9)

* New article `vignette("solomon-methods")` labels each analysis as a
  historical procedure, contemporary recommendation, published Solomon
  proposal, or solomonR extension, and states its estimand, assumptions,
  limitations, and sources. Historical claims were checked against the
  sources, including Campbell and Stanley (1963), Braver and Braver (1988),
  and Sawilowsky et al. (1994).

## Method comparison (#7)

* New `compare_solomon_methods()` fits the unified GLM, maximum likelihood,
  classic Tests A-F and H, and SEM to the same data and aligns their estimates
  of the four Solomon contrasts with each method's pretest adjustment,
  variance assumption, reference distribution, and interval. Analyses that do
  not estimate a raw-scale contrast (permutation tests, Test I, latent SEM,
  Hedges' g) are listed separately with the reason (Lin, 2013; Lundberg et al.,
  2021).

## Sensitization equivalence testing (#4)

* New `equivalence_solomon()` performs a two one-sided tests (TOST)
  equivalence test for a Solomon contrast, by default the Pretest x Treatment
  sensitization contrast (Schuirmann, 1987; Lakens, 2017). It uses the fitted
  model's reference distribution, reports both one-sided tests, the
  1 - 2 alpha interval, and the test against zero, and classifies the result
  as equivalent, trivially small, different, or inconclusive. Bounds have no
  default and are documented as a prespecified smallest effect size of
  interest.

## Design validation and missingness (#5, #6)

* New `validate_solomon()` checks input lengths, 0/1 coding, the presence of
  all four cells, and observed outcomes per cell, and returns every problem
  as an error, warning, or note instead of stopping at the first.
* New `check_solomon_missing()` separates structurally absent pretests in the
  unpretested groups from incidental pretest and posttest missingness,
  unexpected pretest scores, and unassigned participants. It reports counts by
  cell and cited guidance for each category (Solomon, 1949; Rubin, 1976;
  Graham et al., 2006; White & Thompson, 2005; Groenwold et al., 2012; Little &
  Rubin, 2019).

## Safeguards

* `fit_solomon_glm()`, `fit_solomon_sem()`, `fit_solomon_sem_latent()`, and
  `check_solomon_assumptions()` require `treat` and `pretested` to be coded 0/1
  (or logical) and inputs to have equal lengths, instead of returning empty or
  rescaled results (#5).
* `fit_solomon_glm()` warns when pretested participants are missing pretest
  scores (they are excluded, not imputed) and when pretest scores are supplied
  for unpretested participants (#5, #6).
* `power_solomon()` now warns that it is experimental. Interim fixes: `sigma`
  applies to every cell, the Stouffer arm uses the historical one-tailed
  direction, disabled metrics are `NA`, and `delta` is documented as the effect
  among unpretested participants. A validated rebuild is planned (#18).
* `perm_solomon(seed = )` restores the global random-number state on exit.

## Output and documentation

* Printed GLM results name the covariance estimator and reference
  distribution. The summary no longer labels conventional standard errors as
  robust, and `print()` returns the fitted object invisibly.
* `plot_perm()` reports the same corrected p-value as `perm_solomon()`.
* `check_solomon_assumptions()` presents its results as descriptive
  diagnostics rather than gates for choosing an analysis, and uses complete
  pretested cases for the slope-homogeneity test.
* Help pages now cite the methodological sources for each procedure (#9).
* Test I is labeled as the Braver & Braver (1988) Stouffer combination in
  output and help pages, and its result gains a `procedure` column.
* `inst/CITATION` uses `bibentry()`.
* Removed the unused internal `spr2_ci()`, which relied on an invalid interval
  transformation, and duplicated print helpers.

## Development

* Set the v0.3.0 scope around graduate students and applied researchers and
  added v0.4.0, v0.5.0, v0.6.0, and v1.0.0 milestones that mirror the roadmap.

* solomonR 0.3.0 and later is licensed under GNU GPL version 3 only (GPL-3).
  Version 0.2.0 and earlier retain their original MIT terms and notices.
* Aligned development citation metadata, license pages, and roadmap/issue links.
* Distinguished completed implementation, deferred work, and research proposals
  in the roadmap without changing analytical behavior.
* Formatted `LICENSE.md` as markdown so the pkgdown license page has a title
  and readable headings; the license wording is unchanged.
* Added `Language: en-US` and a spelling word list (`inst/WORDLIST`) so
  documentation spell checks report only genuine errors.

* Began development toward solomonR 0.3.0.
* Added R-universe distribution and stable-release installation instructions.
* Reconciled the development roadmap with functionality delivered in v0.2.0.

# solomonR 0.2.0

## Historical Solomon workflow

* Added the full historical Tests A-I workflow in `fit_solomon_classic()`.
* Added Tests B and C for treatment simple effects.
* Added Test D as the equal-weighted treatment effect across pretest conditions.
* Added ANCOVA, gain-score, repeated-measures, and posttest-only historical analyses.
* Added historical decision-path reporting with `[PATH]` markers.
* Added directional Stouffer Test I for historical teaching and replication.
* Added explicit cautions regarding later Type I error critiques of conditional
  Solomon testing procedures.

## Modern observed-variable analysis

* Expanded `fit_solomon_glm()` as the primary unified Solomon analysis.
* Corrected the equal-weighted average treatment effect to
  `treat + 0.5 * treat:pretested`.
* Added Solomon-specific contrasts for:
  * average treatment effect,
  * pretest-by-treatment sensitization,
  * treatment effect among pretested participants,
  * treatment effect among unpretested participants.
* Preserved all four Solomon groups when pretest scores are structurally absent
  in the unpretested groups.
* Added HC3 heteroskedasticity-robust covariance estimation.
* Added CR2 cluster-robust covariance estimation.
* Added Wald-based partial R-squared reporting for Gaussian models.

## Randomization inference

* Rebuilt `perm_solomon()` as a stratified treatment-label permutation test.
* Treatment labels are permuted within pretest-assignment strata.
* Added HC3-studentized permutation statistics.
* Added finite Monte Carlo p-value correction.
* Added concise printing for permutation objects.
* Added permutation-distribution plotting with `plot_perm()`.

## Maximum likelihood

* Added `fit_solomon_ml()` for full-information maximum-likelihood analysis.
* Structural absence of pretest scores in the unpretested groups is handled as
  part of the Solomon design rather than as ordinary missing baseline data.
* Added Solomon-specific ATE, sensitization, and simple treatment effects.
* Added automated validation against corresponding component regressions.

## Structural equation models

* Added observed-variable Solomon SEM with `fit_solomon_sem()`.
* Added an ANCOVA-style SEM pathway for the pretested groups.
* Added latent-variable Solomon SEM with `fit_solomon_sem_latent()`.
* Added measurement-invariance controls for latent outcomes.
* Latent Solomon mean contrasts now require scalar measurement invariance.
* Added explicit model-convergence checks.
* Clarified that global fit indices are not diagnostic for the saturated
  four-group observed mean model.

## Diagnostics and effect sizes

* Added Brown-Forsythe variance checks.
* Added cell-level normality summaries.
* Added ANCOVA slope-homogeneity diagnostics.
* Added Hedges' g for the posttest-only comparison.
* Corrected and clarified Wald-based partial R-squared reporting.

## Documentation and reliability

* Expanded automated tests for historical, GLM, permutation, ML, SEM, printing,
  and plotting workflows.
* Added known-estimand tests for Solomon contrasts.
* Rebuilt the classic and modern analysis vignettes.
* Expanded the package README with installation, method-selection, and
  interpretation guidance.
* Added APA-oriented rounding and concise public-facing output.
* Package checks pass with zero errors, warnings, or notes.
