# Changelog

## solomonR (development version)

### Sample-size planning ([\#24](https://github.com/JUhalt/solomonR/issues/24))

- New
  [`plan_solomon()`](https://juhalt.github.io/solomonR/reference/plan_solomon.md)
  finds the smallest Solomon design, at a chosen allocation across the
  four cells, whose power for each estimand reaches a target. Analytic
  planning uses the normal-theory benchmarks validated in
  [\#18](https://github.com/JUhalt/solomonR/issues/18) and returns the
  exact minimum. `method = "simulation"` plans for the package’s own GLM
  test through
  [`power_solomon()`](https://juhalt.github.io/solomonR/reference/power_solomon.md),
  which matters with small cells, where HC3 standard errors are
  conservative.
- The help page shows that the sensitization contrast always has four
  times the sampling variance of the average treatment effect, so
  detecting sensitization as large as the average effect needs about
  four times as many participants.

### Power simulation rebuilt and validated ([\#18](https://github.com/JUhalt/solomonR/issues/18))

- [`power_solomon()`](https://juhalt.github.io/solomonR/reference/power_solomon.md)
  is rebuilt. It reports the rejection rate of each Solomon test with
  its Monte Carlo standard error, names the estimand, the test, and the
  true effect behind every row, and gains `alpha` and `seed` arguments.
  The global random number state is restored after use.
- The simulator draws a latent baseline that only pretested participants
  observe, so structural pretest absence is part of the design, and
  `sigma` applies to every cell.
- The experimental warning is removed. A pre-specified simulation study
  of 126 scenarios and 315,000 replications, with the protocol and its
  amendment posted on
  [\#18](https://github.com/JUhalt/solomonR/issues/18) before any
  results were examined, found exact agreement with the analytic
  benchmark for the 2x2 ANOVA interaction, nominal size for Test I under
  the complete null, exact scale invariance, and no fit failures.
- Rejection rates from the unified GLM are conservative with small
  cells, following the HC3 standard errors used by default: Type I error
  averaged 0.041 with 10 participants per cell and 0.049 with 100. The
  help page and the new article “Validating power_solomon()” report the
  findings.
- Test I rows report `NA` rather than zero power when
  `stouffer = FALSE`.

## solomonR 0.3.0

### Inference corrections

- `fit_solomon_glm(robust = "CR2")` now uses Satterthwaite degrees of
  freedom for coefficient and contrast tests (Pustejovsky & Tipton,
  2018), matching clubSandwich. The previous normal-reference CR2 tests
  over-rejected with few clusters
  ([\#14](https://github.com/JUhalt/solomonR/issues/14)). Coefficient
  and contrast tables gain a `df` column (`Inf` for normal-reference
  tests).
- `fit_solomon_glm(robust = "CR2")` no longer fails when outcomes are
  missing; the clustering variable is aligned with the rows used by the
  model ([\#14](https://github.com/JUhalt/solomonR/issues/14)).
- [`perm_solomon()`](https://juhalt.github.io/solomonR/reference/perm_solomon.md)
  refuses fits that include a clustering variable, because permuting
  individuals is not a valid randomization test when treatment was
  assigned to clusters
  ([\#15](https://github.com/JUhalt/solomonR/issues/15)). Cluster-level
  randomization inference is planned
  ([\#19](https://github.com/JUhalt/solomonR/issues/19)).
- [`fit_solomon_sem_latent()`](https://juhalt.github.io/solomonR/reference/fit_solomon_sem_latent.md)
  fixes a reference latent mean (U0; P0 in the latent ANCOVA) so the
  latent mean structure is identified. Solomon contrasts are unchanged;
  model degrees of freedom and fit indices are now correct
  ([\#16](https://github.com/JUhalt/solomonR/issues/16)).

### Confidence intervals and reference distributions ([\#8](https://github.com/JUhalt/solomonR/issues/8))

- [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md)
  now defaults to `robust = "HC3"`, following Long & Ervin
  2000. and Hayes & Cai (2007).
- Tests and intervals share one reference distribution: t with residual
  degrees of freedom for Gaussian GLMs (conventional and HC3), the
  normal distribution for binomial and Poisson models, and Satterthwaite
  t for CR2.
- Effect summaries from
  [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md),
  [`fit_solomon_classic()`](https://juhalt.github.io/solomonR/reference/fit_solomon_classic.md)
  (Tests A-H),
  [`fit_solomon_ml()`](https://juhalt.github.io/solomonR/reference/fit_solomon_ml.md),
  [`fit_solomon_sem()`](https://juhalt.github.io/solomonR/reference/fit_solomon_sem.md),
  and
  [`fit_solomon_sem_latent()`](https://juhalt.github.io/solomonR/reference/fit_solomon_sem_latent.md)
  gain `conf.low` and `conf.high` columns and a `conf_level` argument.
  ML and SEM intervals are large-sample Wald intervals.
- The Wald partial R-squared has a noncentral F confidence interval for
  conventional Gaussian fits (Steiger, 2004); none is reported under
  robust covariance.
- The Groups 3-4 Hedges’ g interval now uses the noncentral t method
  (Cumming & Finch, 2001; Kelley, 2007) instead of a normal
  approximation.
- Printed output shows each interval with its confidence level.

### Teaching data and getting-started guide ([\#17](https://github.com/JUhalt/solomonR/issues/17), [\#21](https://github.com/JUhalt/solomonR/issues/21))

- New `solomon_example` data set: a simulated Solomon four-group study
  with 30 participants per group, generated by a documented script with
  a seed and parameters fixed in advance, so users can compare estimates
  with known true effects. It is now the example in the README, help
  pages, and vignettes.
- `solomon_demo` is kept as a second example. Its documentation now
  describes the negative pretest-posttest correlation and chance pretest
  imbalance that make unadjusted, ANCOVA, and gain-score estimates
  disagree.
- New article
  [`vignette("getting-started")`](https://juhalt.github.io/solomonR/articles/getting-started.md)
  follows `solomon_example` from design checks through the historical
  workflow, the recommended analysis, an equivalence test for
  sensitization, and a method comparison to a reporting example, with an
  annotated reading list.

### Small-sample inference for maximum likelihood ([\#22](https://github.com/JUhalt/solomonR/issues/22))

- [`fit_solomon_ml()`](https://juhalt.github.io/solomonR/reference/fit_solomon_ml.md)
  gains `inference = c("wald", "satterthwaite")`. The default keeps van
  Engelenburg’s (1999) Wald inference. The small-sample option keeps the
  maximum-likelihood point estimates but uses unbiased residual
  variances within each pretest condition, t tests within a condition,
  and Welch-Satterthwaite degrees of freedom for contrasts that combine
  conditions (Satterthwaite, 1946; Welch, 1947). Coefficient and effect
  tables gain a `df` column.
- When the smallest cell has fewer than 40 participants and `inference`
  is not supplied,
  [`fit_solomon_ml()`](https://juhalt.github.io/solomonR/reference/fit_solomon_ml.md)
  warns with class `solomonR_small_sample_warning`, and printed output
  notes the small cells. Supplying `inference = "wald"` keeps the
  default without the warning.
- [`compare_solomon_methods()`](https://juhalt.github.io/solomonR/reference/compare_solomon_methods.md)
  reports both inference options for maximum likelihood.

### Maximum-likelihood validation ([\#10](https://github.com/JUhalt/solomonR/issues/10), [\#22](https://github.com/JUhalt/solomonR/issues/22))

- A pre-specified simulation study (Morris, White & Crowther, 2019),
  extended for [\#22](https://github.com/JUhalt/solomonR/issues/22) to
  84 scenarios with 2,000 replications each, validated estimand recovery
  by
  [`fit_solomon_ml()`](https://juhalt.github.io/solomonR/reference/fit_solomon_ml.md).
  Its default Wald intervals were too narrow in small samples (mean
  coverage 0.893 with 6 participants per cell and 0.936 with 20;
  sensitization Type I error 0.099 with 6 per cell). The small-sample
  option was calibrated at every cell size studied (mean coverage 0.949
  to 0.950; Type I error 0.050 to 0.053). The unified GLM with HC3 was
  conservative at 10 or fewer per cell and close to nominal from 20.
- The warning threshold of 40 participants per cell comes from a rule
  posted on [\#22](https://github.com/JUhalt/solomonR/issues/22) before
  the extended results were examined.
- Results, script, and scenario definitions are in the article
  “Validating fit_solomon_ml()”, and the help pages report the findings.

### Method guide ([\#9](https://github.com/JUhalt/solomonR/issues/9))

- New article
  [`vignette("solomon-methods")`](https://juhalt.github.io/solomonR/articles/solomon-methods.md)
  labels each analysis as a historical procedure, contemporary
  recommendation, published Solomon proposal, or solomonR extension, and
  states its estimand, assumptions, limitations, and sources. Historical
  claims were checked against the sources, including Campbell and
  Stanley (1963), Braver and Braver (1988), and Sawilowsky et
  al. (1994).

### Method comparison ([\#7](https://github.com/JUhalt/solomonR/issues/7))

- New
  [`compare_solomon_methods()`](https://juhalt.github.io/solomonR/reference/compare_solomon_methods.md)
  fits the unified GLM, maximum likelihood, classic Tests A-F and H, and
  SEM to the same data and aligns their estimates of the four Solomon
  contrasts with each method’s pretest adjustment, variance assumption,
  reference distribution, and interval. Analyses that do not estimate a
  raw-scale contrast (permutation tests, Test I, latent SEM, Hedges’ g)
  are listed separately with the reason (Lin, 2013; Lundberg et al.,
  2021).

### Sensitization equivalence testing ([\#4](https://github.com/JUhalt/solomonR/issues/4))

- New
  [`equivalence_solomon()`](https://juhalt.github.io/solomonR/reference/equivalence_solomon.md)
  performs a two one-sided tests (TOST) equivalence test for a Solomon
  contrast, by default the Pretest x Treatment sensitization contrast
  (Schuirmann, 1987; Lakens, 2017). It uses the fitted model’s reference
  distribution, reports both one-sided tests, the 1 - 2 alpha interval,
  and the test against zero, and classifies the result as equivalent,
  trivially small, different, or inconclusive. Bounds have no default
  and are documented as a prespecified smallest effect size of interest.

### Design validation and missingness ([\#5](https://github.com/JUhalt/solomonR/issues/5), [\#6](https://github.com/JUhalt/solomonR/issues/6))

- New
  [`validate_solomon()`](https://juhalt.github.io/solomonR/reference/validate_solomon.md)
  checks input lengths, 0/1 coding, the presence of all four cells, and
  observed outcomes per cell, and returns every problem as an error,
  warning, or note instead of stopping at the first.
- New
  [`check_solomon_missing()`](https://juhalt.github.io/solomonR/reference/check_solomon_missing.md)
  separates structurally absent pretests in the unpretested groups from
  incidental pretest and posttest missingness, unexpected pretest
  scores, and unassigned participants. It reports counts by cell and
  cited guidance for each category (Solomon, 1949; Rubin, 1976; Graham
  et al., 2006; White & Thompson, 2005; Groenwold et al., 2012; Little &
  Rubin, 2019).

### Safeguards

- [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md),
  [`fit_solomon_sem()`](https://juhalt.github.io/solomonR/reference/fit_solomon_sem.md),
  [`fit_solomon_sem_latent()`](https://juhalt.github.io/solomonR/reference/fit_solomon_sem_latent.md),
  and
  [`check_solomon_assumptions()`](https://juhalt.github.io/solomonR/reference/check_solomon_assumptions.md)
  require `treat` and `pretested` to be coded 0/1 (or logical) and
  inputs to have equal lengths, instead of returning empty or rescaled
  results ([\#5](https://github.com/JUhalt/solomonR/issues/5)).
- [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md)
  warns when pretested participants are missing pretest scores (they are
  excluded, not imputed) and when pretest scores are supplied for
  unpretested participants
  ([\#5](https://github.com/JUhalt/solomonR/issues/5),
  [\#6](https://github.com/JUhalt/solomonR/issues/6)).
- [`power_solomon()`](https://juhalt.github.io/solomonR/reference/power_solomon.md)
  now warns that it is experimental. Interim fixes: `sigma` applies to
  every cell, the Stouffer arm uses the historical one-tailed direction,
  disabled metrics are `NA`, and `delta` is documented as the effect
  among unpretested participants. A validated rebuild is planned
  ([\#18](https://github.com/JUhalt/solomonR/issues/18)).
- `perm_solomon(seed = )` restores the global random-number state on
  exit.

### Output and documentation

- Printed GLM results name the covariance estimator and reference
  distribution. The summary no longer labels conventional standard
  errors as robust, and [`print()`](https://rdrr.io/r/base/print.html)
  returns the fitted object invisibly.
- [`plot_perm()`](https://juhalt.github.io/solomonR/reference/plot_perm.md)
  reports the same corrected p-value as
  [`perm_solomon()`](https://juhalt.github.io/solomonR/reference/perm_solomon.md).
- [`check_solomon_assumptions()`](https://juhalt.github.io/solomonR/reference/check_solomon_assumptions.md)
  presents its results as descriptive diagnostics rather than gates for
  choosing an analysis, and uses complete pretested cases for the
  slope-homogeneity test.
- Help pages now cite the methodological sources for each procedure
  ([\#9](https://github.com/JUhalt/solomonR/issues/9)).
- Test I is labeled as the Braver & Braver (1988) Stouffer combination
  in output and help pages, and its result gains a `procedure` column.
- `inst/CITATION` uses
  [`bibentry()`](https://rdrr.io/r/utils/bibentry.html).
- Removed the unused internal `spr2_ci()`, which relied on an invalid
  interval transformation, and duplicated print helpers.

### Development

- Set the v0.3.0 scope around graduate students and applied researchers
  and added v0.4.0, v0.5.0, v0.6.0, and v1.0.0 milestones that mirror
  the roadmap.

- solomonR 0.3.0 and later is licensed under GNU GPL version 3 only
  (GPL-3). Version 0.2.0 and earlier retain their original MIT terms and
  notices.

- Aligned development citation metadata, license pages, and
  roadmap/issue links.

- Distinguished completed implementation, deferred work, and research
  proposals in the roadmap without changing analytical behavior.

- Formatted `LICENSE.md` as markdown so the pkgdown license page has a
  title and readable headings; the license wording is unchanged.

- Added `Language: en-US` and a spelling word list (`inst/WORDLIST`) so
  documentation spell checks report only genuine errors.

- Began development toward solomonR 0.3.0.

- Added R-universe distribution and stable-release installation
  instructions.

- Reconciled the development roadmap with functionality delivered in
  v0.2.0.

## solomonR 0.2.0

### Historical Solomon workflow

- Added the full historical Tests A-I workflow in
  [`fit_solomon_classic()`](https://juhalt.github.io/solomonR/reference/fit_solomon_classic.md).
- Added Tests B and C for treatment simple effects.
- Added Test D as the equal-weighted treatment effect across pretest
  conditions.
- Added ANCOVA, gain-score, repeated-measures, and posttest-only
  historical analyses.
- Added historical decision-path reporting with `[PATH]` markers.
- Added directional Stouffer Test I for historical teaching and
  replication.
- Added explicit cautions regarding later Type I error critiques of
  conditional Solomon testing procedures.

### Modern observed-variable analysis

- Expanded
  [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md)
  as the primary unified Solomon analysis.
- Corrected the equal-weighted average treatment effect to
  `treat + 0.5 * treat:pretested`.
- Added Solomon-specific contrasts for:
  - average treatment effect,
  - pretest-by-treatment sensitization,
  - treatment effect among pretested participants,
  - treatment effect among unpretested participants.
- Preserved all four Solomon groups when pretest scores are structurally
  absent in the unpretested groups.
- Added HC3 heteroskedasticity-robust covariance estimation.
- Added CR2 cluster-robust covariance estimation.
- Added Wald-based partial R-squared reporting for Gaussian models.

### Randomization inference

- Rebuilt
  [`perm_solomon()`](https://juhalt.github.io/solomonR/reference/perm_solomon.md)
  as a stratified treatment-label permutation test.
- Treatment labels are permuted within pretest-assignment strata.
- Added HC3-studentized permutation statistics.
- Added finite Monte Carlo p-value correction.
- Added concise printing for permutation objects.
- Added permutation-distribution plotting with
  [`plot_perm()`](https://juhalt.github.io/solomonR/reference/plot_perm.md).

### Maximum likelihood

- Added
  [`fit_solomon_ml()`](https://juhalt.github.io/solomonR/reference/fit_solomon_ml.md)
  for full-information maximum-likelihood analysis.
- Structural absence of pretest scores in the unpretested groups is
  handled as part of the Solomon design rather than as ordinary missing
  baseline data.
- Added Solomon-specific ATE, sensitization, and simple treatment
  effects.
- Added automated validation against corresponding component
  regressions.

### Structural equation models

- Added observed-variable Solomon SEM with
  [`fit_solomon_sem()`](https://juhalt.github.io/solomonR/reference/fit_solomon_sem.md).
- Added an ANCOVA-style SEM pathway for the pretested groups.
- Added latent-variable Solomon SEM with
  [`fit_solomon_sem_latent()`](https://juhalt.github.io/solomonR/reference/fit_solomon_sem_latent.md).
- Added measurement-invariance controls for latent outcomes.
- Latent Solomon mean contrasts now require scalar measurement
  invariance.
- Added explicit model-convergence checks.
- Clarified that global fit indices are not diagnostic for the saturated
  four-group observed mean model.

### Diagnostics and effect sizes

- Added Brown-Forsythe variance checks.
- Added cell-level normality summaries.
- Added ANCOVA slope-homogeneity diagnostics.
- Added Hedges’ g for the posttest-only comparison.
- Corrected and clarified Wald-based partial R-squared reporting.

### Documentation and reliability

- Expanded automated tests for historical, GLM, permutation, ML, SEM,
  printing, and plotting workflows.
- Added known-estimand tests for Solomon contrasts.
- Rebuilt the classic and modern analysis vignettes.
- Expanded the package README with installation, method-selection, and
  interpretation guidance.
- Added APA-oriented rounding and concise public-facing output.
- Package checks pass with zero errors, warnings, or notes.
