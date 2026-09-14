# solomonR 0.2.0.9000

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

* Changed current development licensing to GNU GPL version 3 only (GPL-3).
  Previously published releases retain their original MIT terms and notices.
* Aligned development citation metadata, license pages, and roadmap/issue links.
* Distinguished completed implementation, deferred work, and research proposals
  in the roadmap without changing analytical behavior.

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
