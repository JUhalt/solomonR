# solomonR 0.2.0.9000

## Development

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
