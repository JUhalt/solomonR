# Package index

## Decide and plan

Sample sizes and power for a Solomon study, before data are collected,
and simulated studies with known effects for teaching and for checking
an analysis plan.

- [`plan_solomon()`](https://juhalt.github.io/solomonR/reference/plan_solomon.md)
  **\[stable\]** : Sample-size planning for Solomon designs
- [`power_solomon()`](https://juhalt.github.io/solomonR/reference/power_solomon.md)
  **\[stable\]** : Power simulation for Solomon designs
- [`plot_power_solomon()`](https://juhalt.github.io/solomonR/reference/plot_power_solomon.md)
  **\[stable\]** : Power curves for a Solomon design
- [`analysis_plan_solomon()`](https://juhalt.github.io/solomonR/reference/analysis_plan_solomon.md)
  **\[experimental\]** : Write an analysis plan for a Solomon four-group
  study
- [`simulate_solomon()`](https://juhalt.github.io/solomonR/reference/simulate_solomon.md)
  **\[stable\]** : Simulate a Solomon four-group study with known
  effects

## Check the design

Design coding, cell structure, clustering, missingness, and baseline
balance, before analysis.

- [`validate_solomon()`](https://juhalt.github.io/solomonR/reference/validate_solomon.md)
  **\[stable\]** : Validate the structure and coding of a Solomon design
- [`check_solomon_missing()`](https://juhalt.github.io/solomonR/reference/check_solomon_missing.md)
  **\[stable\]** : Distinguish structural and incidental missingness in
  a Solomon design
- [`check_solomon_assumptions()`](https://juhalt.github.io/solomonR/reference/check_solomon_assumptions.md)
  **\[stable\]** : Descriptive assumption diagnostics for Solomon
  analyses
- [`baseline_solomon()`](https://juhalt.github.io/solomonR/reference/baseline_solomon.md)
  **\[stable\]** : Baseline comparison of the pretested arms

## Analyze: recommended methods

One model for all the groups, with robust, likelihood-based, and
randomization-based inference. fit_solomon_glm() also analyzes designs
with several treatments.

- [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md)
  **\[stable\]** : Fit the unified GLM for a Solomon Four-Group design
- [`fit_solomon_ml()`](https://juhalt.github.io/solomonR/reference/fit_solomon_ml.md)
  **\[stable\]** : Full-information ML analysis for a Solomon Four-Group
  Design
- [`marginal_solomon()`](https://juhalt.github.io/solomonR/reference/marginal_solomon.md)
  **\[stable\]** : Marginal Solomon contrasts for binary and count
  outcomes
- [`perm_solomon()`](https://juhalt.github.io/solomonR/reference/perm_solomon.md)
  **\[stable\]** : Permutation test for a Solomon contrast
- [`equivalence_solomon()`](https://juhalt.github.io/solomonR/reference/equivalence_solomon.md)
  **\[stable\]** : Equivalence test for a Solomon contrast
- [`compare_solomon_methods()`](https://juhalt.github.io/solomonR/reference/compare_solomon_methods.md)
  **\[stable\]** : Compare Solomon analyses and their estimands

## Analyze: missing posttests

Multiple imputation of missing posttests under missing at random, with
delta-adjusted sensitivity analyses and a tipping-point summary.

- [`fit_solomon_mi()`](https://juhalt.github.io/solomonR/reference/fit_solomon_mi.md)
  **\[experimental\]** : Solomon analysis with multiply imputed
  posttests
- [`tipping_point_solomon()`](https://juhalt.github.io/solomonR/reference/tipping_point_solomon.md)
  **\[experimental\]** : Tipping-point analysis for missing posttests

## Analyze: longitudinal designs

Several posttest occasions, analyzed with a mixed model for repeated
measures that is valid when participants drop out at random.

- [`fit_solomon_mmrm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_mmrm.md)
  **\[experimental\]** : Mixed model for repeated measures in a
  longitudinal Solomon design

## Analyze: latent variables

Observed- and latent-variable structural equation models, and
measurement invariance across the Solomon groups.

- [`fit_solomon_sem()`](https://juhalt.github.io/solomonR/reference/fit_solomon_sem.md)
  **\[experimental\]** : SEM analysis for Solomon Four-Group designs
  (mean-structure; optional ANCOVA)
- [`fit_solomon_sem_latent()`](https://juhalt.github.io/solomonR/reference/fit_solomon_sem_latent.md)
  **\[experimental\]** : Latent SEM for Solomon Four-Group designs (POST
  means; optional latent ANCOVA)
- [`invariance_solomon()`](https://juhalt.github.io/solomonR/reference/invariance_solomon.md)
  **\[experimental\]** : Measurement invariance across the four Solomon
  groups

## Report

APA 7 results text with the references for the methods used.

- [`report_solomon()`](https://juhalt.github.io/solomonR/reference/report_solomon.md)
  **\[stable\]** : APA 7 results text for a Solomon analysis

## Synthesize

Reanalyze published studies from summary statistics and compute effect
sizes for meta-analysis.

- [`solomon_from_summary()`](https://juhalt.github.io/solomonR/reference/solomon_from_summary.md)
  **\[stable\]** : Solomon analysis from summary statistics
- [`solomon_effect_sizes()`](https://juhalt.github.io/solomonR/reference/solomon_effect_sizes.md)
  **\[stable\]** : Solomon effect sizes for meta-analysis

## History: historical procedures

Procedures retained for teaching and for reproducing historical
analyses, labeled by author and year. They are not the recommended
analyses.

- [`fit_solomon_1949()`](https://juhalt.github.io/solomonR/reference/fit_solomon_1949.md)
  **\[stable\]** : Solomon's (1949) improvement-score analysis
- [`fit_solomon_classic()`](https://juhalt.github.io/solomonR/reference/fit_solomon_classic.md)
  **\[stable\]** : Historical Solomon Four-Group Analysis
- [`plot_classic_flow()`](https://juhalt.github.io/solomonR/reference/plot_classic_flow.md)
  **\[stable\]** : Historical Solomon decision path
- [`fit_solomon_steyn()`](https://juhalt.github.io/solomonR/reference/fit_solomon_steyn.md)
  **\[stable\]** : Steyn's (2009) analysis of the extended Solomon
  design
- [`fisher_solomon()`](https://juhalt.github.io/solomonR/reference/fisher_solomon.md)
  **\[stable\]** : Historical categorical analysis of a binary Solomon
  outcome
- [`stouffer_solomon()`](https://juhalt.github.io/solomonR/reference/stouffer_solomon.md)
  **\[stable\]** : Stouffer's Z combiner (Walton Braver & Braver, 1988,
  Test I)
- [`p_to_z()`](https://juhalt.github.io/solomonR/reference/p_to_z.md)
  **\[stable\]** : Convert p-value to Z (one-tailed) for Stouffer's
  method

## Plotting

- [`plot_solomon_means()`](https://juhalt.github.io/solomonR/reference/plot_solomon_means.md)
  **\[stable\]** : Plot Solomon posttest cell means
- [`plot_tipping_point()`](https://juhalt.github.io/solomonR/reference/plot_tipping_point.md)
  **\[experimental\]** : Plot a tipping-point analysis
- [`plot_perm()`](https://juhalt.github.io/solomonR/reference/plot_perm.md)
  **\[stable\]** : Plot permutation distribution for a Solomon contrast
- [`plot_solomon_effects()`](https://juhalt.github.io/solomonR/reference/plot_solomon_effects.md)
  **\[stable\]** : Forest plot of the Solomon contrasts
- [`plot_sensitization()`](https://juhalt.github.io/solomonR/reference/plot_sensitization.md)
  **\[stable\]** : Pretest sensitization figure
- [`plot_solomon_design()`](https://juhalt.github.io/solomonR/reference/plot_solomon_design.md)
  **\[stable\]** : Schematic of a Solomon design
- [`plot_solomon_change()`](https://juhalt.github.io/solomonR/reference/plot_solomon_change.md)
  **\[stable\]** : Pretest-to-posttest change in a Solomon design

## Data

- [`solomon_example`](https://juhalt.github.io/solomonR/reference/solomon_example.md)
  : Simulated Solomon four-group study (primary example)
- [`solomon_demo`](https://juhalt.github.io/solomonR/reference/solomon_demo.md)
  : Demo Solomon four-group data set (second example)
- [`elkarkri2025a`](https://juhalt.github.io/solomonR/reference/elkarkri2025a.md)
  : Pretest and posttest statistics from El Karkri et al. (2025a)
- [`jordaan2014`](https://juhalt.github.io/solomonR/reference/jordaan2014.md)
  : Group statistics from Jordaan's (2014) study with three posttest
  occasions
- [`kvalem1996`](https://juhalt.github.io/solomonR/reference/kvalem1996.md)
  : Condom use in the Solomon study of Kvalem et al. (1996)
- [`lana1959`](https://juhalt.github.io/solomonR/reference/lana1959.md)
  : Posttest statistics from Lana's (1959) attitude experiment
- [`mai2020`](https://juhalt.github.io/solomonR/reference/mai2020.md) :
  Transfer-intervention data from Mai et al. (2020)
- [`solomon1949`](https://juhalt.github.io/solomonR/reference/solomon1949.md)
  : Solomon's (1949) spelling experiment
- [`steyn2005`](https://juhalt.github.io/solomonR/reference/steyn2005.md)
  : Group statistics from Steyn's (2005) eight-group study
- [`waltonbraver1988`](https://juhalt.github.io/solomonR/reference/waltonbraver1988.md)
  : The worked example of Walton Braver and Braver (1988)

## Deprecated

Former names that still work, with a warning, through v1.x. Each help
page names the replacement.

- [`plot_solomon()`](https://juhalt.github.io/solomonR/reference/plot_solomon.md)
  **\[deprecated\]** : Plot Solomon posttest cell means in base graphics
  (deprecated)
- [`plot_solomon_gg()`](https://juhalt.github.io/solomonR/reference/plot_solomon_gg.md)
  **\[deprecated\]** : Plot Solomon posttest cell means with ggplot2
  (deprecated)

## Package

An overview of the package, including the lifecycle stages marked by the
badge beside each function: stable, experimental, or deprecated.

- [`solomonR`](https://juhalt.github.io/solomonR/reference/solomonR.md)
  : solomonR: Analyze Solomon Four-Group Designs
