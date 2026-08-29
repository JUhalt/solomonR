# solomonR Roadmap

`solomonR` aims to provide a comprehensive R ecosystem for the
Solomon four-group design: preserving its methodological history,
supporting teaching and replication, and providing defensible modern
analysis, visualization, design-planning, and reporting tools.

The roadmap distinguishes historically important methods from
contemporary recommendations and from new extensions introduced by
`solomonR`.

---

## Guiding principles

1. **Preserve the history**
   - Reproduce the major historical approaches faithfully.
   - Identify which procedures are historical rather than currently recommended.
   - Preserve methodological debates, including the Braver & Braver /
     Sawilowsky et al. meta-analytic controversy.

2. **Teach the design**
   - Make the logic of the Solomon four-group design visible.
   - Provide readable summaries, decision pathways, diagrams, and worked examples.
   - Explain what each analysis estimates and why.

3. **Support modern inference**
   - Prefer clearly defined estimands and unified models.
   - Support robust and randomization-based inference.
   - Add Solomon-specific maximum-likelihood modeling.

4. **Extend the design**
   - Support latent-variable SEM, clustered/longitudinal designs,
     generalized outcomes, and other modern extensions where justified.

5. **Validate before expanding**
   - New methods require mathematical or simulation-based tests.
   - Package checks should remain clean.
   - Historical fidelity and statistical correctness take priority over feature count.

---

# Release plan

## v0.2.0 - Core Solomon toolkit

**Goal:** First public, usable release.

### Modern observed-variable analysis
- [x] Unified GLM
- [x] HC3 robust inference
- [x] CR2 cluster-robust covariance option
- [x] Treatment effect among pretested participants
- [x] Treatment effect among unpretested participants
- [x] Pretest x Treatment sensitization contrast
- [x] Equal-weighted average treatment effect
- [x] Known-estimand unit tests

### Randomization / permutation inference
- [x] Stratified treatment-label permutation
- [x] Studentized permutation statistic
- [x] Finite-Monte-Carlo p-value correction
- [x] Permutation distribution plot
- [x] Unit tests

### Historical / teaching analysis
- [ ] Reconstruct full historical Solomon Test A-I workflow
- [ ] Test A: Pretest x Treatment
- [ ] Tests B/C: Simple-effect follow-ups
- [ ] Test D: Treatment main effect
- [ ] Test E: ANCOVA
- [ ] Test F: Gain-score analysis
- [ ] Test G: Repeated-measures alternative
- [ ] Test H: Posttest-only comparison
- [ ] Test I: Stouffer meta-analytic combination
- [ ] Clearly label historically proposed vs currently recommended procedures
- [ ] Correct directional Stouffer implementation
- [ ] Document later Type I error critiques

### Diagnostics and effect sizes
- [x] Brown-Forsythe checks
- [x] ANCOVA slope-homogeneity check
- [x] Cell normality summaries
- [x] Hedges g for posttest-only comparison
- [ ] Correct / rename Wald-derived R-squared effect size
- [ ] Add confidence intervals consistently

### SEM
- [x] Observed mean-structure SEM
- [x] Latent POST SEM
- [x] Optional latent PRE -> POST model in pretested groups
- [x] Measurement-invariance options
- [ ] Add SEM unit tests
- [ ] Clarify established methodology vs solomonR extensions

### Package quality
- [x] testthat infrastructure
- [x] pkgdown infrastructure
- [x] bundled demo data
- [x] print / summary methods
- [ ] `R CMD check`: 0 errors, 0 warnings, 0 notes
- [ ] Public GitHub repository
- [ ] GitHub release
- [ ] Public pkgdown site
- [ ] Installation instructions

---

## v0.3.0 - Modern Solomon methods

**Goal:** Make solomonR a modern analysis toolkit rather than only
an implementation of historical workflows.

- [ ] `fit_solomon_ml()`
  - Solomon-specific maximum-likelihood regression
  - Based on van Engelenburg (1999)
  - Structural pretest missingness handled explicitly
  - Known-result / simulation validation

- [ ] Sensitization equivalence testing
  - CI-based / TOST-style inference
  - User-specified smallest effect size of interest
  - Distinguish "no significant sensitization" from evidence of negligible sensitization

- [ ] `validate_solomon()`
  - Verify all four cells
  - Inspect group sizes
  - Validate coding
  - Detect expected structural pretest missingness
  - Detect unexpected missingness
  - Identify sparse / empty cells

- [ ] `check_solomon_missing()`
  - Distinguish structural from incidental missingness
  - Prevent inappropriate imputation of deliberately absent pretests

- [ ] `compare_solomon_methods()`
  - Side-by-side classic, GLM, ML, permutation, and SEM results
  - Explicitly identify differing estimands

---

## v0.4.0 - Design planning and visualization

**Goal:** Make the design easy to understand visually and useful before
data collection begins.

### Visualization
- [ ] `plot_solomon_design()`
  - Solomon four-group design schematic
  - Optional cell N / means

- [ ] `plot_sensitization()`
  - Treatment x Pretest interaction visualization
  - Direct graphical representation of sensitization

- [ ] `plot_solomon_effects()`
  - Forest plot of ATE, simple treatment effects, and sensitization

- [ ] `plot_solomon_change()`
  - Pre/post trajectories for pretested groups

- [ ] `plot_classic_flow()`
  - Historical decision tree
  - Optionally highlight the path taken by a fitted dataset

- [ ] `plot_power_solomon()`
  - Power curves / surfaces across N, effect size, rho, and sensitization

### Design planning
- [ ] Rebuild and validate `power_solomon()`
- [ ] Correct data-generating mechanism
- [ ] Validate Type I error and power
- [ ] `plan_solomon()`
  - Required sample size
  - Four-cell allocation
  - Power for multiple Solomon estimands
- [ ] Explore unequal allocation strategies

---

## v0.5.0 - Extended outcomes and designs

- [ ] Binary outcomes
  - Logistic models
  - Risk difference
  - Risk ratio
  - Odds ratio

- [ ] Count outcomes
  - Poisson / negative-binomial models
  - Rate ratios

- [ ] Mixed / multilevel Solomon models
  - Clustered assignment
  - Classrooms / schools / sites

- [ ] Longitudinal Solomon models
  - Repeated follow-ups
  - Treatment x Pretest x Time

- [ ] Modified Solomon designs
  - More than two treatment conditions

- [ ] Quasi-experimental Solomon designs
  - Explicitly distinguish causal interpretation from randomized designs

---

## v0.6.0 - Teaching and reporting

- [ ] Complete introductory vignette
- [ ] Historical-analysis vignette
- [ ] Modern-analysis vignette
- [ ] SEM / latent-variable vignette
- [ ] Power / planning vignette
- [ ] Single worked social-psychology example across methods
- [ ] APA-style reporting helper
- [ ] Method-selection guide
- [ ] Historical timeline / decision-tree documentation
- [ ] Expanded pkgdown site

---

# v1.0.0 - Stable comprehensive release

Candidate requirements:

- Stable public API
- Historical workflow validated
- GLM and ML estimands validated
- Randomization inference validated
- Power simulation validated
- Core visualizations complete
- SEM pathway documented and tested
- Strong automated test suite
- Clean R CMD check
- Cross-platform CI checks
- Public documentation website
- CRAN-ready package
- Reproducible manuscript examples and simulations
- Companion methodological/software manuscript

---

## Publication milestone

Development of the companion methodological/software manuscript will
continue alongside package development.

### v0.8.x - Methodological feature freeze
- Freeze major analyses intended for the v1.0 paper
- Finalize simulation conditions
- Finalize worked example
- Validate all primary estimands and inferential procedures
- Begin full manuscript drafting

### v0.9.x - Release candidate and manuscript freeze
- Stable candidate API
- Complete reproducible simulations
- Complete manuscript tables and figures
- Public preprint when appropriate
- CRAN pre-submission checks

### v1.0.0 - Stable package + manuscript submission
- Stable public release
- Submit to CRAN
- Archive release and reproducible materials
- Submit companion peer-reviewed methodological/software paper

Publication acceptance is not required before the v1.0 release.

# Later / exploratory

These are intentionally outside the initial release path.

- Bayesian Solomon modeling
- Rank-based / nonparametric unified methods
- Multivariate outcomes
- Latent change-score Solomon models
- Missing-data sensitivity analyses
- Optimal allocation algorithms
- Shiny teaching / design application
