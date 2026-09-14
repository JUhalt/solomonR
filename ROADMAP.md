# solomonR Roadmap

`solomonR` aims to provide a comprehensive R ecosystem for the
Solomon four-group design: preserving its methodological history,
supporting teaching and replication, and providing defensible modern
analysis, visualization, design-planning, and reporting tools.

The roadmap distinguishes historically important methods from
contemporary recommendations and from new extensions introduced by
`solomonR`.

The primary audience is graduate students (Master's and doctoral) and
applied researchers who need to plan, analyze, interpret, and report a
Solomon four-group study. Every procedure should be traceable to published
methodology or clearly labeled as a `solomonR`-specific extension.

**Distribution:** releases are published on
[R-universe](https://juhalt.r-universe.dev/solomonR). CRAN submission is
planned for v1.0.0.

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

## v0.2.0 - Core Solomon toolkit (completed)

**Status:** Released August 30, 2026

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
- [x] Reconstruct full historical Solomon Test A-I workflow
- [x] Test A: Pretest x Treatment
- [x] Tests B/C: Simple-effect follow-ups
- [x] Test D: Treatment main effect
- [x] Test E: ANCOVA
- [x] Test F: Gain-score analysis
- [x] Test G: Repeated-measures alternative
- [x] Test H: Posttest-only comparison
- [x] Test I: Stouffer meta-analytic combination
- [x] Clearly label historically proposed vs currently recommended procedures
- [x] Correct directional Stouffer implementation
- [x] Document later Type I error critiques

### Diagnostics and effect sizes
- [x] Brown-Forsythe checks
- [x] ANCOVA slope-homogeneity check
- [x] Cell normality summaries
- [x] Hedges g for posttest-only comparison
- [x] Correct / rename Wald-derived R-squared effect size
Broader confidence-interval consistency was deferred to v0.3; see [#8](https://github.com/JUhalt/solomonR/issues/8).

### SEM
- [x] Observed mean-structure SEM
- [x] Latent POST SEM
- [x] Optional latent PRE -> POST model in pretested groups
- [x] Measurement-invariance options
- [x] Add SEM unit tests
Expanded method-positioning documentation was deferred to v0.3; see [#9](https://github.com/JUhalt/solomonR/issues/9).

### Package quality
- [x] testthat infrastructure
- [x] pkgdown infrastructure
- [x] bundled demo data
- [x] print / summary methods
- [x] `R CMD check`: 0 errors, 0 warnings, 0 notes
- [x] Public GitHub repository
- [x] GitHub release
- [x] Public pkgdown site
- [x] Installation instructions
- [x] R-universe distribution

---

## v0.3.0 - Modern Solomon methods

**Status:** Active development — [v0.3.0 milestone](https://github.com/JUhalt/solomonR/milestone/1).

**Goal:** Make solomonR a modern analysis toolkit rather than only
an implementation of historical workflows.

Already delivered in v0.2.0: `fit_solomon_ml()`. Broader validation remains open in [#10](https://github.com/JUhalt/solomonR/issues/10).
  - Solomon-specific maximum-likelihood regression
  - Based on van Engelenburg (1999)
  - Structural pretest missingness handled explicitly
  - Existing known-result validation; broader simulation evidence is tracked separately.

- [ ] Sensitization equivalence testing — [#4](https://github.com/JUhalt/solomonR/issues/4)
  - CI-based / TOST-style inference
  - User-specified smallest effect size of interest
  - Distinguish "no significant sensitization" from evidence of negligible sensitization

- [x] `validate_solomon()` — [#5](https://github.com/JUhalt/solomonR/issues/5)
  - Verify all four cells
  - Inspect group sizes
  - Validate coding
  - Detect expected structural pretest missingness
  - Detect unexpected missingness
  - Identify sparse / empty cells

- [x] `check_solomon_missing()` — [#6](https://github.com/JUhalt/solomonR/issues/6)
  - Distinguish structural from incidental missingness
  - Prevent inappropriate imputation of deliberately absent pretests

- [ ] `compare_solomon_methods()` — [#7](https://github.com/JUhalt/solomonR/issues/7)
  - Side-by-side classic, GLM, ML, permutation, and SEM results
  - Explicitly identify differing estimands

- [x] Add confidence intervals consistently across release-defining effect summaries, with literature-based reference distributions and an HC3 default — [#8](https://github.com/JUhalt/solomonR/issues/8).

- [ ] Expand documentation distinguishing established Solomon methodology,
  contemporary recommendations, and solomonR-specific extensions — [#9](https://github.com/JUhalt/solomonR/issues/9).

- [ ] Broaden ML simulation validation — [#10](https://github.com/JUhalt/solomonR/issues/10).

- [ ] Reconcile release documentation, licensing and distribution — [#12](https://github.com/JUhalt/solomonR/issues/12).

### Correctness (September 2026 review)

- [x] Satterthwaite small-sample tests for CR2 contrasts — [#14](https://github.com/JUhalt/solomonR/issues/14).
- [x] Refuse individual-level permutation of clustered fits — [#15](https://github.com/JUhalt/solomonR/issues/15).
- [x] Identify the latent Solomon mean structure — [#16](https://github.com/JUhalt/solomonR/issues/16).
- [x] Interim `power_solomon()` safeguards (experimental warning, corrected simulator details); the validated rebuild remains in v0.4.0 — [#18](https://github.com/JUhalt/solomonR/issues/18).

### Documentation for the intended audience

- [ ] Getting-started guide for graduate students and applied researchers, moved forward from v0.6.0 — [#17](https://github.com/JUhalt/solomonR/issues/17).

---

## v0.4.0 - Design planning and visualization

**Milestone:** [v0.4.0](https://github.com/JUhalt/solomonR/milestone/2)

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
- [ ] Rebuild and validate `power_solomon()` — [#18](https://github.com/JUhalt/solomonR/issues/18)
- [ ] Correct data-generating mechanism
- [ ] Validate Type I error and power
- [ ] `plan_solomon()`
  - Required sample size
  - Four-cell allocation
  - Power for multiple Solomon estimands
- [ ] Explore unequal allocation strategies

---

## v0.5.0 - Extended outcomes and designs

**Milestone:** [v0.5.0](https://github.com/JUhalt/solomonR/milestone/3)

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
  - Cluster-level randomization inference — [#19](https://github.com/JUhalt/solomonR/issues/19)

- [ ] Longitudinal Solomon models
  - Repeated follow-ups
  - Treatment x Pretest x Time

- [ ] Modified Solomon designs
  - More than two treatment conditions

- [ ] Quasi-experimental Solomon designs
  - Explicitly distinguish causal interpretation from randomized designs

---

## v0.6.0 - Teaching and reporting

**Milestone:** [v0.6.0](https://github.com/JUhalt/solomonR/milestone/4)

- [ ] Complete introductory vignette — the getting-started guide moved to v0.3.0 ([#17](https://github.com/JUhalt/solomonR/issues/17)); extend it here as needed.
- [ ] Expand the existing historical-analysis vignette; define additional coverage when scoped.
- [ ] Expand the existing modern-analysis vignette; define additional coverage when scoped.
- [ ] SEM / latent-variable vignette
- [ ] Power / planning vignette
- [ ] Single worked social-psychology example across methods
- [ ] APA-style reporting helper
- [ ] Method-selection guide — initial guidance ships with v0.3.0 ([#9](https://github.com/JUhalt/solomonR/issues/9), [#17](https://github.com/JUhalt/solomonR/issues/17))
- [ ] Historical timeline / decision-tree documentation
- [ ] Expanded pkgdown site

---

# v1.0.0 - Stable comprehensive release

**Milestone:** [v1.0.0](https://github.com/JUhalt/solomonR/milestone/5). Distribution remains R-universe until this release.

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

These are intentionally outside the initial release path. Proposals become
release commitments only when their scope and acceptance criteria are agreed
in a linked issue and assigned to a milestone.

- [Proposal: reproducible simulation benchmark reports with Monte Carlo uncertainty (#11)](https://github.com/JUhalt/solomonR/issues/11).
  This research-informed reporting extension would build on existing ML-validation
  and method-comparison work; no additional release scope is committed.

- Bayesian Solomon modeling
- Rank-based / nonparametric unified methods
- Multivariate outcomes
- Latent change-score Solomon models
- Missing-data sensitivity analyses
- Optimal allocation algorithms
- Shiny teaching / design application
