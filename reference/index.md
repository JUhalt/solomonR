# Package index

## Modern Solomon analysis

Unified regression, likelihood, and randomization-based methods.

- [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md)
  : Fit the unified GLM for a Solomon Four-Group design
- [`fit_solomon_ml()`](https://juhalt.github.io/solomonR/reference/fit_solomon_ml.md)
  : Full-information ML analysis for a Solomon Four-Group Design
- [`perm_solomon()`](https://juhalt.github.io/solomonR/reference/perm_solomon.md)
  : Permutation test for a Solomon contrast
- [`equivalence_solomon()`](https://juhalt.github.io/solomonR/reference/equivalence_solomon.md)
  : Equivalence test for a Solomon contrast
- [`compare_solomon_methods()`](https://juhalt.github.io/solomonR/reference/compare_solomon_methods.md)
  : Compare Solomon analyses and their estimands

## Historical Solomon analysis

Classical procedures retained for teaching and methodological
replication.

- [`fit_solomon_classic()`](https://juhalt.github.io/solomonR/reference/fit_solomon_classic.md)
  : Historical Solomon Four-Group Analysis
- [`stouffer_solomon()`](https://juhalt.github.io/solomonR/reference/stouffer_solomon.md)
  : Stouffer's Z combiner (Braver & Braver, 1988, Test I)
- [`p_to_z()`](https://juhalt.github.io/solomonR/reference/p_to_z.md) :
  Convert p-value to Z (one-tailed) for Stouffer's method

## Structural equation models

Observed- and latent-variable SEM approaches to the Solomon design.

- [`fit_solomon_sem()`](https://juhalt.github.io/solomonR/reference/fit_solomon_sem.md)
  : SEM analysis for Solomon Four-Group designs (mean-structure;
  optional ANCOVA)
- [`fit_solomon_sem_latent()`](https://juhalt.github.io/solomonR/reference/fit_solomon_sem_latent.md)
  : Latent SEM for Solomon Four-Group designs (POST means; optional
  latent ANCOVA)

## Design validation and diagnostics

Check design coding, cell structure, and missingness before analysis.

- [`validate_solomon()`](https://juhalt.github.io/solomonR/reference/validate_solomon.md)
  : Validate the structure and coding of a Solomon four-group design
- [`check_solomon_missing()`](https://juhalt.github.io/solomonR/reference/check_solomon_missing.md)
  : Distinguish structural and incidental missingness in a Solomon
  design
- [`check_solomon_assumptions()`](https://juhalt.github.io/solomonR/reference/check_solomon_assumptions.md)
  : Descriptive assumption diagnostics for Solomon analyses
- [`power_solomon()`](https://juhalt.github.io/solomonR/reference/power_solomon.md)
  : Power simulation for Solomon designs (experimental)

## Plotting

- [`plot_solomon()`](https://juhalt.github.io/solomonR/reference/plot_solomon.md)
  : Plot Solomon Posttest Cell Means
- [`plot_solomon_gg()`](https://juhalt.github.io/solomonR/reference/plot_solomon_gg.md)
  : Plot Solomon Posttest Cell Means with ggplot2
- [`plot_perm()`](https://juhalt.github.io/solomonR/reference/plot_perm.md)
  : Plot permutation distribution for a Solomon contrast

## Data

- [`solomon_example`](https://juhalt.github.io/solomonR/reference/solomon_example.md)
  : Simulated Solomon four-group study (primary example)
- [`solomon_demo`](https://juhalt.github.io/solomonR/reference/solomon_demo.md)
  : Demo Solomon four-group data set (second example)

## Package

- [`solomonR`](https://juhalt.github.io/solomonR/reference/solomonR.md)
  : solomonR: Analyze Solomon Four-Group Designs
