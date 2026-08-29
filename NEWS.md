# solomonR 0.2.0

## New analysis methods

* Added observed-variable SEM support.
* Added latent-variable Solomon SEM with measurement-invariance options.
* Added robust unified GLM contrasts for:
  * equal-weighted average treatment effect,
  * pretest-by-treatment sensitization,
  * treatment effect among pretested participants,
  * treatment effect among unpretested participants.

## Randomization inference

* Rebuilt `perm_solomon()` as a stratified permutation test.
* Treatment labels are permuted within pretest-assignment strata.
* Added studentized permutation statistics.
* Added finite-Monte-Carlo p-value correction.
* Added permutation-distribution plotting.

## Reliability and testing

* Added known-estimand tests for Solomon GLM contrasts.
* Added permutation-inference tests.
* Improved package documentation and dependency declarations.
* Cleaned package-check issues.

## In development

* Full historical Test A-I workflow.
* Solomon-specific maximum-likelihood regression.
* Expanded visualization and design-planning tools.
