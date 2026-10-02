# solomonR: Analyze Solomon Four-Group Designs

Provides tools for analyzing, teaching, and studying the Solomon
four-group experimental design. The package includes the historical
Tests A-I workflow, unified generalized linear models with robust
covariance estimation, stratified randomization inference,
Solomon-specific full-information maximum likelihood, and observed- and
latent-variable structural equation models.

## Details

The package distinguishes historically important procedures from
contemporary recommendations and emphasizes explicit Solomon-specific
estimands, structural pretest missingness, reproducible diagnostics, and
modern inference.

Key historical references include Huck & Sandler (1973), Walton Braver &
Braver (1988), the Sawilowsky and Markman methodological exchanges, and
van Engelenburg (1999).

## Lifecycle

Each exported function's help page, and the reference index of the
package website, shows the function's lifecycle stage, in the stages of
the lifecycle package (Henry & Wickham, 2026):

- **Stable.** The interface is settled. Any change goes through
  deprecation: the old name keeps working, with a warning, through v1.x.

- **Experimental.** The function is tested, but its interface or
  defaults may change.
  [`fit_solomon_sem()`](https://juhalt.github.io/solomonR/reference/fit_solomon_sem.md)
  has no simulation study yet at Solomon sample sizes. The issue \#55
  study found no measurement-invariance criterion that holds its
  false-rejection rate in Solomon-sized groups, which affects
  [`fit_solomon_sem_latent()`](https://juhalt.github.io/solomonR/reference/fit_solomon_sem_latent.md)
  and
  [`invariance_solomon()`](https://juhalt.github.io/solomonR/reference/invariance_solomon.md).

- **Deprecated.** A former name that still works, with a warning. Its
  help page names the replacement.

solomonR 0.9.0 settled the interface (issue \#83). The outcome arguments
are `y_post` and `y_pre` everywhere, a fitted model is passed as `fit`,
and the planning functions share one argument order: `n` (or `power`),
`delta`, `sens`, `rho`, `sigma`, `alpha`. Each function that takes data
vectors also takes an optional `data` data frame.

## References

Campbell, D. T., & Stanley, J. C. (1966). *Experimental and
quasi-experimental designs for research*. Rand McNally. (Original work
published 1963)

Henry, L., & Wickham, H. (2026). *lifecycle: Manage the life cycle of
your package functions* (Version 1.0.5) \[R package\].
https://doi.org/10.32614/CRAN.package.lifecycle

Huck, S. W., & Sandler, H. M. (1973). A note on the Solomon 4-group
design: Appropriate statistical analyses. *The Journal of Experimental
Education, 42*(2), 54–55. https://doi.org/10.1080/00220973.1973.11011460

Sawilowsky, S. S., Kelley, D. L., Blair, R. C., & Markman, B. S. (1994).
Meta-analysis and the Solomon four-group design. *The Journal of
Experimental Education, 62*(4), 361–376.
https://doi.org/10.1080/00220973.1994.9944140

Solomon, R. L. (1949). An extension of control group design.
*Psychological Bulletin, 46*(2), 137–150.
https://doi.org/10.1037/h0062958

van Engelenburg, G. (1999). *Statistical analysis for the Solomon
four-group design* (Research Report 99-06). University of Twente.

Walton Braver, M. C., & Braver, S. L. (1988). Statistical treatment of
the Solomon four-group design: A meta-analytic approach. *Psychological
Bulletin, 104*(1), 150–154. https://doi.org/10.1037/0033-2909.104.1.150
