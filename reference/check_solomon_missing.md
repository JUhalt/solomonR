# Distinguish structural and incidental missingness in a Solomon design

**\[stable\]** Classifies missing values in Solomon design data and
explains the supported response to each kind. Pretest scores are
*structurally absent* for participants assigned to the unpretested
groups: withholding the pretest is the experimental manipulation
(Solomon, 1949), so those values must never be imputed. Other missing
values are *incidental* and are handled according to the missing-data
literature (Rubin, 1976; Little & Rubin, 2019).

## Usage

``` r
check_solomon_missing(
  y_post,
  treat,
  pretested,
  y_pre = NULL,
  control = NULL,
  data = NULL
)
```

## Arguments

- y_post:

  Numeric posttest scores.

- treat:

  Treatment indicator coded 0/1 (or logical); or a factor or character
  vector of conditions, with the control named by `control`.

- pretested:

  Pretest indicator coded 0/1 (or logical).

- y_pre:

  Optional numeric pretest scores, missing by design for unpretested
  participants.

- control:

  The control condition when `treat` is a factor or character vector
  with more than two conditions, a Solomon N-group design; see
  [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md).
  With two conditions the result is the same as with a 0/1 `treat`.

- data:

  Optional data frame. When supplied, the other data arguments are
  looked up in it first, as bare column names (`y_post = post`) or as
  strings (`y_post = "post"`).

## Value

An object of class `solomon_missing` with `by_cell` (counts by Solomon
cell), `counts` (totals by category), `pattern` (`"none"`,
`"structural"`, `"incidental"`, or `"mixed"`), and `guidance` (the
interpretation, supported response, and sources for each category
present). For a design with several treatments, `by_cell` has one row
for each of the 2(k + 1) groups, its `treat` column holds the condition,
and `conditions` names the control and the treatments.

## Details

Categories:

- **Structural pretest absence**: unpretested participants have no
  pretest. Unlike planned missing-data designs, in which unmeasured
  values exist and can be imputed (Graham et al., 2006), an imputed
  pretest here would describe a measurement that never occurred.

- **Incidental pretest missingness**: pretested participants without a
  pretest score. solomonR analyses currently use complete cases and
  warn. Because treatment is randomized, deterministic mean imputation
  of the pretest (without using treatment or outcome), with a
  missingness indicator when pretests may not be missing completely at
  random, retains these participants without biasing the treatment
  effect (White & Thompson, 2005; Groenwold et al., 2012).

- **Incidental posttest missingness**: missing outcomes. Complete-case
  analysis is unbiased when missingness is unrelated to the outcome
  given the variables in the model (Little & Rubin, 2019); attrition
  that differs across the groups should be reported.

- **Unexpected pretest scores**: pretest values recorded for unpretested
  participants, which usually indicate a coding or assignment error.
  They are ignored by solomonR analyses.

- **Unassigned participants**: missing treatment or pretest assignment.

Pretest categories are `NA` when `y_pre` is not supplied.

Designs with several treatments: give `treat` as a factor or character
vector of conditions and name the control with `control`. The counts are
then given for each of the 2(k + 1) groups of a design with k treatments
(Steyn, 2009): each treatment and the control, with and without a
pretest.

Not supported: this function does not test the missingness mechanism,
perform imputation, or provide sensitivity analyses for outcome
missingness that depends on unobserved values.

## References

Graham, J. W., Taylor, B. J., Olchowski, A. E., & Cumsille, P. E.
(2006). Planned missing data designs in psychological research.
*Psychological Methods, 11*(4), 323–343.
https://doi.org/10.1037/1082-989X.11.4.323

Groenwold, R. H. H., White, I. R., Donders, A. R. T., Carpenter, J. R.,
Altman, D. G., & Moons, K. G. M. (2012). Missing covariate data in
clinical research: When and when not to use the missing-indicator method
for analysis. *Canadian Medical Association Journal, 184*(11),
1265–1269. https://doi.org/10.1503/cmaj.110977

Little, R. J. A., & Rubin, D. B. (2019). *Statistical analysis with
missing data* (3rd ed.). Wiley. https://doi.org/10.1002/9781119482260

Mai, N. N., Takahashi, Y., & Oo, M. M. (2020). Testing the effectiveness
of transfer interventions using Solomon four-group designs. *Education
Sciences, 10*(4), Article 92. https://doi.org/10.3390/educsci10040092

Rubin, D. B. (1976). Inference and missing data. *Biometrika, 63*(3),
581–592. https://doi.org/10.1093/biomet/63.3.581

Solomon, R. L. (1949). An extension of control group design.
*Psychological Bulletin, 46*(2), 137–150.
https://doi.org/10.1037/h0062958

Steyn, R. (2009). Re-designing the Solomon four-group: Can we improve on
this exemplary model? *Design Principles and Practices: An International
Journal, 3*(1), 383–394.
https://doi.org/10.18848/1833-1874/CGP/v03i01/37588

White, I. R., & Thompson, S. G. (2005). Adjusting for partially missing
baseline measurements in randomized trials. *Statistics in Medicine,
24*(7), 993–1007. https://doi.org/10.1002/sim.1981

## See also

[`validate_solomon()`](https://juhalt.github.io/solomonR/reference/validate_solomon.md)

## Examples

``` r
data(solomon_example)

# Structural absence only: pretests are absent by design in Groups 3 and 4.
with(solomon_example, check_solomon_missing(y_post, treat, pretested, y_pre))
#> Solomon missingness check
#> Pattern: structural pretest absence only (expected in a Solomon design)
#> 
#>  Group                   Cell  n Post missing Pre absent (design) Pre missing
#>      1   Pretested, treatment 30            0                   0           0
#>      2     Pretested, control 30            0                   0           0
#>      3 Unpretested, treatment 30            0                  30           0
#>      4   Unpretested, control 30            0                  30           0
#>  Pre unexpected
#>               0
#>               0
#>               0
#>               0
#> 
#> Structural pretest absence (n = 60)
#>   Participants assigned to the unpretested groups were never pretested; the
#>   absence of a pretest is the experimental manipulation.
#>   Response: Do not impute. Use analyses that respect the design, such as
#>     fit_solomon_glm(), fit_solomon_ml(), fit_solomon_classic(), or SEM.
#>     Unlike planned missing-data designs, where unmeasured values exist and
#>     can be imputed, an imputed pretest here would describe a measurement that
#>     never occurred.
#>   Sources: Solomon (1949); Graham et al. (2006)

# Mixed: add a lost pretest and a missing posttest.
d <- solomon_example
d$y_pre[which(d$pretested == 1)[1]] <- NA
d$y_post[which(d$pretested == 0)[1]] <- NA
with(d, check_solomon_missing(y_post, treat, pretested, y_pre))
#> Solomon missingness check
#> Pattern: structural pretest absence and incidental missingness
#> 
#>  Group                   Cell  n Post missing Pre absent (design) Pre missing
#>      1   Pretested, treatment 30            0                   0           1
#>      2     Pretested, control 30            0                   0           0
#>      3 Unpretested, treatment 30            1                  30           0
#>      4   Unpretested, control 30            0                  30           0
#>  Pre unexpected
#>               0
#>               0
#>               0
#>               0
#> 
#> Structural pretest absence (n = 60)
#>   Participants assigned to the unpretested groups were never pretested; the
#>   absence of a pretest is the experimental manipulation.
#>   Response: Do not impute. Use analyses that respect the design, such as
#>     fit_solomon_glm(), fit_solomon_ml(), fit_solomon_classic(), or SEM.
#>     Unlike planned missing-data designs, where unmeasured values exist and
#>     can be imputed, an imputed pretest here would describe a measurement that
#>     never occurred.
#>   Sources: Solomon (1949); Graham et al. (2006)
#> 
#> Incidental pretest missingness (n = 1)
#>   Participants assigned to be pretested have no pretest score. Determine
#>   whether the pretest was administered but the score was lost, or never
#>   administered, which is a departure from the assigned pretest condition.
#>   Response: solomonR analyses currently use complete cases and warn. Because
#>     treatment is randomized, deterministic mean imputation of the pretest
#>     (without using treatment or outcome), with a missingness indicator when
#>     pretests may not be missing completely at random, retains these
#>     participants without biasing the treatment effect. If the pretest was
#>     never administered, report the departure and analyze participants as
#>     assigned.
#>   Sources: White & Thompson (2005); Groenwold et al. (2012)
#> 
#> Incidental posttest missingness (n = 1)
#>   Posttest scores are missing, so these participants are excluded from
#>   complete-case analyses.
#>   Response: Complete-case analysis is unbiased when missingness is unrelated
#>     to the outcome given the variables in the model. Report missingness by
#>     cell, because attrition that differs across the four groups can undermine
#>     the randomized comparisons, and consider sensitivity analyses if
#>     missingness may depend on the unobserved outcome.
#>   Sources: Rubin (1976); Little & Rubin (2019)

# A six-group design: two treatments and a control (Mai et al., 2020).
check_solomon_missing(post_behavior, condition, pretested, pre_behavior,
                      control = "Control", data = mai2020)
#> Solomon missingness check
#> Solomon N-group design: two treatments (RP, GS) and a control (Control), six groups
#> Pattern: structural pretest absence and incidental missingness
#> 
#>  Group                 Cell  n Post missing Pre absent (design) Pre missing
#>      1        Pretested, RP 35           11                   0           0
#>      2        Pretested, GS 33           10                   0           0
#>      3   Pretested, Control 50           23                   0           0
#>      4      Unpretested, RP 31            9                  31           0
#>      5      Unpretested, GS 27           12                  27           0
#>      6 Unpretested, Control 35           13                  35           0
#>  Pre unexpected
#>               0
#>               0
#>               0
#>               0
#>               0
#>               0
#> 
#> Structural pretest absence (n = 93)
#>   Participants assigned to the unpretested groups were never pretested; the
#>   absence of a pretest is the experimental manipulation.
#>   Response: Do not impute. Use analyses that respect the design, such as
#>     fit_solomon_glm(), fit_solomon_ml(), fit_solomon_classic(), or SEM.
#>     Unlike planned missing-data designs, where unmeasured values exist and
#>     can be imputed, an imputed pretest here would describe a measurement that
#>     never occurred.
#>   Sources: Solomon (1949); Graham et al. (2006)
#> 
#> Incidental posttest missingness (n = 78)
#>   Posttest scores are missing, so these participants are excluded from
#>   complete-case analyses.
#>   Response: Complete-case analysis is unbiased when missingness is unrelated
#>     to the outcome given the variables in the model. Report missingness by
#>     cell, because attrition that differs across the six groups can undermine
#>     the randomized comparisons, and consider sensitivity analyses if
#>     missingness may depend on the unobserved outcome.
#>   Sources: Rubin (1976); Little & Rubin (2019)
```
