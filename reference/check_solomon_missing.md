# Distinguish structural and incidental missingness in a Solomon design

Classifies missing values in Solomon four-group data and explains the
supported response to each kind. Pretest scores are *structurally
absent* for participants assigned to the unpretested groups: withholding
the pretest is the experimental manipulation (Solomon, 1949), so those
values must never be imputed. Other missing values are *incidental* and
are handled according to the missing-data literature (Rubin, 1976;
Little & Rubin, 2019).

## Usage

``` r
check_solomon_missing(y_post, treat, pretested, y_pre = NULL)
```

## Arguments

- y_post:

  Numeric posttest scores.

- treat:

  Treatment indicator coded 0/1 (or logical).

- pretested:

  Pretest indicator coded 0/1 (or logical).

- y_pre:

  Optional numeric pretest scores, missing by design for unpretested
  participants.

## Value

An object of class `solomon_missing` with `by_cell` (counts by Solomon
cell), `counts` (totals by category), `pattern` (`"none"`,
`"structural"`, `"incidental"`, or `"mixed"`), and `guidance` (the
interpretation, supported response, and sources for each category
present).

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
  that differs across the four groups should be reported.

- **Unexpected pretest scores**: pretest values recorded for unpretested
  participants, which usually indicate a coding or assignment error.
  They are ignored by solomonR analyses.

- **Unassigned participants**: missing treatment or pretest assignment.

Pretest categories are `NA` when `y_pre` is not supplied.

Not supported: this function does not test the missingness mechanism,
perform imputation, or provide sensitivity analyses for outcome
missingness that depends on unobserved values.

## References

Graham, J. W., Taylor, B. J., Olchowski, A. E., & Cumsille, P. E.
(2006). Planned missing data designs in psychological research.
*Psychological Methods, 11*(4), 323-343.

Groenwold, R. H. H., White, I. R., Donders, A. R. T., Carpenter, J. R.,
Altman, D. G., & Moons, K. G. M. (2012). Missing covariate data in
clinical research: When and when not to use the missing-indicator method
for analysis. *Canadian Medical Association Journal, 184*(11),
1265-1269.

Little, R. J. A., & Rubin, D. B. (2019). *Statistical analysis with
missing data* (3rd ed.). Wiley.

Rubin, D. B. (1976). Inference and missing data. *Biometrika, 63*(3),
581-592.

Solomon, R. L. (1949). An extension of control group design.
*Psychological Bulletin, 46*(2), 137-150.

White, I. R., & Thompson, S. G. (2005). Adjusting for partially missing
baseline measurements in randomized trials. *Statistics in Medicine,
24*(7), 993-1007.

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
```
