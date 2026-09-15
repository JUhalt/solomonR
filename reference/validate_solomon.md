# Validate the structure and coding of a Solomon four-group design

Checks that data can support a Solomon four-group analysis before a
model is fitted: equal input lengths, 0/1 coding of the design
indicators, all four cells present, enough observed outcomes per cell,
and the distinction between structurally absent and incidentally missing
values (see
[`check_solomon_missing()`](https://juhalt.github.io/solomonR/reference/check_solomon_missing.md)).
Problems are returned as a table of issues rather than stopping at the
first one, so every problem is reported at once.

## Usage

``` r
validate_solomon(y_post, treat, pretested, y_pre = NULL, min_cell_n = 2)
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

- min_cell_n:

  Minimum number of observed posttest scores required in each cell (at
  least 2).

## Value

An object of class `solomon_validation` with `valid` (`TRUE` when no
errors were found), `issues` (severity, check, and message), `cells`
(counts by cell), and `missing` (the
[`check_solomon_missing()`](https://juhalt.github.io/solomonR/reference/check_solomon_missing.md)
result).

## Details

Accepted coding: `treat` and `pretested` must be numeric 0/1 or logical.
Factors and character codes are rejected so that group membership is
never inferred from level order. The Solomon design requires all four
cells (Solomon, 1949): pretested treatment, pretested control,
unpretested treatment, and unpretested control.

Severity:

- **error**: the data cannot support a Solomon analysis as supplied
  (unequal lengths, invalid coding, an empty cell, fewer than
  `min_cell_n` observed posttest scores in a cell, or no observed
  pretests among pretested participants when `y_pre` is supplied).

- **warning**: analyses can run, but some participants are excluded or
  values are inconsistent with the design (unassigned participants,
  incidental missingness, pretest scores recorded for unpretested
  participants).

- **note**: information for reporting, such as the range of cell sizes.

`min_cell_n` is a technical minimum, not a sample-size recommendation:
at least two observations per cell are needed to estimate within-cell
variability, and HC3 standard errors can be undefined for a cell with a
single observation. Plan cell sizes with a power analysis.

## References

Solomon, R. L. (1949). An extension of control group design.
*Psychological Bulletin, 46*(2), 137-150.

## See also

[`check_solomon_missing()`](https://juhalt.github.io/solomonR/reference/check_solomon_missing.md)

## Examples

``` r
data(solomon_example)

# A valid design.
with(solomon_example, validate_solomon(y_post, treat, pretested, y_pre))
#> Solomon design validation: no errors found
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
#> [NOTE] Cell sizes range from 30 to 30.

# An empty cell makes the design invalid.
no_group_3 <- solomon_example[!(solomon_example$pretested == 0 & solomon_example$treat == 1), ]
with(no_group_3, validate_solomon(y_post, treat, pretested, y_pre))
#> Solomon design validation: errors found
#> 
#>  Group                   Cell  n Post missing Pre absent (design) Pre missing
#>      1   Pretested, treatment 30            0                   0           0
#>      2     Pretested, control 30            0                   0           0
#>      3 Unpretested, treatment  0            0                   0           0
#>      4   Unpretested, control 30            0                  30           0
#>  Pre unexpected
#>               0
#>               0
#>               0
#>               0
#> 
#> [ERROR] The Solomon design requires all four cells; no participants are in:
#>   Unpretested, treatment.
```
