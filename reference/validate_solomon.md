# Validate the structure and coding of a Solomon four-group design

**\[stable\]** Checks that data can support a Solomon four-group
analysis before a model is fitted: equal input lengths, 0/1 coding of
the design indicators, all four cells present, enough observed outcomes
per cell, and the distinction between structurally absent and
incidentally missing values (see
[`check_solomon_missing()`](https://juhalt.github.io/solomonR/reference/check_solomon_missing.md)).
Problems are returned as a table of issues rather than stopping at the
first one, so every problem is reported at once.

## Usage

``` r
validate_solomon(
  y_post,
  treat,
  pretested,
  y_pre = NULL,
  min_cell_n = 2,
  cluster = NULL,
  data = NULL
)
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

- cluster:

  Optional cluster identifier (for example, class, school, or site), one
  value per participant.

- data:

  Optional data frame. When supplied, the other data arguments are
  looked up in it first, as bare column names (`y_post = post`) or as
  strings (`y_post = "post"`).

## Value

An object of class `solomon_validation` with `valid` (`TRUE` when no
errors were found), `issues` (severity, check, and message), `cells`
(counts by cell, with clusters and cluster sizes when `cluster` is
supplied), and `missing` (the
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

Clustered designs: when participants belong to classes, schools, or
sites, supply `cluster`. The number of clusters and the range of cluster
sizes are then reported for each cell, together with whether treatment
and pretesting were assigned to whole clusters or varied within them.

A cell that contains a single cluster is an error. Class and condition
are then completely confounded, as when each Solomon condition is one
intact class (El Karkri et al., 2025a): a difference between the classes
cannot be separated from the treatment, pretest, or sensitization
effect, and no between-cluster variability can be estimated.
Cluster-robust inference (Pustejovsky & Tipton, 2018) needs several
clusters in every cell. For example, Kvalem et al. (1996) randomized 124
school classes to the four Solomon conditions. When whole clusters are
randomized, a cell with two or three clusters is a warning: Hayes and
Moulton (2017, p. 128) regard four clusters per arm as an absolute
minimum.

## References

El Karkri, M., Quesada, A., & Romero-Ariza, M. (2025a). The dual impact
of pretest sensitisation and the cognitive acceleration through science
education programme in the Solomon four-group design. *Brain Sciences,
16*(1), Article 64. https://doi.org/10.3390/brainsci16010064

Hayes, R. J., & Moulton, L. H. (2017). *Cluster randomised trials* (2nd
ed.). Chapman and Hall/CRC. https://doi.org/10.4324/9781315370286

Kvalem, I. L., Sundet, J. M., Rivø, K. I., Eilertsen, D. E., &
Bakketeig, L. S. (1996). The effect of sex education on adolescents' use
of condoms: Applying the Solomon four-group design. *Health Education
Quarterly, 23*(1), 34–47. https://doi.org/10.1177/109019819602300103

Pustejovsky, J. E., & Tipton, E. (2018). Small-sample methods for
cluster-robust variance estimation and hypothesis testing in fixed
effects models. *Journal of Business & Economic Statistics, 36*(4),
672–683. https://doi.org/10.1080/07350015.2016.1247004

Solomon, R. L. (1949). An extension of control group design.
*Psychological Bulletin, 46*(2), 137–150.
https://doi.org/10.1037/h0062958

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

# One intact class per condition confounds class with condition.
one_class <- with(solomon_example, 2 * pretested + treat)
with(solomon_example, validate_solomon(y_post, treat, pretested, y_pre,
                                       cluster = one_class))
#> Solomon design validation: errors found
#> 
#>  Group                   Cell  n Post missing Pre absent (design) Pre missing
#>      1   Pretested, treatment 30            0                   0           0
#>      2     Pretested, control 30            0                   0           0
#>      3 Unpretested, treatment 30            0                  30           0
#>      4   Unpretested, control 30            0                  30           0
#>  Pre unexpected Clusters
#>               0        1
#>               0        1
#>               0        1
#>               0        1
#> 
#> [ERROR] A single cluster makes up each of these cells: Pretested, treatment;
#>   Pretested, control; Unpretested, treatment; Unpretested, control. Cluster
#>   and condition are completely confounded, so a difference between clusters
#>   cannot be separated from the Solomon effects and cluster-robust inference
#>   is unavailable.
#> [NOTE] Cell sizes range from 30 to 30.
#> [NOTE] 4 clusters; 1 to 1 per cell, with 30 to 30 participants per cluster.
#>   Treatment is constant within every cluster (assigned to whole clusters).
#>   Pretesting is constant within every cluster (assigned to whole clusters).
```
