# Demo Solomon four-group data set (second example)

A synthetic data set with 25 observations in each of the four
combinations of treatment and pretest assignment. It is kept as a second
example because two of its features make different analyses of the same
effect disagree:

## Format

A data frame with 100 rows and 4 variables:

- y_post:

  Numeric posttest score.

- treat:

  Treatment indicator: 0 = control, 1 = treatment.

- pretested:

  Pretest indicator: 0 = not pretested, 1 = pretested.

- y_pre:

  Numeric pretest score. Structurally missing for participants assigned
  to the unpretested groups.

## Source

Synthetic data generated for examples and package documentation.

## Details

- among pretested participants, the pretest and posttest are negatively
  correlated (r = -0.14); and

- the pretested treatment group's mean pretest score is about 4.5 points
  higher than the pretested control group's.

As a result, unadjusted, ANCOVA, and gain-score estimates of the
treatment effect among pretested participants differ markedly (see
[`compare_solomon_methods()`](https://juhalt.github.io/solomonR/reference/compare_solomon_methods.md)
and
[`vignette("getting-started")`](https://juhalt.github.io/solomonR/articles/getting-started.md)).
The parameters used to generate these data were not recorded, so their
true effects are unknown. For the primary example with documented true
values, see
[solomon_example](https://juhalt.github.io/solomonR/reference/solomon_example.md).

## See also

[solomon_example](https://juhalt.github.io/solomonR/reference/solomon_example.md)
