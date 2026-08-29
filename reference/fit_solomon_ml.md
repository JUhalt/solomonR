# Full-information ML analysis for a Solomon Four-Group Design

Fits the maximum-likelihood regression model described by van
Engelenburg (1999). Pretest information is incorporated for the
pretested groups while structurally missing pretests in the unpretested
groups are handled through a separate residual variance.

## Usage

``` r
fit_solomon_ml(
  y_post,
  treat,
  pretested,
  y_pre,
  weights = c("equal"),
  control = list()
)
```

## Arguments

- y_post:

  Numeric posttest scores.

- treat:

  Treatment indicator coded 0 = control and 1 = treatment.

- pretested:

  Pretest indicator coded 0 = unpretested and 1 = pretested.

- y_pre:

  Numeric pretest scores. These should be missing by design for
  participants assigned to the unpretested groups.

- weights:

  Character. How to define the average treatment effect across pretest
  conditions. Currently `"equal"` gives equal weight to the pretested
  and unpretested treatment effects.

- control:

  Optional list passed to
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html).

## Value

An object of class `solomon_ml`.

## Details

The model estimates the treatment effect, pretest effect, Treatment x
Pretest interaction, pretest-posttest slope, and separate residual
standard deviations for pretested and unpretested participants.

## References

van Engelenburg, G. (1999). Statistical analysis for the Solomon
four-group design. University of Twente Research Report 99-06.
