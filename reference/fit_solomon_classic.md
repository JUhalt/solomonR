# Historical Solomon Four-Group Analysis

Implements the historical Test A-I framework associated with analysis of
the Solomon four-group design. All tests are calculated when possible,
while the historical decision pathway is stored separately.

## Usage

``` r
fit_solomon_classic(
  y_post,
  treat,
  pretested,
  y_pre,
  alpha = 0.05,
  pretested_test = c("ancova", "gain", "repeated"),
  combine_with_stouffer = TRUE,
  stouffer_direction = c("greater", "less")
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

  Numeric pretest scores. Values should be missing for participants
  assigned to the unpretested groups.

- alpha:

  Significance level used to reconstruct the historical decision
  pathway. Default is 0.05.

- pretested_test:

  Which historical pretested-group analysis should be followed in the
  decision pathway: `"ancova"`, `"gain"`, or `"repeated"`. All three are
  still calculated and returned.

- combine_with_stouffer:

  Logical. If `TRUE`, include historical Test I in the decision pathway
  when earlier treatment tests are nonsignificant.

- stouffer_direction:

  Direction of the historical one-tailed treatment hypothesis used for
  Test I: `"greater"` or `"less"`.

## Value

An object of class `solomon_classic`. The `tests` component contains
Tests A-I, while `path` records the historical decision sequence for the
observed data.

## Details

The historical sequence includes the four-group posttest factorial
model, simple treatment effects, ANCOVA, gain-score analysis, the
equivalent two-wave repeated-measures interaction, the posttest-only
comparison, and the optional Stouffer meta-analytic combination.

Test I is included for historical replication and teaching. Later work
raised concerns about experiment-wise Type I error when the procedure is
used conditionally. Its presence in this function should not be
interpreted as a general contemporary recommendation.
