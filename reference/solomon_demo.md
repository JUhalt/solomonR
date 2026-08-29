# Demo Solomon Four-Group Dataset

A synthetic dataset illustrating a Solomon four-group experimental
design. The dataset contains 25 observations in each of the four
combinations of treatment and pretest assignment.

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
