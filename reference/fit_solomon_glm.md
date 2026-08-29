# Fit the unified GLM for a Solomon Four-Group design

Fit the unified GLM for a Solomon Four-Group design

## Usage

``` r
fit_solomon_glm(
  y,
  treat,
  pretested,
  pretest_score = NULL,
  covariates = NULL,
  robust = c("none", "HC3", "CR2"),
  cluster = NULL,
  family = stats::gaussian()
)
```

## Arguments

- y:

  numeric posttest vector

- treat:

  0/1 indicator (1 = treatment)

- pretested:

  0/1 indicator (1 = group received pretest)

- pretest_score:

  numeric vector for those pretested; NA for others

- covariates:

  optional data.frame of additional covariates

- robust:

  character: "none","HC3","CR2" (cluster-robust via clubSandwich if
  `cluster` supplied)

- cluster:

  optional clustering id (e.g., class/site)

- family:

  model family (default gaussian())

## Value

a list with model, tidy tables, and predefined contrasts (ATE,
interaction, simple effects)
