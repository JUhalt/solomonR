# Historical categorical analysis of a binary Solomon outcome

**\[stable\]** Reproduces the categorical path described by El Karkri et
al. (2025b) for qualitative outcomes: the treatment comparison is tested
separately among pretested and unpretested participants with Fisher's
exact test (Pearson's chi-square test is also reported), and pretest
sensitization is "judged to exist if a significant effect is observed
for pre-tested groups but not for non-pre-tested groups" (p. 7).

## Usage

``` r
fisher_solomon(
  y_post,
  treat,
  pretested,
  alpha = 0.05,
  data = NULL,
  y = deprecated()
)
```

## Arguments

- y_post:

  Binary posttest outcome coded 0/1 (or logical).

- treat:

  Treatment indicator coded 0/1 (or logical). Designs with several
  treatments are not supported; see
  [`fit_solomon_glm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_glm.md).

- pretested:

  Pretest indicator coded 0/1 (or logical).

- alpha:

  Significance level for the historical rule. Default 0.05.

- data:

  Optional data frame. When supplied, the other data arguments are
  looked up in it first, as bare column names (`y_post = post`) or as
  strings (`y_post = "post"`).

- y:

  **\[deprecated\]** Use `y_post`.

## Value

An object of class `solomon_fisher` with `tests` (one row per pretest
condition and for both combined: counts, proportions, the uncorrected
Pearson chi-square, and Fisher's exact p-value) and `sensitization` (the
historical rule's verdict).

## Details

This is provided for teaching and replication. The rule compares
significance, not effects: when the treatment effect is the same in both
pretest conditions, it still declares sensitization whenever the
pretested comparison reaches significance and the unpretested one does
not. As Gelman and Stern (2006) put it, "even large changes in
significance levels can correspond to small, nonsignificant changes in
the underlying quantities" (p. 328). A test of sensitization compares
the effects themselves; see
[`marginal_solomon()`](https://juhalt.github.io/solomonR/reference/marginal_solomon.md).

## References

El Karkri, M., Quesada, A., & Romero-Ariza, M. (2025b). Methodological
aspects of the Solomon four-group design: Detecting pre-test
sensitisation and analysing qualitative and quantitative variables in
education research. *Review of Education, 13*(1), Article e70050.
https://doi.org/10.1002/rev3.70050

Gelman, A., & Stern, H. (2006). The difference between "significant" and
"not significant" is not itself statistically significant. *The American
Statistician, 60*(4), 328–331. https://doi.org/10.1198/000313006X152649

## See also

[`marginal_solomon()`](https://juhalt.github.io/solomonR/reference/marginal_solomon.md)

## Examples

``` r
# Kvalem et al. (1996): condom use at most recent intercourse, 6 months.
kvalem <- data.frame(
  treat = c(1, 1, 0, 0), pretested = c(1, 0, 1, 0),
  events = c(51, 21, 76, 69), n = c(73, 49, 148, 133)
)
d <- kvalem[rep(1:4, kvalem$n), c("treat", "pretested")]
d$y <- unlist(lapply(1:4, function(i) {
  rep(c(1, 0), c(kvalem$events[i], kvalem$n[i] - kvalem$events[i]))
}))
fisher_solomon(y, treat, pretested, data = d)
#> Historical categorical Solomon analysis (El Karkri et al., 2025b)
#> 
#>    Condition Events (T) n (T) Events (C) n (C) Chi-square Fisher p
#>    Pretested         51    73         76   148       6.85   0.0095
#>  Unpretested         21    49         69   133       1.17   0.3180
#>     Combined         72   122        145   281       1.88   0.1922
#> 
#> Historical rule (significant among pretested, not among unpretested, alpha = 0.05): sensitization declared
#> Caution: this rule compares significance, not effects, and is not a test of
#> the Pretest x Treatment interaction. See marginal_solomon().
```
