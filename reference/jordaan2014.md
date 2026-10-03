# Group statistics from Jordaan's (2014) study with three posttest occasions

The sample sizes, means, and standard deviations of a Solomon four-group
design measured on three posttest occasions: an evaluation of a life
skills program for young adult male offenders (Jordaan, 2014, Table 7.3,
p. 112). The outcomes are the three subscales of the Coping Strategy
Indicator. Numbers reported in the publication are reused with citation.

## Usage

``` r
jordaan2014
```

## Format

A data frame with 42 rows, one per subscale, group, and occasion, and 8
variables:

- subscale:

  `Social support`, `Problem solving`, or `Avoidance`.

- group:

  The group's number in the thesis: 1 = program with pretest, 2 =
  program without pretest, 3 = control with pretest, 4 = control without
  pretest.

- treat:

  1 for the program groups, 0 for the control groups.

- pretested:

  1 if the group was pretested, 0 if not.

- occasion:

  `Pretest` (groups 1 and 3 only), `Posttest` (after the program),
  `Follow-up 1` (3 months later), or `Follow-up 2` (6 months later).

- n:

  Sample size.

- mean, sd:

  Mean and standard deviation.

## Source

Jordaan, J. (2014). *The development and evaluation of a life skills
programme for young adult prisoners* \[Doctoral thesis, University of
the Free State\]. KovsieScholar. https://hdl.handle.net/11660/832

## Details

**Design.** 120 offenders aged 21 to 25 with long sentences were
selected by systematic random sampling in a maximum-security
correctional center and assigned at random to the program or the control
condition. Half of each condition was assigned at random to be pretested
(pp. 86–87). The program ran for six months; the control group followed
the center's normal daily activities (p. 87). All groups were tested
after the program and again 3 and 6 months later (p. 98).

**Attrition.** Transfers removed 17 offenders from the program groups
and 7 from the control groups, leaving 96 (p. 109). Those 96 completed
every posttest: every analysis in the thesis has 92 error degrees of
freedom. Because more were lost from the program groups, the groups
analyzed are not guaranteed to be comparable, even though they were
randomized.

**Measures.** Each subscale of the Coping Strategy Indicator has 11
items scored from 1 (not at all) to 3 (a lot), so scores range from 11
to 33. High scores on problem solving and on seeking social support, and
low scores on avoidance, indicate better coping (p. 92).

**Known results.** These statistics reproduce, within rounding, the 2 x
2 analyses of variance of the posttests on each occasion (Tables
7.4–7.19, pp. 113–127). The Pretest x Treatment interaction F values
are:

- **Social support:** 9.678 (p = .002) after the program, 0.266 at 3
  months, and 2.306 at 6 months.

- **Problem solving:** 0.819, 0.563, and 5.556 (p = .021).

- **Avoidance:** 0.373, 0.688, and 0.327.

The thesis followed the decision sequence of Walton Braver and Braver
(1988) on each occasion separately
([`fit_solomon_classic()`](https://juhalt.github.io/solomonR/reference/fit_solomon_classic.md)).
The correlations between occasions are not reported, so analyses of
change across occasions cannot be reproduced; see the article "Worked
Example: Repeated Posttests".

## References

Walton Braver, M. C., & Braver, S. L. (1988). Statistical treatment of
the Solomon four-group design: A meta-analytic approach. *Psychological
Bulletin, 104*(1), 150–154. https://doi.org/10.1037/0033-2909.104.1.150

## See also

[`solomon_from_summary()`](https://juhalt.github.io/solomonR/reference/solomon_from_summary.md),
[`fit_solomon_mmrm()`](https://juhalt.github.io/solomonR/reference/fit_solomon_mmrm.md),
[steyn2005](https://juhalt.github.io/solomonR/reference/steyn2005.md)

## Examples

``` r
# Problem solving, 6 months after the program (Table 7.13, p. 121).
ps6 <- subset(jordaan2014, subscale == "Problem solving" & occasion == "Follow-up 2")
with(ps6, solomon_from_summary(n, mean, sd, treat = treat, pretested = pretested))
#> Solomon analysis from summary statistics
#> ----------------------------------------
#> Pooled error variance: 11.657 on 92 df (equal variances assumed)
#> 
#> Two-way ANOVA on the posttest (Type III sums of squares)
#>   Treatment            SS =    0.115  df = 1  F = 0.01  p = 0.921
#>   Pretest              SS =    0.059  df = 1  F = 0.01  p = 0.944
#>   Treatment x Pretest  SS =   64.539  df = 1  F = 5.54  p = 0.021
#>   Error                SS = 1072.442  df = 92
#> 
#> Contrasts with 95% confidence intervals
#>   Test A: Pretest x Treatment               3.320 [0.518, 6.122], t(92) = 2.35, p = 0.021
#>   Test B: Treatment | pretested             1.590 [-0.455, 3.635], t(92) = 1.54, p = 0.126
#>   Test C: Treatment | unpretested          -1.730 [-3.646, 0.186], t(92) = -1.79, p = 0.076
#>   Test D: ATE (avg over pretest)           -0.070 [-1.471, 1.331], t(92) = -0.10, p = 0.921
#>   Pretest main effect                      -0.050 [-1.451, 1.351], t(92) = -0.07, p = 0.944
```
