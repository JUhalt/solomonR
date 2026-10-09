# a partial-invariance fit is reported with its freed parameter and fit details

    Code
      print_masked(r)
    Output
      The design was a Solomon four-group design (Solomon, 1949), with 150, 150,
      150, and 150 participants analyzed in the pretested treatment, pretested
      control, unpretested treatment, and unpretested control groups, respectively.
      
      Latent posttest means were compared across the four Solomon groups in a
      multiple-group structural equation model in which four posttest indicators
      (y1, y2, y3, and y4) measured one latent posttest, fitted with lavaan
      (Version <version>; Rosseel, 2012) by maximum likelihood with Huber-White
      robust standard errors and a scaled test statistic asymptotically equal to
      the Yuan-Bentler statistic (MLR; Yuan & Bentler, 2000), using
      full-information maximum likelihood for missing indicator values. The
      loadings and intercepts of the indicators were constrained to be equal across
      the four groups, the scalar measurement invariance that latent mean
      comparisons require (Meredith, 1993; Vandenberg & Lance, 2000), except the
      intercept of y4, which was estimated separately in each group (partial
      invariance; Byrne et al., 1989). For identification, the latent posttest mean
      of the unpretested control group was fixed at 0, and the latent variance at 1
      in the pretested treated group; the contrasts are therefore in standard
      deviations of the latent posttest in that group. Measurement invariance was
      tested by comparing configural, metric, and scalar models in sequence, with
      the same parameters freed (Vandenberg & Lance, 2000); the scaled chi-square
      difference test (Satorra & Bentler, 2001) and the change in fit (Chen, 2007)
      both supported partial scalar invariance. The contrasts between latent means
      were tested with Wald z tests.
      
      The four-group model gave a scaled χ²(23) = 34.52, p = .058, CFI = .991,
      RMSEA = .059, and SRMR = .036; the CFI and RMSEA are computed from the
      unscaled maximum likelihood χ²(23) = 34.95.
      The average treatment effect across pretest conditions was 0.36, 95% CI
      [0.20, 0.52], z = 4.38, p < .001. The Pretest x Treatment interaction
      (pretest sensitization) was 0.01, 95% CI [-0.30, 0.32], z = 0.07, p = .945.
      The treatment effect among pretested participants was 0.37, 95% CI [0.13,
      0.60], z = 3.05, p = .002. The treatment effect among unpretested
      participants was 0.36, 95% CI [0.14, 0.57], z = 3.28, p = .001. The pretest
      effect on the latent posttest (pretested minus unpretested participants) was
      -0.06 among control participants, 95% CI [-0.28, 0.15], z = -0.59, p = .556,
      and -0.05 among treated participants, 95% CI [-0.28, 0.17], z = -0.46, p =
      .645; their average, the pretest main effect, was -0.06, 95% CI [-0.21,
      0.10], z = -0.74, p = .459.
      
      References
      
      Byrne, B. M., Shavelson, R. J., & Muthén, B. (1989). Testing for the
          equivalence of factor covariance and mean structures: The issue of
          partial measurement invariance. Psychological Bulletin, 105(3), 456–466.
          https://doi.org/10.1037/0033-2909.105.3.456
      
      Chen, F. F. (2007). Sensitivity of goodness of fit indexes to lack of
          measurement invariance. Structural Equation Modeling: A Multidisciplinary
          Journal, 14(3), 464–504. https://doi.org/10.1080/10705510701301834
      
      Meredith, W. (1993). Measurement invariance, factor analysis and factorial
          invariance. Psychometrika, 58(4), 525–543.
          https://doi.org/10.1007/BF02294825
      
      Rosseel, Y. (2012). lavaan: An R package for structural equation modeling.
          Journal of Statistical Software, 48(2), 1–36.
          https://doi.org/10.18637/jss.v048.i02
      
      Satorra, A., & Bentler, P. M. (2001). A scaled difference chi-square test
          statistic for moment structure analysis. Psychometrika, 66(4), 507–514.
          https://doi.org/10.1007/BF02296192
      
      Solomon, R. L. (1949). An extension of control group design. Psychological
          Bulletin, 46(2), 137–150. https://doi.org/10.1037/h0062958
      
      Vandenberg, R. J., & Lance, C. E. (2000). A review and synthesis of the
          measurement invariance literature: Suggestions, practices, and
          recommendations for organizational research. Organizational Research
          Methods, 3(1), 4–70. https://doi.org/10.1177/109442810031002
      
      Yuan, K.-H., & Bentler, P. M. (2000). Three likelihood-based methods for mean
          and covariance structure analysis with nonnormal missing data.
          Sociological Methodology, 30(1), 165–200.
          https://doi.org/10.1111/0081-1750.00078
      

# a fit whose invariance check fails says so instead of claiming scalar invariance

    Code
      print_masked(r)
    Output
      The design was a Solomon four-group design (Solomon, 1949), with 150, 150,
      150, and 150 participants analyzed in the pretested treatment, pretested
      control, unpretested treatment, and unpretested control groups, respectively.
      
      Latent posttest means were compared across the four Solomon groups in a
      multiple-group structural equation model in which four posttest indicators
      (y1, y2, y3, and y4) measured one latent posttest, fitted with lavaan
      (Version <version>; Rosseel, 2012) by maximum likelihood with Huber-White
      robust standard errors and a scaled test statistic asymptotically equal to
      the Yuan-Bentler statistic (MLR; Yuan & Bentler, 2000), using
      full-information maximum likelihood for missing indicator values. The
      loadings and intercepts of the indicators were constrained to be equal across
      the four groups, the scalar measurement invariance that latent mean
      comparisons require (Meredith, 1993; Vandenberg & Lance, 2000). For
      identification, the latent posttest mean of the unpretested control group was
      fixed at 0, and the latent variance at 1 in the pretested treated group; the
      contrasts are therefore in standard deviations of the latent posttest in that
      group. Measurement invariance was tested by comparing configural, metric, and
      scalar models in sequence (Vandenberg & Lance, 2000): the scaled chi-square
      difference test (Satorra & Bentler, 2001) and the change in fit (Chen, 2007)
      both supported only metric invariance, so neither criterion supported the
      scalar invariance that the latent mean contrasts assume. The contrasts
      between latent means were tested with Wald z tests.
      
      The four-group model gave a scaled χ²(26) = 173.45, p < .001, CFI = .883,
      RMSEA = .195, and SRMR = .104; the CFI and RMSEA are computed from the
      unscaled maximum likelihood χ²(26) = 174.54.
      The average treatment effect across pretest conditions was 0.31, 95% CI
      [0.15, 0.46], z = 3.80, p < .001. The Pretest x Treatment interaction
      (pretest sensitization) was 0.16, 95% CI [-0.15, 0.47], z = 0.99, p = .320.
      The treatment effect among pretested participants was 0.38, 95% CI [0.15,
      0.62], z = 3.26, p = .001. The treatment effect among unpretested
      participants was 0.23, 95% CI [0.01, 0.44], z = 2.10, p = .036. The pretest
      effect on the latent posttest (pretested minus unpretested participants) was
      -0.19 among control participants, 95% CI [-0.41, 0.03], z = -1.72, p = .085,
      and -0.03 among treated participants, 95% CI [-0.25, 0.19], z = -0.28, p =
      .778; their average, the pretest main effect, was -0.11, 95% CI [-0.27,
      0.04], z = -1.41, p = .158.
      
      References
      
      Chen, F. F. (2007). Sensitivity of goodness of fit indexes to lack of
          measurement invariance. Structural Equation Modeling: A Multidisciplinary
          Journal, 14(3), 464–504. https://doi.org/10.1080/10705510701301834
      
      Meredith, W. (1993). Measurement invariance, factor analysis and factorial
          invariance. Psychometrika, 58(4), 525–543.
          https://doi.org/10.1007/BF02294825
      
      Rosseel, Y. (2012). lavaan: An R package for structural equation modeling.
          Journal of Statistical Software, 48(2), 1–36.
          https://doi.org/10.18637/jss.v048.i02
      
      Satorra, A., & Bentler, P. M. (2001). A scaled difference chi-square test
          statistic for moment structure analysis. Psychometrika, 66(4), 507–514.
          https://doi.org/10.1007/BF02296192
      
      Solomon, R. L. (1949). An extension of control group design. Psychological
          Bulletin, 46(2), 137–150. https://doi.org/10.1037/h0062958
      
      Vandenberg, R. J., & Lance, C. E. (2000). A review and synthesis of the
          measurement invariance literature: Suggestions, practices, and
          recommendations for organizational research. Organizational Research
          Methods, 3(1), 4–70. https://doi.org/10.1177/109442810031002
      
      Yuan, K.-H., & Bentler, P. M. (2000). Three likelihood-based methods for mean
          and covariance structure analysis with nonnormal missing data.
          Sociological Methodology, 30(1), 165–200.
          https://doi.org/10.1111/0081-1750.00078
      

