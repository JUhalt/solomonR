# reports for the teaching data are stable

    Code
      print(report_solomon(fit_solomon_glm(d$y_post, d$treat, d$pretested, d$y_pre),
      design = list(randomized = c(32, 31, 30, 33), prespecified = TRUE, measurement = "the same questionnaire in all four groups.")))
    Output
      The design was a Solomon four-group design (Solomon, 1949), with 30, 30, 30,
      and 30 participants analyzed in the pretested treatment, pretested control,
      unpretested treatment, and unpretested control groups, respectively. Of 32,
      31, 30, and 33 participants assigned to these groups, 30, 30, 30, and 30 were
      analyzed (attrition of 6.2%, 3.2%, 0.0%, and 9.1%, respectively). The
      analysis of pretest sensitization was pre-specified. Measurement: the same
      questionnaire in all four groups.
      
      Posttest outcomes were analyzed with a linear model containing treatment,
      pretesting, and their interaction, adjusting for the pretest score among
      pretested participants (Lin, 2013), with HC3 heteroskedasticity-consistent
      standard errors (MacKinnon & White, 1985; Long & Ervin, 2000).
      
      The average treatment effect across pretest conditions was 2.66, 95% CI
      [-0.47, 5.80], t(115) = 1.68, p = .095.
      The Pretest x Treatment interaction (pretest sensitization) was -1.94, 95% CI
      [-8.21, 4.33], t(115) = -0.61, p = .541.
      The treatment effect among pretested participants was 1.69, 95% CI [-2.76,
      6.15], t(115) = 0.75, p = .453.
      The treatment effect among unpretested participants was 3.63, 95% CI [-0.78,
      8.05], t(115) = 1.63, p = .106.
      
      References
      
      Lin, W. (2013). Agnostic notes on regression adjustments to experimental
          data: Reexamining Freedman's critique. The Annals of Applied Statistics,
          7(1), 295–318. https://doi.org/10.1214/12-AOAS583
      
      Long, J. S., & Ervin, L. H. (2000). Using heteroscedasticity consistent
          standard errors in the linear regression model. The American
          Statistician, 54(3), 217–224.
          https://doi.org/10.1080/00031305.2000.10474549
      
      MacKinnon, J. G., & White, H. (1985). Some heteroskedasticity-consistent
          covariance matrix estimators with improved finite sample properties.
          Journal of Econometrics, 29(3), 305–325.
          https://doi.org/10.1016/0304-4076(85)90158-7
      
      Solomon, R. L. (1949). An extension of control group design. Psychological
          Bulletin, 46(2), 137–150. https://doi.org/10.1037/h0062958
      

---

    Code
      print(report_solomon(fit_solomon_classic(d$y_post, d$treat, d$pretested, d$
      y_pre, flow = "1990")))
    Output
      The design was a Solomon four-group design (Solomon, 1949), with 30, 30, 30,
      and 30 participants analyzed in the pretested treatment, pretested control,
      unpretested treatment, and unpretested control groups, respectively.
      
      The historical Tests A-I sequence was followed (Walton Braver & Braver,
      1988), as amended by Braver and Walton Braver (1990).
      
      The Pretest x Treatment interaction (Test A) was -1.90, F(1, 116) = 0.29, p =
      .589.
      The treatment main effect (Test D) was 2.68, F(1, 116) = 2.34, p = .129.
      The ANCOVA treatment effect in the pretested groups (Test E) was 1.69, F(1,
      57) = 0.59, p = .444.
      The posttest-only treatment effect (Test H) was 3.63, t(58) = 1.66, p = .103.
      Test I, the Stouffer combination (Stouffer et al., 1949; Walton Braver &
      Braver, 1988), gave z = 1.69, p = .090. The experiment-wise Type I error of
      this sequence exceeds its nominal level (Sawilowsky et al., 1994).
      The 1990 sequence ended at Test I, which the 1990 amendment regards as the
      most definitive test.
      In the history/maturation check (Campbell & Stanley, 1963/1966; Mai et al.,
      2020), the unpretested control posttest differed from the treated-group
      pretest by 1.47, t(58) = 0.54, p = .590, and from the control-group pretest
      by 1.53, t(58) = 0.61, p = .543.
      
      References
      
      Braver, S. L., & Walton Braver, M. C. (1990). Meta-analysis for Solomon
          four-group designs reconsidered: A reply to Sawilowsky and Markman.
          Perceptual and Motor Skills, 71(1), 321–322.
          https://doi.org/10.2466/pms.1990.71.1.321
      
      Campbell, D. T., & Stanley, J. C. (1966). Experimental and quasi-experimental
          designs for research. Rand McNally. (Original work published 1963)
      
      Huck, S. W., & Sandler, H. M. (1973). A note on the Solomon 4-group design:
          Appropriate statistical analyses. The Journal of Experimental Education,
          42(2), 54–55. https://doi.org/10.1080/00220973.1973.11011460
      
      Mai, N. N., Takahashi, Y., & Oo, M. M. (2020). Testing the effectiveness of
          transfer interventions using Solomon four-group designs. Education
          Sciences, 10(4), Article 92. https://doi.org/10.3390/educsci10040092
      
      Sawilowsky, S. S., Kelley, D. L., Blair, R. C., & Markman, B. S. (1994).
          Meta-analysis and the Solomon four-group design. The Journal of
          Experimental Education, 62(4), 361–376.
          https://doi.org/10.1080/00220973.1994.9944140
      
      Solomon, R. L. (1949). An extension of control group design. Psychological
          Bulletin, 46(2), 137–150. https://doi.org/10.1037/h0062958
      
      Stouffer, S. A., Suchman, E. A., DeVinney, L. C., Star, S. A., & Williams, R.
          M., Jr. (1949). The American soldier: Adjustment during army life (Vol.
          1). Princeton University Press.
      
      Walton Braver, M. C., & Braver, S. L. (1988). Statistical treatment of the
          Solomon four-group design: A meta-analytic approach. Psychological
          Bulletin, 104(1), 150–154. https://doi.org/10.1037/0033-2909.104.1.150
      

---

    Code
      print(report_solomon(solomon_from_summary(c(9, 25, 17, 37), c(10.94, 7.8, 8.94,
        9.35), c(2.26, 2.29, 1.98, 2.11))))
    Output
      The design was a Solomon four-group design (Solomon, 1949), with 9, 25, 17,
      and 37 participants analyzed in the pretested treatment, pretested control,
      unpretested treatment, and unpretested control groups, respectively.
      
      The published cell statistics were reanalyzed with a two-way analysis of
      variance of the posttest (Type III sums of squares), the model behind Tests
      A-D (Walton Braver & Braver, 1988).
      
      The Pretest x Treatment interaction was F(1, 84) = 11.46, p = .001.
      The treatment main effect was F(1, 84) = 6.78, p = .011, and the pretest main
      effect F(1, 84) = 0.18, p = .669.
      The treatment effect among pretested participants was 3.14, 95% CI [1.47,
      4.81], t(84) = 3.75, p < .001.
      The treatment effect among unpretested participants was -0.41, 95% CI [-1.67,
      0.85], t(84) = -0.65, p = .518.
      
      References
      
      Solomon, R. L. (1949). An extension of control group design. Psychological
          Bulletin, 46(2), 137–150. https://doi.org/10.1037/h0062958
      
      Walton Braver, M. C., & Braver, S. L. (1988). Statistical treatment of the
          Solomon four-group design: A meta-analytic approach. Psychological
          Bulletin, 104(1), 150–154. https://doi.org/10.1037/0033-2909.104.1.150
      

