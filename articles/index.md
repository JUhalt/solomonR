# Articles

### Decide

- [Should I Use a Solomon
  Design?](https://juhalt.github.io/solomonR/articles/should-i-use-a-solomon-design.md):

  When pretest sensitization is a real threat, what the Solomon
  four-group design costs, and what the alternatives are.

### Plan

- [Planning a Solomon
  Study](https://juhalt.github.io/solomonR/articles/planning.md):

  Choosing sample sizes for a Solomon four-group study with
  plan_solomon(), power_solomon(), and plot_power_solomon(), using
  planning values from the published literature.

### Analyze

- [Getting Started: Analyzing a Solomon Four-Group
  Study](https://juhalt.github.io/solomonR/articles/getting-started.md):

- [Modern Analysis of the Solomon Four-Group
  Design](https://juhalt.github.io/solomonR/articles/glm-solomon.md):

- [Solomon Methods: History, Recommendations, and
  Extensions](https://juhalt.github.io/solomonR/articles/solomon-methods.md):

- [Structural Equation Models for Solomon
  Designs](https://juhalt.github.io/solomonR/articles/sem.md):

  Observed- and latent-variable SEM for the Solomon four-group design:
  the four-group mean structure, the pretested-groups ANCOVA, latent
  outcomes, and measurement invariance across the Solomon groups.

- [Missing and Repeated
  Posttests](https://juhalt.github.io/solomonR/articles/missing-and-repeated.md):

  How to analyze a Solomon study when some posttests are missing, and
  when the posttest is measured on several occasions and participants
  drop out.

- [Worked Example: A Published Solomon
  Study](https://juhalt.github.io/solomonR/articles/worked-example.md):

  The full analysis path, from checking the design to drafting the
  report, on the published individual data of Mai et al. (2020), checked
  at each step against the published results.

### Report

- [Reporting a Solomon
  Study](https://juhalt.github.io/solomonR/articles/reporting.md):

  What to report from a Solomon four-group study, and how
  report_solomon() drafts the results, the design statement, and the
  references.

### Synthesize

- [Reanalysis and
  Synthesis](https://juhalt.github.io/solomonR/articles/synthesis.md):

  Reanalyzing a published Solomon study from its summary statistics, and
  computing effect sizes for a meta-analysis.

### Teach

- [Teaching with
  solomonR](https://juhalt.github.io/solomonR/articles/teaching.md):

  Six lessons with exercises and solutions for teaching the Solomon
  four-group design: what each group contributes, sensitization, the
  design’s history, the historical test sequence, absence of evidence,
  and real studies.

### History

- [A History of the Solomon Design and Its
  Analysis](https://juhalt.github.io/solomonR/articles/history.md):

  From Solomon (1949) to the MERIT recommendations: how the four-group
  design was proposed, how its analysis was debated, and which
  historical procedures solomonR reproduces.

- [The Classic Solomon Four-Group
  Analysis](https://juhalt.github.io/solomonR/articles/classic-solomon.md):

### Validation evidence

- [Validation
  Evidence](https://juhalt.github.io/solomonR/articles/validation-evidence.md):

- [Validating fit_solomon_ml(): A Simulation
  Study](https://juhalt.github.io/solomonR/articles/ml-validation.md):

- [Validating power_solomon(): A Simulation
  Study](https://juhalt.github.io/solomonR/articles/power-validation.md):

- [Checking plan_solomon(): Re-simulating Planned
  Designs](https://juhalt.github.io/solomonR/articles/plan-validation.md):

- [Binary Outcomes: Validating
  marginal_solomon()](https://juhalt.github.io/solomonR/articles/binary-validation.md):

  Simulation study of the marginal Solomon contrasts for binary
  outcomes, the conditional logistic contrast, and the historical
  categorical rule.

- [Count Outcomes: Validating Poisson Fits and
  marginal_solomon()](https://juhalt.github.io/solomonR/articles/count-validation.md):

  Simulation study of the Solomon contrasts for count outcomes, with and
  without overdispersion, from fit_solomon_glm(family = poisson()) and
  marginal_solomon().

- [Clustered Designs: Validating Cluster-Level Randomization
  Inference](https://juhalt.github.io/solomonR/articles/cluster-validation.md):

  Simulation study of perm_solomon() for cluster-randomized Solomon
  designs, with raw and studentized cluster-level statistics, compared
  with CR2 tests from fit_solomon_glm().

- [Clustered Designs: Validating Marginal Risk
  Contrasts](https://juhalt.github.io/solomonR/articles/cluster-marginal-validation.md):

  Simulation study of marginal_solomon() on cluster-robust (CR2)
  logistic fits of cluster-randomized Solomon designs, with a
  cluster-level comparator.

- [Historical Tests: Replicating the Published Error
  Rates](https://juhalt.github.io/solomonR/articles/classic-validation.md):

  Replication of the published Type I error rates of the Walton Braver
  and Braver sequence of Tests A-I, its 1995 revision, and Sawilowsky’s
  alpha allocations, with fit_solomon_classic().

- [Latent Contrasts: Validating the Invariance
  Check](https://juhalt.github.io/solomonR/articles/invariance-validation.md):

  Simulation study of which measurement-invariance criterion should
  govern latent Solomon mean contrasts, how often each criterion falsely
  rejects or detects noninvariance in Solomon-sized groups, and how
  noninvariance biases the sensitization contrast.

- [Missing Posttests: Validating the Sensitivity
  Analysis](https://juhalt.github.io/solomonR/articles/mi-validation.md):

  Simulation study of fit_solomon_mi(): multiple imputation of missing
  posttests under missing at random and under known departures from it.

- [Longitudinal Designs: Validating the Repeated-Measures
  Analysis](https://juhalt.github.io/solomonR/articles/mmrm-validation.md):

  Simulation study of fit_solomon_mmrm(): the Solomon contrasts at each
  of three posttest occasions under dropout that is missing at random.

### References

- [References and the Solomon
  Literature](https://juhalt.github.io/solomonR/articles/references.md):

  Every source solomonR draws on, in APA Style (7th ed.), with a note on
  what each Solomon-design work contributes and where the package uses
  it.

- [Coverage of the Published
  Methodology](https://juhalt.github.io/solomonR/articles/coverage.md):

  Every published work on the Solomon design that solomonR draws on,
  what it contributes to the analysis of the design, and where the
  package implements it, or why it does not yet.

- [How to Cite solomonR and the Methods It
  Implements](https://juhalt.github.io/solomonR/articles/citing.md):

  Cite the package and the published methods behind each analysis.
