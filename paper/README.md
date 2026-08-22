# Technical report plan

Working title: **Fast Weighted Dependence Measures and Independence Tests**.

The report should center on three contributions that are not adequately
recorded elsewhere:

1. the finite-sample representation and fast algorithm for weighted
   Hoeffding's dependence measure;
2. the asymptotic framework used for weighted independence tests; and
3. the weighted extension of Chatterjee's correlation, including its exact
   conditional path-statistic mean and first-order null variance.

Pearson correlation, Spearman's rho, Kendall's tau, and Blomqvist's beta
provide context and a common notation, but should receive only a concise
treatment. The report is about the statistical and algorithmic ideas, not a
second user manual for the C++ or R packages.

## Proposed outline

### 1. Introduction

- Motivate observation weights in dependence measurement and testing.
- Explain why naively replacing counts by weights is insufficient for
  higher-order statistics and null inference.
- State the three main contributions precisely.
- Mention the C++ library and R interface briefly as reproducible
  implementations, not as separate contributions.

### 2. Weighted empirical framework

- Define the observations, nonnegative case weights, normalized masses, and
  Kish effective sample size.
- Separate frequency-weight interpretations from deterministic or
  covariate-dependent weights.
- Define weighted marginal and joint ranks, including conventions for ties.
- State the asymptotic regime explicitly, especially the required diffuseness
  of the normalized weights and whether weights may depend on the predictor.

### 3. Standard weighted dependence measures

- Give compact definitions of weighted Pearson correlation, Spearman's rho,
  Kendall's tau-b, and Blomqvist's beta.
- Show reduction to the ordinary sample measures under uniform weights.
- Record invariance to positive rescaling, handling of zero weights, and tie
  conventions.
- Summarize computation and test transformations in one table; move lengthy
  tie-correction formulas to an appendix.

### 4. Weighted Hoeffding's D and Chatterjee's xi

- Begin with the population target and the unweighted sample statistic.
- Derive the weighted distinct-index sums rather than presenting a heuristic
  plug-in estimator.
- Express the statistic through marginal ranks, bivariate ranks, powers of the
  weights, and elementary symmetric weight sums.
- Present the bivariate-rank algorithm based on sorting and weighted inversion
  counting.
- Prove equivalence to the direct combinatorial definition and establish
  time and memory complexity.
- Discuss ties, zero weights, numerical behavior, and the unweighted special
  case.
- Define the directional weighted Chatterjee estimator after sorting by the
  predictor.
- Motivate the base-point edge weights and weighted-rank denominator.
- Explain response-independent predictor-tie breaking, scale invariance,
  zero-mass equivalence, directionality, and response ties.

### 5. Asymptotic independence tests

- State the null hypothesis, sampling assumptions, and admissible weight
  sequences.
- Develop the limiting argument in a common weighted empirical-process or
  projection framework.
- Explain when Kish effective sample size supplies the correct scaling and
  where additional weight functionals or tie corrections enter.
- Derive the normal limits and transformations for Pearson, Spearman,
  Kendall, and Blomqvist.
- Treat Hoeffding's degenerate null separately, including the connection to
  the Blum--Kiefer--Rosenblatt limit and the approximation used for p-values.
- Derive the weighted Chatterjee path statistic's conditional null mean,
  projection variance, and central limit theorem for a continuous response.
- Distinguish proved asymptotic results from finite-sample approximations used
  by the software.

### 6. Numerical experiments

- Compare each optimized estimator with a small direct implementation.
- Benchmark the weighted Hoeffding algorithm against the direct combinatorial
  calculation over increasing sample sizes.
- Check null size across sample size, weight concentration, ties, and several
  weight-generating mechanisms.
- Validate the weighted Chatterjee null variance and compare power under
  directional and nonmonotone alternatives.
- Keep the design targeted at the theoretical claims rather than producing a
  broad software benchmark.

### Appendices

- Combinatorial derivation for weighted Hoeffding's measure.
- Proofs for the weighted independence-test asymptotics.
- Conditional-null calculations for weighted Chatterjee's correlation.
- Standard-measure formulas, tie corrections, and additional simulations.

## Decisions before drafting

- Fix the precise interpretation and admissible stochastic dependence of the
  weights in each theorem.
- Reconstruct every mathematical claim from the derivation; implementation
  behavior alone is not evidence for a theorem.
- Perform a focused literature review before making novelty claims, especially
  for weighted Hoeffding statistics and weighted versions of Chatterjee's
  correlation.
- Decide authorship and attribution for the Chatterjee work before filling the
  title page.
- Choose a target venue or technical-report series only after the main results
  and proof lengths are clear.

## Building

The numerical results require R package `wdm` version 0.3.0 or newer.
From the repository root, reproduce the committed tables and run metadata
with:

```sh
Rscript paper/simulation/run-study.R
```

The default study is single-threaded and is designed to finish in less than
five minutes.  For a quick pipeline check, reduce the replications with the
`WDM_NULL_REPS` and `WDM_POWER_REPS` environment variables.

From the repository root, run:

```sh
latexmk -cd -pdf paper/main.tex
```
