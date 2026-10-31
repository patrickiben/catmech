# Exploratory simulation plan

Written before executing the revised package benchmark, 2026-10-31.

## Aim and estimand

Estimate the coefficient of the continuous true laboratory value in the linear
analysis `outcome ~ lab_value + age + sex + inflammation`. The data-generating
coefficient is 0.8. These simulations assess this one downstream estimand. They do
not establish recovery of individual latent values or broad clinical validity.
The units of the coefficient are outcome units per laboratory unit.

## Data and methods

Use 300 independent records per replicate, 100 replicates per scenario and 10
imputations by default. Seed each data set separately and use the same generated
data across methods. Include outcome in every imputation predictor set. Exclude
`lab_true` and `corrupted` from every imputation input; only the oracle analysis
can use the hidden truth. The thresholds are `(-Inf, 5.5, 7, Inf)`, chosen to avoid
the original demonstration's nearly empty lowest category.

| Scenario | Missing rate | Flag rate | Fraction corrupted among flags | Category error | Error distribution |
|---|---:|---:|---:|---:|---|
| missing_only | 0.20 | 0 | 1 | 0 | normal |
| perfect_flags | 0.20 | 0.10 | 1 | 0 | normal |
| false_positive_flags | 0.20 | 0.10 | 0.50 | 0 | normal |
| category_error | 0.20 | 0.10 | 1 | 0.10 | normal |
| heavy_tails | 0.20 | 0.10 | 1 | 0 | heavy_tailed |

Methods are catmech; ordinary mice PMM retaining observed flagged values; mice PMM
with every flagged value treated as missing; observed-only regression retaining
flagged values; observed-only regression excluding flagged values; and oracle
regression on complete true values. The two PMM baselines receive the same
covariates, outcome and category information. Ordinary PMM additionally receives
the quality flag. The trusted PMM baseline directly tests whether any benefit
comes from interval constraints beyond simply excluding suspect observations.
Both use standard mice type-1 matching and five donors.

For catmech, use `q = 1` except in `false_positive_flags`, where `q = 0.5`.
This is a sensitivity replacement probability, not a fitted posterior probability
of corruption. Even when its numeric value equals the generative fraction of
corrupt flags, conditioning on a compatible observed value and exact category can
change the appropriate corruption probability. This scenario is a stress test,
not a claim that the sensitivity mixture reproduces the data-generating posterior.
Use `conflict_policy = "reconstruct"` explicitly. The category-error scenario
violates the package's exact-category assumption and probes consequences of that
violation; it is not within the supported inferential scope. Do not silently
change the specification to make a scenario perform better.

Only the continuous target is incomplete in this benchmark. One mice iteration
therefore draws the target from fully observed predictors without an iterative
chain convergence problem. The package may internally use more iterations; report
its actual settings in the saved environment and source. Ten imputations remain a
small exploratory choice, not a recommendation for final analysis.

## Performance measures and uncertainty

For each method and scenario report the requested and successful replicate counts,
failure proportion, bias, RMSE, mean estimated standard error, mean interval
width and 95% interval coverage. Log all warnings and failures by replicate, including unsupported donor
configurations. Performance measures based on successful runs are explicitly
conditional on success; a failing method must not appear competitive by silently
dropping failed runs.

Monte Carlo standard errors (MCSEs), with R successful replicates:

- Bias: `sd(estimate - truth) / sqrt(R)`.
- RMSE: `sd((estimate - truth)^2) / (2 * RMSE * sqrt(R))` (delta method).
- Mean standard error: `sd(standard_error) / sqrt(R)`.
- Mean interval width: `sd(upper - lower) / sqrt(R)`.
- Coverage: `sqrt(coverage * (1 - coverage) / R)`.
- Failure proportion: the same binomial formula, using all requested replicates.

Pair each catmech result with each comparator on the same replicate. Report mean
paired differences in error, squared error and coverage with `sd(difference) /
sqrt(R_paired)` MCSEs and the number of complete pairs. Negative squared-error
differences favor catmech. Do not use a sub-two-MCSE gap as evidence of superiority.

At 95% coverage, 100 replicates give about 0.022 coverage MCSE. This can detect gross
problems but cannot establish nominal coverage to within a percentage point. About
1,900 independent replicates would be needed for a target MCSE of 0.005 at 95%
coverage. The default run is deliberately labeled exploratory.

## Success, failure and transfer boundary

Software success means reproducible execution, finite completed target values,
zero category violations in supported exact-category runs, intact original inputs,
and explicit accounting of failures. No superiority claim is prespecified.
Numerical feasibility or zero interval violations alone does not validate Rubin
pooling. Coverage substantially below 0.95, beyond its Monte Carlo uncertainty,
weakens an inferential claim for that scenario. Good results here apply only to
these data-generating processes and sample size; they do not establish robustness
to dependent observations, multiple incomplete predictors, informative selection
of trusted donors, category misclassification, or arbitrary analysis models.

## Reproduction

From the repository root, after installing development dependencies:

```sh
Rscript tools/benchmark.R 100 10
```

Optional third argument is a new or empty output directory. Existing results are
never overwritten. Results contain replicate-level
estimates and failures, aggregate measures, paired comparisons, warnings, session
information and source hashes captured before execution. The script rejects any
source or plan edits made during execution. Use a fresh directory for a confirmatory run after
changing code. The archived v0.9 coefficient from one seed is not comparative
evidence and is not carried forward as a result of this benchmark.

## Reporting additions after the first exploratory run

The first 100-replicate run completed before these reporting additions. Mean
interval width with its MCSE and 95% Wilson intervals for coverage and failure
probabilities were added to improve interpretation, especially when a plug-in
binomial MCSE is zero. Source hashes were moved to the start of execution, with
a consistency check at completion. A machine-generated results report was added.
No scenario, estimand, seed schedule, method or performance tuning changed as a
result of the first run. Any further implementation fixes require rerunning the
same plan; they must not be selected to produce a favorable result.
