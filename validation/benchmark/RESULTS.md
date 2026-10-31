# Exploratory benchmark results

100 data sets per scenario, 300 records per data set, 10 imputations. The true analysis coefficient is 0.8. All results below are conditional on a successful fit; failures are retained in replicates.csv.

100 fits emitted recorded warnings; 0 fits failed. In missing_only, ordinary PMM logs removal of the constant all-No flag; this does not change its predictor information.

Each parenthesized value is a Monte Carlo standard error. Coverage also includes a 95% Wilson interval for its simulation probability. A zero plug-in MCSE at 0% or 100% is not certainty.

| Scenario | Method | Successful/requested | Bias (MCSE) | RMSE (MCSE) | Coverage (MCSE); Wilson interval | Mean interval width (MCSE) |
|---|---|---:|---:|---:|---:|---:|
| category_error | catmech | 100/100 | -0.151 (0.007) | 0.167 (0.007) | 0.37 (0.048); [0.282, 0.468] | 0.252 (0.002) |
| false_positive_flags | catmech | 100/100 | -0.311 (0.010) | 0.326 (0.009) | 0.43 (0.050); [0.337, 0.528] | 0.560 (0.017) |
| heavy_tails | catmech | 100/100 | 0.013 (0.007) | 0.072 (0.006) | 0.95 (0.022); [0.888, 0.978] | 0.267 (0.005) |
| missing_only | catmech | 100/100 | 0.001 (0.006) | 0.061 (0.006) | 0.93 (0.026); [0.863, 0.966] | 0.237 (0.002) |
| perfect_flags | catmech | 100/100 | -0.002 (0.006) | 0.056 (0.004) | 0.97 (0.017); [0.915, 0.990] | 0.240 (0.002) |
| category_error | naive_observed | 100/100 | -0.722 (0.003) | 0.723 (0.003) | 0.00 (0.000); [0.000, 0.037] | 0.102 (0.001) |
| false_positive_flags | naive_observed | 100/100 | -0.661 (0.004) | 0.662 (0.004) | 0.00 (0.000); [0.000, 0.037] | 0.133 (0.002) |
| heavy_tails | naive_observed | 100/100 | -0.722 (0.006) | 0.724 (0.005) | 0.00 (0.000); [0.000, 0.037] | 0.100 (0.001) |
| missing_only | naive_observed | 100/100 | 0.003 (0.007) | 0.065 (0.005) | 0.95 (0.022); [0.888, 0.978] | 0.257 (0.002) |
| perfect_flags | naive_observed | 100/100 | -0.719 (0.003) | 0.720 (0.003) | 0.00 (0.000); [0.000, 0.037] | 0.101 (0.001) |
| category_error | oracle | 100/100 | 0.004 (0.006) | 0.061 (0.004) | 0.96 (0.020); [0.902, 0.984] | 0.229 (0.001) |
| false_positive_flags | oracle | 100/100 | 0.002 (0.006) | 0.056 (0.004) | 0.96 (0.020); [0.902, 0.984] | 0.231 (0.001) |
| heavy_tails | oracle | 100/100 | 0.011 (0.006) | 0.059 (0.004) | 0.96 (0.020); [0.902, 0.984] | 0.238 (0.004) |
| missing_only | oracle | 100/100 | 0.003 (0.006) | 0.057 (0.004) | 0.95 (0.022); [0.888, 0.978] | 0.229 (0.001) |
| perfect_flags | oracle | 100/100 | 0.002 (0.005) | 0.054 (0.004) | 0.96 (0.020); [0.902, 0.984] | 0.229 (0.001) |
| category_error | ordinary_pmm | 100/100 | -0.707 (0.003) | 0.707 (0.003) | 0.00 (0.000); [0.000, 0.037] | 0.100 (0.001) |
| false_positive_flags | ordinary_pmm | 100/100 | -0.638 (0.005) | 0.640 (0.005) | 0.00 (0.000); [0.000, 0.037] | 0.131 (0.002) |
| heavy_tails | ordinary_pmm | 100/100 | -0.708 (0.005) | 0.711 (0.005) | 0.00 (0.000); [0.000, 0.037] | 0.098 (0.001) |
| missing_only | ordinary_pmm | 100/100 | 0.001 (0.006) | 0.061 (0.005) | 0.94 (0.024); [0.875, 0.972] | 0.239 (0.002) |
| perfect_flags | ordinary_pmm | 100/100 | -0.704 (0.003) | 0.705 (0.003) | 0.00 (0.000); [0.000, 0.037] | 0.098 (0.001) |
| category_error | trusted_observed | 100/100 | 0.001 (0.007) | 0.069 (0.005) | 0.95 (0.022); [0.888, 0.978] | 0.271 (0.002) |
| false_positive_flags | trusted_observed | 100/100 | 0.003 (0.007) | 0.065 (0.004) | 0.97 (0.017); [0.915, 0.990] | 0.274 (0.002) |
| heavy_tails | trusted_observed | 100/100 | 0.012 (0.007) | 0.073 (0.005) | 0.94 (0.024); [0.875, 0.972] | 0.284 (0.006) |
| missing_only | trusted_observed | 100/100 | 0.003 (0.007) | 0.065 (0.005) | 0.95 (0.022); [0.888, 0.978] | 0.257 (0.002) |
| perfect_flags | trusted_observed | 100/100 | 0.005 (0.007) | 0.067 (0.005) | 0.95 (0.022); [0.888, 0.978] | 0.272 (0.002) |
| category_error | trusted_pmm | 100/100 | 0.002 (0.006) | 0.063 (0.004) | 0.98 (0.014); [0.930, 0.994] | 0.254 (0.002) |
| false_positive_flags | trusted_pmm | 100/100 | 0.002 (0.006) | 0.060 (0.005) | 0.94 (0.024); [0.875, 0.972] | 0.247 (0.002) |
| heavy_tails | trusted_pmm | 100/100 | 0.023 (0.007) | 0.074 (0.005) | 0.96 (0.020); [0.902, 0.984] | 0.266 (0.005) |
| missing_only | trusted_pmm | 100/100 | 0.001 (0.006) | 0.061 (0.005) | 0.94 (0.024); [0.875, 0.972] | 0.239 (0.002) |
| perfect_flags | trusted_pmm | 100/100 | -0.001 (0.006) | 0.057 (0.004) | 0.95 (0.022); [0.888, 0.978] | 0.243 (0.002) |

## Paired comparison with the strongest matched PMM baseline

The baseline treats all flagged measurements as missing and supplies the same covariates, outcome and category to standard mice PMM. Negative squared-error differences favor catmech. Differences smaller than two MCSEs do not support a superiority claim.

| Scenario | Pairs | Mean squared-error difference | MCSE | Coverage difference | MCSE |
|---|---:|---:|---:|---:|---:|
| missing_only | 100 | 0.00009 | 0.00011 | -0.01 | 0.017 |
| perfect_flags | 100 | -0.00017 | 0.00011 | 0.02 | 0.014 |
| false_positive_flags | 100 | 0.10236 | 0.00608 | -0.51 | 0.052 |
| category_error | 100 | 0.02374 | 0.00226 | -0.61 | 0.053 |
| heavy_tails | 100 | -0.00026 | 0.00024 | -0.01 | 0.017 |

## Interpretation and boundary

The perfect-flags scenario uses q = 1. The uncertain-flags scenario uses q = 0.5, equal to the generative corruption fraction but not a calibrated conditional posterior. The highest category has no finite upper bound: a severely corrupted value may remain category-compatible and be retained by the sensitivity mixture. Exact category consistency alone therefore cannot validate inference.

The category-error scenario deliberately violates the exact-category assumption. Forcing incorrect categories can damage both bias and coverage. The heavy-tail scenario changes the latent residual distribution to a standardized Student t with three degrees of freedom.

These comparisons test a specific coefficient and five synthetic settings. They neither prove proper multiple imputation nor validate confidence intervals for new data or analysis models. Inspect the trusted PMM comparison before attributing a gain over naive observed-data analysis to a new method.

## Artifacts

replicates.csv records every fit, warning and failure; summary.csv includes all MCSEs, failure rates and Wilson intervals; paired_comparisons.csv includes every comparator; scenarios.csv records the data-generating settings; source-hashes.csv and session-info.txt identify the code and environment. The script rejects source edits during a run.
