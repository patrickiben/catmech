# Coverage evaluation: confirmatory

The table reports conditional coverage only among returned intervals. Refusal, nonconvergence and other errors remain in the requested denominator. Bias-centered coverage is a retrospective simulation diagnostic using hidden truth, not a corrected confidence interval.

| Scenario | Method | Intervals/requested | Refused/MCMC/error | Bias (MCSE) | RMSE | Coverage [95% Wilson] | Mean width | Bias-centered diagnostic |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| synthetic/category_error | legacy | 1000/1000 | 0/0/0 | -0.154 (0.002) | 0.170 | 0.346 [0.317, 0.376] | 0.251 | 0.928 |
| synthetic/false_positive_flags | legacy | 1000/1000 | 0/0/0 | -0.296 (0.003) | 0.315 | 0.397 [0.367, 0.428] | 0.513 | 0.945 |
| synthetic/heavy_tails | legacy | 1000/1000 | 0/0/0 | 0.006 (0.002) | 0.073 | 0.952 [0.937, 0.964] | 0.275 | 0.944 |
| synthetic/missing_only | legacy | 1000/1000 | 0/0/0 | -0.003 (0.002) | 0.061 | 0.950 [0.935, 0.962] | 0.237 | 0.947 |
| synthetic/null_effect | legacy | 1000/1000 | 0/0/0 | 0.002 (0.001) | 0.042 | 0.983 [0.973, 0.989] | 0.200 | 0.983 |
| synthetic/perfect_flags | legacy | 1000/1000 | 0/0/0 | -0.002 (0.002) | 0.062 | 0.953 [0.938, 0.964] | 0.240 | 0.954 |
| synthetic/category_error | model_exact | 0/1000 | 1000/0/0 | NA (NA) | NA | NA [NA, NA] | NA | NA |
| synthetic/false_positive_flags | model_exact | 250/250 | 0/0/0 | 0.007 (0.004) | 0.065 | 0.932 [0.894, 0.957] | 0.241 | 0.932 |
| synthetic/heavy_tails | model_exact | 249/250 | 0/1/0 | 0.015 (0.005) | 0.083 | 0.932 [0.893, 0.957] | 0.266 | 0.924 |
| synthetic/missing_only | model_exact | 250/250 | 0/0/0 | -0.001 (0.004) | 0.059 | 0.976 [0.949, 0.989] | 0.238 | 0.976 |
| synthetic/null_effect | model_exact | 250/250 | 0/0/0 | -0.001 (0.004) | 0.058 | 0.960 [0.928, 0.978] | 0.239 | 0.960 |
| synthetic/perfect_flags | model_exact | 250/250 | 0/0/0 | 0.002 (0.004) | 0.055 | 0.984 [0.960, 0.994] | 0.240 | 0.988 |
| synthetic/category_error | model_ignore | 1000/1000 | 0/0/0 | -0.002 (0.002) | 0.063 | 0.952 [0.937, 0.964] | 0.255 | 0.954 |
| synthetic/false_positive_flags | model_ignore | 1000/1000 | 0/0/0 | -0.003 (0.002) | 0.064 | 0.952 [0.937, 0.964] | 0.255 | 0.949 |
| synthetic/heavy_tails | model_ignore | 1000/1000 | 0/0/0 | 0.009 (0.003) | 0.086 | 0.925 [0.907, 0.940] | 0.274 | 0.926 |
| synthetic/missing_only | model_ignore | 1000/1000 | 0/0/0 | -0.002 (0.002) | 0.062 | 0.950 [0.935, 0.962] | 0.245 | 0.950 |
| synthetic/null_effect | model_ignore | 1000/1000 | 0/0/0 | 0.003 (0.002) | 0.068 | 0.944 [0.928, 0.957] | 0.274 | 0.946 |
| synthetic/perfect_flags | model_ignore | 1000/1000 | 0/0/0 | -0.000 (0.002) | 0.064 | 0.963 [0.949, 0.973] | 0.254 | 0.963 |
| synthetic/category_error | model_noisy | 250/250 | 0/0/0 | -0.001 (0.004) | 0.065 | 0.932 [0.894, 0.957] | 0.247 | 0.936 |
| synthetic/category_error | model_wrong_noise | 250/250 | 0/0/0 | -0.020 (0.004) | 0.069 | 0.920 [0.880, 0.948] | 0.248 | 0.920 |
| synthetic/category_error | naive_observed | 1000/1000 | 0/0/0 | -0.721 (0.001) | 0.721 | 0.000 [0.000, 0.004] | 0.101 | 0.916 |
| synthetic/false_positive_flags | naive_observed | 1000/1000 | 0/0/0 | -0.655 (0.002) | 0.657 | 0.000 [0.000, 0.004] | 0.133 | 0.830 |
| synthetic/heavy_tails | naive_observed | 1000/1000 | 0/0/0 | -0.723 (0.001) | 0.725 | 0.000 [0.000, 0.004] | 0.101 | 0.859 |
| synthetic/missing_only | naive_observed | 1000/1000 | 0/0/0 | -0.000 (0.002) | 0.065 | 0.943 [0.927, 0.956] | 0.256 | 0.943 |
| synthetic/null_effect | naive_observed | 1000/1000 | 0/0/0 | 0.001 (0.001) | 0.027 | 0.949 [0.934, 0.961] | 0.108 | 0.950 |
| synthetic/perfect_flags | naive_observed | 1000/1000 | 0/0/0 | -0.721 (0.001) | 0.721 | 0.000 [0.000, 0.004] | 0.102 | 0.915 |
| synthetic/category_error | oracle | 1000/1000 | 0/0/0 | -0.000 (0.002) | 0.057 | 0.945 [0.929, 0.958] | 0.229 | 0.945 |
| synthetic/false_positive_flags | oracle | 1000/1000 | 0/0/0 | -0.001 (0.002) | 0.057 | 0.955 [0.940, 0.966] | 0.229 | 0.956 |
| synthetic/heavy_tails | oracle | 1000/1000 | 0/0/0 | 0.004 (0.002) | 0.061 | 0.944 [0.928, 0.957] | 0.241 | 0.938 |
| synthetic/missing_only | oracle | 1000/1000 | 0/0/0 | -0.001 (0.002) | 0.058 | 0.959 [0.945, 0.970] | 0.229 | 0.958 |
| synthetic/null_effect | oracle | 1000/1000 | 0/0/0 | 0.002 (0.002) | 0.058 | 0.955 [0.940, 0.966] | 0.229 | 0.952 |
| synthetic/perfect_flags | oracle | 1000/1000 | 0/0/0 | 0.001 (0.002) | 0.058 | 0.952 [0.937, 0.964] | 0.229 | 0.952 |
| synthetic/category_error | trusted_normal | 1000/1000 | 0/0/0 | -0.002 (0.002) | 0.063 | 0.952 [0.937, 0.964] | 0.255 | 0.952 |
| synthetic/false_positive_flags | trusted_normal | 1000/1000 | 0/0/0 | -0.003 (0.002) | 0.064 | 0.953 [0.938, 0.964] | 0.255 | 0.952 |
| synthetic/heavy_tails | trusted_normal | 1000/1000 | 0/0/0 | 0.009 (0.003) | 0.086 | 0.921 [0.903, 0.936] | 0.274 | 0.923 |
| synthetic/missing_only | trusted_normal | 1000/1000 | 0/0/0 | -0.002 (0.002) | 0.062 | 0.951 [0.936, 0.963] | 0.245 | 0.951 |
| synthetic/null_effect | trusted_normal | 1000/1000 | 0/0/0 | 0.003 (0.002) | 0.068 | 0.941 [0.925, 0.954] | 0.273 | 0.940 |
| synthetic/perfect_flags | trusted_normal | 1000/1000 | 0/0/0 | -0.000 (0.002) | 0.064 | 0.959 [0.945, 0.970] | 0.254 | 0.959 |
| synthetic/category_error | trusted_observed | 1000/1000 | 0/0/0 | -0.002 (0.002) | 0.067 | 0.957 [0.943, 0.968] | 0.271 | 0.956 |
| synthetic/false_positive_flags | trusted_observed | 1000/1000 | 0/0/0 | -0.001 (0.002) | 0.068 | 0.950 [0.935, 0.962] | 0.272 | 0.951 |
| synthetic/heavy_tails | trusted_observed | 1000/1000 | 0/0/0 | 0.004 (0.002) | 0.073 | 0.943 [0.927, 0.956] | 0.290 | 0.942 |
| synthetic/missing_only | trusted_observed | 1000/1000 | 0/0/0 | -0.000 (0.002) | 0.065 | 0.943 [0.927, 0.956] | 0.256 | 0.943 |
| synthetic/null_effect | trusted_observed | 1000/1000 | 0/0/0 | 0.002 (0.002) | 0.068 | 0.946 [0.930, 0.958] | 0.272 | 0.944 |
| synthetic/perfect_flags | trusted_observed | 1000/1000 | 0/0/0 | 0.000 (0.002) | 0.067 | 0.958 [0.944, 0.969] | 0.271 | 0.958 |
| synthetic/category_error | trusted_pmm | 1000/1000 | 0/0/0 | -0.004 (0.002) | 0.063 | 0.952 [0.937, 0.964] | 0.253 | 0.948 |
| synthetic/false_positive_flags | trusted_pmm | 1000/1000 | 0/0/0 | -0.002 (0.002) | 0.063 | 0.947 [0.931, 0.959] | 0.244 | 0.948 |
| synthetic/heavy_tails | trusted_pmm | 1000/1000 | 0/0/0 | 0.015 (0.002) | 0.074 | 0.947 [0.931, 0.959] | 0.274 | 0.942 |
| synthetic/missing_only | trusted_pmm | 1000/1000 | 0/0/0 | -0.002 (0.002) | 0.062 | 0.951 [0.936, 0.963] | 0.239 | 0.946 |
| synthetic/null_effect | trusted_pmm | 1000/1000 | 0/0/0 | 0.003 (0.002) | 0.060 | 0.953 [0.938, 0.964] | 0.240 | 0.951 |
| synthetic/perfect_flags | trusted_pmm | 1000/1000 | 0/0/0 | -0.000 (0.002) | 0.063 | 0.954 [0.939, 0.965] | 0.243 | 0.952 |

All MCSEs, interval yields, status rates and paired comparisons are in the adjacent CSV files. A model is not automatically validated by completing this script. Gibbs results use a smaller prespecified subset; primary and stress settings have different interpretation boundaries. Consult COVERAGE_PLAN.md before making a claim.

The direct ignore-category model is an implementation of established normal multiple imputation after masking quality flags. Exact/noisy category models additionally assume the declared category likelihood. Missed corruption, inaccurate category likelihoods, model misspecification and incompatible analysis can defeat inference.

source-hashes.csv verifies the files read before execution against the end of this run. settings.json records the fixed computational settings. The checkpoint was removed only after successful source-integrity verification and complete result generation.
