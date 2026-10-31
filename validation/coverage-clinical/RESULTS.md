# Coverage evaluation: clinical

The table reports conditional coverage only among returned intervals. Refusal, nonconvergence and other errors remain in the requested denominator. Bias-centered coverage is a retrospective simulation diagnostic using hidden truth, not a corrected confidence interval.

| Scenario | Method | Intervals/requested | Refused/MCMC/error | Bias (MCSE) | RMSE | Coverage [95% Wilson] | Mean width | Bias-centered diagnostic |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| clinical/clean | model_ignore | 250/250 | 0/0/0 | 0.006 (0.003) | 0.045 | 0.948 [0.913, 0.969] | 0.173 | 0.960 |
| clinical/flagged_errors | model_ignore | 250/250 | 0/0/0 | 0.005 (0.003) | 0.045 | 0.952 [0.918, 0.972] | 0.180 | 0.944 |
| clinical/missed_flags | model_ignore | 250/250 | 0/0/0 | -0.073 (0.003) | 0.088 | 0.644 [0.583, 0.701] | 0.175 | 0.916 |
| clinical/rounded_records | model_ignore | 250/250 | 0/0/0 | 0.005 (0.003) | 0.044 | 0.948 [0.913, 0.969] | 0.172 | 0.960 |
| clinical/site_shift | model_ignore | 250/250 | 0/0/0 | -0.134 (0.003) | 0.143 | 0.148 [0.109, 0.197] | 0.173 | 0.916 |
| clinical/clean | naive_observed | 250/250 | 0/0/0 | 0.006 (0.003) | 0.047 | 0.944 [0.908, 0.966] | 0.181 | 0.948 |
| clinical/flagged_errors | naive_observed | 250/250 | 0/0/0 | -0.147 (0.003) | 0.156 | 0.096 [0.065, 0.139] | 0.172 | 0.932 |
| clinical/missed_flags | naive_observed | 250/250 | 0/0/0 | -0.147 (0.003) | 0.156 | 0.096 [0.065, 0.139] | 0.172 | 0.932 |
| clinical/rounded_records | naive_observed | 250/250 | 0/0/0 | 0.006 (0.003) | 0.047 | 0.948 [0.913, 0.969] | 0.181 | 0.948 |
| clinical/site_shift | naive_observed | 250/250 | 0/0/0 | -0.147 (0.003) | 0.156 | 0.096 [0.065, 0.139] | 0.172 | 0.932 |
| clinical/clean | oracle | 250/250 | 0/0/0 | 0.005 (0.003) | 0.041 | 0.960 [0.928, 0.978] | 0.161 | 0.968 |
| clinical/flagged_errors | oracle | 250/250 | 0/0/0 | 0.005 (0.003) | 0.041 | 0.960 [0.928, 0.978] | 0.161 | 0.968 |
| clinical/missed_flags | oracle | 250/250 | 0/0/0 | 0.005 (0.003) | 0.041 | 0.960 [0.928, 0.978] | 0.161 | 0.968 |
| clinical/rounded_records | oracle | 250/250 | 0/0/0 | 0.005 (0.003) | 0.041 | 0.960 [0.928, 0.978] | 0.161 | 0.968 |
| clinical/site_shift | oracle | 250/250 | 0/0/0 | 0.005 (0.003) | 0.041 | 0.960 [0.928, 0.978] | 0.161 | 0.968 |
| clinical/clean | trusted_normal | 250/250 | 0/0/0 | 0.006 (0.003) | 0.045 | 0.948 [0.913, 0.969] | 0.173 | 0.960 |
| clinical/flagged_errors | trusted_normal | 250/250 | 0/0/0 | 0.005 (0.003) | 0.045 | 0.952 [0.918, 0.972] | 0.180 | 0.944 |
| clinical/missed_flags | trusted_normal | 250/250 | 0/0/0 | -0.073 (0.003) | 0.088 | 0.636 [0.575, 0.693] | 0.176 | 0.916 |
| clinical/rounded_records | trusted_normal | 250/250 | 0/0/0 | 0.005 (0.003) | 0.044 | 0.948 [0.913, 0.969] | 0.172 | 0.960 |
| clinical/site_shift | trusted_normal | 250/250 | 0/0/0 | -0.134 (0.003) | 0.143 | 0.136 [0.099, 0.184] | 0.174 | 0.920 |
| clinical/clean | trusted_observed | 250/250 | 0/0/0 | 0.006 (0.003) | 0.047 | 0.944 [0.908, 0.966] | 0.181 | 0.948 |
| clinical/flagged_errors | trusted_observed | 250/250 | 0/0/0 | 0.006 (0.003) | 0.049 | 0.948 [0.913, 0.969] | 0.194 | 0.952 |
| clinical/missed_flags | trusted_observed | 250/250 | 0/0/0 | -0.072 (0.003) | 0.088 | 0.648 [0.587, 0.705] | 0.186 | 0.932 |
| clinical/rounded_records | trusted_observed | 250/250 | 0/0/0 | 0.006 (0.003) | 0.047 | 0.948 [0.913, 0.969] | 0.181 | 0.948 |
| clinical/site_shift | trusted_observed | 250/250 | 0/0/0 | -0.133 (0.003) | 0.143 | 0.224 [0.177, 0.280] | 0.184 | 0.912 |
| clinical/clean | trusted_pmm | 250/250 | 0/0/0 | 0.008 (0.003) | 0.045 | 0.944 [0.908, 0.966] | 0.169 | 0.948 |
| clinical/flagged_errors | trusted_pmm | 250/250 | 0/0/0 | 0.011 (0.003) | 0.046 | 0.936 [0.899, 0.960] | 0.173 | 0.948 |
| clinical/missed_flags | trusted_pmm | 250/250 | 0/0/0 | -0.068 (0.003) | 0.083 | 0.672 [0.612, 0.727] | 0.173 | 0.932 |
| clinical/rounded_records | trusted_pmm | 250/250 | 0/0/0 | 0.008 (0.003) | 0.045 | 0.944 [0.908, 0.966] | 0.169 | 0.944 |
| clinical/site_shift | trusted_pmm | 250/250 | 0/0/0 | -0.130 (0.003) | 0.138 | 0.184 [0.141, 0.237] | 0.171 | 0.920 |

All MCSEs, interval yields, status rates and paired comparisons are in the adjacent CSV files. A model is not automatically validated by completing this script. Gibbs results use a smaller prespecified subset; primary and stress settings have different interpretation boundaries. Consult COVERAGE_PLAN.md before making a claim.

The direct ignore-category model is an implementation of established normal multiple imputation after masking quality flags. Exact/noisy category models additionally assume the declared category likelihood. Missed corruption, inaccurate category likelihoods, model misspecification and incompatible analysis can defeat inference.

source-hashes.csv verifies the files read before execution against the end of this run. settings.json records the fixed computational settings. The checkpoint was removed only after successful source-integrity verification and complete result generation.
