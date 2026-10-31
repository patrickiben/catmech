# v0.3 extension evaluation: reference

Prespecified settings, assumptions and interpretation boundaries are in ../EXTENSION_PLAN.md. Coverage is conditional on returned intervals; request-status.csv retains failures. Same-reference results across recording scenarios reuse cohorts and reference subsets, so they are not independent replications.

| Setting | Method | Intervals/requested | Bias (MCSE) | Coverage (MCSE) [95% Wilson] | Width (MCSE) |
|---|---|---:|---:|---:|---:|
| reference/clean | gold_only | 1000/1000 | 0.005 (0.003) | 0.953 (0.007) [0.938, 0.964] | 0.330 (0.001) |
| reference/flagged_errors | gold_only | 1000/1000 | 0.005 (0.003) | 0.953 (0.007) [0.938, 0.964] | 0.330 (0.001) |
| reference/missed_flags | gold_only | 1000/1000 | 0.005 (0.003) | 0.953 (0.007) [0.938, 0.964] | 0.330 (0.001) |
| reference/rounded_records | gold_only | 1000/1000 | 0.005 (0.003) | 0.953 (0.007) [0.938, 0.964] | 0.330 (0.001) |
| reference/site_shift | gold_only | 1000/1000 | 0.005 (0.003) | 0.953 (0.007) [0.938, 0.964] | 0.330 (0.001) |
| reference/clean | mice_reference_normal | 1000/1000 | -0.003 (0.002) | 0.947 (0.007) [0.931, 0.959] | 0.279 (0.001) |
| reference/flagged_errors | mice_reference_normal | 1000/1000 | -0.003 (0.002) | 0.947 (0.007) [0.931, 0.959] | 0.279 (0.001) |
| reference/missed_flags | mice_reference_normal | 1000/1000 | -0.003 (0.002) | 0.947 (0.007) [0.931, 0.959] | 0.279 (0.001) |
| reference/rounded_records | mice_reference_normal | 1000/1000 | -0.003 (0.002) | 0.947 (0.007) [0.931, 0.959] | 0.279 (0.001) |
| reference/site_shift | mice_reference_normal | 1000/1000 | -0.003 (0.002) | 0.947 (0.007) [0.931, 0.959] | 0.279 (0.001) |
| reference/clean | mice_reference_pmm | 1000/1000 | 0.007 (0.002) | 0.851 (0.011) [0.828, 0.872] | 0.212 (0.001) |
| reference/flagged_errors | mice_reference_pmm | 1000/1000 | 0.007 (0.002) | 0.851 (0.011) [0.828, 0.872] | 0.212 (0.001) |
| reference/missed_flags | mice_reference_pmm | 1000/1000 | 0.007 (0.002) | 0.851 (0.011) [0.828, 0.872] | 0.212 (0.001) |
| reference/rounded_records | mice_reference_pmm | 1000/1000 | 0.007 (0.002) | 0.851 (0.011) [0.828, 0.872] | 0.212 (0.001) |
| reference/site_shift | mice_reference_pmm | 1000/1000 | 0.007 (0.002) | 0.851 (0.011) [0.828, 0.872] | 0.212 (0.001) |
| reference/clean | oracle | 1000/1000 | 0.002 (0.001) | 0.952 (0.007) [0.937, 0.964] | 0.161 (0.000) |
| reference/flagged_errors | oracle | 1000/1000 | 0.002 (0.001) | 0.952 (0.007) [0.937, 0.964] | 0.161 (0.000) |
| reference/missed_flags | oracle | 1000/1000 | 0.002 (0.001) | 0.952 (0.007) [0.937, 0.964] | 0.161 (0.000) |
| reference/rounded_records | oracle | 1000/1000 | 0.002 (0.001) | 0.952 (0.007) [0.937, 0.964] | 0.161 (0.000) |
| reference/site_shift | oracle | 1000/1000 | 0.002 (0.001) | 0.952 (0.007) [0.937, 0.964] | 0.161 (0.000) |
| reference/clean | reference_normal | 1000/1000 | -0.003 (0.002) | 0.952 (0.007) [0.937, 0.964] | 0.278 (0.001) |
| reference/flagged_errors | reference_normal | 1000/1000 | -0.003 (0.002) | 0.952 (0.007) [0.937, 0.964] | 0.278 (0.001) |
| reference/missed_flags | reference_normal | 1000/1000 | -0.003 (0.002) | 0.952 (0.007) [0.937, 0.964] | 0.278 (0.001) |
| reference/rounded_records | reference_normal | 1000/1000 | -0.003 (0.002) | 0.952 (0.007) [0.937, 0.964] | 0.278 (0.001) |
| reference/site_shift | reference_normal | 1000/1000 | -0.003 (0.002) | 0.952 (0.007) [0.937, 0.964] | 0.278 (0.001) |
| reference/clean | reference_pmm | 1000/1000 | 0.007 (0.002) | 0.851 (0.011) [0.828, 0.872] | 0.212 (0.001) |
| reference/flagged_errors | reference_pmm | 1000/1000 | 0.007 (0.002) | 0.851 (0.011) [0.828, 0.872] | 0.212 (0.001) |
| reference/missed_flags | reference_pmm | 1000/1000 | 0.007 (0.002) | 0.851 (0.011) [0.828, 0.872] | 0.212 (0.001) |
| reference/rounded_records | reference_pmm | 1000/1000 | 0.007 (0.002) | 0.851 (0.011) [0.828, 0.872] | 0.212 (0.001) |
| reference/site_shift | reference_pmm | 1000/1000 | 0.007 (0.002) | 0.851 (0.011) [0.828, 0.872] | 0.212 (0.001) |
| reference/clean | v02_mask_normal | 1000/1000 | 0.001 (0.001) | 0.964 (0.006) [0.951, 0.974] | 0.172 (0.000) |
| reference/flagged_errors | v02_mask_normal | 1000/1000 | 0.001 (0.001) | 0.964 (0.006) [0.951, 0.974] | 0.180 (0.000) |
| reference/missed_flags | v02_mask_normal | 1000/1000 | -0.077 (0.001) | 0.581 (0.016) [0.550, 0.611] | 0.175 (0.000) |
| reference/rounded_records | v02_mask_normal | 1000/1000 | 0.000 (0.001) | 0.965 (0.006) [0.952, 0.975] | 0.172 (0.000) |
| reference/site_shift | v02_mask_normal | 1000/1000 | -0.139 (0.001) | 0.127 (0.011) [0.108, 0.149] | 0.173 (0.000) |

Normal/PMM comparisons must hold category information constant. Reference workflows obtain additional measured truth; they do not identify missed errors without that information. Reference-only and independent mice comparators use the same gold measurements. Nonignorable selection, noisy gold and sparse donor support are deliberate failure probes.

Summary, paired contrasts, bias decomposition and all requests are in adjacent CSVs. Bias-centered coverage is retrospective and uses hidden truth; it is not a real-data correction. Source hashes match before and after execution. No result is made valid merely by this script completing.
