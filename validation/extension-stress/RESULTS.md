# v0.3 extension evaluation: stress

Prespecified settings, assumptions and interpretation boundaries are in ../EXTENSION_PLAN.md. Coverage is conditional on returned intervals; request-status.csv retains failures. Same-reference results across recording scenarios reuse cohorts and reference subsets, so they are not independent replications.

| Setting | Method | Intervals/requested | Bias (MCSE) | Coverage (MCSE) [95% Wilson] | Width (MCSE) |
|---|---|---:|---:|---:|---:|
| support/missed_flags_gold10 | gold_only | 250/250 | 0.006 (0.009) | 0.956 (0.013) [0.923, 0.975] | 0.567 (0.006) |
| support/missed_flags_gold50 | gold_only | 250/250 | 0.001 (0.004) | 0.960 (0.012) [0.928, 0.978] | 0.231 (0.001) |
| stress/noisy_reference | gold_only | 500/500 | -0.104 (0.004) | 0.772 (0.019) [0.733, 0.807] | 0.321 (0.001) |
| stress/nonignorable_selection | gold_only | 1000/1000 | 0.268 (0.003) | 0.136 (0.011) [0.116, 0.159] | 0.339 (0.001) |
| support/site_shift_gold10 | gold_only | 250/250 | 0.006 (0.009) | 0.956 (0.013) [0.923, 0.975] | 0.567 (0.006) |
| support/site_shift_gold50 | gold_only | 250/250 | 0.001 (0.004) | 0.960 (0.012) [0.928, 0.978] | 0.231 (0.001) |
| support/missed_flags_gold10 | mice_reference_normal | 250/250 | -0.020 (0.007) | 0.960 (0.012) [0.928, 0.978] | 0.467 (0.008) |
| support/missed_flags_gold50 | mice_reference_normal | 250/250 | -0.002 (0.003) | 0.960 (0.012) [0.928, 0.978] | 0.205 (0.001) |
| stress/noisy_reference | mice_reference_normal | 500/500 | -0.112 (0.003) | 0.674 (0.021) [0.632, 0.714] | 0.276 (0.002) |
| stress/nonignorable_selection | mice_reference_normal | 1000/1000 | 0.264 (0.002) | 0.056 (0.007) [0.043, 0.072] | 0.281 (0.001) |
| support/site_shift_gold10 | mice_reference_normal | 250/250 | -0.020 (0.007) | 0.960 (0.012) [0.928, 0.978] | 0.467 (0.008) |
| support/site_shift_gold50 | mice_reference_normal | 250/250 | -0.002 (0.003) | 0.960 (0.012) [0.928, 0.978] | 0.205 (0.001) |
| support/missed_flags_gold10 | mice_reference_pmm | 250/250 | 0.003 (0.008) | 0.748 (0.027) [0.691, 0.798] | 0.300 (0.007) |
| support/missed_flags_gold50 | mice_reference_pmm | 250/250 | 0.002 (0.003) | 0.908 (0.018) [0.866, 0.938] | 0.188 (0.001) |
| stress/noisy_reference | mice_reference_pmm | 500/500 | -0.106 (0.003) | 0.512 (0.022) [0.468, 0.556] | 0.214 (0.001) |
| stress/nonignorable_selection | mice_reference_pmm | 1000/1000 | 0.281 (0.002) | 0.014 (0.004) [0.008, 0.023] | 0.211 (0.001) |
| support/site_shift_gold10 | mice_reference_pmm | 250/250 | 0.003 (0.008) | 0.748 (0.027) [0.691, 0.798] | 0.300 (0.007) |
| support/site_shift_gold50 | mice_reference_pmm | 250/250 | 0.002 (0.003) | 0.908 (0.018) [0.866, 0.938] | 0.188 (0.001) |
| support/missed_flags_gold10 | oracle | 250/250 | 0.002 (0.002) | 0.940 (0.015) [0.903, 0.963] | 0.161 (0.000) |
| support/missed_flags_gold50 | oracle | 250/250 | 0.002 (0.002) | 0.940 (0.015) [0.903, 0.963] | 0.161 (0.000) |
| stress/noisy_reference | oracle | 500/500 | 0.003 (0.002) | 0.954 (0.009) [0.932, 0.969] | 0.161 (0.000) |
| stress/nonignorable_selection | oracle | 1000/1000 | 0.002 (0.001) | 0.952 (0.007) [0.937, 0.964] | 0.161 (0.000) |
| support/site_shift_gold10 | oracle | 250/250 | 0.002 (0.002) | 0.940 (0.015) [0.903, 0.963] | 0.161 (0.000) |
| support/site_shift_gold50 | oracle | 250/250 | 0.002 (0.002) | 0.940 (0.015) [0.903, 0.963] | 0.161 (0.000) |
| support/missed_flags_gold10 | reference_normal | 250/250 | -0.020 (0.007) | 0.960 (0.012) [0.928, 0.978] | 0.471 (0.008) |
| support/missed_flags_gold50 | reference_normal | 250/250 | -0.002 (0.003) | 0.960 (0.012) [0.928, 0.978] | 0.206 (0.001) |
| stress/noisy_reference | reference_normal | 500/500 | -0.112 (0.003) | 0.678 (0.021) [0.636, 0.717] | 0.275 (0.002) |
| stress/nonignorable_selection | reference_normal | 1000/1000 | 0.264 (0.002) | 0.055 (0.007) [0.042, 0.071] | 0.281 (0.001) |
| support/site_shift_gold10 | reference_normal | 250/250 | -0.020 (0.007) | 0.960 (0.012) [0.928, 0.978] | 0.471 (0.008) |
| support/site_shift_gold50 | reference_normal | 250/250 | -0.002 (0.003) | 0.960 (0.012) [0.928, 0.978] | 0.206 (0.001) |
| support/missed_flags_gold10 | reference_pmm | 250/250 | 0.003 (0.008) | 0.748 (0.027) [0.691, 0.798] | 0.300 (0.007) |
| support/missed_flags_gold50 | reference_pmm | 250/250 | 0.002 (0.003) | 0.908 (0.018) [0.866, 0.938] | 0.188 (0.001) |
| stress/noisy_reference | reference_pmm | 500/500 | -0.106 (0.003) | 0.512 (0.022) [0.468, 0.556] | 0.214 (0.001) |
| stress/nonignorable_selection | reference_pmm | 1000/1000 | 0.281 (0.002) | 0.014 (0.004) [0.008, 0.023] | 0.211 (0.001) |
| support/site_shift_gold10 | reference_pmm | 250/250 | 0.003 (0.008) | 0.748 (0.027) [0.691, 0.798] | 0.300 (0.007) |
| support/site_shift_gold50 | reference_pmm | 250/250 | 0.002 (0.003) | 0.908 (0.018) [0.866, 0.938] | 0.188 (0.001) |
| support/missed_flags_gold10 | v02_mask_normal | 250/250 | -0.080 (0.003) | 0.580 (0.031) [0.518, 0.640] | 0.174 (0.001) |
| support/missed_flags_gold50 | v02_mask_normal | 250/250 | -0.080 (0.003) | 0.580 (0.031) [0.518, 0.640] | 0.174 (0.001) |
| stress/noisy_reference | v02_mask_normal | 500/500 | -0.077 (0.002) | 0.596 (0.022) [0.552, 0.638] | 0.174 (0.000) |
| stress/nonignorable_selection | v02_mask_normal | 1000/1000 | 0.001 (0.001) | 0.964 (0.006) [0.951, 0.974] | 0.172 (0.000) |
| support/site_shift_gold10 | v02_mask_normal | 250/250 | -0.143 (0.003) | 0.088 (0.018) [0.059, 0.130] | 0.172 (0.001) |
| support/site_shift_gold50 | v02_mask_normal | 250/250 | -0.143 (0.003) | 0.088 (0.018) [0.059, 0.130] | 0.172 (0.001) |

Normal/PMM comparisons must hold category information constant. Reference workflows obtain additional measured truth; they do not identify missed errors without that information. Reference-only and independent mice comparators use the same gold measurements. Nonignorable selection, noisy gold and sparse donor support are deliberate failure probes.

Summary, paired contrasts, bias decomposition and all requests are in adjacent CSVs. Bias-centered coverage is retrospective and uses hidden truth; it is not a real-data correction. Source hashes match before and after execution. No result is made valid merely by this script completing.
