# v0.3 extension evaluation: controls

Prespecified settings, assumptions and interpretation boundaries are in ../EXTENSION_PLAN.md. Coverage is conditional on returned intervals; request-status.csv retains failures. Same-reference results across recording scenarios reuse cohorts and reference subsets, so they are not independent replications.

| Setting | Method | Intervals/requested | Bias (MCSE) | Coverage (MCSE) [95% Wilson] | Width (MCSE) |
|---|---|---:|---:|---:|---:|
| control/all_reference | reference_normal | 20/20 | 0.010 (0.009) | 0.900 (0.067) [0.699, 0.972] | 0.162 (0.001) |
| control/no_reference | reference_normal | 0/20 | NA (NA) | NA (NA) [NA, NA] | NA (NA) |
| control/all_reference | reference_pmm | 20/20 | 0.010 (0.009) | 0.900 (0.067) [0.699, 0.972] | 0.162 (0.001) |
| control/no_reference | reference_pmm | 0/20 | NA (NA) | NA (NA) [NA, NA] | NA (NA) |

Normal/PMM comparisons must hold category information constant. Reference workflows obtain additional measured truth; they do not identify missed errors without that information. Reference-only and independent mice comparators use the same gold measurements. Nonignorable selection, noisy gold and sparse donor support are deliberate failure probes.

Summary, paired contrasts, bias decomposition and all requests are in adjacent CSVs. Bias-centered coverage is retrospective and uses hidden truth; it is not a real-data correction. Source hashes match before and after execution. No result is made valid merely by this script completing.
