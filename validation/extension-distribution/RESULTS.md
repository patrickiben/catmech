# v0.3 extension evaluation: distribution

Prespecified settings, assumptions and interpretation boundaries are in ../EXTENSION_PLAN.md. Coverage is conditional on returned intervals; request-status.csv retains failures. Same-reference results across recording scenarios reuse cohorts and reference subsets, so they are not independent replications.

| Setting | Method | Intervals/requested | Bias (MCSE) | Coverage (MCSE) [95% Wilson] | Width (MCSE) |
|---|---|---:|---:|---:|---:|
| distribution/category_error | mice_normal_cat | 1000/1000 | 0.001 (0.002) | 0.955 (0.007) [0.940, 0.966] | 0.250 (0.001) |
| distribution/false_positive_flags | mice_normal_cat | 1000/1000 | -0.002 (0.002) | 0.944 (0.007) [0.928, 0.957] | 0.244 (0.001) |
| distribution/heavy_category_error | mice_normal_cat | 1000/1000 | -0.001 (0.003) | 0.918 (0.009) [0.899, 0.933] | 0.270 (0.001) |
| distribution/heavy_tails | mice_normal_cat | 1000/1000 | 0.006 (0.003) | 0.908 (0.009) [0.888, 0.924] | 0.265 (0.001) |
| distribution/missing_only | mice_normal_cat | 1000/1000 | 0.000 (0.002) | 0.944 (0.007) [0.928, 0.957] | 0.239 (0.000) |
| distribution/null_effect | mice_normal_cat | 1000/1000 | -0.000 (0.002) | 0.949 (0.007) [0.934, 0.961] | 0.242 (0.001) |
| distribution/perfect_flags | mice_normal_cat | 1000/1000 | -0.002 (0.002) | 0.948 (0.007) [0.932, 0.960] | 0.245 (0.001) |
| distribution/strong_mar | mice_normal_cat | 1000/1000 | 0.002 (0.002) | 0.947 (0.007) [0.931, 0.959] | 0.273 (0.001) |
| distribution/category_error | mice_normal_no_cat | 1000/1000 | 0.002 (0.002) | 0.954 (0.007) [0.939, 0.965] | 0.254 (0.001) |
| distribution/false_positive_flags | mice_normal_no_cat | 1000/1000 | -0.003 (0.002) | 0.954 (0.007) [0.939, 0.965] | 0.255 (0.001) |
| distribution/heavy_category_error | mice_normal_no_cat | 1000/1000 | 0.000 (0.003) | 0.910 (0.009) [0.891, 0.926] | 0.272 (0.001) |
| distribution/heavy_tails | mice_normal_no_cat | 1000/1000 | 0.006 (0.003) | 0.917 (0.009) [0.898, 0.933] | 0.273 (0.002) |
| distribution/missing_only | mice_normal_no_cat | 1000/1000 | 0.000 (0.002) | 0.940 (0.008) [0.924, 0.953] | 0.246 (0.001) |
| distribution/null_effect | mice_normal_no_cat | 1000/1000 | -0.001 (0.002) | 0.950 (0.007) [0.935, 0.962] | 0.273 (0.001) |
| distribution/perfect_flags | mice_normal_no_cat | 1000/1000 | -0.003 (0.002) | 0.960 (0.006) [0.946, 0.970] | 0.255 (0.001) |
| distribution/strong_mar | mice_normal_no_cat | 1000/1000 | -0.005 (0.002) | 0.947 (0.007) [0.931, 0.959] | 0.294 (0.001) |
| distribution/category_error | mice_pmm_cat | 1000/1000 | 0.001 (0.002) | 0.952 (0.007) [0.937, 0.964] | 0.252 (0.001) |
| distribution/false_positive_flags | mice_pmm_cat | 1000/1000 | -0.002 (0.002) | 0.945 (0.007) [0.929, 0.958] | 0.244 (0.001) |
| distribution/heavy_category_error | mice_pmm_cat | 1000/1000 | 0.006 (0.002) | 0.930 (0.008) [0.912, 0.944] | 0.279 (0.002) |
| distribution/heavy_tails | mice_pmm_cat | 1000/1000 | 0.014 (0.002) | 0.926 (0.008) [0.908, 0.941] | 0.272 (0.002) |
| distribution/missing_only | mice_pmm_cat | 1000/1000 | 0.001 (0.002) | 0.947 (0.007) [0.931, 0.959] | 0.239 (0.000) |
| distribution/null_effect | mice_pmm_cat | 1000/1000 | 0.000 (0.002) | 0.946 (0.007) [0.930, 0.958] | 0.240 (0.000) |
| distribution/perfect_flags | mice_pmm_cat | 1000/1000 | -0.002 (0.002) | 0.943 (0.007) [0.927, 0.956] | 0.244 (0.001) |
| distribution/strong_mar | mice_pmm_cat | 1000/1000 | -0.012 (0.002) | 0.930 (0.008) [0.912, 0.944] | 0.275 (0.001) |
| distribution/category_error | mice_pmm_no_cat | 1000/1000 | 0.003 (0.002) | 0.941 (0.007) [0.925, 0.954] | 0.249 (0.001) |
| distribution/false_positive_flags | mice_pmm_no_cat | 1000/1000 | -0.002 (0.002) | 0.945 (0.007) [0.929, 0.958] | 0.249 (0.001) |
| distribution/heavy_category_error | mice_pmm_no_cat | 1000/1000 | 0.004 (0.002) | 0.928 (0.008) [0.910, 0.942] | 0.278 (0.002) |
| distribution/heavy_tails | mice_pmm_no_cat | 1000/1000 | 0.009 (0.002) | 0.940 (0.008) [0.924, 0.953] | 0.278 (0.002) |
| distribution/missing_only | mice_pmm_no_cat | 1000/1000 | 0.001 (0.002) | 0.938 (0.008) [0.921, 0.951] | 0.243 (0.001) |
| distribution/null_effect | mice_pmm_no_cat | 1000/1000 | -0.001 (0.002) | 0.950 (0.007) [0.935, 0.962] | 0.272 (0.001) |
| distribution/perfect_flags | mice_pmm_no_cat | 1000/1000 | -0.002 (0.002) | 0.945 (0.007) [0.929, 0.958] | 0.249 (0.001) |
| distribution/strong_mar | mice_pmm_no_cat | 1000/1000 | -0.015 (0.003) | 0.911 (0.009) [0.892, 0.927] | 0.280 (0.001) |
| distribution/category_error | normal_cat | 1000/1000 | 0.001 (0.002) | 0.956 (0.006) [0.941, 0.967] | 0.250 (0.001) |
| distribution/false_positive_flags | normal_cat | 1000/1000 | -0.002 (0.002) | 0.945 (0.007) [0.929, 0.958] | 0.244 (0.001) |
| distribution/heavy_category_error | normal_cat | 1000/1000 | -0.001 (0.003) | 0.912 (0.009) [0.893, 0.928] | 0.270 (0.001) |
| distribution/heavy_tails | normal_cat | 1000/1000 | 0.006 (0.003) | 0.911 (0.009) [0.892, 0.927] | 0.265 (0.001) |
| distribution/missing_only | normal_cat | 1000/1000 | 0.000 (0.002) | 0.945 (0.007) [0.929, 0.958] | 0.239 (0.000) |
| distribution/null_effect | normal_cat | 1000/1000 | -0.000 (0.002) | 0.949 (0.007) [0.934, 0.961] | 0.242 (0.001) |
| distribution/perfect_flags | normal_cat | 1000/1000 | -0.002 (0.002) | 0.949 (0.007) [0.934, 0.961] | 0.244 (0.001) |
| distribution/strong_mar | normal_cat | 1000/1000 | 0.002 (0.002) | 0.943 (0.007) [0.927, 0.956] | 0.272 (0.001) |
| distribution/category_error | normal_no_cat | 1000/1000 | 0.002 (0.002) | 0.954 (0.007) [0.939, 0.965] | 0.254 (0.001) |
| distribution/false_positive_flags | normal_no_cat | 1000/1000 | -0.003 (0.002) | 0.953 (0.007) [0.938, 0.964] | 0.255 (0.001) |
| distribution/heavy_category_error | normal_no_cat | 1000/1000 | 0.000 (0.003) | 0.913 (0.009) [0.894, 0.929] | 0.273 (0.001) |
| distribution/heavy_tails | normal_no_cat | 1000/1000 | 0.006 (0.003) | 0.917 (0.009) [0.898, 0.933] | 0.273 (0.002) |
| distribution/missing_only | normal_no_cat | 1000/1000 | 0.000 (0.002) | 0.940 (0.008) [0.924, 0.953] | 0.246 (0.001) |
| distribution/null_effect | normal_no_cat | 1000/1000 | -0.001 (0.002) | 0.955 (0.007) [0.940, 0.966] | 0.273 (0.001) |
| distribution/perfect_flags | normal_no_cat | 1000/1000 | -0.003 (0.002) | 0.958 (0.006) [0.944, 0.969] | 0.255 (0.001) |
| distribution/strong_mar | normal_no_cat | 1000/1000 | -0.005 (0.002) | 0.949 (0.007) [0.934, 0.961] | 0.294 (0.001) |
| distribution/category_error | oracle | 1000/1000 | 0.004 (0.002) | 0.954 (0.007) [0.939, 0.965] | 0.229 (0.000) |
| distribution/false_positive_flags | oracle | 1000/1000 | -0.001 (0.002) | 0.957 (0.006) [0.943, 0.968] | 0.230 (0.000) |
| distribution/heavy_category_error | oracle | 1000/1000 | -0.001 (0.002) | 0.946 (0.007) [0.930, 0.958] | 0.241 (0.001) |
| distribution/heavy_tails | oracle | 1000/1000 | 0.001 (0.002) | 0.949 (0.007) [0.934, 0.961] | 0.240 (0.001) |
| distribution/missing_only | oracle | 1000/1000 | 0.001 (0.002) | 0.949 (0.007) [0.934, 0.961] | 0.229 (0.000) |
| distribution/null_effect | oracle | 1000/1000 | 0.001 (0.002) | 0.947 (0.007) [0.931, 0.959] | 0.229 (0.000) |
| distribution/perfect_flags | oracle | 1000/1000 | -0.001 (0.002) | 0.941 (0.007) [0.925, 0.954] | 0.229 (0.000) |
| distribution/strong_mar | oracle | 1000/1000 | -0.002 (0.002) | 0.953 (0.007) [0.938, 0.964] | 0.229 (0.000) |
| distribution/category_error | pmm_cat | 1000/1000 | 0.001 (0.002) | 0.952 (0.007) [0.937, 0.964] | 0.252 (0.001) |
| distribution/false_positive_flags | pmm_cat | 1000/1000 | -0.002 (0.002) | 0.945 (0.007) [0.929, 0.958] | 0.244 (0.001) |
| distribution/heavy_category_error | pmm_cat | 1000/1000 | 0.006 (0.002) | 0.930 (0.008) [0.912, 0.944] | 0.279 (0.002) |
| distribution/heavy_tails | pmm_cat | 1000/1000 | 0.014 (0.002) | 0.926 (0.008) [0.908, 0.941] | 0.272 (0.002) |
| distribution/missing_only | pmm_cat | 1000/1000 | 0.001 (0.002) | 0.947 (0.007) [0.931, 0.959] | 0.239 (0.000) |
| distribution/null_effect | pmm_cat | 1000/1000 | 0.000 (0.002) | 0.946 (0.007) [0.930, 0.958] | 0.240 (0.000) |
| distribution/perfect_flags | pmm_cat | 1000/1000 | -0.002 (0.002) | 0.943 (0.007) [0.927, 0.956] | 0.244 (0.001) |
| distribution/strong_mar | pmm_cat | 1000/1000 | -0.012 (0.002) | 0.930 (0.008) [0.912, 0.944] | 0.275 (0.001) |
| distribution/category_error | pmm_no_cat | 1000/1000 | 0.003 (0.002) | 0.941 (0.007) [0.925, 0.954] | 0.249 (0.001) |
| distribution/false_positive_flags | pmm_no_cat | 1000/1000 | -0.002 (0.002) | 0.945 (0.007) [0.929, 0.958] | 0.249 (0.001) |
| distribution/heavy_category_error | pmm_no_cat | 1000/1000 | 0.004 (0.002) | 0.928 (0.008) [0.910, 0.942] | 0.278 (0.002) |
| distribution/heavy_tails | pmm_no_cat | 1000/1000 | 0.009 (0.002) | 0.940 (0.008) [0.924, 0.953] | 0.278 (0.002) |
| distribution/missing_only | pmm_no_cat | 1000/1000 | 0.001 (0.002) | 0.938 (0.008) [0.921, 0.951] | 0.243 (0.001) |
| distribution/null_effect | pmm_no_cat | 1000/1000 | -0.001 (0.002) | 0.950 (0.007) [0.935, 0.962] | 0.272 (0.001) |
| distribution/perfect_flags | pmm_no_cat | 1000/1000 | -0.002 (0.002) | 0.945 (0.007) [0.929, 0.958] | 0.249 (0.001) |
| distribution/strong_mar | pmm_no_cat | 1000/1000 | -0.015 (0.003) | 0.911 (0.009) [0.892, 0.927] | 0.280 (0.001) |
| distribution/category_error | trusted_observed | 1000/1000 | 0.005 (0.002) | 0.958 (0.006) [0.944, 0.969] | 0.271 (0.001) |
| distribution/false_positive_flags | trusted_observed | 1000/1000 | -0.000 (0.002) | 0.947 (0.007) [0.931, 0.959] | 0.272 (0.001) |
| distribution/heavy_category_error | trusted_observed | 1000/1000 | -0.002 (0.002) | 0.943 (0.007) [0.927, 0.956] | 0.288 (0.001) |
| distribution/heavy_tails | trusted_observed | 1000/1000 | 0.001 (0.002) | 0.946 (0.007) [0.930, 0.958] | 0.288 (0.001) |
| distribution/missing_only | trusted_observed | 1000/1000 | 0.002 (0.002) | 0.935 (0.008) [0.918, 0.949] | 0.257 (0.001) |
| distribution/null_effect | trusted_observed | 1000/1000 | -0.001 (0.002) | 0.951 (0.007) [0.936, 0.963] | 0.272 (0.001) |
| distribution/perfect_flags | trusted_observed | 1000/1000 | -0.004 (0.002) | 0.955 (0.007) [0.940, 0.966] | 0.271 (0.001) |
| distribution/strong_mar | trusted_observed | 1000/1000 | -0.001 (0.003) | 0.949 (0.007) [0.934, 0.961] | 0.329 (0.001) |

Normal/PMM comparisons must hold category information constant. Reference workflows obtain additional measured truth; they do not identify missed errors without that information. Reference-only and independent mice comparators use the same gold measurements. Nonignorable selection, noisy gold and sparse donor support are deliberate failure probes.

Summary, paired contrasts, bias decomposition and all requests are in adjacent CSVs. Bias-centered coverage is retrospective and uses hidden truth; it is not a real-data correction. Source hashes match before and after execution. No result is made valid merely by this script completing.
