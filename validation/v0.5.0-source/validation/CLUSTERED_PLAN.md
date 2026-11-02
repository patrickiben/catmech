# Prospective test: imputation for clustered audits

Written 2026-11-01 and committed before the option exists in the package, so version history records the order. Amendments, if any, go in the section at the end.

## Why

Real audits cluster by site, and `catmech_reference(use_recorded = TRUE)` treats rows as independent. In the [real-audit test](REAL_AUDIT_SYNTHESIS.md) its clustered coverage at a 25% facility audit ran 8.6 to 10.7 points under a recount's. A scratch version that resampled audited facilities in each imputation reached about 91% on 200 cohorts, but that was post hoc, on the analyses that suggested it. This test decides whether the package should offer and recommend the idea.

The idea is not new. It combines the bootstrap at the imputation step (von Hippel and Bartlett, 2021) with the clustered bootstrap and cluster-robust inference (Cameron and Miller, 2015). No new statistical method is claimed.

## What is being tested

Two additions, and existing calls must stay numerically identical:

- `catmech_reference(..., use_recorded = TRUE, cluster = "site", cluster_strata = NULL)`. For each imputation, audited clusters are resampled with replacement, within `cluster_strata` if given. The recorded-value regression is fitted to that resample, and unaudited rows get its prediction plus normal noise with its residual variance. Records with no recorded value use the same steps without that predictor. Fewer than five audited clusters, or a resample that keeps losing design support, is refused.
- `catmech_pool_lm(fit, formula, cluster = "site")`. A CR1 cluster-robust variance within each completed dataset, pooled by Rubin's rules with complete-data degrees of freedom G - 1. A fit made with `cluster` uses that column by default.

## Held-out real data

The Rwanda recount (`validation/real-audit/rwanda-pool.csv`), 69 health centers, facility-months, using items that have not been analyzed by any method. Only their disagreement counts have been seen, with hospitals included, during the file inspection. They share facilities and months with the analyses that suggested the idea, so this is held out by item and not by facility.

| Analysis | Recorded | Reference | Covariate, as reported |
|---|---|---|---|
| `lbw` | `lbw_notpremature_h + lbw_premature_h` | `lbw_r` | `hf_deliveries_h` |
| `resus` | `newborns_nocryres_succ_h` | `nocry_res_success_r` | `hf_deliveries_h` |
| `stillbirths` | `sbmacerated_h + sbfresh_h` | `sbmacerated_r + sbfresh_r` | `hf_deliveries_h` |
| `deliveries` | `hf_deliveries_h` | `hf_deliveries_r` | `ancnewreg_h` |
| `anc1_early` | `anc1vis1trim_h` | `anc1vis_16weeks_r` | `ancnewreg_h` |

Each is `target ~ covariate + catchment`, on rows with all values observed. The target recounts one item with the covariate as reported, as in the earlier test. `anc1_early`'s two sides are defined slightly differently, as the earlier plan notes. Cohorts and audits follow that plan: facilities resampled within catchment, then a stratified random audit of facilities at 25% and 50%, with `cluster_strata` set to catchment. 500 cohorts per cell.

Seeds: cohort 860,100,000 + index x 1,000,000 + replicate. Selection 870,100,000 and imputation 880,000,000, each + index x 1,000,000 + fraction index x 100,000 + replicate.

## Simulations with known truth

Clusters i with three rows each: `x_ij = s_i + noise`, with skewed cluster sizes `s_i`, and `T_ij = b0 + b1 x_ij + u_i + e_ij`. Recorded values equal the target, or add an error. The audit is a simple random sample of clusters, and `cluster_strata` is not used.

| Factor | Levels |
|---|---|
| Clusters | 30, 70, 150 |
| Intraclass correlation of the target residual | 0, 0.2, 0.5 |
| Errors | `iid` (40% of rows, additive normal), `persistent` (cluster-level error rate with mean 40%), `tracking` (`persistent`, with an error mean that grows with `x`, which biases the recorded slope) |
| Audit fraction | 25%, 50% |

That is 54 cells at 500 cohorts each. Two control cells have one row per cluster and 210 clusters, at each fraction, where the option should change little. The target is the true `b1`.

Seeds: cohort 900,100,000 + cell index x 100,000 + replicate. Imputation 910,000,000 + cell index x 100,000 + replicate.

## Methods and measures

`oracle` (reference on every row), `recorded`, `linear_mi` (the existing option) and `linear_mi_clustered` (the new one). All use 20 imputations where they impute, and all intervals come from `catmech_pool_lm(cluster = ...)`, so only the imputation differs.

Measures are clustered coverage of the target, mean clustered width, returned intervals, recovery RMSE (distance from the same cohort's oracle) and, for real data, RMSE against the pool slope. Ratios get a 95% percentile interval from 2,000 paired resamples.

The four earlier analyses are rerun with the study's own cohort seeds as a labeled development record. They enter no rule.

## Decision rule

At facility-months at a 25% audit, on the five held-out analyses:

- **Coverage holds** if `linear_mi_clustered` is closer to the oracle's coverage than `linear_mi` in at least four of five, and within 5 points of the oracle in at least four of five.
- **Accuracy holds** if the upper limit of its recovery RMSE ratio against `linear_mi` is at most 1.05 in at least four of five.

In the simulations, at 25% and 50%:

- **Simulation holds** if coverage is at least 92% and within 3 points of the oracle's in at least 80% of cells with 70 or more clusters.

The verdict:

- **GO** if all three hold. The documentation then recommends the option for clustered audits.
- **PARTIAL** if coverage holds and either other condition fails. The option stays, documented as experimental and with the failure stated.
- **NO-GO** if coverage fails. The option stays labeled not validated and is not recommended, as with the audit planner.

The verdict is reported as it falls, with no change of thresholds after viewing.

## Predictions fixed in advance

1. On the held-out analyses at 25%, `linear_mi_clustered` is closer to the oracle's coverage than `linear_mi` in at least four of five.
2. Its mean clustered width is 15% to 60% wider than `linear_mi`'s in all five.
3. Its recovery RMSE is no worse than `linear_mi`'s by more than 5% in at least four of five.
4. In the simulations with 30 clusters and a 25% audit, about eight audited clusters, its coverage falls below 92% in most cells.
5. In the one-row-per-cluster controls, its width is within 5% of `linear_mi`'s.
6. The gap between the two methods shrinks from 25% to 50% on the held-out analyses.

Each prediction is scored in the write-up, including those that fail.

## Process

An engineering pilot, with separate seeds and small counts, is read for request statuses and timing only. If it shows an engineering fault that needs a change, the change is listed below before the study runs. Package sources are frozen before the study, with hashes checked before and after, and existing workflows are rechecked against frozen sources. Results are not retuned after viewing. A change requires a new versioned run.

## Sources

von Hippel PT, Bartlett JW. Maximum likelihood multiple imputation: faster imputations and consistent standard errors without posterior draws. *Statistical Science* 2021. [doi:10.1214/20-STS793](https://doi.org/10.1214/20-STS793). Cameron AC, Miller DL. A practitioner's guide to cluster-robust inference. *Journal of Human Resources* 2015. [doi:10.3368/jhr.50.2.317](https://doi.org/10.3368/jhr.50.2.317).

## Amendments

None yet.
