# Prospective test: replacing the resampling in clustered imputation

Written 2026-11-01 and committed before either candidate is in the package or has run on any real item, so version history records the order. Amendments, if any, go in the section at the end.

## Why

`catmech_reference(cluster = , cluster_strata = )` resamples audited clusters within strata. A bootstrap mean within a stratum of n audited clusters has (n - 1)/n of the usual variance, and a stratum with one audited cluster contributes none. The [first test](CLUSTERED_PLAN.md), which came out [PARTIAL](CLUSTERED_SYNTHESIS.md), did not anticipate this, and a 25% audit there holds 17 clusters in seven catchments (1, 2, 2, 2, 3, 3 and 4). A review also found that the first plan's simulations barely contain the clustering the option is for. This test picks a replacement for the resampling and decides whether it earns a recommendation.

The textbook correction, a rescaled bootstrap (Rao and Wu, 1988), is a comparator here and not a candidate, because of the scratch results listed next.

## What has already been seen

This test is held out by method, not by data.

- **All nine Rwanda items of the first test.** Results for linear imputation, the first option, leaving out strata, resampling rows and a scratch rescaled bootstrap, at both audit sizes, on the same cohorts and seeds used below. Both candidates were designed knowing those results, including the cost on rare items.
- **A scratch run on simulated audits, on other seeds, which led to the choice of candidates.** This is coverage at a 25% audit, 500 cohorts per setting, 69 clusters in seven strata. [Scripts and output](clustered-fix-scratch/) are kept.

| Setting | Oracle | Existing | First option | No strata | Rescaled | Cluster-robust normal draw | Jackknife draw | Random intercept |
|---|---|---|---|---|---|---|---|---|
| `shift` | 93.0 | 85.6 | 93.6 | 96.4 | 97.2 | 85.2 | 94.4 | 96.6 |
| `under` | 94.0 | 86.4 | 93.6 | 97.0 | 98.0 | 84.6 | 94.8 | 93.8 |
| `rare` | 82.2 | 88.4 | 88.2 | 91.2 | 95.8 | 85.0 | 89.4 | 90.2 |
| `row_only` | 90.4 | 95.4 | 95.4 | 97.8 | 98.6 | 89.8 | 95.2 | 96.8 |

  The rescaled bootstrap overshoots, and its recovery RMSE ratio against the existing option ran 1.076 to 1.302. The cluster-robust normal draw was dropped because of this table. In that run the jackknife draw had no rule for a stratum with one audited cluster.
- **Simulation settings.** `shift`, `under`, `rare` and `row_only`, in the `rwanda` design, were used in that look. `shift_low`, `outlier`, the `pairs` and `none` designs, and the one-row controls are new.
- **The two unused items.** Their agreement counts and means, during the file inspection and again while planning: 32 of 207 rows disagree on `ancnewreg` and 49 of 207 on `livebirths`.
- **Residual correlation on very large draws of each setting,** used to set its parameters. No method was run on those draws.

Neither candidate has run on a real item.

## The candidates

Both replace `catmech_reference_cluster_fill()`, which the option calls for records with a recorded value and for those without. Both center each imputation's coefficients on the audit's own fit and never refit on resampled clusters, the likely source of the first option's cost on rare items. Both refuse fewer than five audited clusters or an audit design without full column rank. Neither is a new statistical method.

**`jackknife`.** Fit the recorded-value regression to all audited rows, giving coefficients b and residual variance s2. In each stratum with at least two audited clusters, leave out each cluster in turn and refit, giving coefficient differences d_hi. The covariance is the stratified jackknife (Rust and Rao, 1996):

V = H / (H - s) x sum over those strata of (n_h - 1) / k_h x sum over i of d_hi d_hi'

Here H is the number of strata, s is the number that cannot contribute, n_h is the stratum's audited clusters, and k_h is how many of them could be left out. A left-out cluster that loses full rank or leaves no residual degrees of freedom is skipped, and k_h falls by one. A stratum with one audited cluster, or none that can be left out, cannot contribute, and the H / (H - s) factor replaces its contribution by the average across strata, as the `survey` package's `"average"` option does for such strata. If no stratum can contribute, the fit is refused. For a stratified mean, this V equals the usual variance with no shortfall. With no `cluster_strata`, it is the delete-one-cluster jackknife (MacKinnon, Nielsen and Webb, 2023).

Each imputation draws coefficients from N(b, V), draws the residual variance as s2 (n - p) / chi-square(n - p), and adds normal noise. Imputed noise stays independent across rows.

**`random_intercept`.** The model is `target = x'beta + c_i + e_ij`, with strata entering only as predictors. Variance components come from moments (Henderson, 1953, method 3). From the fit with cluster indicators, SSE_w has df_w degrees of freedom. The rest of the pooled fit's residual sum of squares is SS_b with df_b. Each imputation draws sigma2 = SSE_w / chi-square(df_w) and tau2 = max(0, df_b x (SS_b / chi-square(df_b) - sigma2) / c), where c is the trace of the cluster incidence matrix projected off the predictors. Given those, it takes the generalized least squares fit and draws coefficients from its covariance. Unaudited clusters get c_i from N(0, tau2), and audited clusters with unaudited rows get c_i from its conditional distribution given their residuals. Noise e_ij comes from N(0, sigma2). With no way to separate the components, as with one row per cluster, tau2 is zero and sigma2 comes from the pooled fit. This is ordinary multilevel imputation (Schafer and Yucel, 2002; Grund, Lüdtke and Robitzsch, 2018) in closed form, with no added dependency.

## Comparators

`oracle` (reference on every row), `linear_mi` (the existing option), `linear_mi_clustered` (the first option, within strata), `no_strata` (the first option without strata) and `rao_wu` (rescaled bootstrap: n_h - 1 clusters drawn in each stratum and weighted n_h / (n_h - 1), a stratum with one cluster keeping it). All use 20 imputations where they impute, and all intervals come from `catmech_pool_lm(cluster = )`, so only the imputation differs.

## Real data

The Rwanda recount (`validation/real-audit/rwanda-pool.csv`), 69 health centers, facility-months. Each analysis is `target ~ covariate + catchment`, on rows with all values observed. Cohorts and audits follow the first plan: facilities resampled within catchment, then a stratified random audit of facilities at 25% and 50%, with `cluster_strata` set to catchment where a method uses it. The target is the pool slope. 500 cohorts per cell.

- **The five items of the first test:** `lbw`, `resus`, `stillbirths`, `deliveries` and `anc1_early`, with the first plan's definitions and seeds, so cohorts and audits are identical.
- **The four development items:** `pnc`, `syphilis`, `iron` and `anc4`, with the first study's seeds.
- **Two unused items, the only ones no method has been run on.** Both are common events with few disagreements, so they check for harm and cannot show a gain. Their seeds follow the first plan's held-out formula with index 6 and 7.

| Analysis | Recorded | Reference | Covariate, as reported |
|---|---|---|---|
| `ancnewreg` | `ancnewreg_h` | `ancnewreg_r` | `hf_deliveries_h` |
| `livebirths` | `livebirths_h` | `livebirths_r` | `ancnewreg_h` |

## Simulations with known truth

Clusters i in strata h, with three rows each: `x_ij = s_i + noise`, with `s_i` lognormal (log mean 3, log SD 0.5) and noise SD 3, and `T_ij = 2 + 0.6 x_ij + a_h + u_i + e_ij`, with `a_h` from N(0, 4), `u_i` from N(0, 20) and `e_ij` from N(0, 80). The target is 0.6.

| Design | Clusters per stratum | Audited at 25% | Audited at 50% | `cluster_strata` |
|---|---|---|---|---|
| `rwanda` | 14, 9, 8, 13, 5, 9, 11 | 17, one stratum with one | 33 | stratum |
| `pairs` | 8 in each of 10 strata | 20, two per stratum | 40, four per stratum | stratum |
| `none` | 70, no strata (`a_h` = 0) | 18 | 35 | none |

Stratified designs use `max(1, round(fraction x size))` in each stratum and include the stratum in the analysis and imputation models.

| Errors | Recorded value | Residual correlation on large draws |
|---|---|---|
| `shift_low` | `T + v_i + n_ij`, `v_i` SD 3, `n_ij` SD 6 | 0.20 |
| `shift` | `T + v_i + n_ij`, `v_i` SD 5, `n_ij` SD 3 | 0.60 |
| `under` | `(1 - d_i) T + n_ij`, `d_i` from Beta(2, 6), `n_ij` SD 2, which biases the recorded slope | 0.48 to 0.49 |
| `outlier` | `T + v_i + n_ij`, `v_i` zero in 90% of clusters and SD 20 in the rest, `n_ij` SD 3 | 0.62 |
| `rare` | Counts: `T_ij` Poisson with mean `(0.1 + 0.02 x_ij + a_h) f_i`, the bracket floored at 0.01, `f_i` from Gamma(0.5, 0.5), `a_h` uniform on 0 to 0.2. The target is 0.02. Recorded equals `T` except that 20% of rows are off by one, floored at zero | 0.04 to 0.05 |
| `row_only` | `T + n_ij`, `n_ij` SD 6, with `T` carrying no cluster effect | 0.00 to 0.01 |

The first four are the **clustered cells**, 24 in all at both fractions, and their correlation is that of the imputation model's residuals, within the 0.20 to 0.67 seen on Rwanda's common items. With `rare` and `row_only`, and two one-row controls (`row_only`, 210 clusters, no strata, 25% and 50%), that makes **38 cells**.

Seeds: cohort 1,300,100,000 + cell index x 100,000 + replicate. Imputation 1,310,000,000 + cell index x 100,000 + replicate. Pilots use 1,400,100,000 and 1,410,000,000 for simulations, and 1,420,100,000, 1,430,100,000 and 1,440,000,000 for real items.

## Measures

Clustered coverage of the target, mean and median clustered width, returned intervals, recovery RMSE (distance from the same cohort's oracle), RMSE against the target, between-imputation variance, degrees of freedom and widths over five times their cell median. Ratios against `linear_mi` get a 95% percentile interval from 2,000 paired resamples.

## Selecting the candidate

Simulations alone select the candidate, before any real result is read:

1. A candidate that refuses more than 1% of its simulated requests is out.
2. The one with the smaller mean absolute distance of coverage from 95% across the 24 clustered cells is selected.
3. If the two are within 0.5 points, the one with the lower mean recovery RMSE ratio against `linear_mi` across all 38 cells is selected.

Only the selected candidate enters the verdict. The other is reported and is never substituted.

## Decision rule

At a 25% audit on real items:

- **Coverage holds** if, on at least eight of the nine items, coverage is at least 90% or at least the oracle's, and at most 98%, with `anc1_early`, the item the first option missed, among the eight.
- **Accuracy holds** if the upper limit of the recovery RMSE ratio against `linear_mi` is at most 1.05 on at least four of the five items. This is the first test's condition, unchanged.
- **Unused items hold** if, on both, coverage is within 5 points of the oracle's and that upper limit is at most 1.05.

In the simulations, at 25% and 50%:

- **Simulation coverage holds** if, in at least 80% of the 24 clustered cells, coverage is at least 92% or at least the oracle's, and at most 98%.
- **Simulation accuracy holds** if that upper limit is at most 1.05 in at least 80% of the 38 cells.

The verdict:

- **GO** if all five hold.
- **PARTIAL** if coverage holds and any other condition fails.
- **NO-GO** if coverage fails.

The first test judged coverage partly by moving closer to the oracle. A review found that this flattered the option, because falling from 96% toward an oracle of 94% counted as a gain, and it would count leaving a healthy item alone as a failure. So coverage is now judged against fixed limits, on all nine items and not four of five. The first test's coverage condition is still scored for comparison and enters no verdict.

The verdict is reported as it falls, with no change of thresholds after viewing.

## What each verdict does

- **GO.** The selected candidate replaces the resampling, as version 0.6.0. The documentation drops "experimental" and recommends the option for audits that sample whole clusters, stating its limits: one health system, items already seen under other methods, simulation settings chosen by the same assistant, and no human review.
- **PARTIAL.** The same replacement, and the option stays experimental with the failure stated.
- **NO-GO.** Package code is unchanged, the flaw stays documented, and nothing further is tried in this round.

If `random_intercept` is selected, `cluster_strata` is removed, because strata then enter only as predictors. Any replacement changes results for calls with `cluster`. Calls without it must stay numerically identical.

## Predictions fixed in advance

1. The simulations select `random_intercept`.
2. `rao_wu` covers above 96% in most clustered cells of the stratified designs at 25%.
3. Both candidates are closer to 95% on average than `linear_mi_clustered` across the 24 clustered cells.
4. On `anc1_early` at 25%, the selected candidate covers at least 90%.
5. On `lbw`, `resus` and `stillbirths` at 25%, its recovery RMSE ratio is below the first option's on all three.
6. On the six common items at 25%, its recovery RMSE ratio is at most 1.02 on all six.
7. In the one-row controls, its mean width is within 5% of `linear_mi`'s.
8. On the two unused items at 25%, its coverage is within 3 points of `linear_mi`'s.
9. On the five items at 25%, its mean width is at most 40% wider than `linear_mi`'s on all five.
10. Across all requests, it has fewer widths over five times their cell median than `linear_mi_clustered`.

Each prediction is scored in the write-up, including those that fail.

## Process

An engineering pilot, with separate seeds and three cohorts per scenario, is read for request statuses and timing only. If it shows an engineering fault that needs a change, the change is listed below before the study runs. Hashes of the package sources, candidates, runner, scorer, this plan and the data are checked before and after the study. As a control, `linear_mi` and `linear_mi_clustered` on the nine items must reproduce the first study's requests, and `no_strata` and `rao_wu` the review's.

The package's code does not change until the study is scored. The selected candidate then moves into `R/reference.R`, and a replay through the package must reproduce every study request for that method. Results are not retuned after viewing. A further change requires a new plan and a new run.

## Sources

Rao JNK, Wu CFJ. Resampling inference with complex survey data. *Journal of the American Statistical Association* 1988. [doi:10.1080/01621459.1988.10478591](https://doi.org/10.1080/01621459.1988.10478591). Rust KF, Rao JNK. Variance estimation for complex surveys using replication techniques. *Statistical Methods in Medical Research* 1996. [doi:10.1177/096228029600500305](https://doi.org/10.1177/096228029600500305). MacKinnon JG, Nielsen MØ, Webb MD. Cluster-robust inference: a guide to empirical practice. *Journal of Econometrics* 2023. [doi:10.1016/j.jeconom.2022.04.001](https://doi.org/10.1016/j.jeconom.2022.04.001). Henderson CR. Estimation of variance and covariance components. *Biometrics* 1953. [doi:10.2307/3001853](https://doi.org/10.2307/3001853). Schafer JL, Yucel RM. Computational strategies for multivariate linear mixed-effects models with missing values. *Journal of Computational and Graphical Statistics* 2002. [doi:10.1198/106186002760180608](https://doi.org/10.1198/106186002760180608). Grund S, Lüdtke O, Robitzsch A. Multiple imputation of missing data for multilevel models: simulations and recommendations. *Organizational Research Methods* 2018. [doi:10.1177/1094428117703686](https://doi.org/10.1177/1094428117703686). The `survey` package's [options page](https://r-survey.r-forge.r-project.org/survey/html/surveyoptions.html) describes the `"average"` setting for a stratum with one primary sampling unit.

## Amendments

None yet.
