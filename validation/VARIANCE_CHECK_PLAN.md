# Plan: scoring the recount under corrected cluster-robust variances

Written 2026-11-01 before this scoring runs. No imputation is involved, and no package code changes in this round.

## Why

In the [replacement test](CLUSTERED_FIX_SYNTHESIS.md), the oracle, a recount of every row, covered only 88.4% on `lbw` and about 83% on simulated rare counts. That analysis uses `catmech_pool_lm(cluster = )`, whose variance is what `clubSandwich` calls CR1S, with t(G - 1). If the analysis itself undercovers, a method can look better by offsetting a variance that is too small, so coverage is a weak measure of imputation. This check separates the two.

## What has already been seen

The oracle's coverage, width and error under the current variance on every cell below, all from earlier studies. A literature search's complete-data findings, which favor corrections, shaped the predictions. Nothing here has been run with a corrected variance.

## Data

The existing cohorts of the [replacement study](CLUSTERED_FIX_PLAN.md), reused through the runner's own functions and seeds. The oracle analyzes the reference value on every row, so the audit fraction does not change a real item's cohort: the 11 real items are scored once, using the 25% cells. The 38 simulated cells are scored as they are, each with its own cohorts. 500 cohorts per cell, 24,500 in all. The analysis is `target ~ covariate + catchment` on real items and the simulation's analysis model in each cell, and the target is the same as in that study.

## Variances

For coefficient k with G clusters, all use ordinary least squares estimates:

- **`cr1s`, the current one.** Factor (G / (G - 1)) ((n - 1) / (n - p)), with t(G - 1).
- **`cr2`.** Bias-reduced linearization with the symmetric inverse square root of I - H_gg, and Bell-McCaffrey Satterthwaite degrees of freedom under a working model of independent, equal-variance errors (Pustejovsky and Tipton, 2018).
- **`cv3`.** Delete-one-cluster jackknife, (G - 1) / G x the sum of squared differences from the full estimate, computed with the closed form using (I - H_gg)^-1, with t(G - 1) (MacKinnon, Nielsen and Webb, 2023). A singular I - H_gg uses a pseudoinverse.
- **`wcr`.** Restricted wild cluster bootstrap with Rademacher weights and 9,999 replications, testing the true target with the `cr1s` t statistic. A cohort covers if p > 0.05, so this method has coverage and no width. The seed is 1,600,000,000 + cell index x 100,000 + replicate, with real items as cells 101 to 111.

## Measures

Coverage of the target, mean width, mean standard error over the empirical standard deviation of the estimate, mean Satterthwaite degrees of freedom and mean largest partial leverage of a cluster, as diagnostics. With 500 cohorts, coverage has a Monte Carlo standard error of about 1 point at 95% and 1.3 at 90%.

## Reading fixed in advance

A variance is a **candidate to try through imputation** if both hold:

- On the 11 real items, coverage is at least 92% on at least nine and at most 98% on all.
- In the 38 simulated cells, coverage is at least 92% and at most 98% in at least 80%.

If more than one qualifies, none is picked here. This reading changes no package code and no earlier verdict. If none qualifies, the shortfall in the analysis is not repaired by these corrections at these sizes, and later comparisons must keep reading coverage against the oracle.

## Predictions fixed in advance

1. `cr2` and `cv3` both raise `lbw` to at least 92%.
2. No variance brings the `rare` cells to a mean of 92%.
3. `cv3` is wider than `cr2`, and `cr2` wider than `cr1s`, on all 11 real items.
4. `wcr` coverage is within 2 points of `cv3`'s on at least nine of the 11 real items.
5. The `cr1s` ratio of mean standard error to empirical standard deviation is below 1 on at least nine of the 11 real items.
6. At least one variance qualifies as a candidate.

## Controls

Estimates and `cr1s` standard errors must reproduce the study's saved oracle rows to within 1e-8. The fast `wcr`, `cr2` and `cv3` formulas must match brute-force refits and explicit matrices on small datasets, and `cr2` must reduce to the HC2 form with one row per cluster. `clubSandwich` and `sandwich` are not installed here, so `cr2` is not checked against them.

## Sources

Pustejovsky JE, Tipton E. *Journal of Business and Economic Statistics* 2018. [doi:10.1080/07350015.2016.1247004](https://doi.org/10.1080/07350015.2016.1247004). MacKinnon JG, Nielsen MØ, Webb MD. *Journal of Econometrics* 2023. [doi:10.1016/j.jeconom.2022.04.001](https://doi.org/10.1016/j.jeconom.2022.04.001).
