# Prospective test: corrected cluster-robust variances through imputation

Written 2026-11-02 and committed before the runner exists, so version history records the order. Amendments, if any, go in the section at the end.

## Why

[With every row recounted](VARIANCE_CHECK_RESULTS.md), the variance in `catmech_pool_lm(cluster = )` undercovers, and CR2, CV3 and a restricted wild cluster bootstrap repair most of that. No source found tests these through multiple imputation; the only one is a blog post with no simulation. This test decides whether `catmech_pool_lm(cluster = )` should change its default variance. The risk it looks for is over-coverage: an imputation variance that is conservative could combine with a corrected analysis variance to give intervals that are too wide.

## What has already been seen

- Coverage of `linear_mi` and `random_intercept` under the current variance, on the 500 cohorts of every existing cell.
- The recount under all four variances on those cohorts. This shaped the rule and predictions.
- The Zambia file's structure and agreement counts. It is not used here and stays available for a later test.

Nobody has computed a corrected variance on imputed data for these cohorts. The real items are held out by method, not by data.

## Imputers and variances

The imputers are `oracle` (reference on every row), `linear_mi` (`catmech_reference(use_recorded = TRUE)`) and `random_intercept` (`catmech_reference(use_recorded = TRUE, cluster = )`), with 20 imputations. Point estimates do not depend on the variance, so accuracy results are unchanged by anything here.

Each completed dataset gets these variances, for coefficient k, pooled by Rubin's rules with Barnard-Rubin degrees of freedom:

- **`cr1s`, the current one,** with complete-data degrees of freedom G - 1.
- **`cr2`,** with complete-data degrees of freedom equal to the mean over imputations of the Bell-McCaffrey Satterthwaite value. These are constant when the imputed variable is the outcome.
- **`cv3`,** with complete-data degrees of freedom G - 1.
- **`wcr`, exploratory.** A restricted wild cluster bootstrap reference distribution for the pooled statistic. In each replication, the same Rademacher weights apply to every imputation's restricted residuals, and the reference statistic is the mean bootstrap numerator over the root of the mean bootstrap `cr1s` variance. The observed statistic is the pooled estimate minus the target, over the root of the Rubin total under `cr1s`. It uses 999 replications, gives coverage and no width, and has no source behind it. It runs only when the imputed variable is the outcome, and it cannot be adopted here.

## Data

- **Real items.** The [replacement study's](CLUSTERED_FIX_PLAN.md) 11 Rwanda items at 25% and 50%, with its cohort, audit and imputation seed formulas. At 25%, replicates run to 2,000 and not 500, which brings the Monte Carlo error near 90% coverage to about 0.7 points. The first 500 are the existing cohorts. At 50% there are 500.
- **Simulated cells, imputed outcome.** That study's 38 cells, 500 cohorts each, with its seeds.
- **Simulated cells, imputed covariate, new.** The clustered settings above have the audited value as the outcome, but the package's main example has it as a predictor. So there are six more cells, in the `rwanda` design with `shift`, `under` and `row_only` errors at 25% and 50%. Given the study's `x`, `T` and strata, the analysis outcome is `Y = 1 + 0.5 T + 0.3 x + b_h + c_i + d_ij`, with `b_h` from N(0, 1), `c_i` from N(0, 5) and `d_ij` from N(0, 20). The analysis is `Y ~ T + x + stratum`, the target is 0.5, and the imputation predictors are `x`, the stratum and `Y`. Seeds: cohort 1,700,100,000 + cell index x 100,000 + replicate, imputation 1,710,000,000 + cell index x 100,000 + replicate, wild bootstrap 1,720,000,000 plus the same.

That is 44 simulated cells. Wild bootstrap seeds for the other cells are 1,730,000,000 + cell index x 100,000 + replicate, with real items as cells 101 to 122.

## Measures

Coverage of the target, mean and median width, complete-data and pooled degrees of freedom, and mean standard error over the empirical standard deviation of the estimate.

## Decision rule

Each of `cr2` and `cv3` is judged against the current variance. It is **eligible** if all of these hold, using `random_intercept` unless stated:

- **R1, real items at 25%.** At least 92% on at least nine of 11 items, and at most 98% on at least 10 of 11.
- **R2, real items at 50%.** At least 92% on at least eight of 11, and at most 98% on at least 10 of 11.
- **R3, closer.** Mean distance from 95% across the 11 items at 25% is at least 0.5 points below the current variance's.
- **S1, simulations.** In at least 80% of the 44 cells, coverage is at most 98% and is at least 92% or at least the oracle's under the same variance.
- **L1, no over-coverage with linear imputation.** With `linear_mi`, at most 98% on at least 10 of 11 items at 25%, and in at least 90% of the 44 cells.
- **W1, width.** With `random_intercept` at 25%, mean width is at most 30% above the current variance's on every one of the 11 items.

**What follows:**

- If only one is eligible, it becomes the default for `catmech_pool_lm(cluster = )`, with the current variance kept as an option, as version 0.7.0.
- If both are eligible, `cr2` becomes the default, unless `cv3`'s mean distance from 95% across the 11 items at 25% is more than 0.5 points lower.
- If neither is eligible, the default stays, and the documentation states the complete-data finding and this result.

Calls without `cluster` stay numerically identical in every case. No outcome changes an earlier verdict. `catmech_reference(cluster = )` stays experimental whatever happens here, because its accuracy conditions do not depend on the variance. Results are reported as they fall, with no change of thresholds after viewing.

The imputed-covariate cells are few and new. If a variance is eligible overall but falls under 92% in most of them, the write-up says so plainly.

## Predictions fixed in advance

1. `cr2` is eligible.
2. `cv3` fails W1 on at least one item.
3. With `random_intercept` at 25%, `syphilis` stays below 92% under `cr2`.
4. `resus` exceeds 98% under `cv3` with `random_intercept` at 25%.
5. Corrected coverage rises less through imputation than on the recount: the mean gain of `cr2` over the current variance across the 11 items at 25% is smaller with `random_intercept` than with `oracle`.
6. Pooled degrees of freedom under `cr2` are below 15 on average on all 11 items at 25%.

## Process

An engineering pilot, with separate seeds and a few cohorts per cell, is read for statuses and timing only. Any change it forces is listed below before the study runs. Hashes of the package, runner, scorer, plan and pool are checked before and after. As controls, on replicates 1 to 500, `linear_mi` and `random_intercept` under the current variance must reproduce the replacement study's estimates, standard errors and degrees of freedom, and the corrected recount must reproduce the earlier scoring. If a variance is adopted, it moves into the package, and a replay must reproduce the study's intervals through `catmech_pool_lm()`.

No independent reviewer is planned for this round, and no human has reviewed it. The same assistant wrote the plan, runner and scoring, so this is not a neutral comparison.

## Sources

Pustejovsky JE, Tipton E. *Journal of Business and Economic Statistics* 2018. [doi:10.1080/07350015.2016.1247004](https://doi.org/10.1080/07350015.2016.1247004). MacKinnon JG, Nielsen MØ, Webb MD. *Journal of Econometrics* 2023. [doi:10.1016/j.jeconom.2022.04.001](https://doi.org/10.1016/j.jeconom.2022.04.001). Barnard J, Rubin DB. Small-sample degrees of freedom with multiple imputation. *Biometrika* 1999. [doi:10.1093/biomet/86.4.948](https://doi.org/10.1093/biomet/86.4.948). Pustejovsky JE. [Blog post on multiple imputation with `clubSandwich`](https://jepusto.netlify.app/posts/MI-with-clubSandwich/index.html), which proposes the mean Satterthwaite value as complete-data degrees of freedom and reports no simulation.

## Amendments

None yet.
