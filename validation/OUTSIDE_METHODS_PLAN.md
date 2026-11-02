# Prospective test: the clustered option against outside methods

Written 2026-11-02 and committed before the runner exists, so version history records the order. Amendments, if any, go in the section at the end.

## Why

Every test of `catmech_reference(cluster = )` so far compares it with the oracle and with catmech's own variants. So nothing yet shows the option has a reason to exist. A literature search made for this package, kept outside the repository, ranked this comparison first. This test decides whether the option stays, and the burden is on the option: if an outside method does about as well, the option is retired.

## What has already been seen

- Coverage, width and accuracy of `oracle`, `linear_mi` and `random_intercept` on these exact cohorts, in the [corrected-variance study](VARIANCE_MI_SYNTHESIS.md). Those saved results are reused here and not rerun.
- Earlier comparisons of recorded-only, audited-only and a home-built raking with linear MI, on [other cohorts and a different audit unit](REAL_AUDIT_SYNTHESIS.md). None involved the random intercept.
- Timing and edge cases of `mice` and `survey` on synthetic fixtures with seeds used nowhere else. Those runs computed no coverage or accuracy.

Nobody has run any outside multilevel imputation, `survey` raking or a model-free difference estimator on these cohorts.

## The option being judged

`random_intercept` is `catmech_reference(use_recorded = TRUE, cluster = )` with 20 imputations, pooled by `catmech_pool_lm(cluster = )` under `cv3`, which is the package's current default.

## Methods compared

All use the analysis model of the cell. Tuning after the pilot is not allowed. Software options stay at their documented defaults except where stated.

- **`recorded`, a baseline.** Ordinary least squares on the recorded value in all rows, with `cv3` and t(G - 1). It uses no audit, and it cannot retire the option.
- **`audited_only`.** Unweighted least squares on the audited rows with the reference value, with `cv3` and t(G_a - 1), where G_a is the number of audited clusters.
- **`difference`, new and model-free.** Write theta(value, rows) for the analysis coefficient. The estimate is theta(recorded, all rows, unweighted) + theta(reference, audited rows, weighted) - theta(recorded, audited rows, weighted), with weights d = N_h / m_h, the clusters in the stratum over the audited clusters in it. Row scores are psi_i = (B x_i)_k e_i, where B is the inverse of the weighted cross-product and e_i the residual of that fit. Sums by cluster give psi^R_g from the reference fit and delta_g = psi^R_g - psi^W_g from the two audited fits. The variance is V1 + V2:
  - V1 = (G / (G - 1)) x the sum over audited clusters of d_g (psi^R_g)^2.
  - V2 = the sum over strata of (1 - f_h) (m_h / (m_h - 1)) x the sum over audited clusters in the stratum of (d_g delta_g minus its stratum mean)^2. A stratum with one audited cluster instead contributes (1 - f_h) (d_g delta_g)^2.

  The interval uses t with max(1, G_a - H) degrees of freedom, H being the number of strata. It needs `V1` and `V2` to match the phase components from `survey` on a fixed dataset before the study runs.
- **`raking`.** `survey::twophase` with clusters at both phases, the strata at phase two, and audited clusters as the subset, then `survey::calibrate(phase = 2, calfun = "raking", formula = ~ if_k)`, then `survey::svyglm()` with the analysis formula. Here `if_k` is the row influence function for the analysis coefficient from the recorded fit on all rows. It is the single auxiliary, because the literature warns against many constraints and because an earlier pilot here failed often with one per coefficient. A stratum with one audited cluster uses `survey.lonely.psu = "adjust"`. The interval uses t with `survey::degf()` degrees of freedom. A failed calibration counts as no interval.
- **`mice_2l_lmer` and `mice_2l_pan`.** `mice::mice()` with `2l.lmer` and `2l.pan`, 20 imputations and `maxit = 1`, for the reference value given the recorded value and the cell's predictors, with a random intercept for the cluster. That is the same model as the option. Each completed dataset goes through the same analysis as the option, `cv3` pooled by `mice::pool.scalar(n = G, k = 1)`. `mice` removes a predictor it finds collinear, by default, and the number of those logged events is recorded. When `lmer` cannot fit, `mice` warns "lmer does not run" and returns unmodeled draws, so any request with that warning counts as no interval.

**Exploratory, accuracy only, and in no rule:** `fitted_value`, least squares on the fitted values from a design-weighted regression of the reference value on the recorded value and predictors. The earlier "design-weighted difference estimator" in the [real-audit study](REAL_AUDIT_SYNTHESIS.md) reduces to this. Its correction term is a weighted residual from a regression containing the same predictors, so it is zero by the normal equations. So that study's label overstated what it tested, and that estimator was never a difference estimator in the usual sense. The write-up will say so.

## Data

The same cohorts as the [corrected-variance study](VARIANCE_MI_PLAN.md), from the same seed formulas:

- **Real items.** 11 Rwanda items at 25% with 2,000 cohorts and at 50% with 500.
- **Simulated cells, audited value as outcome.** 38 cells, 500 cohorts each.
- **Simulated cells, audited value as predictor.** Six cells, 500 cohorts each. `catmech_pool_lm()` warns against the option in this use, so these enter no rule.

That is 49,500 requests for each of six methods, 297,000 in all. Imputation seeds are `mice_2l_lmer` 1,950,000,000 and `mice_2l_pan` 1,970,000,000, plus cell number x 100,000 plus replicate, with simulated cells numbered by index, real items 101 to 122, and predictor cells 201 to 206. The other methods draw nothing after the cohort.

## Measures

Point estimate, standard error, degrees of freedom, interval width, coverage of the target, RMSE against the target, bias, availability and seconds per request. Accuracy ratios use cohorts where both methods gave an interval. Speed is reported and enters no rule.

## Decision rule

The rule applies to each of `audited_only`, `difference`, `raking`, `mice_2l_lmer` and `mice_2l_pan`, called C. The comparison is always with `random_intercept`. C is **no worse** if all of these hold:

- **N1, accuracy on real items at 25%.** The geometric mean over 11 items of RMSE(C) / RMSE(option) is at most 1.05.
- **N2, accuracy on real items at 50%.** The same, at most 1.05.
- **N3, coverage on real items at 25%.** At least 92% on at least as many items as the option, at most 98% on at least 10 of 11, and mean distance from 95% at most 0.25 points above the option's.
- **N4, coverage on real items at 50%.** At least 92% on at least the option's count minus one, at most 98% on at least 10 of 11, and mean distance from 95% at most 0.50 points above the option's.
- **N5, simulations.** Over the 36 audited-outcome cells other than the two `control_one_row` cells, the geometric mean RMSE ratio is at most 1.05. The number of cells at most 98% and either at least 92% or at least the oracle's under `cv3` is at least the option's count minus two.
- **N6, availability.** An interval on at least 99% of requests in each of real items at 25%, real items at 50%, and those 36 simulated cells.

`control_one_row` has one row per cluster, which does not identify a random intercept, so the two cells are reported and excluded from N5 and N6 for every method.

The margins favor the comparators on purpose. A new experimental option has to do more than tie maintained software or simpler estimators.

**What follows:**

- **RETIRE** if any C is no worse. `catmech_reference(cluster = )` then gets a deprecation warning in version 0.8.0, naming the passing method with the smallest N1 ratio, with ties broken by mean distance in N3. Documentation changes to match. The code stays, so that earlier results can still be reproduced.
- **KEEP** otherwise. The option stays experimental, since the earlier PARTIAL verdict stands. The documentation reports this comparison.

Results are reported as they fall, with no change of thresholds after viewing. No outcome changes an earlier verdict.

**Description only, no rule:**

- How often `recorded` has lower RMSE than every audit-using method. If it is often, a partial audit can cost more than it gives.
- In the six predictor cells, how each method compares with `linear_mi` under `cv3`, which is the current advice.
- Rare counts, by cell.
- Seconds per request.

## Predictions fixed in advance

1. The outcome is RETIRE, through `mice_2l_lmer`.
2. `raking` fails N1, with a ratio above 1.05.
3. `difference` reaches 92% on at least as many real items at 25% as the option, and fails N1.
4. `audited_only` has an N1 ratio above 1.20.
5. The mean distance from 95% of `mice_2l_pan` on real items at 25% is within 0.50 points of `mice_2l_lmer`'s.
6. `recorded` has lower RMSE than every audit-using method on at least three of 11 real items at 25%.

## Process

An engineering pilot, with separate seeds and a few cohorts per cell, is read for statuses and timing only. Any change it forces is listed below before the study runs. Hashes of the package, runner, scorer, plan, pool and reused results are checked before and after.

As checks:

- On a fixed dataset, `difference` must match `survey`'s point estimate through `svyglm` and its phase one and phase two variance components, and `twophase` weights must equal d = N_h / m_h.
- `fitted_value` must have a correction term of zero, to numerical accuracy.
- As a control, on the first 10 cohorts of every cell, `random_intercept` run afresh through the package must reproduce the saved estimate, standard error and degrees of freedom, and each saved target must match. Those runs also give the option's seconds per request.

The results reused from the earlier study were computed before the variance moved into the package and reproduced through it by [its replay](variance-mi-replay/README.md).

No independent reviewer is planned for this round, and no human has reviewed it. The same assistant built the option, wrote this plan and the scoring, and chose how each comparator is set up. So this is not a neutral comparison, and that matters most here, where the comparators' recipes decide the result. Each is a reasonable default use. None is the best possible use of that method.

## Sources

The [`survey` documentation for `twophase`](https://r-survey.r-forge.r-project.org/survey/html/twophase.html). Breslow NE, Lumley T, Ballantyne CM, Chambless LE, Kulich M. Improved Horvitz-Thompson estimation of model parameters from two-phase stratified samples: applications in epidemiology. *Statistics in Biosciences* 2009. [doi:10.1007/s12561-009-9001-6](https://doi.org/10.1007/s12561-009-9001-6). Han K, Shaw PA, Lumley T. Combining multiple imputation with raking of weights: an efficient and robust approach in the setting of nearly true models. *Statistics in Medicine* 2021. [doi:10.1002/sim.9210](https://doi.org/10.1002/sim.9210). van Buuren S, Groothuis-Oudshoorn K. mice: multivariate imputation by chained equations in R. *Journal of Statistical Software* 2011. [doi:10.18637/jss.v045.i03](https://doi.org/10.18637/jss.v045.i03). Audigier V, White IR, Jolani S, and others. Multiple imputation for multilevel data with continuous and binary variables. *Statistical Science* 2018. [doi:10.1214/18-sts646](https://doi.org/10.1214/18-sts646).

## Amendments
