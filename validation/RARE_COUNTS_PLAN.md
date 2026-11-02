# Prospective test: what is going wrong with rare counts

Written 2026-11-02 and committed before the runner exists, so version history records the order. Amendments, if any, go in the section at the end.

## Why

The documents say rare counts are unsolved. On six simulated rare-count cells, even a full recount averages 87.5% under `cv3`, so no imputation can reach 95% there by being more accurate. But the real rare items do not show this. Under `cv3` with the option at a 25% audit, `lbw`, `resus` and `stillbirths` cover 94.0%, 96.7% and 95.1%. The simulated cells have a cluster effect drawn from a gamma distribution with shape 0.5, which is very heavy-tailed. So the problem may be heavy-tailed cluster effects and not rarity. This test separates the two, tries other analyses and count-aware imputation, and tests a warning. "No fix" is an accepted outcome.

**What has been seen.** Coverage of the recount and of `linear_mi` and `random_intercept` on the six rare cells and three rare real items. Nothing has been computed on the new cells, candidates or diagnostic.

## Part 1: the analysis, with every row recounted

No imputation. Cells cross three factors, with three rows per cluster and the existing `rwanda`, `pairs` and `none` designs and the existing covariate: **frailty**, a cluster multiplier that is 1, or gamma with shape and rate 2, or gamma with shape and rate 0.5; and **level**, the existing rare mean, or ten times it. Counts are Poisson with mean `pmax(.01, .1 + .02 x + stratum level) x frailty x scale`, so the target slope is 0.02 x scale. That is 18 cells with 1,000 cohorts each. The analysis is `y ~ x + stratum`, with no strata in `none`.

Reference distributions for the slope:

- `cv3`, the current default, `cr2` with its Satterthwaite degrees of freedom, and the restricted wild cluster bootstrap `wcr`, with 999 replications, as in the earlier [recount check](VARIANCE_CHECK_PLAN.md).
- `boot_t`, new here: a pairs cluster bootstrap-t, resampling clusters within strata, 399 replications, studentized by the `cr1s` variance. It is a standard way to handle skewed statistics.
- `poisson`: a log-link Poisson regression with a cluster-robust variance and t(G - 1). **This answers a different question,** a log slope, whose target for each cell comes from one fit on 200,000 simulated clusters. The linear mean means that log model is misspecified here, so the target is a best-fitting value.

**The diagnostic.** D is the largest share of the `cv3` variance for the slope that comes from one cluster: the largest squared delete-one change over the sum of squared changes.

## Part 2: count-aware imputation

On the six rare cells and on `lbw`, `resus` and `stillbirths` at 25% and 50%, using the same cohorts, seeds and audits as the [corrected-variance study](VARIANCE_MI_PLAN.md), with its saved `oracle`, `linear_mi` and `random_intercept` results reused. Candidates, each with 20 imputations, pooled with `cv3`:

- `nb_mi`: a negative binomial regression of the reference value on `log1p(recorded)` and the predictors, fitted on audited rows by `MASS::glm.nb`. Each imputation draws coefficients from their normal approximation and log theta from its own, then draws counts. It falls back to Poisson when the fit fails. It assumes no cluster effect.
- `pmm_recorded`: `mice` predictive mean matching given the recorded value and predictors, with `eps = 0` so that the recorded value stays. It keeps values on the observed support.

These are established tools and not new methods. Requests have a 60-second limit, and one over it counts as no interval.

## Part 3: false alarms for the diagnostic

D on the recount, 500 cohorts each, for the 32 existing simulated cells that are not rare and the 11 real items at 25%.

## Decision rules

**P1, what the shortfall is.** Under `cv3` on the recount, let M0 be the mean coverage over the three rare cells with no frailty and M5 over the three rare cells with shape 0.5.

- **HEAVY TAILS** if M0 is at least 93 and M5 at most 90.
- **RARITY** if M0 is below 92.
- **MIXED** otherwise.

The documents then state the finding in place of "rare counts are unsolved".

**P2, a better analysis.** Each of `cr2`, `wcr` and `boot_t` is **eligible** if its mean over the three rare shape 0.5 cells is at least 92, it is at most 98 in all 18 cells, and its mean over the other 15 cells is no more than 0.5 points below `cv3`'s. `poisson` is judged the same way but reported as a different question. No package code changes under this rule, because these do not yet pass through multiple imputation. An eligible method is named in the documents for complete-data analysis.

**P3, count-aware imputation.** On the set of six rare cells and three real items at both audits, 12 in all, a candidate is **eligible** if its mean distance from 95% is at least 0.5 points below `linear_mi`'s, the geometric mean of its RMSE over `linear_mi`'s is at most 1.00, it is at most 98% on at least 10 of 12, and it gives intervals on at least 99% of requests. An eligible candidate is named in the documents, and any package option needs its own later plan.

**P4, a warning.** Flag a request when D is above 0.5, meaning one cluster supplies more than half of the variance. The warning is **ADOPTED** if, pooled over the 18 cells in Part 1, coverage under `cv3` in flagged requests is at least 5 points below unflagged, at least 25% of requests in the rare shape 0.5 cells are flagged, and at most 5% are flagged across Part 3. Then `catmech_pool_lm(cluster = )` warns when any coefficient has D above 0.5, computed on the first completed dataset, naming the coefficients. Numbers do not change. If it is not adopted, nothing is added.

Results are reported as they fall, with no change of thresholds after viewing.

## Predictions fixed in advance

1. P1 is HEAVY TAILS.
2. No method in P2 is eligible, with `boot_t` the closest.
3. Neither candidate in P3 is eligible.
4. The warning in P4 is adopted.

## Limits known in advance

These cells use one covariate pattern and three rows per cluster. The gamma frailty is an assumption, and shape 0.5 may be more extreme than real audits, which is partly what this tests. The real rare items are only three, all from one recount, and all have been seen. The threshold for D is a judgment call, set before seeing D. The same assistant wrote the plan, runner and scoring, and no reviewer has read it.

## Seeds and process

Part 1 cohorts use 1,210,000,000 + cell x 100,000 + replicate, with bootstrap seeds 1,220,000,000 and wild seeds 1,230,000,000 plus the same, and targets for the log slope 1,235,000,000 + cell. Part 2 imputation uses 1,240,000,000 for `nb_mi` and 1,250,000,000 for `pmm_recorded`, plus cell x 100,000 + replicate, with real items numbered from 101. An engineering pilot on separate seeds is read for statuses and timing only. Any change it forces is listed below before the study. Hashes are checked before and after.

## Amendments
