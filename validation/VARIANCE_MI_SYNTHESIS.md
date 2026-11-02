# Test of corrected variances through imputation: results and interpretation

Under the [rule fixed in advance](VARIANCE_MI_PLAN.md) the outcome is **CV3**. The delete-one-cluster jackknife met every condition, narrowly, and CR2 failed. As the plan requires, `cv3` becomes the default for `catmech_pool_lm(cluster = )` in version 0.7.0, with the earlier variance kept as `cluster_variance = "cr1s"`. Calls without `cluster` are unchanged.

The test also found something it was not designed to decide. **When the audited value is a predictor in the analysis, `catmech_reference(cluster = )` biased the slope toward zero**, by up to 16% in simulations. The rule did not turn on this, but it matters more than the choice of variance. `catmech_pool_lm()` now warns about that use, and the README's example no longer shows it. `catmech_reference(cluster = )` stays experimental, and earlier verdicts stand.

Everything labeled post hoc below came after the results and enters no rule. There were **148,500 requests: 148,486 returned intervals and 14 refusals**. All came at 25%: 12 on `iron`, in six cohorts, where the audited rows had no full-rank design, and two on `stillbirths`, in one cohort with no residual variance. The runner labeled one of those an error, and the plan's amendment explains why it is a refusal. The study took about 11 minutes. [Scored output](variance-mi-decisions.txt) and [post-hoc numbers](variance-mi-posthoc.txt) come from separate scripts.

## Decision rule

Conditions use `random_intercept` unless stated. The 11 real items have 2,000 cohorts at 25% and 500 at 50%.

| Condition | `cr2` | `cv3` | Needed |
|---|---|---|---|
| R1: real items at 25%, at least 92%, at most 98% | 10 and 10 of 11 | 9 and 11 of 11 | 9 and 10 |
| R2: real items at 50%, at least 92%, at most 98% | 10 and 10 of 11 | 10 and 11 of 11 | 8 and 10 |
| R3: mean distance from 95% at 25%, against 3.18 now | 1.49 | 1.47 | At least 0.5 lower |
| S1: simulated cells within limits | **23 of 44** | 36 of 44 | 36 |
| L1: linear MI at most 98%, real items and cells | 10 of 11, 41 of 44 | 11 of 11, 44 of 44 | 10, 40 |
| W1: largest width increase at 25% | **+95%** | +15% | At most +30% |
| Eligible | **No** | Yes | |

**This was a narrow pass.** `cv3` met R1 exactly, at nine of 11 with `iron` at 91.9%, and met S1 exactly, at 36 of 44. Coverage near 92% with 2,000 cohorts has a Monte Carlo standard error of 0.61 points, and with 500 cohorts 1.21. A little noise the other way would have kept the old default. Against that, the complete-data evidence for leaving the old variance is strong: under it, the recount reached 92% on only seven of 11 items.

## What changed with the variance

Coverage in percent, at 25%, on 2,000 cohorts. The current variance here is the one that was the default before this change.

| Item | Recount, current | Recount, `cv3` | Random intercept, current | Random intercept, `cv3` | Linear MI, current | Linear MI, `cv3` |
|---|---|---|---|---|---|---|
| `lbw` | 86.5 | 92.0 | 92.5 | 94.0 | 90.3 | 92.5 |
| `resus` | 94.3 | 96.7 | 96.2 | 96.7 | 95.6 | 96.0 |
| `stillbirths` | 94.3 | 96.8 | 93.7 | 95.1 | 93.1 | 94.3 |
| `deliveries` | 90.1 | 92.6 | 90.8 | 93.6 | 87.6 | 91.5 |
| `anc1_early` | 93.9 | 95.5 | 91.6 | 92.8 | 79.7 | 81.6 |
| `pnc` | 94.5 | 97.0 | 93.1 | 94.4 | 84.4 | 86.1 |
| `syphilis` | 93.3 | 96.0 | 87.5 | 90.0 | 81.7 | 84.2 |
| `iron` | 91.8 | 94.3 | 90.0 | 91.9 | 79.0 | 81.0 |
| `anc4` | 95.8 | 97.5 | 93.3 | 94.8 | 89.5 | 91.8 |
| `ancnewreg` | 93.3 | 95.7 | 91.8 | 94.2 | 85.5 | 88.8 |
| `livebirths` | 91.3 | 94.5 | 92.0 | 94.8 | 89.8 | 92.8 |

- **Real items came closer to 95%, at modest cost.** With the random intercept, mean distance went from 3.18 to 1.47 points, and widths rose 1% to 15%.
- **`cv3` does not rescue linear imputation.** Without `cluster`, mean distance goes from 8.18 to 6.04 points, and coverage is still 81.0% to 86.1% on the worst items. What remains there is an imputation shortfall.
- **Simulations pay for it in over-coverage.** With the random intercept, clustered cells average 95.3% under the old variance and 96.6% under `cv3`, and cells within limits fall from 40 of 44 to 36. Together with the real-item gains, this suggests the option's between-imputation variance can be conservative when simulated residuals fit its model well. Post hoc, and not tested.
- **The earlier PARTIAL verdict is not changed.** Point estimates do not depend on the variance, so the accuracy conditions fail just as before.
- **More cohorts moved earlier numbers.** The recount on `lbw` is 86.5% with 2,000 cohorts, against 88.4% with the first 500. That is a reminder of how noisy the 500-cohort rules were.

## Why CR2 failed

CR2 itself did well on the recount, so the failure lies in combining it with imputation. Its small Satterthwaite degrees of freedom enter Barnard-Rubin as the complete-data value, following the one source found on this, a blog post with no simulation. Pooled degrees of freedom then ran as low as 4.3 on average, and 19.01% of `random_intercept` requests fell under 5, against 0.22% for `cv3`. Widths rose by up to 95% with the random intercept and 195% with linear MI. This tests that one recipe. It does not show that CR2 cannot work through imputation.

## The exploratory wild bootstrap

The reference distribution invented for this test did poorly. With the random intercept at 25%, mean distance from 95% was 5.18 points, with `iron` at 78.2% and `resus` at 88.4%, although it does as well as the others on the recount. Its numerator ignores imputation error, which probably makes the reference too narrow. It should not be used, and no package option was added.

## When the audited value is a predictor

This is post hoc as a reading. The six cells were planned, but no rule or prediction addressed bias. The analysis is `Y ~ T + x + stratum`, with T imputed and a target of 0.5.

| Errors | Audit | Linear MI, bias % | Random intercept, bias % | Linear MI, coverage with `cv3` | Random intercept, coverage with `cv3` | RMSE ratio |
|---|---|---|---|---|---|---|
| `shift` | 25% | -3.0 | -16.1 | 95.4 | 84.4 | 1.654 |
| `shift` | 50% | -1.1 | -10.1 | 95.4 | 86.0 | 1.495 |
| `under` | 25% | -2.1 | -7.8 | 94.4 | 92.4 | 1.211 |
| `under` | 50% | -0.4 | -3.8 | 96.0 | 95.4 | 1.111 |
| `row_only` | 25% | -0.0 | -2.6 | 96.2 | 96.8 | 1.016 |
| `row_only` | 50% | -0.5 | -1.2 | 97.2 | 97.4 | 0.995 |

Worse with errors that shift whole clusters, as expected if the cause is what it appears to be: an unaudited cluster gets a new intercept drawn independently of its outcomes, which acts like added noise in a predictor and attenuates the slope. That explanation has not been tested. Every earlier real-data test of this option had the audited value as the outcome, so none could show this, and the README's example used the option as a predictor until now. This rests on six simulated cells in one design, and a real audit may behave differently.

**Advice now:** use `catmech_reference(cluster = )` only when the audited value is the analysis outcome. When it is a predictor, fit without `cluster` and pool with `catmech_pool_lm(cluster = )`.

## Predictions, scored

| # | Prediction | Result |
|---|---|---|
| 1 | `cr2` is eligible | **Failed**, on S1 and W1 |
| 2 | `cv3` fails W1 on at least one item | **Failed**, largest +15% |
| 3 | `syphilis` below 92% under `cr2` at 25% | **Held**, 91.9% |
| 4 | `resus` above 98% under `cv3` at 25% | **Failed**, 96.7% |
| 5 | `cr2` gains less through imputation than on the recount | **Failed**: +3.11 against +2.06 points |
| 6 | Pooled df under `cr2` below 15 on all 11 items | **Failed**, 4.3 to 19.1 |

One of six held. Most of the forecasts favored the wrong variance, which is some evidence that the rule was not fitted to the expected answer, though not proof, since it was written after the recount results.

## What this changes

- **`catmech_pool_lm(cluster = )` uses `cv3` by default.** Intervals with `cluster` get somewhat wider. Use `cluster_variance = "cr1s"` to reproduce earlier results.
- **`catmech_pool_lm()` warns** when a fit made with `cluster` has the audited value as a predictor. This changes no numbers.
- **Linear imputation without `cluster` still undercovers** on the items where it did before.
- **Rare counts remain unsolved.** On simulated rare counts, the recount averages 87.5% under `cv3`.

## Limits

The real items are held out by method, not by data, and it is one health system. The same assistant wrote the plan, runner, scoring and this write-up, and no reviewer, AI or human, has read this round. So this is not a neutral comparison. The rule's conditions were written after the recount results, and `cv3` passed two of them exactly. It has not been compared with `clubSandwich`, which is not installed here, though `cv3` matches brute-force refits. The study's control found a small runner fault, described with a status label in the plan's amendments: degrees of freedom differ by under 0.01 in a few requests, with no coverage flips. Nothing here uses the Zambia file, which stays available for a later test.

On process: the plan was committed before the runner existed. An engineering pilot was read for statuses and timing only. Hashes were checked before and after the run. As controls, 90,000 requests reproduce the replacement study under the old variance, and 24,500 reproduce the earlier scoring of the recount. After scoring, `cv3` moved into the package, and [the replay](variance-mi-replay/README.md) reproduces the study through `catmech_pool_lm()`. Commit order in a local repository is not third-party proof of timing.

## Sources

See the [plan](VARIANCE_MI_PLAN.md) for references.
