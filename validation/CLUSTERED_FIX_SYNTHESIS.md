# Test of replacements for the clustered resampling: results and interpretation

Under the [rule fixed in advance](CLUSTERED_FIX_PLAN.md) the verdict is **PARTIAL**. Simulations selected random-intercept imputation over a stratified jackknife draw. Its coverage conditions and the check on unused items held, and both accuracy conditions failed. As the plan requires, it replaces the resampling in `catmech_reference(cluster = )`, `cluster_strata` is removed, and the option stays experimental. It is not recommended as a default.

The verdict rests on the plan alone. Everything labeled post hoc below came after the results and enters no rule.

What this round achieved:

- **The flaw is gone from the code.** Nothing is resampled, so the within-stratum shortfall cannot arise. Strata now enter as predictors.
- **Coverage rose where it was worst, but barely moved on average.** On `anc1_early`, the item the first option missed, it went from 88.8% to 91.0%, against 79.6% for linear imputation without `cluster`. `syphilis` fell from 93.0% to 89.0%, the one item outside the limits. Post hoc, mean distance from 95% across nine items went from 3.25 to 3.07 points, against 7.92 without `cluster` and 2.40 for a recount of every row.
- **Error improved.** Against linear imputation without `cluster`, the recovery RMSE ratio was 0.748 to 0.949 on the six common items, clearly below 1 on each. On the three rare items it was 1.039 to 1.054, where the first option's was 1.050 to 1.137. That remaining cost is why the accuracy condition failed.
- **Intervals got narrower and steadier.** They were narrower than the first option's on seven of nine items. Widths over five times their cell median fell from 24 to 5 of 29,998 requests, and the smallest degrees of freedom rose from 0.79 to 2.05.
- **The two items no method had touched did well.** On `ancnewreg`, coverage rose from 85.8% to 91.4% (oracle 92.4%), with a ratio of 0.774. On `livebirths`, it rose from 90.2% to 92.2% (oracle 91.8%), with a ratio of 0.921. They share facilities with the other items.

Higher coverage was available at a price. The jackknife, which the selection rule passed over, covered 91.4% to 97.6% on the nine items, but its intervals were up to 169% wider than those without `cluster`, its error showed no gain, and it had 119 very wide requests.

There were **210,000 requests: 209,988 returned intervals, 12 refusals and no other errors**. All 12 came on development iron at 25%, in cohorts 257 and 355, where the audited rows gave no full-rank design for any imputing method. The study took about 11 minutes at 14 workers. [Scored output](clustered-fix-decisions.txt) comes from a script frozen before the run, and post-hoc numbers come from [a separate table](clustered-fix-posthoc.txt).

## Selection

Simulations alone selected the candidate, and neither candidate refused a simulated request.

| Method | Mean distance from 95%, 24 clustered cells, points | Mean recovery RMSE ratio, 38 cells |
|---|---|---|
| `jackknife` | 1.31 | 1.013 |
| `random_intercept` | 1.41 | 0.994 |
| First option, within strata | 1.43 | 1.004 |
| `rao_wu` | 1.80 | 1.199 |
| First option, no strata | 1.92 | 1.010 |
| Existing option | 4.72 | |

The two candidates were within 0.5 points, so the lower mean ratio decided. This was a close call settled by the tie-break.

## Real items at a 25% audit

Coverage is clustered coverage of the pool slope, in percent, on 500 cohorts. Ratios are against linear imputation without `cluster` (the existing option), and a ratio below 1 favors the method. The first five items are the first test's, the next four are its development items, and the last two are the unused items.

| Item | Oracle | Existing | First option | Replacement | Ratio, first option | Ratio, replacement | Mean width: existing, first option, replacement |
|---|---|---|---|---|---|---|---|
| `lbw` | 88.4 | 92.0 | 89.8 | 93.6 | 1.137 | 1.054 [1.023, 1.087] | 0.055, 0.056, 0.060 |
| `resus` | 94.0 | 96.0 | 93.6 | 96.4 | 1.104 | 1.039 [1.000, 1.079] | 0.050, 0.059, 0.057 |
| `stillbirths` | 93.4 | 94.4 | 94.0 | 93.8 | 1.050 | 1.041 [1.004, 1.078] | 0.013, 0.014, 0.015 |
| `deliveries` | 91.6 | 87.0 | 91.4 | 90.4 | 0.903 | 0.902 [0.866, 0.936] | 0.234, 0.249, 0.240 |
| `anc1_early` | 95.6 | 79.6 | 88.8 | 91.0 | 0.850 | 0.748 [0.703, 0.797] | 0.257, 0.337, 0.289 |
| `pnc` | 95.2 | 85.4 | 93.6 | 92.8 | 0.894 | 0.945 [0.902, 0.984] | 0.237, 0.307, 0.272 |
| `syphilis` | 93.2 | 84.6 | 93.0 | 89.0 | 0.868 | 0.949 [0.910, 0.988] | 0.762, 0.847, 0.804 |
| `iron` | 90.2 | 79.5 | 90.8 | 91.0 | 0.855 | 0.808 [0.751, 0.868] | 0.480, 0.678, 0.518 |
| `anc4` | 96.6 | 87.2 | 90.8 | 92.2 | 0.970 | 0.948 [0.920, 0.979] | 0.313, 0.362, 0.347 |
| `ancnewreg` | 92.4 | 85.8 | 91.8 | 91.4 | 0.906 | 0.774 [0.726, 0.827] | 0.450, 0.547, 0.482 |
| `livebirths` | 91.8 | 90.2 | 94.2 | 92.2 | 0.928 | 0.921 [0.887, 0.957] | 0.261, 0.282, 0.274 |

## Decision rule

| Condition | Result | Needed | Holds? |
|---|---|---|---|
| Coverage: at least 90% or the oracle's, at most 98%, on the nine items at 25%, with `anc1_early` among them | 8 of 9 (`syphilis` missed) | 8 of 9 | Yes |
| Accuracy: upper limit of the recovery ratio at most 1.05 on the first test's five items | 2 of 5 | 4 of 5 | **No** |
| Unused items: within 5 points of the oracle, and that upper limit at most 1.05 | 2 of 2 | Both | Yes |
| Simulation coverage: at least 92% or the oracle's, at most 98%, in the clustered cells | 23 of 24 (`rwanda_shift` at 25%, 98.4%) | 80% | Yes |
| Simulation accuracy: that upper limit at most 1.05 across all cells | 22 of 38 | 80% | **No** |

**The first test's own coverage condition would not have held.** The replacement is closer to the oracle than the existing option on only three of five items, because on `lbw` and `resus` it rose further above an oracle that was already lower. Under that condition this would be NO-GO. The plan replaced it with fixed limits before any real result with these candidates, for the reason it gives, but readers should know the difference matters here. The first option scores 5 of 5 on the old condition and 8 of 9 on the new one.

## Reading the result

**Against the first option, the replacement is better on error and steadiness, and about even on coverage.** Post hoc, at 25% it had the lower recovery ratio on six of nine items and the higher on `pnc` and `syphilis`, with `deliveries` level. It was closer to 95% on four items and farther on four, with `resus` level. At a 50% audit the coverage gain is larger: mean distance from 95% across the nine items was 3.07 points, against 3.87 for the first option and 5.71 for the existing one.

**The jackknife got higher coverage from much wider intervals.** On `iron` its mean width was 1.292, against 0.518 for the replacement and 0.480 for the existing option. On `anc1_early` it was 0.510 against 0.289. Its recovery ratio on the nine items ran 0.985 to 1.057, so its error showed no gain. Some of its coverage above the oracle's may offset the cluster-robust variance's own undercoverage, since the oracle reached only 88.4% on `lbw`. It is not in the package, and its code stays in [the candidates file](../tools/clustered-fix-candidates.R).

**The textbook correction overshoots.** The rescaled bootstrap covered above 96% in seven of eight clustered cells of the stratified designs at 25%, up to 99.0%. Its recovery ratio reached 3.088 in one cell and averaged 1.199. On the nine real items it came closest to 95% (1.11 points), but with error up to 1.446 on `lbw`.

**How much the shortfall itself mattered depends on the evidence.** On real items, leaving out strata, rescaling and the jackknife all took `anc1_early` from 88.8% to between 93.8% and 94.4%, which supports the earlier reading that the shortfall drove much of that miss. In the simulations it mattered little: the first option's mean distance from 95% was 1.43 points, close to both candidates, and in the stratified designs at 25% it averaged 94.4% and 95.5%. Post hoc and untested, refits on few resampled clusters may add instability that offsets the shortfall at small audits. Its between-imputation variance, relative to its own error, was below the jackknife's in all 12 clustered cells at 50%, but in only 6 of 12 at 25%.

**The replacement does not fix `syphilis`.** It covered 89.0% at 25% and 86.8% at 50%, against an oracle's 93.2%, and the first option and the jackknife did better at both sizes. Every method but the oracle stays under 91% on that item at 50%.

**Both failed conditions use the upper limit of a ratio.** At 500 cohorts, that limit often sits about 0.04 above the point ratio, as in `none_under` at 25%: 1.016 [0.983, 1.051]. Post hoc, the point ratio was at most 1.05 in 35 of 38 simulated cells, against 22 on the upper limit that the rule uses. The first option's counts were 33 and 30, so by the rule's measure the replacement did worse on this condition than the option it replaces. The rule stands as written. The cost on rare counts is also real and not just noise: the lower limits on the three rare items are 1.000 or above, and `rwanda_rare` at 25% has a ratio of 1.090 [1.045, 1.134].

**These simulations do contain the clustering.** The first plan's did not. Here the existing option covers 86.0% to 92.6% across the clustered cells and is under 92% in 18 of 24, which is the shortfall the option is meant to address. Means by setting, in percent:

| Setting | Audit | Oracle | Existing | First option | No strata | `rao_wu` | Jackknife | Replacement |
|---|---|---|---|---|---|---|---|---|
| Clustered, `rwanda` design | 25% | 92.7 | 89.7 | 94.4 | 96.6 | 98.4 | 94.8 | 96.4 |
| Clustered, `rwanda` design | 50% | 90.7 | 90.3 | 92.1 | 93.1 | 93.3 | 92.5 | 94.4 |
| Clustered, `pairs` design | 25% | 91.4 | 88.1 | 95.5 | 98.1 | 96.5 | 94.9 | 95.4 |
| Clustered, `pairs` design | 50% | 92.6 | 92.1 | 94.0 | 93.7 | 95.3 | 94.5 | 95.3 |
| Clustered, no strata | 25% | 92.5 | 89.8 | 93.7 | 93.7 | 93.6 | 93.5 | 94.9 |
| Clustered, no strata | 50% | 92.6 | 91.7 | 93.3 | 93.3 | 93.5 | 93.2 | 95.3 |
| Rare counts | 25% | 83.1 | 85.5 | 86.1 | 88.6 | 90.6 | 87.6 | 87.4 |
| Rare counts | 50% | 83.3 | 82.7 | 82.7 | 82.5 | 83.3 | 82.8 | 83.1 |
| Row-only errors | 25% | 91.9 | 96.0 | 95.9 | 97.3 | 96.7 | 95.2 | 97.2 |
| Row-only errors | 50% | 93.9 | 94.9 | 94.5 | 95.2 | 95.4 | 94.5 | 96.0 |
| One-row controls | 25% | 92.8 | 94.4 | 93.6 | 93.6 | 93.4 | 94.0 | 94.4 |
| One-row controls | 50% | 93.8 | 95.0 | 95.0 | 95.0 | 94.8 | 94.6 | 95.0 |

On rare counts the oracle itself covers about 83%, because the cluster-robust variance undercovers on heavy-tailed counts. No imputation option repairs that, and `catmech_pool_lm(cluster = )` is unchanged.

## Predictions, scored

| # | Prediction | Result |
|---|---|---|
| 1 | The simulations select `random_intercept` | **Held**, by the tie-break |
| 2 | `rao_wu` above 96% in most clustered cells of the stratified designs at 25% | **Held**, 7 of 8 |
| 3 | Both candidates closer to 95% on average than the first option | **Held**, narrowly: 1.31, 1.41 and 1.43 |
| 4 | `anc1_early` at 25% at least 90% | **Held**, 91.0% |
| 5 | Rare items' ratios below the first option's on all three | **Held**: 1.054 against 1.137, 1.039 against 1.104, 1.041 against 1.050 |
| 6 | Common items' ratios at most 1.02 on all six | **Held**, 0.748 to 0.949 |
| 7 | One-row controls' mean width within 5% | **Held**, unchanged |
| 8 | Unused items' coverage within 3 points of the existing option | **Failed**: `ancnewreg` rose by 5.6 points, and `livebirths` by 2.0 |
| 9 | Mean width at most 40% wider on the five items | **Held**: +8%, +14%, +10%, +3%, +12% |
| 10 | Fewer very wide requests than the first option | **Held**, 5 against 24 |

The simulation predictions lean on the scratch run the plan discloses, so their holding is weak evidence. The plan also said the unused items could not show a gain. It was wrong: `ancnewreg` shows one.

## At a 50% audit

| Item | Oracle | Existing | First option | Replacement | Ratio, first option | Ratio, replacement |
|---|---|---|---|---|---|---|
| `lbw` | 88.4 | 88.0 | 87.0 | 88.8 | 0.999 | 0.998 [0.976, 1.022] |
| `resus` | 94.0 | 95.8 | 96.8 | 96.6 | 1.045 | 1.049 [1.011, 1.091] |
| `stillbirths` | 93.4 | 94.0 | 93.6 | 94.2 | 1.028 | 1.007 [0.976, 1.042] |
| `deliveries` | 91.6 | 88.8 | 91.0 | 91.2 | 0.964 | 0.953 [0.924, 0.979] |
| `anc1_early` | 95.6 | 88.0 | 92.4 | 94.2 | 0.976 | 0.877 [0.817, 0.934] |
| `pnc` | 95.2 | 86.0 | 91.4 | 92.8 | 0.939 | 0.909 [0.876, 0.945] |
| `syphilis` | 93.2 | 85.2 | 88.0 | 86.8 | 0.984 | 1.019 [0.995, 1.044] |
| `iron` | 90.2 | 84.6 | 89.6 | 92.8 | 0.964 | 0.957 [0.911, 1.002] |
| `anc4` | 96.6 | 95.2 | 96.0 | 96.8 | 0.978 | 0.936 [0.902, 0.969] |
| `ancnewreg` | 92.4 | 88.4 | 91.0 | 92.2 | 0.946 | 0.849 [0.799, 0.902] |
| `livebirths` | 91.8 | 89.2 | 91.0 | 91.6 | 0.973 | 0.994 [0.966, 1.021] |

## What this changes

- **`catmech_reference(cluster = )` is now random-intercept imputation.** Each imputation draws the variance components, then coefficients given them, then a cluster intercept and row noise. Results for calls with `cluster` differ from the earlier version. Calls without it are numerically unchanged.
- **`cluster_strata` is gone.** If the audit is stratified, put the strata in `predictors`.
- **It stays experimental.** Treat it as a check, not a default. Expect some added error on rare counts, and read the `df` column.
- **It assumes one normal intercept per cluster.** It does not model cluster-specific slopes or non-constant variance. The `under` setting gives such errors, and its coverage held there, but that is one setting.
- **`catmech_pool_lm(cluster = )` is unchanged,** with the limits in the [first write-up](CLUSTERED_SYNTHESIS.md).

## Limits

This test is held out by method, not by data. Both candidates were designed knowing the earlier results on these items, including the cost on rare items, and the scratch run the plan discloses shaped the choice of candidates. The same assistant chose the simulation settings, and four of the six error settings were used in that scratch run. The two unused items are new to every method but share facilities and months with the others, and the plan chose them because they were all that remained. It is one health system, and `anc1_early`'s two sides are defined slightly differently, as the first plan notes. With 500 cohorts, the Monte Carlo error on a paired coverage difference is about one point.

On process: the plan was committed before either candidate was in the package or had run on a real item, and no amendment was needed. An engineering pilot on separate seeds (1,260 requests, all returned) was read for statuses and timing only. Hashes of the sources, candidates, runner, scorer, plan and pool were checked before and after the run. As controls, 27,000 requests reproduce the first study, and 18,000 reproduce the review's reruns. After scoring, the selected candidate moved into the package, and all 30,000 of its study requests replay through the package to within 1e-10. Commit order in a local repository is not third-party proof of timing. This round was written and scored by an AI assistant. Unlike the first test, it has had no independent review, and no human has reviewed it.

## Sources

See the [plan](CLUSTERED_FIX_PLAN.md) for references.
