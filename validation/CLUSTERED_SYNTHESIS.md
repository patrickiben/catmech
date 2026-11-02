# Test of imputation for clustered audits: results and interpretation

The option tested here was later replaced, and this write-up records the first test. See the [second test](CLUSTERED_FIX_SYNTHESIS.md).

Under the [rule fixed in advance](CLUSTERED_PLAN.md) the verdict is **PARTIAL**. Resampling audited facilities in each imputation passed the coverage and simulation conditions, but it failed the accuracy condition on three of five held-out Rwanda items. As the plan requires, `catmech_reference(cluster = )` stays in the package as an experimental option with that failure stated. It is not recommended as a default.

The verdict rests on the plan alone. Everything labeled post hoc below came after the results and enters no rule.

What holds:

- **It helps where the existing option falls short of a recount.** On the two held-out items where `catmech_reference(use_recorded = TRUE)` covered 4.6 and 16.0 points below a recount, clustered coverage rose by 4.4 and 9.2 points, and error fell. The four development slopes moved the same way.
- **It costs where the existing option was already fine.** On the other three held-out items, coverage fell toward the oracle's from above and error rose by 14%, 10% and 5%, the last within Monte Carlo noise.
- **Post hoc, two features separate those cases.** The three items where it cost accuracy are rare events, whose residuals show no within-facility correlation. The six where it helped are common events, whose residuals do. Rarity may be a proxy, because within items the gain follows that correlation. Nine items from one health system cannot settle this.
- **Both the gain and the cost come mainly from resampling whole clusters.** In an unplanned rerun, resampling single rows produced little of either.
- **Resampling within strata understates variance when few clusters per stratum are audited,** as they were here. The plan did not anticipate this. In the same rerun, leaving out strata raised coverage on all nine items at both audit sizes.
- **The plan's simulations do not test what the option is for.** Their residuals are barely clustered once the recorded value is in the model.

`catmech_pool_lm(cluster = )` is a separate matter. It is the standard CR1 variance, it matches a hand calculation in the tests, and every method here used it, so this verdict does not judge it. Use it on clustered data, knowing that it undercovered here with 30 clusters (a mean of 90.4% on complete data, with 14 of 18 cells under 92%) and on some items with 69 (`lbw`, 88.4%).

There were **148,000 requests: 147,996 returned intervals, 4 refusals and no other errors**. All 4 refusals came on development iron at 25%, in cohorts 257 and 355, where both linear MI versions found no full-rank design. The run took 753 seconds at 14 workers. [Scored output](clustered-decisions.txt) comes from a separate script. The prospective numbers below come from it or from the [study's summary](clustered-study/summary.csv), and the post-hoc numbers come from [a separate table](clustered-posthoc.txt).

## Held-out items at a 25% facility audit

Coverage is clustered coverage of the pool slope, in percent, on 500 cohorts. The difference is paired, clustered minus existing. Ratios are clustered over existing, with 95% paired percentile limits. A ratio below 1 favors the clustered option.

| Item | Oracle | Existing | Clustered | Difference (MCSE) | Mean width, existing to clustered | Recovery RMSE ratio | RMSE ratio against pool slope |
|---|---|---|---|---|---|---|---|
| `lbw` | 88.4 | 92.0 | 89.8 | -2.2 (1.2) | 0.055 to 0.056 | 1.137 [1.104, 1.176] | 1.043 [1.019, 1.069] |
| `resus` | 94.0 | 96.0 | 93.6 | -2.4 (1.2) | 0.050 to 0.059 | 1.104 [1.058, 1.159] | 1.101 [1.057, 1.145] |
| `stillbirths` | 93.4 | 94.4 | 94.0 | -0.4 (0.9) | 0.013 to 0.014 | 1.050 [0.995, 1.116] | 1.036 [0.995, 1.084] |
| `deliveries` | 91.6 | 87.0 | 91.4 | +4.4 (1.0) | 0.234 to 0.249 | 0.903 [0.870, 0.934] | 0.931 [0.912, 0.949] |
| `anc1_early` | 95.6 | 79.6 | 88.8 | +9.2 (1.5) | 0.257 to 0.337 | 0.850 [0.827, 0.875] | 0.886 [0.864, 0.909] |

## Decision rule

| Condition | Result | Needed | Holds? |
|---|---|---|---|
| Coverage: closer to the oracle than the existing option | 5 of 5 | 4 of 5 | Yes |
| Coverage: within 5 points of the oracle | 4 of 5 (`anc1_early` missed) | 4 of 5 | Yes |
| Accuracy: upper limit of the recovery ratio at most 1.05 | 2 of 5 | 4 of 5 | **No** |
| Simulation: coverage at least 92% and within 3 points of the oracle, 70 or more clusters | 33 of 36 cells, 92% | 80% | Yes |

The three simulation cells that missed were 70 clusters, intraclass correlation 0.2, `iid`, at 50% (90.4%); 70 clusters, 0.5, `persistent`, at 25% (91.8%); and 150 clusters, 0.2, `persistent`, at 50% (91.2%). The oracle was also under 92% in all three (90.8%, 90.6%, 91.4%).

Two of these conditions separate the options less than the table suggests. The existing option also comes within 5 points of the oracle on four of five held-out items, with the same miss on `anc1_early`, and it meets the simulation condition in 28 of 36 cells, one short of 80%.

## Reading the result

**The count of 5 of 5 on coverage flatters the option.** On `lbw`, `resus` and `stillbirths`, the existing option already covered at or above the oracle, and the clustered one came down toward it. That satisfies the rule but rescues nothing, and these are the items where it lost accuracy. Measured against 95% rather than the oracle, a comparison made after seeing the results, the clustered option is closer on only `deliveries` and `anc1_early`.

**The gain is real where it appears.** Across all nine items at 25%, the worst clustered coverage was 88.8%, against 79.5% for the existing option. The oracle's own worst was 88.4%.

**The development slopes agree, but they suggested the idea.** They enter no rule.

| Development slope, 25% | Oracle | Existing | Clustered | Recovery RMSE ratio |
|---|---|---|---|---|
| PNC | 95.2 | 85.4 | 93.6 | 0.894 [0.856, 0.933] |
| Syphilis | 93.2 | 84.6 | 93.0 | 0.868 [0.835, 0.902] |
| Iron | 90.2 | 79.5 | 90.8 | 0.855 [0.808, 0.904] |
| ANC 4th visits | 96.6 | 87.2 | 90.8 | 0.970 [0.943, 0.999] |

**The simulations pass, but they do not test what the option is for.** The recorded value carries the cluster effect, so once it is in the model the residuals are barely clustered: their within-cluster correlation is about 0.00 to 0.11 across settings, against 0.20 to 0.67 on Rwanda's six common items. That is a flaw in the plan's design. On average the existing option is within a point of the oracle or above it in these cells, so they leave little to fix.

| Clusters | Cells | Oracle | Existing | Clustered |
|---|---|---|---|---|
| 30 | 18 | 90.4 | 92.4 | 93.0 |
| 70 | 18 | 93.0 | 93.1 | 93.5 |
| 150 | 18 | 93.8 | 93.2 | 93.5 |

They do show a small cost without any rare counts. With `iid` errors, intraclass correlation 0.5 and a 25% audit, the ratios are 1.072, 1.037 and 1.041, the only three of 54 cells clearly above 1. Six cells are clearly below 1, four of them with `tracking` errors.

**Very small audits can give absurd intervals.** One simulation cell has a mean width of 82,870: with 30 clusters at 25%, about eight audited clusters, one request had 0.16 degrees of freedom and a width of about 41 million. Nearly all the variance in such a request lies between imputations, and that collapses the Barnard-Rubin degrees of freedom. That cell's recovery RMSE barely moves (0.0867 against 0.0860), but such requests are themselves less accurate: across the nine 30-cluster cells at 25%, those with widths over five times their cell median have a recovery RMSE of 0.205, against 0.086 for the rest. There were 31 such requests of 28,000 simulated for the clustered option and 13 for the existing one, all in those nine cells. Median widths there were 3% to 7% wider. There were none in 5,000 held-out requests for either option, but 7 and 3 of 3,998 development requests, at about 17 audited clusters. On iron at 25%, 18.5% of clustered requests had under 5 degrees of freedom, against 6.4% for the existing option. Read the `df` column, and summarize repeated runs with a median. No code changed in response, because the sources were frozen for this study.

## Post hoc: what a recheck and an independent review found

None of this enters the rule or changes the verdict. The first version of this write-up was rechecked by the assistant that wrote it and then read by one independent AI reviewer. The unplanned runs below reuse the study's cohorts, audits and seeds, and change only the resampling. As a control, rerunning the study's own option reproduced all 8,998 of its returned requests exactly. See [the table](clustered-posthoc.txt) and [the records](clustered-review/variants.csv).

### Which items gained

| Item | Mean recorded count per facility-month | Zeros in the recount | Within-facility correlation of residuals | Coverage difference at 25% | Recovery RMSE ratio at 25% | Ratio at 50% |
|---|---|---|---|---|---|---|
| `stillbirths` | 0.2 | 86% | +0.01 | -0.4 | 1.050 | 1.028 |
| `resus` | 0.8 | 59% | -0.05 | -2.4 | 1.104 | 1.045 |
| `lbw` | 1.2 | 38% | -0.16 | -2.2 | 1.137 | 0.999 |
| ANC 4th visits | 19.8 | 2% | +0.37 | +3.6 | 0.970 | 0.978 |
| `anc1_early` | 22.4 | 0% | +0.67 | +9.2 | 0.850 | 0.976 |
| PNC | 29.0 | 0% | +0.37 | +8.2 | 0.894 | 0.939 |
| `deliveries` | 29.8 | 0% | +0.66 | +4.4 | 0.903 | 0.964 |
| Syphilis | 32.7 | 21% | +0.20 | +8.4 | 0.868 | 0.984 |
| Iron | 40.1 | 12% | +0.64 | +11.2 | 0.855 | 0.964 |

The correlation is that of the imputation model's residuals (recount on recorded value, covariate and catchment) on the whole pool.

Both how common the event is and that correlation separate the groups. Rarity cannot vary within an item, but the gain follows the correlation measured on each audit's own rows. Pooled over the six common items, coverage gained 11.5 points where it was at least .15 and 4.6 where it was lower, with a similar pattern at cutoffs from 0 to .3. On `deliveries` the gain is +19.2 against +0.7, with ratios of 0.62 and 1.04. `anc1_early` is the exception (+9.0 against +10.9).

Treat both as leads. The .15 cutoff was chosen after the results. Without three facilities, the pool correlation for `deliveries` falls from 0.66 to -0.01. All three rare items use deliveries as the covariate, and five of the six common items use antenatal registrations. Four of those six are also the slopes that suggested the idea. Errors that grow with facility size do not separate the groups (the correlation of absolute residual with the covariate is 0.14 to 0.47 in the first group and 0.25 to 0.48 in the second), and neither does interval width (`resus` widened by 19% and lost accuracy, `deliveries` widened by 7% and gained).

### Clusters, rows and strata

Coverage and recovery RMSE ratio against the existing option, at 25%. "Rescaled" is a scratch version and not a package feature.

| Item | Existing | Clusters within strata (the study) | Rows within strata | Clusters, no strata | Rescaled | Oracle | Ratio: rows | Ratio: no strata | Ratio: rescaled |
|---|---|---|---|---|---|---|---|---|---|
| `stillbirths` | 94.4 | 94.0 | 95.0 | 96.6 | 96.4 | 93.4 | 0.986 | 1.014 | 1.254 |
| `resus` | 96.0 | 93.6 | 95.2 | 96.4 | 96.6 | 94.0 | 0.984 | 1.172 | 1.375 |
| `lbw` | 92.0 | 89.8 | 91.8 | 92.8 | 95.0 | 88.4 | 0.983 | 1.389 | 1.446 |
| ANC 4th visits | 87.2 | 90.8 | 90.0 | 93.4 | 94.6 | 96.6 | 0.996 | 0.941 | 1.030 |
| `anc1_early` | 79.6 | 88.8 | 77.8 | 94.4 | 94.2 | 95.6 | 1.007 | 0.814 | 0.833 |
| PNC | 85.4 | 93.6 | 88.6 | 95.2 | 95.4 | 95.2 | 0.988 | 0.943 | 1.001 |
| `deliveries` | 87.0 | 91.4 | 87.8 | 92.2 | 91.8 | 91.6 | 0.979 | 0.886 | 0.963 |
| Syphilis | 84.6 | 93.0 | 87.0 | 94.6 | 95.8 | 93.2 | 0.954 | 0.836 | 0.901 |
| Iron | 79.5 | 90.8 | 78.9 | 95.2 | 93.6 | 90.2 | 1.002 | 0.888 | 0.920 |

**Rows.** Resampling single rows within strata changed coverage on the six common items by -1.8 to +3.2 points, against +3.6 to +11.2 for whole clusters. Its recovery RMSE ratios run from 0.954 to 1.007, and on the rare items from 0.983 to 0.986, against 1.050 to 1.137 for clusters. So both the gain and the cost come mainly from resampling clusters.

**Strata.** The study resampled within catchment, as the plan says. A bootstrap mean within a stratum of n audited clusters has (n - 1)/n of the usual variance, and a stratum with one audited cluster contributes none. A 25% audit here has 17 audited clusters in seven catchments (1, 2, 2, 2, 3, 3 and 4), so the shortfall is large, and the documentation did not warn of it. Leaving out strata raised coverage on all nine items at both audit sizes, including `anc1_early` from 88.8% to 94.4%, and iron from 90.8% to 95.2%. It also raised error further on `lbw` and `resus`. A scratch rescaled bootstrap in the style of Rao and Wu (1988) covered 91.8% to 96.6% on all nine, with median between-imputation variance 1.54 to 2.15 times the study's, and it too cost more accuracy on the rare items. The 6.8-point miss on `anc1_early` may therefore come largely from this flaw and not from the idea. A rescaled version needs its own planned test.

### Against a larger audit

Doubling the audit from 25% to 50% cut the existing option's recovery RMSE by 34% to 56% on all nine items, against a change of -15% to +14% from the clustered option. On error a larger audit matters more.

On coverage the option can match or beat doubling on some items, with no more auditing, but it does so by widening intervals (+7% to +41% on the six common items), while doubling mostly narrows them (-33% to +1%). Option at 25% minus doubling, the paired difference is clearly positive on PNC, syphilis and iron (+7.6, +7.8 and +6.0, MCSE about 2), unclear on `deliveries` and `anc1_early` (+2.6 and +0.8), and negative on ANC 4th visits (-4.4). On `anc1_early` the option reached 88.8% at a width of 0.337, and doubling reached 88.0% at 0.205.

### As the audit grows

The option's advantage on error fades at 50%, with ratios of 0.964 to 1.045 on the held-out items. Its coverage gains do not all fade: +4.4, +5.4 and +5.0 points on `anc1_early`, PNC and iron. The existing option stays 7.6 to 9.2 points under the oracle on `anc1_early`, PNC and syphilis, and on syphilis the clustered option also falls, to 88.0%.

### Running both

Within an item, the existing option covered less in the third of audits where the clustered interval was widest relative to it, on 8 of 9 items. On syphilis it covered 69.9% in that third and 91.9% in the rest. This gives some support to running both and treating a much wider clustered interval as a warning, but no threshold has been tested.

## Predictions, scored

| # | Prediction | Result |
|---|---|---|
| 1 | Closer to the oracle's coverage in at least 4 of 5 | **Held**, 5 of 5, with the caveat above |
| 2 | Mean width 15% to 60% wider in all five | **Failed**: +1%, +19%, +9%, +7%, +31% |
| 3 | Recovery RMSE no worse by more than 5% in at least 4 of 5 | **Failed**, 2 of 5 on the upper limit the rule uses and 3 of 5 on point estimates |
| 4 | With 30 clusters at 25%, coverage below 92% in most cells | **Failed**, 0 of 9; the range was 92.8% to 96.8%, so eight audited clusters did better on coverage than predicted |
| 5 | Width within 5% in the one-row controls | **Held**: -0.9% and -0.5% |
| 6 | The gap shrinks from 25% to 50% | **Held in 4 of 5**, tied on `stillbirths`: 2.2 to 1.0, 2.4 to 1.0, 0.4 to 0.4, 4.4 to 2.2, 9.2 to 4.4 |

## At a 50% audit

| Item | Coverage difference (MCSE) | Recovery RMSE ratio |
|---|---|---|
| `lbw` | -1.0 (0.8) | 0.999 [0.974, 1.028] |
| `resus` | +1.0 (0.7) | 1.045 [1.009, 1.084] |
| `stillbirths` | -0.4 (0.7) | 1.028 [1.001, 1.058] |
| `deliveries` | +2.2 (0.8) | 0.964 [0.948, 0.978] |
| `anc1_early` | +4.4 (1.0) | 0.976 [0.957, 0.994] |

## What this changes

- **Use `catmech_pool_lm(cluster = )` on clustered data,** with the limits stated at the top. This is unaffected by the verdict.
- **Treat `catmech_reference(cluster = )` as experimental.** The trade is better coverage where the existing option falls short, against up to 14% more error elsewhere. Be wary of it on counts averaging about one or fewer, though that reading is post hoc and nothing was tested between 1.2 and 19.8. Running both and reporting any disagreement is a reasonable precaution, with the support above and no planned test.
- **Be careful with `cluster_strata`.** With few audited clusters per stratum, it understates variance. Leaving it out raised coverage here, but the test never compared the two, and on rare items it cost accuracy.
- **A larger audit matters more for error.** The clustered option can matter as much for coverage, at the price of wider intervals.
- **The earlier wording was too strong.** The [real-audit synthesis](REAL_AUDIT_SYNTHESIS.md) said clustered audits need clustered imputation. It now points here.
- **A next version needs a new plan.** It should test a rescaled bootstrap, simulations with clustered residuals, and whether within-facility correlation on audited rows shows when the option helps.

## Limits

The held-out items are held out by item and not by facility: they share 69 health centers and their months with the analyses that suggested the idea, so a facility-level oddity could affect all nine. It is one health system, and Rwanda's items were chosen from what the file offers, not sampled. `anc1_early`'s two sides are defined slightly differently, as the plan notes, and 85% of its rows disagree, so its gain may reflect the definition as much as recording errors. With 500 cohorts, the Monte Carlo error on a paired coverage difference is about one point.

On process: the plan was committed before the option existed, and no amendment was needed. An engineering pilot on separate seeds (4,144 requests, all returned) was read for statuses and timing, plus a check that its development rows matched the earlier study, with a largest difference of about 1e-15. Hashes of the sources were checked before and after the run, and existing workflows matched their frozen versions. Commit order in a local repository is not third-party proof of timing. This round was written and scored by an AI assistant. AI review is not independent human peer review, and no human has reviewed it.

## What the recheck and review changed

The reviewer re-derived the results and found no error in the verdict, in the code for `cluster`, or in the comparisons. It and the assistant's own recheck corrected these statements in the first version:

- "A single audit cannot show which case applies." Both features above can be read from an audit, though neither is tested.
- "A larger audit matters more than this choice." That holds for error and not always for coverage.
- The split was attributed to how common the event is, with "the cause unknown." Within-facility correlation separates the groups equally well.
- "This test cannot say how much of the gain comes from clustering." An unplanned run on the same cohorts answers much of that.
- The simulations "could not have shown the cost." They show a small one.
- "Point estimates are unaffected" by very wide intervals. Cell RMSE barely moves, but those requests are less accurate.
- "Use it on any clustered data." It undercovers with few clusters.
- The advantage "fades as the audit grows." That holds for error only.
- "Mostly zeros" as the warning sign. `lbw`, the largest cost, has 38% zeros in the recount.
- The stratification shortfall and the simulations' lack of clustered residuals were not disclosed.

## Sources

The idea combines the bootstrap at the imputation step (von Hippel and Bartlett, 2021) with the clustered bootstrap and cluster-robust inference (Cameron and Miller, 2015). See the [plan](CLUSTERED_PLAN.md) for references. The rescaled bootstrap follows Rao JNK, Wu CFJ. Resampling inference with complex survey data. *Journal of the American Statistical Association* 1988. [doi:10.1080/01621459.1988.10478591](https://doi.org/10.1080/01621459.1988.10478591).
