# Test of imputation for clustered audits: results and interpretation

Under the [rule fixed in advance](CLUSTERED_PLAN.md) the verdict is **PARTIAL**. Resampling audited facilities in each imputation passed the coverage and simulation conditions, but it failed the accuracy condition on three of five held-out Rwanda items. As the plan requires, `catmech_reference(cluster = )` stays in the package as an experimental option with that failure stated. It is not recommended as a default.

What holds:

- **It helps where the existing option undercovers.** On the two held-out items where `catmech_reference(use_recorded = TRUE)` fell well short of a recount, clustered coverage rose by 4.4 and 9.2 points, and error fell. The four development slopes moved the same way.
- **It costs where the existing option was already fine.** On the other three held-out items, coverage fell toward the oracle's from above and error rose by 14%, 10% and 5%, the last within Monte Carlo noise.
- **After the test, the split turned out to follow how common the event was.** The three items where it cost accuracy are rare events, averaging 1.2 or fewer recorded per facility-month, with 38% to 86% zeros in the recount. The six where it helped average 19.8 to 40.1, with at most 21% zeros. That was seen only after the results, on nine items from one health system, so it is a lead and not a rule. Interval width does not sort them: `resus` widened by 19% and lost accuracy, while `deliveries` widened by 7% and gained.
- **Its advantage fades as the audit grows.** At a 50% audit, the error ratios run from 0.964 to 1.045.

`catmech_pool_lm(cluster = )` is a separate matter. It is the standard CR1 variance, it matches a hand calculation in the tests, and every method here used it, so this verdict does not judge it. Use it on any clustered data.

There were **148,000 requests: 147,996 returned intervals, 4 refusals and no other errors**. All 4 refusals came on development iron at 25%, in cohorts 257 and 355, where both linear MI versions found no full-rank design. The run took 753 seconds at 14 workers. [Scored output](clustered-decisions.txt) comes from a separate script, and every number below comes from it or from the [study's summary](clustered-study/summary.csv).

## Held-out items at a 25% facility audit

Coverage is clustered coverage of the same cohort-level target, in percent, on 500 cohorts. The difference is paired, clustered minus existing. Ratios are clustered over existing, with 95% paired percentile limits. A ratio below 1 favors the clustered option.

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

**The simulations show that the option does little harm, not that it helps.** On average the existing option is within a point of the oracle or above it in them, so they do not reproduce the Rwanda shortfall and leave nothing to fix. They also have no rare-event counts, so they could not have shown the cost either.

| Clusters | Cells | Oracle | Existing | Clustered |
|---|---|---|---|---|
| 30 | 18 | 90.4 | 92.4 | 93.0 |
| 70 | 18 | 93.0 | 93.1 | 93.5 |
| 150 | 18 | 93.8 | 93.2 | 93.5 |

**This test cannot say how much of the gain comes from clustering.** The simulated data have clusters but no shortfall, and the Rwanda data have a shortfall along with [unequal variance](REAL_AUDIT_SYNTHESIS.md#what-limits-the-intervals). A bootstrap also absorbs unequal variance, so resampling single rows might have helped too. At quarter totals, with one row per facility, the earlier test found linear MI covering 84.9% to 97.0%, which suggests clustering plays a part without settling it. A row-level comparison was not in the plan and was not run. The option's name states what it does, not why it helped.

**Very small audits can give absurd intervals.** One simulation cell has a mean width of 82,870: with 30 clusters at 25%, about eight audited clusters, one request had 0.16 degrees of freedom and a width of about 41 million. Nearly all the variance in such a request lies between imputations, and that collapses the Barnard-Rubin degrees of freedom. Point estimates are unaffected, with recovery RMSE of 0.0867 against 0.0860 in that cell. Widths over five times their cell median occurred in 31 of 28,000 simulated requests for the clustered option and 13 of 28,000 for the existing one, all in the nine 30-cluster cells at 25%. There were none in 5,000 held-out requests for either. Median widths in those nine cells were 3% to 7% wider. Read the `df` column, and summarize repeated runs with a median. No code changed in response, because the sources were frozen for this study.

## Seen after the results: how common the event is

This section is post hoc and enters no rule. It came from checking an earlier version of this write-up, which said that nothing tells a user in advance which case applies.

| Item | Mean recorded count per facility-month | Zeros in the recount | Coverage at 25%, existing and oracle | Coverage difference at 25% | Recovery RMSE ratio at 25% | Ratio at 50% |
|---|---|---|---|---|---|---|
| `stillbirths` | 0.2 | 86% | 94.4, 93.4 | -0.4 | 1.050 | 1.028 |
| `resus` | 0.8 | 59% | 96.0, 94.0 | -2.4 | 1.104 | 1.045 |
| `lbw` | 1.2 | 38% | 92.0, 88.4 | -2.2 | 1.137 | 0.999 |
| ANC 4th visits | 19.8 | 2% | 87.2, 96.6 | +3.6 | 0.970 | 0.978 |
| `anc1_early` | 22.4 | 0% | 79.6, 95.6 | +9.2 | 0.850 | 0.976 |
| PNC | 29.0 | 0% | 85.4, 95.2 | +8.2 | 0.894 | 0.939 |
| `deliveries` | 29.8 | 0% | 87.0, 91.6 | +4.4 | 0.903 | 0.964 |
| Syphilis | 32.7 | 21% | 84.6, 93.2 | +8.4 | 0.868 | 0.984 |
| Iron | 40.1 | 12% | 79.5, 90.2 | +11.2 | 0.855 | 0.964 |

A user can see the count in the recorded column before any audit. Treat this as a lead to test: it covers nine items from one health system, it was noticed after the results, and four of the six items on the favorable side are the ones that suggested the idea. The cause is unknown. Errors that grow with facility size do not separate the groups, since the correlation of absolute residual with the covariate is 0.14 to 0.47 in the first group and 0.25 to 0.48 in the second. See [the table](clustered-posthoc.txt).

**Against a larger audit.** Doubling the audit from 25% to 50% cut the existing option's recovery RMSE by 34% to 56% on all nine items, against a change of -15% to +14% from the clustered option, so on error a larger audit matters more. On coverage it can be the reverse. On the six common items, the clustered option gained more coverage than doubling the audit on five, all but ANC 4th visits, and it needs no more auditing.

## Predictions, scored

| # | Prediction | Result |
|---|---|---|
| 1 | Closer to the oracle's coverage in at least 4 of 5 | **Held**, 5 of 5, with the caveat above |
| 2 | Mean width 15% to 60% wider in all five | **Failed**: +1%, +19%, +9%, +7%, +31% |
| 3 | Recovery RMSE no worse by more than 5% in at least 4 of 5 | **Failed**, 2 of 5 |
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

- **Use `catmech_pool_lm(cluster = )` on clustered data.** This is unaffected by the verdict.
- **Treat `catmech_reference(cluster = )` as experimental.** The trade is better coverage where the existing option falls short, against up to 14% more error elsewhere. On these items the losses came on rare events, so be wary of it on counts that are mostly zeros, though that reading is post hoc. Running both and reporting any disagreement is a reasonable precaution, but this study did not test it.
- **The earlier wording was too strong.** The [real-audit synthesis](REAL_AUDIT_SYNTHESIS.md) said clustered audits need clustered imputation. It now points here.
- **A larger audit matters more for error, and the clustered option can matter more for coverage.** See above; this too is post hoc.

## Limits

The held-out items are held out by item and not by facility: they share 69 health centers and their months with the analyses that suggested the idea, so a facility-level oddity could affect all nine. It is one health system, and Rwanda's items were chosen from what the file offers, not sampled. `anc1_early`'s two sides are defined slightly differently, as the plan notes, and 85% of its rows disagree, so its gain may reflect the definition as much as recording errors. With 500 cohorts, the Monte Carlo error on a paired coverage difference is about one point.

On process: the plan was committed before the option existed, and no amendment was needed. An engineering pilot on separate seeds (4,144 requests, all returned) was read for statuses and timing, plus a check that its development rows matched the earlier study, with a largest difference of about 1e-15. Hashes of the sources were checked before and after the run, and existing workflows matched their frozen versions. Commit order in a local repository is not third-party proof of timing. This round was written and scored by an AI assistant and has had no independent review, human or AI, unlike the [real-audit test](REAL_AUDIT_SYNTHESIS.md). When the same assistant rechecked its own write-up, it found two statements too strong: that nothing tells a user in advance which case applies, and that a larger audit matters more than this choice. Both are corrected here.

## Sources

The idea combines the bootstrap at the imputation step (von Hippel and Bartlett, 2021) with the clustered bootstrap and cluster-robust inference (Cameron and Miller, 2015). See the [plan](CLUSTERED_PLAN.md) for references.
