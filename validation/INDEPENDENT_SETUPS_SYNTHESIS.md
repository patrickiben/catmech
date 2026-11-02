# Test of comparators set up by fresh agents: results and interpretation

Under the [rule fixed in advance](INDEPENDENT_SETUPS_PLAN.md), the same rule as the [earlier comparison](OUTSIDE_METHODS_SYNTHESIS.md), the outcome is **KEEP**. No setup was no worse than `catmech_reference(cluster = )` on all six conditions.

**I do not think that is the message. `ind_mice` looks better than the option, and it missed on a technicality.** It reached 92% on all 11 real items at 25%, against the option's 9, with accuracy within about 1%. It failed only because it covered above 98% on `resus` and `stillbirths`, including 98.05% on one where the limit was 98. With the audited value as a predictor, where the option's slope is biased by -16.1%, it had +1.4%, with coverage of 94.6% against 84.4%. These real items were never shown to the agent. The rule's result stands and thresholds do not move, so this does not retire the option under the rule. But the documents no longer imply that the option is the best available.

There were 148,500 requests with 2 giving no interval, plus 54,000 on the second health system. The study took about 15 minutes. A fresh rerun of 198 real-item requests [reproduced the saved results](independent-setups/replay-check.txt). [Scored output](independent-setups-decisions.txt) and [post-hoc numbers](independent-setups-posthoc.txt) come from separate scripts. The agents' [code](../tools/independent-setups/) is as written, with [their reports](independent-setups/agent-reports.json) saved. Each reported keeping to the rule against opening the repository, which is an instruction and was not enforced.

## Decision rule

The option reaches 92% on 9 of 11 real items at 25%, with mean distance from 95% of 1.47, and 10 of 11 at 50%, with 1.25, and has 30 of 36 simulated cells within limits.

| Setup | Accuracy, 25% (N1, at most 1.05) | Accuracy, 50% (N2, at most 1.05) | Real, 25%: at 92 or more, at most 98, distance (N3) | Real, 50% (N4) | Simulations, ratio and cells (N5) | Available (N6) | No worse |
|---|---|---|---|---|---|---|---|
| `ind_mice` | 1.012 | 1.009 | 11, **9**, 1.68 | 11, **9**, 1.62 | 0.960, 33 | Pass | No |
| `ind_survey` | **1.141** | **1.106** | 11, 10, **2.02** | 11, **8**, 2.00 | 0.998, 36 | Pass | No |
| `ind_design` | **1.165** | **1.079** | 11, 10, **1.75** | 11, **6**, **2.40** | 1.030, 35 | Pass | No |

Bold marks what failed. Each condition is a set of checks, and N3 and N4 also need at most 98% on at least 10 of 11.

## What the agents did differently

These are their own choices, from [their reports](independent-setups/agent-reports.json). All three are more elaborate than a default use of a package, and all three use heuristics that they say they tuned on the generator.

- **`ind_mice`** uses `mice` only as a wrapper around its own method. It imputes the error and not the value, with separate within-cluster and between-cluster slopes, a jackknife analysis variance and its own degrees of freedom. When the audited value is a predictor, it shrinks the imputation variance by a heuristic factor, which is a deliberate departure from proper multiple imputation. Separate within-cluster and between-cluster slopes are the visible difference from the option, and the statistics reviewer in the [overall check](OVERALL_QC.md) pointed to a similar contrast. Whether they are what helps has not been tested on their own.
- **`ind_survey`** does not use `twophase`. It computes a linear calibration estimator on partialled-out values, with a jackknife variance and a count-based floor.
- **`ind_design`** is a difference estimator with an estimated coefficient on the recorded value, leverage corrections, and a jackknife variance for counts.

## What the results show

- **All three repair coverage on real items.** Each reaches 92% on all 11 at both audit sizes. The earlier recipes for `raking` and `difference` reached none at 25%. So the [earlier comparison's](OUTSIDE_METHODS_SYNTHESIS.md) low coverage for design-based methods mostly reflected my recipes, as it suspected.
- **Their accuracy on real items does not catch the option's, apart from `ind_mice`.** `ind_survey` and `ind_design` are 14% and 16% worse at 25%, and their intervals are wider, by +15% to +91% and +50% to +189%. A model pays where its assumptions roughly hold.
- **`ind_mice` over-covers somewhat on rare real items**, with RMSE well below the option's there: 0.711 on `resus`, 0.789 on `stillbirths`. On those items the recount itself covers 96.65% to 97.80%. Its intervals run from -15% to +41% of the option's width at 25%. It is less accurate on `iron` (1.458) and `pnc` (1.239).
- **The predictor problem is largely solved by two setups in these cells.** With `shift` errors at 25%, bias was +1.4% for `ind_mice` and +1.2% for `ind_survey`, against -16.1% for the option and -3.0% for linear MI, with RMSE 0.951 and 0.984 of linear MI's. These six cells were available to the agents through the generator, so they may flatter.
- **The six simulated rare cells improve, at a price.** Mean coverage was 94.17% for `ind_survey`, 92.73% for `ind_design` and 91.53% for `ind_mice`, against 88.77% for the option and 87.50% for a recount, with intervals 20%, 31% and 16% wider than the option's. These cells were also in the agents' generator, and they tuned on them. Covering above the recount means these intervals are wider than the recount's own. It does not show better estimates.
- **On the second health system, whole-cluster audits, description only:** `ind_survey` has a mean distance from 95% of 5.83 at `district_33` and 8.38 at `district_50`, against 10.12 and 10.59 for the option. `ind_mice` is close to the option, at 9.70 and 10.34. These items were held from the agents too, and the [same limits on single districts](SECOND_SYSTEM_SYNTHESIS.md) apply.

## Predictions, scored

| # | Prediction | Result |
|---|---|---|
| 1 | The outcome is KEEP | **Held** |
| 2 | `ind_mice` has an N1 ratio at most 1.05, and fails N3 | **Held** |
| 3 | `ind_survey` and `ind_design` fail N1 | **Held** |
| 4 | `ind_survey` reaches 92% on more real items at 25% than the earlier recipe | **Held**, all 11 |

## What this changes

- **No package code changes under this rule.** `catmech_reference(cluster = )` stays, and stays experimental.
- **The documents now say** that it is not the best available, and point to these setups as research code. They have not been reviewed and are not package features.
- **The obvious next step is its own planned test:** whether separate within-cluster and between-cluster slopes improve the option, on data not yet seen. The real items here can no longer serve as a held-out test, because I have now seen these results.

## Limits

This reduces the problem of the option's author setting up the comparators, and does not remove it. The agents are the same model as the author, and they got a brief written by the author, with three warnings from earlier work. They tuned on the generator behind the simulated cells, so only the real items and the second health system are out of sample for them, and those items have shaped earlier choices about the option. Three agents are not a survey of practice, and a statistician might set these up differently or find errors in them. No human has reviewed any of this. The rule's 98% limit, set before these results, decided the outcome for `ind_mice`. Coverage near 98% with 500 cohorts has a Monte Carlo error of about 0.6 points, so the N4 failure is near the limit, on the same two items.

The scorer's control had a lookup fault, fixed before any results were read, as the plan's amendment records.
