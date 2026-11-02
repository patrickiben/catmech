# Test on a second health system: results and interpretation

Under the [rule fixed in advance](SECOND_SYSTEM_PLAN.md) the outcome is **MIXED**. So the Rwanda result neither replicated nor was contradicted. `catmech_reference(cluster = )` stays experimental, with no deprecation under this rule.

**The plain reading is that the option added little here.** Its accuracy matched linear imputation's: the geometric mean RMSE ratio was 1.005 and 0.988 in the two cells the rule uses. In Rwanda, the same comparison gave lower error on common items. Its coverage was closer to 95% by about 1 point in `clinic_25`, enough to support it there, and by 0.23 points in `district_50`, which is not enough. The support in `clinic_25` is fragile: with the both-zero rows kept, the gap shrinks to about 0.1. Errors on this file barely cluster by district, with an intraclass correlation of -0.063 to 0.201, so there was little for a random intercept to capture. That is consistent with the option helping only where errors cluster, though one file cannot establish it.

There were 270,000 requests and 7 without an interval. The study took about 47 minutes. [Scored output](second-system-decisions.txt) and [post-hoc numbers](second-system-posthoc.txt) come from separate scripts. Everything labeled post hoc came after the results and enters no rule.

## Decision rule

| Cell | D, option | D, linear MI | R | Available | Reading |
|---|---|---|---|---|---|
| `clinic_25` | 7.98 | 8.99 | 1.005 | 99.99% | Supports |
| `district_50` | 10.59 | 10.82 | 0.988 | 100.00% | Neither |

D is the mean distance from 95% over nine items, and R the option's RMSE over linear MI's. Support needed D at least 0.5 lower and R at most 1.02.

## Why the coverage numbers look so low

This part is post hoc. Even a full recount averages only about 12 points from 95%, covering 64.0% on `deliveries`, 68.1% on `fpmethods` and 65.5% on `exposedbabies`. Those figures do not mean a method failed in the usual way. On these three items, one district almost decides the slope:

| Item | Clinics in that district, of pool | Pool slope | Without that district | Recount coverage with it in the cohort | Without it |
|---|---|---|---|---|---|
| `deliveries` | 4 of 105 | -75.7 | 32.6 | 100.0 | 0.0 |
| `fpmethods` | 27 of 116 | -73.7 | 141.7 | 98.8 | 10.6 |
| `exposedbabies` | 4 of 78 | -9.3 | 5.4 | 99.5 | 0.0 |

Each cohort resamples the 18 districts, so about a third of cohorts leave out any given one, and those cohorts cannot estimate a slope set by what they lack. No method can get these right from the data in such a cohort, which is why linear MI and the option also fall from about 85% to about 16% to 34% on the two small-district items. So these items say more about the data and this kind of resampling test than about the methods, and they dominate D. **The same design could flatter or punish a method for reasons unrelated to imputation.** A cluster-robust interval also cannot warn about a cluster that is not in the sample. This is a limit of the test, and it was not foreseen in the plan beyond a general warning about 18 unequal clusters.

On the other six items, the recount covers about 90% to 93%, and imputation methods run higher, often at 95% to 99%.

## Results by method

Mean distance from 95%, items at 92% or more of nine, and geometric mean RMSE against linear MI:

| Audit | Recount | Recorded | Audited rows only | Linear MI | Option | `mice_2l_lmer` |
|---|---|---|---|---|---|---|
| `clinic_25` | 12.09, 3, 0.690 | 11.56, 4, 1.102 | 15.72, 4, 1.403 | 8.99, 7, 1.000 | 7.98, 7, 1.005 | 8.87, 7, 1.003 |
| `clinic_50` | 12.09, 3, 0.820 | 11.56, 4, 1.310 | 15.21, 3, 1.158 | 10.68, 6, 1.000 | 10.60, 6, 0.995 | 10.94, 6, 0.999 |
| `district_33` | 12.09, 3, 0.772 | 11.56, 4, 1.232 | 18.18, 0, 1.438 | 10.19, 6, 1.000 | 10.12, 6, 0.987 | 10.44, 6, 0.989 |
| `district_50` | 12.09, 3, 0.853 | 11.56, 4, 1.363 | 17.43, 0, 1.219 | 10.82, 6, 1.000 | 10.59, 6, 0.988 | 10.58, 6, 0.988 |

- **Linear MI's undercoverage from Rwanda did not appear in the same way.** Apart from the three items above, it mostly covered at or above 94%, and often above 98%. In Rwanda it fell to about 81% on items where errors clustered. So the warning about linear MI on clustered audits rests on one system, where residual errors clustered more on those items, measured there by facility and not by district.
- **Imputation paid against the audited rows alone.** On most items audited rows alone had far higher error.
- **The recorded value won less often than in Rwanda.** It had lower error than linear MI on 3 of 9 items at `clinic_25` and on 1 of 9 at `clinic_50`.
- **`mice_2l_lmer`, with the recorded value kept, again matched the option on accuracy** and on these items matched its coverage more closely than in Rwanda.
- **A warning sign to read.** Linear MI's mean width at `district_33` is inflated by a few absurd intervals: one is 3.8e+09 wide against a median of 459, with pooled degrees of freedom of 0.176. Pooled degrees of freedom were under 1 on 0.23% of requests there. These come from Barnard-Rubin degrees of freedom when audited clusters are few and between-imputation variance dominates. The README says to read the `df` column, and this shows why.

The sensitivity run keeping both-zero rows, on 500 cohorts, moved D by less than 2 points. It left the option within about 0.1 points of linear MI in three audits and slightly behind it at `district_33` (9.80 against 9.69), so the small coverage lead on this file is not robust ([scored output](second-system-decisions.txt)).

## Predictions, scored

| # | Prediction | Result |
|---|---|---|
| 1 | The outcome is MIXED | **Held** |
| 2 | The oracle reaches 92% on at most six of nine items at `district_50` | **Held**, three |
| 3 | The option gives an interval on at least 95% of requests in every audit | **Held** |
| 4 | `recorded` has a lower RMSE than `linear_mi` on at least five of nine items at `clinic_25` | **Failed**, three |
| 5 | Linear MI reaches 92% on at least seven of nine items at `clinic_50` | **Failed**, six |

## What this changes

- **The clustered evidence is no longer from one health system,** though this second one is thin and adds little either way.
- **The option's accuracy gain did not replicate,** so the documents now say its gain on Rwanda may depend on how strongly errors cluster.
- **No package code changes under this rule.**

## Limits

There are 18 unequal clusters and one quarter, so a clinic has one row per item, and clustering exists only by district. Items within a clinic are correlated. Facilities were selected with probability proportional to size, and these analyses ignore that. Both-zero rows were dropped by a judgment call, and the record's two README files were not read. The structure and agreement counts had been seen before the plan, and the plan did not foresee how far single districts would drive the results. The same assistant wrote the plan, runner, scoring and this write-up, and no reviewer has read it.

## Sources

See the [plan](SECOND_SYSTEM_PLAN.md) and [the data note](second-system/README.md).
