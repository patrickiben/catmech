# Audit planner test: results and interpretation

**Verdict: NO-GO** under the rule [fixed before the run](AUDIT_PLAN_PLAN.md). The rule required a median error of at most 10% at every target fraction in both datasets, and the 5% target failed with both pilots in both datasets:
11.3% on ClinicalTrials.gov and 13.9% on NHANES with the 10% pilot, and 16.8% and 15.5% with the 25% pilot. Every other target was within 10%, the largest being 7.8%. The rule has no exception for that pattern, so `catmech_audit_plan()` should not be described as validated.

This is a test on two public datasets. It is not a clinical audit, and these are new cohorts from the same pools as the [frozen real-data study](REAL_DATA_SYNTHESIS.md), whose tables had already been seen by hand, as the plan discloses.

The run covered 1,500 cohorts per real dataset and 300 per synthetic cell: 14,400 target fits and 7,200 plans. **236 target fits and 19 plans were refused**, all because a category had no row with both values in a small audit. There were no other errors. It took 878 seconds at 14 workers against my estimate of about 12 minutes, on a machine shared with other jobs. Plain Python recomputed every summary, the verdict and both pool values of `k`.

Sources are in [`audit-plan-study/`](audit-plan-study/). Each table names its file.

## The rule

Median absolute error of predicted against actual width, with a 95% order-statistic interval and the 90th percentile. From `summary.csv`. The limit is 10%.

| Dataset, pilot | Target 5% (30 records) | Target 10% (60) | Target 25% (150) | Target 50% (300) | Rule |
|---|---:|---:|---:|---:|---|
| ClinicalTrials.gov, 10% (60) | **11.3%** (10.2 to 12.4), p90 38.8% | | 2.4% (2.2 to 2.6), p90 8.5% | 1.5% (1.4 to 1.6), p90 4.8% | Fails at 5% |
| NHANES, 10% (60) | **13.9%** (13.1 to 14.8), p90 35.3% | | 3.4% (3.2 to 3.7), p90 9.6% | 2.8% (2.6 to 2.9), p90 7.2% | Fails at 5% |
| ClinicalTrials.gov, 25% (150) | **16.8%** (15.3 to 17.8), p90 43.8% | 6.6% (6.2 to 7.1), p90 20.6% | | 0.9% (0.8 to 1.0), p90 2.6% | Fails at 5% |
| NHANES, 25% (150) | **15.5%** (14.7 to 16.4), p90 37.5% | 7.8% (7.3 to 8.3), p90 19.6% | | 1.6% (1.5 to 1.7), p90 3.8% | Fails at 5% |

The closest failure was ClinicalTrials.gov with the 10% pilot, and its interval lies above the limit. Monte Carlo noise does not explain it.

Refusals at the 5% target: 100 target fits on ClinicalTrials.gov and 117 on NHANES, of 1,500 each, because a rare category had no usable row among 30 records. Those cohorts have no actual width and left the cell. At the 10% target the counts were 9 and 10, and the same number of plans from 60-record pilots were refused.
Refused plans count as infinite error wherever a target fit existed, per the plan. In 11 more ClinicalTrials.gov cohorts, all 30 rows of the 5% audit agreed, so linear MI kept the recorded values with no imputation variance (`targets.csv`, `recorded_fit`).

## Why the 5% target fails

The predictions are nearly centered there, with signed median errors of -3.6%, -3.9%, +1.6% and -1.9% (`summary.csv`). The problem is spread, because the width of a single 30-record audit is itself hard to predict.

| Median error of a reference predictor (`summary.csv`) | ClinicalTrials.gov | NHANES |
|---|---:|---:|
| Noise floor at 5%: knows the cell's median actual width, and nothing about the cohort | 17.4% | 16.1% |
| Noise floor at 10%, 25%, 50% | 9.1%, 5.0%, 4.1% | 9.3%, 5.1%, 4.4% |
| Retest at 5%: same audit, another imputation seed, first 200 cohorts | 3.7% | 9.3% |
| Retest at 10%, 25%, 50% | 2.3%, 1.2%, 0.5% | 5.1%, 2.1%, 0.9% |

On NHANES, imputation noise alone at `m = 20` nearly uses the whole 10% allowance. The planner beat the floor with the 10% pilot, but that pilot contains the target's 30 records, an overlap a real plan would not have. With the 25% pilot, where the overlap is smaller, the planner came within 0.6 points of the floor in both datasets (16.8% against 17.4%, and 15.5% against 16.1%).

Better estimates of `k` would not close that gap by themselves. I expected this from the development split (prediction 1), and the rule was not mine to change. Whether a narrower claim is worth testing is a decision for Patrick and his reviewer. Such a test would need new data, because these results have now been seen.

## Other descriptive results

**Interval coverage of the actual ratio**, meaning actual width over the same cohort's full-audit width (`summary.csv`). The `k` interval and the predictive range are both 95% nominal, and neither reached it. The `k` interval is not meant to cover one audit's width, as this table confirms. Coverage by the predictive range was better but still fell short, most of all when predicting large audits from a small pilot.

| Dataset, pilot, target | `k` interval | Predictive range |
|---|---:|---:|
| ClinicalTrials.gov, 60, then 5%, 25%, 50% | 60.3%, 62.7%, 54.3% | 83.1%, 69.5%, 60.3% |
| ClinicalTrials.gov, 150, then 5%, 10%, 50% | 46.7%, 67.9%, 67.0% | 83.0%, 92.3%, 77.3% |
| NHANES, 60, then 5%, 25%, 50% | 48.5%, 60.6%, 53.4% | 83.9%, 84.2%, 74.1% |
| NHANES, 150, then 5%, 10%, 50% | 35.2%, 44.1%, 47.7% | 82.5%, 87.9%, 76.7% |

**Picking the fraction.** The planner chose the same smallest fraction as the truth for reaching within 25% of the full-audit width in 59% (60-record pilot) and 48% (150) of cohorts on ClinicalTrials.gov, and 57% and 52% on NHANES (`fraction-picks.csv`). The truth is noisy here too. On ClinicalTrials.gov it was 5% in 52% of cohorts, 10% in 28%, 25% in 19% and 50% in 1.5%.

**Pilot estimates against the pool** (`k-summary.csv`, `pool-k.csv`). The pool values are the same calculation with every row audited.

| Dataset | Pool `k` | Pilot | Median | 5th to 95th percentile | `k` interval contains pool `k` |
|---|---:|---:|---:|---:|---:|
| ClinicalTrials.gov | 0.0418 | 60 | 0.0197 | 0.0012 to 0.131 | 53.7% |
| ClinicalTrials.gov | | 150 | 0.0285 | 0.0067 to 0.127 | 71.1% |
| NHANES | 0.0698 | 60 | 0.0598 | 0.0313 to 0.135 | 71.3% |
| NHANES | | 150 | 0.0658 | 0.0392 to 0.110 | 74.3% |

Small pilots estimate `k` poorly, especially with exact-or-different errors: a median of 8 disagreeing rows on ClinicalTrials.gov at 60 records, so a few rows decide `k`. Their median falls below the pool value because the distribution is skewed, and most pilots miss some of the rare large errors. The interval for `k` undercovers at both sizes, in both datasets. Moderate-fraction widths still came out well because `sqrt(1 + 3k)` is not sensitive to `k`.

The direct answer is more sensitive. By the formula, the pool values call for 42 of 600 records (7.0%) on ClinicalTrials.gov and 67 (11.2%) on NHANES to come within 25% of the full-audit width, against 384 and 379 when discarding recorded values (`pool-k.csv`; exact values 41.50, 66.25, 383.00 and 378.20, rounded up because fewer records would miss the target). Pilots gave a median of 21 with a 5th to 95th percentile range of 2 to 114 on ClinicalTrials.gov at 60 records, 29 (8 to 111) at 150, 58 (32 to 117) on NHANES at 60, and 63 (40 to 99) at 150 (`plans.csv`, `answer_records`).

**The two forms of `k`.** The theory's ratio has two readings when residual variance is not constant, and this study separates them (`frozen-comparison.csv`, with the frozen values from `real-data-study/summary.csv`).

| Dataset, frozen audit | Implied by frozen widths | Working-model `k`, pool | Implied by frozen RMSE | Empirical `classical`, pool |
|---|---:|---:|---:|---:|
| ClinicalTrials.gov 5%, 10%, 25%, 50% | 0.049, 0.042, 0.044, 0.044 | 0.042 | 0.040, 0.023, 0.027, 0.010 | 0.028 |
| NHANES 10%, 25%, 50% | 0.073, 0.068, 0.070 | 0.070 | 0.089, 0.103, 0.092 | 0.088 |

The working-model value matches what the intervals report. The empirical value is closer to how much the true error grew, in six of seven cells, though the RMSE ratios are noisy, particularly in the 200-cohort cells. So on NHANES the reported intervals grow less than the true error does, and on ClinicalTrials.gov they grow more. This fits the frozen study's coverage but was not designed as a test of it.

**Discarding recorded values.** The formula understates these widths at small audits, because the analysis there rests on few references (`frozen-comparison.csv`).

| Frozen audit | Frozen width over oracle | Formula with pool `discard_k` | Understated by |
|---|---:|---:|---:|
| ClinicalTrials.gov 5%, 10%, 25%, 50% | 5.51, 3.61, 2.14, 1.46 | 4.46, 3.15, 1.99, 1.41 | 19.0%, 12.7%, 6.7%, 3.6% |
| NHANES 10%, 25%, 50% | 3.57, 2.12, 1.42 | 3.10, 1.97, 1.40 | 13.1%, 7.2%, 1.8% |

The direction of the comparison is unaffected: keeping the recorded value stays far ahead.

## Boundary cells

Synthetic, 300 cohorts each, descriptive only (`summary.csv`, `k-summary.csv`). The last line is from `targets.csv`, with the known slope of 0.8. **It was not in the plan.**

| | `missed_flags` | `error_gross_and_small` |
|---|---:|---:|
| Median error, 60-record pilot: targets 5%, 25%, 50% | 21.5%, 8.7%, 6.0% | 13.7%, 16.4%, 19.4% |
| Median error, 150-record pilot: targets 5%, 10%, 50% | 23.9%, 14.3%, 3.7% | **741.9%**, 20.5%, 7.1% |
| Largest 90th percentile | 73.8% | 1,474.7% |
| Median `k`, 60 and 150 records | 0.307, 0.299 | 0.506, 0.526 |
| `gross_errors` warning | 0 of 600 pilots | 600 of 600 pilots |
| Linear MI coverage of the slope at 5%, 10%, 25%, 50% | 79.7%, 87.7%, 91.7%, 93.3% | 32.0%, 53.7%, 83.0%, 90.7% |

The gross cell failed as expected, and the warning caught every pilot. The median ratio of actual to full-audit width at the 5% target was 0.316, so linear MI there reported far too much precision while coverage collapsed. A likely mechanism, not separately tested, is that a small audit without a gross error imputes very large values on unaudited gross records, and those high-leverage rows shrink the interval.
Width accuracy in that cell means little, since the interval itself is wrong.

`missed_flags` has larger errors, so `k` is about 0.3, and the planner is less accurate there than on real data. It is outside the 10% limit at the 5% and 10% targets.

## Predictions, scored

| Prediction | Result |
|---|---|
| 1. NO-GO, because the 5% target stays above 10% with both pilots in at least one dataset | **Held**, in both datasets. |
| 2. Every other real-data cell at or below 10% | **Held.** The largest was 7.8%. |
| 3. Pool `k` from 0.03 to 0.05 on ClinicalTrials.gov and 0.06 to 0.08 on NHANES, close to values implied by frozen widths | **Held**: 0.0418 and 0.0698, against 0.042 to 0.049 and 0.068 to 0.073. |
| 4. `classical` below `k` on ClinicalTrials.gov and above on NHANES, and tracks RMSE inflation better | **Held**: 0.028 against 0.042, and 0.088 against 0.070, and closer in six of seven cells. |
| 5. The `k` interval contains pool `k` in under 90% of pilots | **Held**, which is a failure of the interval: 53.7% to 74.3%. |
| 6. The predictive range contains the actual ratio more often than the `k` interval in every real-data cell | **Held** in all 12 cells, though it remained below nominal. |
| 7. Gross warning on at least 95% of gross-cell pilots and at most 1% in real and `missed_flags` cells | **Failed in part.** 100% in the gross cell, but 19 of 1,500 (1.27%) on ClinicalTrials.gov at 150 records. It was 9 of 1,491 (0.60%) at 60, and none on NHANES or `missed_flags`. |
| 8. The gross cell gives the largest errors | **Held**: a median of 741.9%. Its 60-record pilot at the 5% target was not unusually bad (13.7%). |
| 9. The discard calculation understates frozen widths at 5% and 10% by at least 5% | **Held**: 12.7% to 19.0%. |

## Failures and limits

- The verdict is NO-GO. That is the result, and nothing was retried or retuned after viewing. `R/audit-plan.R` and the plan match their hashes in `source-hashes.csv`, and the frozen copies are in `v0.4.1-source/`.
- A 60-record pilot does not determine `k` well on exact-or-different data, and the bootstrap interval for it undercovers. The advice to plan on the upper limit lowers the risk but does not remove it.
- The planner does not predict one audit's width to within 10% at 30 records, and its predictive range falls short of nominal.
- The planner cannot see rare errors a pilot missed. On these data, that shows up as a pilot `k` below the pool `k`.
- Two public datasets at one cohort size, with `m = 20`. Neither dataset is a trial's source-data verification, and the ClinicalTrials.gov reference is sponsor-entered. Resampling a pool says nothing about transport.
- The development split shares its source with the evaluation split, and the frozen tables were seen before this work began. This test is prospective for the planner's code, and not fully blind for the pools.
- Choices made on development data, listed in the plan, could still be tuned to those two datasets.

## What this changes

- Do not present the planner as validated. It is a documented calculation with a recorded NO-GO.
- The descriptive pattern is that widths at 10% or more were predicted within 7.8% (median) and 20.6% (90th percentile), and 5% was not. This pattern was not the prespecified test.
- The planning conclusion from the frozen study holds, and the pool values make it concrete: keeping recorded values cuts the audit needed to get within 25% of a full audit from about 380 of 600 records to about 40 to 70 on these data.
