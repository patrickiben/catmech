# Prospective test of the audit planner

Written 2026-11-01 before any run of the planner on the evaluation split. The benchmark records this file's hash
before it starts and again when it ends. It is not edited afterward.

## Question and method

`catmech_audit_plan()` answers "how many records should we audit?" for linear MI with the recorded value,
`catmech_reference(use_recorded = TRUE)`. From a random pilot it estimates `k` and predicts interval width at other audit
fractions:

    width(f) = width_full * sqrt(1 + k * (1 - f) / f)

This is classical two-phase theory, not a new result. The test asks whether it is accurate enough to plan with.

## What had already been seen

These looks came before this test and cannot be undone, so they are stated here:

- The frozen study tables were already seen by hand. `validation/REAL_DATA_SYNTHESIS.md` and
  `validation/real-data-study/summary.csv` report mean widths on the same evaluation pools. The steady `k` of about 0.04 to 0.05
  on ClinicalTrials.gov and about 0.07 on NHANES was worked out by hand from those tables and motivated this work.
- While building, I also compared those tables' RMSE values with development-split results by hand. That suggested the
  empirical form of `k` tracks sampling variance rather than interval width (prediction 4 below).
- The cohorts below are new draws from the same evaluation pools as that study, so this test is not independent of those pools.

No row of the evaluation split was read by the planner code or by any development script. Everything below labeled
development came from the development split, using scratch scripts and `validation/audit-plan-pilot`.

## What development decided

Each choice was made on the development split and is fixed now:

| Choice | Reason on development data |
|---|---|
| `k` evaluated under the working models, in closed form, rather than the empirical ratio of influence-function variances | Rubin's rules give a model-based interval, and this version tracked actual widths with less noise. Pool values: ClinicalTrials.gov 0.036 against 0.024 for the empirical ratio; NHANES 0.067 against 0.101. The empirical version stays in the output as `classical`. |
| Anchor on the pilot fit's within-imputation variance, not its total width | Total width carries the noise of the between-imputation variance. With a 10% pilot and 25% target on NHANES, median error was 3.3% against 6.5%. |
| Plain formula, with no finite-sample terms for linear MI | Adding them did not lower the error. |
| Bias-corrected and accelerated interval, with `n - 2q` in resampled variances | Percentile intervals sat low. Coverage of the pool `k` at 60 records rose from 50% to 63% on ClinicalTrials.gov and from 66% to 83% on NHANES, still below nominal. |
| Predictive range added beside the `k` interval | The `k` interval describes an expected width, so it often misses the width of one audit. |
| Gross-error warning above 5 robust SDs (pilot discrepancy) or 10 robust SDs (recorded value) | On 40 to 280 development cohorts per cell, it was raised on every pilot in the gross cell and on none in the real or `missed_flags` cells. |

## Datasets and analyses

Frozen pools in `validation/real-data/`, evaluation split, as in `validation/REAL_DATA_PLAN.md`:

- ClinicalTrials.gov: `outcome ~ lab_value + phase + sponsor + log_sites + start_year`, with `lab_value` the log recorded enrollment
  and the reference the log participant-flow Started total. Term: `lab_value`.
- NHANES 2017-18: `outcome ~ lab_value + age + sex + race + height`, with `lab_value` self-reported weight and the reference measured
  weight. Term: `lab_value`.

Predictors for imputation are the analysis covariates and the outcome. Survey weights are ignored, as in that study.

## The decision rule

This rule was fixed by Patrick's reviewer. I did not write it and may not alter it.

> Target: in each cohort of 600, drawn with replacement from the evaluation pool, draw nested random audits. Fit linear MI with
> use_recorded = TRUE at target fractions 5%, 10%, 25% and 50% and get its actual interval width. Predict those widths from a pilot
> audit nested inside the same cohort, using catmech_audit_plan alone, with pilots of 10% (60 records) and 25% (150 records). Exclude
> the case where the target fraction equals the pilot. The error for one cohort is |predicted - actual| / actual.
>
> GO: with the 10% pilot, the median absolute error is at most 10% at every target fraction, in both datasets.
>
> PARTIAL: that fails, but it holds with the 25% pilot in both datasets.
>
> NO-GO: otherwise.

## How the rule is applied

- Each cohort gets one random ordering of its 600 records. The audits are its first 30, 60, 150 and 300 records, so every audit
  contains the smaller ones.
- Actual width is the upper limit minus the lower limit for the term from `catmech_pool_lm()`, with `m = 20`.
- Predicted width is `table$width` from `catmech_audit_plan(fractions = c(.05, .10, .25, .50), bootstrap = 500, m = 20)`, given the
  cohort with references on the pilot rows only. Its warnings do not exclude a cohort. All are recorded.
- With the 10% pilot the targets are 5%, 25% and 50%. With the 25% pilot they are 5%, 10% and 50%.
- If the target fit fails or refuses in a cohort, no actual width exists. That cohort leaves that cell and is counted.
- If the planner fails or refuses where the target fit succeeded, that cohort's error is infinite, so a refusal can only raise a median.
- The verdict uses ClinicalTrials.gov and NHANES only. It is not retried.

## Settings

1,500 cohorts per real dataset and 300 per synthetic cell. The development pilot ran 640 cohorts in 123 seconds at 14 workers,
so this should take about 12 minutes. If timing had said otherwise, I would have lowered the counts here before the run.

Seeds, each plus a cell offset (ClinicalTrials.gov 0, NHANES 2,000,000, `missed_flags` 4,000,000,
`error_gross_and_small` 6,000,000) and a replicate term:

| Stream | Evaluation base | Development base | Replicate term |
|---|---:|---:|---:|
| Cohort | 600,100,000 | 700,100,000 | r |
| Audit order | 610,100,000 | 710,100,000 | r |
| Target imputation | 620,000,000 | 720,000,000 | 2r |
| Plan, 10% pilot | 630,000,000 | 730,000,000 | 4r |
| Plan, 25% pilot | 640,000,000 | 740,000,000 | 4r |
| Target retest | 650,000,000 | 750,000,000 | 2r |

Synthetic fixtures use master 430271 (evaluation) or 430272 (development), and stream bases 670,100,000, 676,100,000, 678,100,000 and
690,000,000 (evaluation) or 770,100,000, 776,100,000, 778,100,000 and 790,000,000 (development). The cohort base is not used for synthetic
cells, which get their cohorts from the fixture.

These are disjoint from every earlier study. The earlier tools use bases at or below 505,100,000 plus replicate terms:
`benchmark.R` (100,000 and up), `coverage-benchmark.R` (7, 20 and 30 million), `extension-benchmark.R` (50 to 90 million),
`contamination-fixtures.R` (170 to 190 million), `real-data-benchmark.R` (200 to 225 and 300 to 325 million), and `merged-audit-pilot.R`
(400 to 405 and 500 to 505 million). Fixture masters used earlier are 390173, 410273 and 20261031.

## Descriptive reports

None of these enters the verdict:

- 90th-percentile error, signed median error and failure counts for each cell.
- How often the `k` interval and the predictive range contain the actual ratio, meaning actual width over the same cohort's full-audit width.
- How often the planner picks the same smallest fraction as the truth for reaching within 25% of the full-audit width, over the grid 5%, 10%, 25%, 50%.
  The truth is the smallest fraction whose fit succeeded with ratio at most 1.25, or none. A refused plan counts as a miss.
- Pool-level `k` from the whole evaluation pool against the pilot estimates, and how often the `k` interval contains it.
- The noise floor: the median error of a predictor that knows each cell's median actual width and nothing about the cohort. A planner has
  no access to this. It shows how much error comes from the width itself.
- The retest error on the first 200 cohorts of each cell: the same audit fitted again with a different imputation seed, predicting the first fit.
  It shows how much comes from imputation noise alone.
- The discard calculation against the frozen study's `reference_normal` widths, and pool-level `k` values against the frozen study's implied values.
- Boundary cells on synthetic data from `tools/contamination-fixtures.R`: `missed_flags`, and `error_gross_and_small`, where failure is expected.

## Predictions fixed in advance

On 280 development cohorts per dataset, the rule gave NO-GO. Median errors were:

| Dataset, pilot | 5% | 10% | 25% | 50% |
|---|---:|---:|---:|---:|
| ClinicalTrials.gov, 10% | 12.0% | | 2.6% | 1.6% |
| ClinicalTrials.gov, 25% | 16.2% | 6.5% | | 0.9% |
| NHANES, 10% | 14.5% | | 3.4% | 3.3% |
| NHANES, 25% | 17.3% | 7.6% | | 1.7% |

So I expect:

1. **Verdict NO-GO**, because the 5% target (30 records) stays above 10% with both pilots in at least one dataset. A 30-record fit is noisy on its
   own: on development data the noise floor there was 16% to 18%, and two fits of the same audit differed by a median of about 5% (ClinicalTrials.gov) and 11% (NHANES).
2. Every real-data cell other than the 5% target has median error at or below 10%.
3. Pool `k` is between 0.03 and 0.05 on ClinicalTrials.gov and between 0.06 and 0.08 on NHANES, close to the values implied by the frozen widths.
4. The empirical `classical` value is below `k` on ClinicalTrials.gov and above it on NHANES, and it tracks the frozen study's RMSE inflation better than `k` does.
5. The `k` interval contains the pool `k` in under 90% of pilots at both sizes in both datasets.
6. The predictive range contains the actual ratio more often than the `k` interval does in every real-data cell.
7. The gross-error warning is raised on at least 95% of pilots in the gross cell and at most 1% in the real and `missed_flags` cells.
8. The gross cell produces the largest errors in the study.
9. The discard calculation understates the frozen `reference_normal` widths at 5% and 10% by at least 5%.

Each is scored in `validation/AUDIT_PLAN_SYNTHESIS.md`, including those that fail.

## Limits

Two public datasets, neither a trial's source-data verification. Good results would not validate clinical use, and resampling a pool says nothing about
transport to other data. The test covers cohorts of 600 only. Results are not retuned after viewing, and any change requires a new versioned run.
