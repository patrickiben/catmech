# Prospective test on real audit data

Written 2026-11-01 before any imputation, raking or contamination model was run on these data.

## What has already been seen

These earlier looks shaped the analyses and predictions below, so this test is not blind to how much the errors matter:

- Agreement counts between recorded and reference values in both datasets.
- Reproduction of the published Rwanda table and the published effect-size counts.
- Plain OLS slopes with reported or recount values on both sides, hospitals included: PNC 0.631 against 0.562, syphilis 0.568 against 0.665, iron 1.020 against 0.954, ANC 4th visit 0.435 against 0.284. Effect-size slopes moved under 0.2 SE.

No correction method has been compared on these data. The pools are too small to split, so the engineering pilot uses the same pools with separate seeds, and only its request statuses and timing are read.

## Why this test

The earlier real-data study used recording differences, not audits. Its reference values were not checks of the recorded values, and its errors were small relative to the target. This test uses two real audits in which someone went back to the source. Nothing is injected.

Rwanda is the main test because its errors are large relative to the target, the case where the documentation says the contamination model is strongest. The effect sizes are a control where little should need correcting.

**Not a clinical source-data verification.** Rwanda holds facility counts, not patient values. The effect sizes come from psychology. Good results here do not validate clinical use.

## Data

**Rwanda data verification** ([PLOS ONE 2020](https://doi.org/10.1371/journal.pone.0235823), S1 Dataset, CC BY). Two trained data collectors recounted 14 maternal and newborn items in facility registers and compared them with HMIS, the national health information system, for 76 facilities over three months in 2017. The analysis uses the 69 health centers, because hospitals are few, far larger and do not report antenatal items.

Limits stated by the source or visible in it: the facilities were selected for a scale-up program and are not a random sample, the recount was done in consultation with facility reporting staff, and a register can itself be wrong. For ANC standard visits, register entries and reporting rules do not always match, so part of that gap may reflect how gestational age was handled rather than a recording slip.

| Analysis | Recorded and reference target | Covariate, as reported | Slope means |
|---|---|---|---|
| `pnc` | Newborn PNC within 24 hours | Deliveries | PNC visits per delivery |
| `syphilis` | Syphilis tests at first ANC visit | ANC new registrations | Tests per registration |
| `iron` | Iron and folic acid courses | ANC new registrations | Courses per registration |
| `anc4` | ANC 4th standard visits | ANC new registrations | 4th visits per registration |

Each analysis is `target ~ covariate + catchment`, using rows where target and covariate are all observed. The covariate stays as reported because this package handles one audited target, so the estimand is the slope a recount of that one item would give. `pnc`, `syphilis` and `iron` are the exact-or-wrong analyses. `anc4` has few exact matches and is reported as a stress case outside the decision on the contamination model.

Each analysis runs at two levels:

- **Facility-months** (primary): up to three rows per facility, which breaks the package's independent-records assumption.
- **Quarter totals** (sensitivity): one row per facility, summed over its three months, for facilities with all three months observed, so records are independent but few.

**Effect sizes** ([Maassen et al. 2020](https://doi.org/10.1371/journal.pone.0233107), [OSF](https://osf.io/7nsmd/), CC0). Coders recomputed 500 effect sizes sampled from 1,978 in 33 psychology meta-analyses. Only `info` 0 (reproducible) and 1 (incorrect) are used, 350 rows with an independent reference. Copied values have no reference, and ambiguous ones took the candidate closest to the reported value.

| Analysis | Rows | Analysis | Agreement tolerance |
|---|---|---|---|
| `effect_smd` | `g` or `d`, meta-analyses with at least two rows | `effect ~ log(n) + meta` | 0.049, the paper's margin |
| `effect_r` | `r` or `z`, meta-analyses with at least two rows | `effect ~ log(n) + meta` | 0.025, the paper's margin |

## Design

Each pool acts as the population. Its coefficient, using reference values on every row, is the target.

1. **Draw a cohort.** Rwanda resamples facilities with replacement within catchment, keeping each facility's rows together. Effect sizes resample rows with replacement within meta-analysis.
2. **Hide references.** Rwanda audits a stratified random sample of facilities within catchment. Effect sizes audit a stratified random sample of rows within meta-analysis. The count in each stratum is `max(1, round(fraction * size))`, and the analyst knows these probabilities.
3. **Run every method** on the same cohort and audit.

Fractions are 25% and 50%. The same cohort is used at both fractions.

| Level | Cells | Cohorts per cell | Agreement tolerance |
|---|---|---|---|
| Rwanda facility-months | 4 analyses, 2 fractions | 500 | 0 |
| Rwanda quarter totals | 4 analyses, 2 fractions | 200 | 0 |
| Effect sizes | 2 analyses, 2 fractions | 500 | As above |

If the pilot predicts more than 45 minutes, these counts may be lowered to no fewer than 200 per cell, before the study run and recorded in `settings.json`. Monte Carlo precision is not the main limit here. With 69 facilities, the pool is.

Seeds: cohort 400,100,000 + analysis index x 1,000,000 + replicate. Selection 410,100,000 and imputation 420,000,000, each + analysis index x 1,000,000 + fraction index x 100,000 + replicate. Pilots use bases 500,100,000, 510,100,000 and 520,000,000.

## Methods

| Method | What it does |
|---|---|
| `oracle` | Reference values on every row: what a full audit of this cohort would give |
| `recorded` | Treats reported values as true, no audit used |
| `audited_only` | OLS on audited rows, discarding reported values elsewhere |
| `linear_mi` | `catmech_reference(use_recorded = TRUE)`, 20 imputations |
| `raking` | Generalized raking on an intercept and the slope's influence function, from a proxy built with the recorded value, using the known probabilities |
| `contamination` | `catmech_contamination()`, one corruption rate, 20 imputations, package defaults |

## Measures

- **Recovery RMSE** (primary): root mean square of the estimate minus the same cohort's `oracle` estimate. It answers the audit question directly: how close does this get to a full audit? It removes the sampling noise every method shares.
- **RMSE and bias against the pool coefficient**, for comparison with the earlier study.
- **Coverage of the pool coefficient**, judged against the oracle's, under two intervals. Default intervals assume independent rows, as in the earlier study. Clustered intervals use a CR1 sandwich by facility with Rubin's rules, or the clustered raking variance. Effect sizes cluster by row.
- **Returned intervals** per request, and refusals.
- **Per-record scores**: AUC and Brier for the contamination model's corruption probabilities among unaudited rows.

Ratios of recovery RMSE get a 95% percentile interval from 2,000 paired resamples of cohorts.

## Decision rules

Rules use facility-months at 25% unless stated.

1. **Does a partial audit pay?** For each analysis, yes if the upper limit of `linear_mi` against `recorded` is below 1.
2. **Does the contamination model add anything?** Using `pnc`, `syphilis` and `iron`:
   - **Supported** if the upper limit of `contamination` against `linear_mi` is below 1 in at least two, and in each of those it returns at least 90% of intervals, its clustered coverage is within 3 points of the oracle's, and its quarter-total ratio is also below 1.
   - **Contradicted** if the lower limit is above 1 in at least two.
   - **Inconclusive** otherwise.
3. **Control.** The effect sizes report what auditing costs when there is little to correct: each method's recovery RMSE against `recorded`.

## Predictions fixed in advance

1. `linear_mi` beats `recorded` at 25% on `anc4`, where errors are systematic, and does not on `pnc`, where errors are balanced in direction.
2. `contamination` has lower recovery RMSE than `linear_mi` at 25% on `pnc`, `syphilis` and `iron`, and does not on `anc4`.
3. Default intervals undercover for every method at facility-months, including the oracle. Clustered intervals are closer to 95%, and `contamination` stays within 3 points of the oracle.
4. On both effect-size analyses at 25%, `recorded` has the lowest recovery RMSE of all non-oracle methods.
5. `raking` has wider intervals than `linear_mi` at 25%.
6. Per-record AUC stays below 0.75 at facility-months, so the model still cannot reliably say which record is wrong.
7. Every method that uses the audit does better at 50% than at 25%.

Each prediction is scored in the write-up, including those that fail.

## Pilot and amendments

Two pilots of 14 cohorts per cell ran with the pilot seeds, about 47 seconds each, predicting about 21 minutes for the study, so cohort counts stay as planned. Only request statuses and timing were read. Two changes followed, both before any estimate was seen:

- **Raking auxiliaries.** Calibrating on every coefficient's influence function, as the earlier study did, failed on 61% of effect-size requests and 30% of quarter-total requests at 25%, because each stratum added a column supported by very few audited rows. Leaving that in would have handicapped raking, so it now calibrates on an intercept and the slope's influence function. That pilot returned 279 of 280 raking requests.
- **Prediction 5.** It also predicted that raking would return fewer intervals than linear MI. The pilot statuses revealed that, so the clause was removed.

The pilot also showed `contamination` convergence refusals at 25%: 5 of 56 at facility-months and 11 of 56 at quarter totals. Package defaults stay in place, so these count against its returned intervals. [Statuses](real-audit-pilot-check.txt) are saved.

## What the result changes

If rule 2 is supported, the documentation can say a real audit showed an advantage on exact-or-wrong errors. If it is contradicted or inconclusive, the claim that the model is strongest when errors are large stays synthetic-only, and the documentation must say so.

Results are not retuned after viewing. A change requires a new versioned run. No package code changes as part of this test.
