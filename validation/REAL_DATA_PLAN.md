# Prospective real-data audit study

Written 2026-11-01 before downloading either full dataset. The only earlier look was an
in-memory feasibility probe of 300 ClinicalTrials.gov records, in which 87.0% of enrollment
counts equaled the participant-flow total. NHANES values had not been read.

## Why this design

No public record-level clinical source-data audit was found, so this is not a clinical audit.
It uses two large public datasets where every record has both a recorded value and a reference
value. Hiding most references mimics an audit, and the whole dataset gives the known answer.
Recording errors are the real ones. Nothing is injected.

Cohorts of 600 are drawn with replacement from the evaluation pool, so the pool acts as the
population and its coefficient, computed with references on every row, is exactly the target.
Each cohort's references are then hidden according to the audit design. The pool regression is a
descriptive projection and may be misspecified. Homoskedastic intervals may then be wrong for
every method, including the full-reference oracle, so coverage is judged relative to the oracle.

## Datasets

Each record is assigned by SHA-256 of its identifier to a development split (hash mod 5 = 0) or
the evaluation split. The development split is used only for engineering checks.

**ClinicalTrials.gov** (API v2): interventional studies with posted results, completion date from
2015-01-01 to 2022-12-31, actual enrollment.

- Recorded: log of the registered enrollment count.
- Reference: log of the first-period participant-flow Started total.
- Analysis: `log_duration ~ log_n + phase + sponsor + log_sites + start_year`, where duration is
  start to primary completion in months. The estimand is the log-log slope of duration on size.
- Groups: phase (early, phase 2, phase 3 or 4, not applicable), sponsor (industry, other).
- Excluded: nonpositive counts, nonpositive duration, no listed location, unparseable dates.

Both values are entered by sponsors, and enrolled and started can legitimately differ. The reference
is the results-section count, not independently verified truth. This tests handling of real
discrepancies with an exact-or-different structure.

**NHANES 2017-18** (CDC; DEMO_J, BMX_J, WHQ_J, BPX_J): adults aged 20 or older, not pregnant, with
measured weight, height and at least one systolic reading.

- Recorded: self-reported weight (WHD020) in kg; codes 7777 and 9999 are missing and stay missing.
- Reference: measured weight (BMXWT).
- Analysis: `sbp ~ weight + age + sex + race + height`, using the mean of available systolic readings.
- The survey design is ignored because the pool is only a finite population for resampling. Nothing
  here estimates a U.S. population quantity.

Self-report is rounding and noise, not exact-or-wrong, so this dataset is the deliberate control.

## Cells

Each dataset has `audit_25pct` (1,000 cohorts), `audit_10pct`, `audit_50pct` and `targeted_audit`
(200 cohorts each). ClinicalTrials.gov also has `audit_05pct` (200 cohorts). The targeted design
selects with probability 0.60 when the recorded value is more than 1.25 SDs from the cohort mean,
otherwise 0.15, and also 0.15 when no value was recorded. The analyst knows these probabilities.

## Methods

All receive the same data: `oracle` (reference on every row), `naive_recorded` (recorded value used
as truth, with normal imputation only where it is missing), `gold_only`, `reference_normal`,
`reference_pmm`, `reference_pmm_abb`, `cole_mi`, `raking`, `contamination`, `contamination_strata`.
Prediction-powered inference is not run separately; raking is the same design-based family and
serves as its representative here.

- `contamination`: one corruption rate. ClinicalTrials.gov `agreement_tol = 0`; NHANES `agreement_tol = 1` kg.
- `contamination_strata`: corruption strata by sponsor (ClinicalTrials.gov) or sex (NHANES).
- `contamination_no_tol`: NHANES `audit_25pct` only, `agreement_tol = 0`. Expected to warn.

`catmech_holdout_check()` runs on the first 40 cohorts of `audit_25pct` for each dataset.
`catmech_recording_audit()` also runs once on each whole evaluation pool, by the listed groups.

Seeds: cohort 200,100,000 + replicate, selection 210,100,000 + replicate, imputation 220,000,000 +
replicate; NHANES adds 5,000,000. Engineering pilots use the development split and bases 300,100,000,
310,100,000 and 320,000,000.

## Predictions fixed in advance

1. ClinicalTrials.gov, 25%: `contamination` coverage within two points of the oracle, and narrower
   intervals than `reference_normal`, `cole_mi` and `raking`.
2. NHANES, 25%: `contamination` has no clear advantage over `cole_mi` or `raking`, because the
   exact-or-error form does not describe self-report. A win would be reported, but not as support
   for that form.
3. In both datasets the errors are small relative to the spread of the true values, so
   `naive_recorded` may be close to the oracle. If so, the honest conclusion is that this audit
   finds little to correct, and the value of the model is limited to losing less precision than
   discarding recorded values.
4. `targeted_audit`: `reference_normal` loses coverage because selection depends on W, while
   `contamination` and `raking` hold it.
5. `reference_pmm_abb` covers better than `reference_pmm` in both datasets.

Each prediction is scored in the write-up, including those that fail.

## Interpretation limits

These are real recording differences in two public datasets, neither a trial's source-data
verification. Good results here do not validate clinical use, and resampling one pool does not
address transport to other datasets. ClinicalTrials.gov changes over time, so the pool is frozen
as a dated CSV with its query and hash. Results are not retuned after viewing. A change requires a
new versioned run.
