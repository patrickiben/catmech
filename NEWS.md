# catmech 0.4.0

- Add `catmech_contamination()`, a recording-error model calibrated on reference rows. It uses recorded values without
  trusting them and gives per-record error probabilities. Options cover strata, agreement tolerance, a declared reference
  SD and a selection shift. It is an established contamination model used with validation data, not a new estimator.
- Add `catmech_reference(use_recorded = TRUE)`, linear imputation given the recorded value. It was best or tied among reference-using
  methods at 25% audits on two real datasets and is the recommended starting point in most cases.
- Add `catmech_recording_audit()`, a model-free summary of what references show, and `catmech_holdout_check()`, which
  scores imputations on hidden references without simulation truth and reports when a model shows no skill.
- Add `catmech_reference(method = "pmm_abb")`, resampling donors for each imputation. Synthetic coverage rose from 86.7% to 95.9%, but it
  overcovered on real data. `method = "pmm"` is unchanged and stays gated.
- Existing workflows are numerically unchanged: all 18 comparison runs match frozen v0.3.1 exactly.

Evidence, including failures, is in `validation/CONTAMINATION_SYNTHESIS.md` and `validation/REAL_DATA_SYNTHESIS.md`. The model failed
under heavy-tailed targets, selection on the true value with errors present and undeclared reference noise, and one primary screen failed
narrowly. On real data, linear MI and raking matched or beat it, and it could not identify which records were wrong. `use_recorded` and the hold-out advice
were added after the frozen studies, so those studies did not run them directly; `use_recorded` was checked against the study's linear MI on
the same real-data cohorts.

# catmech 0.3.1

- Require explicit sensitivity opt-in for pooled inference from partial-reference
  PMM after held-out studies found substantial undercoverage. Full-reference
  datasets retain ordinary complete-data inference. No imputation estimates change.
- Reject multivariate-response and coefficient-free analysis formulas with clear
  errors instead of silently returning `NULL`.
- Preserve the exact v0.3.0 numerical study sources and verify successful scalar
  analysis results are unchanged by these guards.

# catmech 0.3.0

- Add `catmech_reference()` for an explicitly measured internal reference subset.
  Preserve supplied reference values and impute every unvalidated target, including
  unflagged recorded values. Normal and PMM options require declared selection
  and reference-accuracy assumptions.
- Add `catmech_pmm()` for standard flag-masking predictive mean matching, with
  complete explicit predictors, original trusted donors and honest ambiguous
  donor-identity auditing. No q retention or forced category intervals apply.
- Distinguish inference labels across Gaussian, predictive-matching, reference
  and legacy sensitivity workflows; refuse unknown workflow engines.
- Correct the fully observed/reference analysis endpoint to ordinary regression
  degrees of freedom. Zero between-imputation variation alone is insufficient.
- Add independent multivariate-t posterior-predictive checks, reference-data
  contracts, and equal-information method comparisons on new held-out streams.
- Preserve v0.2 evidence against its frozen original sources.

These additions require actual reference information or appropriate donor/model
support. They do not detect latent truth from flags, establish novel estimators,
or guarantee robustness to selection, noisy references or misspecification.

# catmech 0.2.0

- Add `catmech_model()` for Gaussian imputation after masking all flagged values.
  Optional exact or fixed-confusion category likelihoods use Gibbs augmentation.
- Preserve trusted numeric values and recorded categories. Reject impossible
  likelihoods, inadequate trusted designs and failed sampler diagnostics.
- Add checked linear-model pooling and finite-scenario sensitivity reporting,
  including omitted planned scenarios and failed fits.
- Retain the legacy q-based PMM sensitivity API without changing its statistical
  semantics. Existing YAML specifications are not converted to the model API.
- Add independent posterior verification, a prespecified held-out study, strong
  normal/PMM comparators and synthetic clinical defect fixtures.
- Preserve historical validation with its original source hashes.

The principal model has stronger inferential foundations, but does not correct
missed contamination, estimate category-error rates or guarantee coverage under
model misspecification. See the validation record for successes and failures.
