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
