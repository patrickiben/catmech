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
