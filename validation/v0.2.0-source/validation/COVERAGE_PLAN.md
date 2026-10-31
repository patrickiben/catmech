# Prospective coverage evaluation for the model workflow

Protocol written 2026-10-31, before running the new workflow. The v0.1 benchmark
remains a historical record. New model changes are evaluated on disjoint data
seeds; the old results are not overwritten or relabeled.

## Question and rationale

Can masking flagged measurements and using a coherent conditional normal model
remove the large biases seen with arbitrary retain/replace probabilities and
forced incorrect categories? Does using a declared category likelihood help when
that likelihood is correct, and fail or refuse when it is wrong?

The original uncertain-flag failure had wide intervals already. This evaluation
therefore requires checking point-estimate bias and interval width as well as
coverage. Widening intervals to hide bias is not the proposed repair.

The estimand is the laboratory coefficient in the prespecified linear outcome
model, 0.8 except in the null control where it is 0. All imputation methods receive
the downstream outcome as a predictor. Hidden true values and corruption labels
are excluded from imputation inputs.

## Connections that change the evaluation

- **A False Floor:** retain the strongest simple PMM and normal-MI alternatives.
  Test both a null control and biased alternatives: null calibration alone cannot
  detect attenuation of a nonzero effect.
- **Residual carryover work:** decompose observed-data error into selection of
  retained observations and residual measurement error using paired synthetic
  oracle analyses. These are diagnostic calculations, not available estimators.
- **Calibrated clinical ML:** keep engineering/pilot and held-out seed streams
  separate. A marginal flag error frequency is not a calibrated conditional
  corruption probability. No synthetic confusion matrix is claimed to have been
  estimated or validated on real clinical data.
- **Anchor Audit:** retain refusals and unresolved MCMC runs in the denominator.
  Assumption-path stability is conditional on the scenarios actually evaluated;
  no finite envelope receives a universal coverage guarantee.

## Fixed data-generating settings

Primary synthetic settings use n=300 and 20 imputations. Thresholds are
`(-Inf, 5.5, 7, Inf)`. The original simulator is retained without retuning its
parameters. `null_effect` uses the uncertain-flag data mechanism but an independently
seeded standard-normal outcome, so the true laboratory coefficient is exactly zero.

| Scenario | Missing rate | Flag rate | Corrupt fraction among flags | Category error | Latent residual | Role |
|---|---:|---:|---:|---:|---|---|
| missing_only | .20 | 0 | 1 | 0 | normal | Primary |
| perfect_flags | .20 | .10 | 1 | 0 | normal | Primary |
| false_positive_flags | .20 | .10 | .50 | 0 | normal | Primary |
| category_error | .20 | .10 | 1 | .10 | normal | Primary |
| null_effect | .20 | .10 | .50 | 0 | normal | Negative control |
| heavy_tails | .20 | .10 | 1 | 0 | standardized t(3) | Model stress |

Run 1,000 independent data sets in each setting for the direct workflows and
baselines. At 95% coverage this yields about .0069 binomial Monte Carlo standard
error (MCSE). Use the first 250 prespecified replicates for applicable Gibbs
models; their approximately .0138 coverage MCSE supports a narrower conclusion.
For the exact-category model in `category_error`, run all 1,000 requests because
trusted-value contradictions should be rejected before MCMC. A refusal is not
an observed confidence interval and must remain visible in requested denominators.

Engineering pilot seeds are `7,000,000 + 100,000 * scenario_index + replicate`.
Held-out seeds are `20,000,000 + 100,000 * scenario_index + replicate`.
Imputation RNG seeds add 1,000,000. No confirmatory seed appeared in the v0.1
benchmark. The five-replicate engineering pilot checks runtime, API integration
and numerical diagnostics, not whether a favorable scientific result can be found.
The planned MCMC settings are four chains, 200 burn-in and 400 retained iterations
per chain. Failed convergence is recorded; it is not silently retried until passing.
Any engineering change after a pilot must be dated before the held-out run.

## Comparators

1. Full-data oracle regression using hidden truth.
2. Naive observed-value regression, retaining corrupted values.
3. Trusted observed-value regression, dropping all flagged values.
4. Trusted PMM: mask all flagged values, use covariates, outcome and category.
5. Trusted normal MI: mask all flags, use covariates and outcome; invoke the
   established `mice.impute.norm` independently of the new implementation.
6. New model with categories ignored: direct normal posterior predictive draws.
7. Legacy sensitivity workflow: q=1 except q=.5 under uncertain flags and null;
   force categorical conflicts as in the original exploratory comparison.
8. New exact-category model for all scenarios on its Gibbs subset, with explicit
   refusals for contradictions among trusted observations.
9. New noisy-category model for `category_error`, with declared matrix H having
   .90 on the diagonal and .05 off diagonal (rows=true bin, columns=observed bin).
10. Wrong-H stress on the same noisy data, using .98 diagonal and .01 off diagonal.

The independent normal-MI comparator uses a complete predictor matrix with no
categories and a standard one-target posterior draw; PMM includes category
indicators. This prevents attributing simple removal of contamination to a new
statistical invention. Category likelihood models mask every flagged W. They do
not estimate the corruption mechanism or confusion matrix from undocumented data.

## Additional clinical-format stress fixtures

Use `tools/lab-fixture.R`, its `heldout` stream, 250 replicates per scenario,
n=600, and the same direct methods except the legacy constraint engine. Scenarios
are clean, flagged_errors, missed_flags, site_shift and rounded_records. The
clinical estimand is the coefficient 0.8 in
`outcome ~ LBSTRESN + age + sex + inflammation + SITEID`.

The fixture's NORMAL interval includes both clinical reference endpoints. The
current category engine uses globally right-closed intervals, so do not silently
apply that convention to these clinical labels. Run the ignore-category and
flag-aware PMM/normal workflows; category-constrained models are outside this
fixture evaluation. `site_shift` is a fixed shared measurement shift at a site,
not a validation of confidence intervals under random cluster-correlated outcomes.
The rounded-record fixture also violates exact numeric measurement assumptions.
These results diagnose transfer limits; they are not clinical validation.

## Performance, failures and diagnostics

Store every request with status `ok`, `refused`, `convergence_failure` or `error`.
Record the maximum split R-hat and minimum bulk and tail ESS supplied by the model, and
keep warning text. Gibbs screening requires R-hat at most 1.05 and both bulk
and tail ESS at least 100 times the number of chains (400 here). These screens
do not prove convergence. The iteration count was raised from the initial design
of 200 to 400 before the first engineering pilot, to support the stronger ESS
screen; no coverage result prompted the change. Report requested count, returned interval count, refusal rate,
convergence-failure rate and other-error rate. Coverage conditional on a returned
interval is distinct from the fraction of all requests returning a covering
interval. The latter is a service-yield measure, not a confidence-level estimate.
Use Wilson intervals for rates, including 0% and 100%.

For successful analyses report bias, RMSE, empirical standard deviation, mean
standard error, mean interval width, nominal 95% coverage and paired differences
against trusted PMM and trusted normal. MCSEs use replicate-level standard errors;
RMSE uses the delta method. Comparisons share the same generated data. Retain
complete-pair counts and do not claim superiority for gaps smaller than two MCSEs.

A retrospective bias-centered diagnostic shifts each successful confidence
interval by that method/scenario's Monte Carlo mean error, then counts coverage.
This distinguishes location error from dispersion error. It uses simulation truth
and is **not** an implementable correction, a recommended confidence interval, or
reported inference. Its uncertainty is estimated by 400 bootstrap resamples of
whole replicate records, recomputing the centering within each resample.

For each generated data set, also fit the full truth, truth restricted to retained
observed unflagged records, and their actual recorded values. Record:

`selection difference = beta(truth, retained rows) - beta(truth, all rows)`

`residual-error difference = beta(recorded, retained rows) - beta(truth, retained rows)`

Their sum equals the trusted observed-data error relative to the sample oracle.
This decomposition says nothing about an unobservable decomposition in real data.

## Interpretation commitments

No requirement forces a method to pass. A coverage Wilson interval wholly within
[.925,.975] in a primary 1,000-replicate setting is evidence of approximate nominal
coverage for that setting; otherwise report the discrepancy or uncertainty.
This tolerance is a development criterion, not an acceptable clinical error rate.
Also inspect bias, width and RMSE; a coverage result alone cannot establish quality.
A calibrated null with biased nonzero effects is a failure of a broader validity
claim. The Gibbs subset and deliberate wrong-model stresses do not establish
universal calibration. Do not choose or recalibrate a method using held-out outcomes.

Reproducibility includes the script, fixture, package-source and metadata hashes,
settings, source-stability checks and session information. Existing output folders
must not be overwritten. A failed or interrupted run is retained as a checkpoint
and must not be relabeled as a completed validation.

## Engineering implementation record

The first engineering pilot executed all requests in 11 seconds with no numerical
or MCMC failures and five expected category-conflict refusals. Its source-stability
check rejected the run because development edits occurred concurrently; that
checkpoint is not treated as verified scientific evidence. No accuracy results
were used to select settings. Before held-out execution, independent replicate
workers were added to reduce runtime. Every data and imputation seed is explicit;
record order is deterministic. The default is four workers (one on Windows),
configurable through `CATMECH_BENCHMARK_WORKERS`. A one-worker/four-worker pilot
comparison must verify equal scientific output before the held-out run.
