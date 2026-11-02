# catmech 0.7.0

- `catmech_pool_lm(cluster = )` now uses a delete-one-cluster jackknife variance (`cv3`) by default, with complete-data degrees of freedom of clusters minus one. Intervals with `cluster` get somewhat wider, and `cluster_variance = "cr1s"` reproduces the earlier default. Calls without `cluster` are numerically unchanged.
- A prospective test, with its plan committed before its runner existed, chose this under a rule fixed in advance. On 11 real items at a 25% audit, with the random-intercept option, mean distance of coverage from 95% went from 3.18 to 1.47 points, and widths rose by up to 15%. It passed narrowly, meeting two conditions exactly. CR2 with Satterthwaite degrees of freedom failed through imputation, with widths up by as much as 95%, although it did well with every row checked. See `validation/VARIANCE_MI_SYNTHESIS.md`.
- `catmech_pool_lm()` now warns when a fit from `catmech_reference(cluster = )` has the audited value as a predictor. In new simulated cells, that use biased the slope toward zero by up to 16%, with coverage as low as 84.4%. Use `cluster` in the fit only when the audited value is the outcome. This was not known before, and the README's earlier example showed the harmful use.
- `catmech_reference(cluster = )` stays experimental, and earlier verdicts stand, because point estimates do not change. No package dependency was added.

# catmech 0.6.0

- `catmech_reference(cluster = )` now adds a random intercept for the cluster, with its variance components drawn in each imputation, instead of resampling audited clusters. This removes the flaw in the earlier version, where resampling within strata understated variance. Results for calls with `cluster` change. `cluster_strata` is removed, so put strata in `predictors`. This is ordinary multilevel imputation in closed form, and not a new method.
- A second prospective test, with its plan committed before either candidate was in the package or had run on a real item, came out **PARTIAL** under a rule fixed in advance. Simulations with clustered residuals selected this candidate over a stratified jackknife draw, by a tie-break on error. On a Rwanda facility recount at a 25% audit, coverage was within its limits on eight of nine items and 89.0% on the ninth. Against the option without `cluster`, error was 5% to 25% lower on six common items and 4% to 5% higher on three rare ones. Both accuracy conditions failed, so the option stays experimental and is not a recommended default. It would also have failed the first test's coverage condition, which the plan replaced before any real result with these candidates. See `validation/CLUSTERED_FIX_SYNTHESIS.md`.
- Calls without `cluster` are numerically unchanged: 18 of 18 runs match frozen v0.3.1, and 12 of 12 match frozen v0.4.0. No package dependency was added.

# catmech 0.5.0

- Add `cluster` and `cluster_strata` to `catmech_reference(use_recorded = TRUE)`, which resample audited clusters in each imputation, and `cluster` to `catmech_pool_lm()`, a CR1 cluster-robust variance pooled with complete-data degrees of freedom of clusters minus one. These combine established ideas and are not a new method.
- A prospective test, with its plan committed before the option existed, came out **PARTIAL** under a rule fixed in advance. On five held-out items of a Rwanda facility recount at a 25% audit, coverage moved closer to the oracle's on all five and came within 5 points on four. The accuracy condition failed, met by two of five. The gain came where linear imputation undercovered, and error rose elsewhere. `catmech_reference(cluster = )` is therefore experimental and is not a recommended default. An independent AI review, after the test, confirmed the verdict and found no bug. It also found that the losses fell on rare events whose residuals show no within-facility correlation, which is untested, and that resampling within strata understates variance when few clusters per stratum are audited. See `validation/CLUSTERED_SYNTHESIS.md`.
- Calls without `cluster` are numerically unchanged: 18 of 18 runs match frozen v0.3.1, and 12 of 12 match frozen v0.4.0. No package dependency was added.
- Credit the real-audit data sources in the README, and state how they were copied in `validation/real-audit/README.md`.

# catmech 0.4.1

- Add `catmech_audit_plan()`, which estimates from a random pilot how much narrower one coefficient's interval would be at each audit fraction, keeping or discarding the recorded value, and how
  many records reach a chosen width. It applies classical two-phase theory, evaluated in closed form under the working models. This is not a new result.
- A prospective test on held-out real data returned **NO-GO** under a rule fixed in advance. Median error in predicting width was within 10% at target audits of 10%, 25% and 50% (7.8% at worst), and 11.3% to 16.8% at a 5% audit. A 30-record audit's width
  is too noisy to predict that closely. The bootstrap interval for `k` undercovered the pool value, and small pilots estimated `k` poorly. The planner should not be described as validated. See `validation/AUDIT_PLAN_SYNTHESIS.md` and `docs/AUDIT_PLANNING.md`.
- Existing workflows are numerically unchanged: 18 of 18 runs match frozen v0.3.1, and 12 of 12 match frozen v0.4.0. No package dependency was added.
- A later test on real audit data, with no package code changed, used a Rwanda health-facility recount and a psychology effect-size audit. The contamination model's advantage over linear imputation was **inconclusive** under a rule fixed in advance: it often had lower point error, but its intervals were too narrow. Linear imputation's coverage held up better. See `validation/REAL_AUDIT_SYNTHESIS.md`.

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
