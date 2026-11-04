# catmech 0.7.2

- Refuse clustered reference imputation when one site needs stochastic fills both
  with and without recorded values. Separate fits would draw independent effects
  for that site and lose shared uncertainty. This case needs a joint model;
  deterministic agreement and supported single-pattern sites remain available.
- Include formula offsets in the audit planner's analytic response. Equivalent
  offset and adjusted-response formulas now give the same planning calculations;
  nonlinear functions of the audited target remain unsupported.
- Identify compound contamination strata by their full tuples, so delimiter
  characters in labels cannot merge distinct strata or bypass calibration support.
  Unambiguous labels retain their existing ordering and numerical behavior.
- Describe clustered variance-component draws as approximate. Balance alone does
  not make them exact with general covariates. Narrow the rare-count interpretation:
  the experiment varied frailty variance, skewness and dependence together, so its
  historical `HEAVY TAILS` decision label does not identify an isolated mechanism.
- Add regression checks for the three defects. Preserve frozen numerical studies
  and their original decision rules; no new coverage guarantee is claimed.

# catmech 0.7.1

An overall check by five independent AI reviewers found cases where the package gave badly wrong results with no signal. These changes add warnings and refusals. Earlier studies' numbers are unchanged, and no estimate or interval changes in a case that still returns one.

- `catmech_reference(use_recorded = TRUE)` now warns when every audited row agrees and unaudited recorded values are kept with no imputation uncertainty, stating the audit size and the error rate it cannot rule out. `catmech_pool_lm()` repeats this and ends the interpretation in `_recorded_kept`. In a simulated check, such intervals did not cover the truth when errors were large and the audit small.
- `catmech_pool_lm(cluster = )` now warns, naming the coefficients, when a single cluster determines a coefficient, where such intervals covered about a third of the time in a simulated check. It also refuses a pooled variance that is only rounding noise.
- `catmech_reference(cluster = )` now refuses when the predictors absorb every audited cluster, as with one audited cluster per stratum, and warns when it cannot separate the two variances. `spec$cluster_uncertainty` now says what was done. The refusal of fewer than five audited clusters now also applies when every audited row agrees. Blank cluster ids are refused, and an exact offset on audited rows is refused as it is without `cluster`.
- `catmech_sensitivity()` gains `cluster` and `cluster_variance`, so fits can be compared under one variance.
- `catmech_recording_audit()` now flags a difference between audited and other rows only when it exceeds sampling noise, and says selection may not be random. Before, it told about 30% of truly random small audits otherwise.
- `catmech_holdout_check()` no longer warns or fails with an internal message in two edge cases.
- Documents: the recommendation table now lists all of linear MI's recorded failures and qualifies raking, the README's examples run as written, the AI disclosure states what the assistant wrote, and data sources and terms are credited. See `validation/OVERALL_QC.md`.
- Documents: the README now opens with how the package works, with a diagram. `docs/learn/` adds a plain-language guide and a browser audit simulator, a teaching model whose engine `docs/learn/qa/compare.js` checks against the package.
- Later planned tests, with no package code changed: on a second health system, the outcome was MIXED, with errors that barely clustered and little added by the option (`validation/SECOND_SYSTEM_SYNTHESIS.md`); the count-data rule returned HEAVY TAILS, with no fix or warning found (`validation/RARE_COUNTS_SYNTHESIS.md`); the v0.7.2 interpretation correction notes that the experiment changed frailty variance and shape together; comparators written by fresh AI agents did not meet the rule, so the option stays, but the best of them looks better than the option on most measures (`validation/INDEPENDENT_SETUPS_SYNTHESIS.md`); separate within-cluster and between-cluster slopes helped when the audited value was a predictor and raised error on small real audits when it was the outcome, so they were not adopted (`validation/WITHIN_BETWEEN_SYNTHESIS.md`).

# catmech 0.7.0

- `catmech_pool_lm(cluster = )` now uses a delete-one-cluster jackknife variance (`cv3`) by default, with complete-data degrees of freedom of clusters minus one. Intervals with `cluster` get somewhat wider, and `cluster_variance = "cr1s"` reproduces the earlier default. Calls without `cluster` are numerically unchanged.
- A prospective test, with its plan committed before its runner was frozen and before the study ran, chose this under a rule fixed in advance. On 11 real items at a 25% audit, with the random-intercept option, mean distance of coverage from 95% went from 3.18 to 1.47 points, and widths rose by up to 15%. It passed narrowly, meeting two conditions exactly. CR2 with Satterthwaite degrees of freedom failed through imputation, with widths up by as much as 95%, although it did well with every row checked. See `validation/VARIANCE_MI_SYNTHESIS.md`.
- `catmech_pool_lm()` now warns when a fit from `catmech_reference(cluster = )` has the audited value as a predictor. In new simulated cells, that use biased the slope toward zero by up to 16%, with coverage as low as 84.4%. Use `cluster` in the fit only when the audited value is the outcome. This was not known before, and the README's earlier example showed the harmful use.
- `catmech_reference(cluster = )` stays experimental, and earlier verdicts stand, because point estimates do not change. No package dependency was added.
- A later comparison with outside methods, with no package code changed, kept `catmech_reference(cluster = )` under a rule fixed in advance that could have retired it. No method was no worse, but much of the gap to `mice` came from its default of dropping a predictor correlated at 0.99 or more with the imputed variable. With that off, in a rerun made after seeing results, `mice` with `2l.lmer` matched the option's accuracy and fell short only on coverage at a 25% audit. See `validation/OUTSIDE_METHODS_SYNTHESIS.md`.

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
