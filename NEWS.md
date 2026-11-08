# catmech 1.0.1

Fixes from a review of 1.0.0 on 2026-11-07; evidence and the scripts that reproduce
it are in `validation/REVIEW_FIXES_2026-11-07.md`. Not submitted to CRAN.

- `catmech_compare_strategies()` now keeps the analysis outcome among the imputation
  predictors it detects. It used to drop it, which biased the linear-imputation and
  contamination-model slopes toward zero whenever the audited value was a predictor:
  in 100 simulated data sets with a true slope of 0.80, its linear-imputation estimate
  averaged 0.63. The `cluster` argument no longer switches on the experimental
  clustered imputation; it now gives every strategy, including the validation subsample
  and the naive analysis, the same cluster-robust variance over the coefficients that
  strategy can estimate. A failed contamination fit now warns with its reason, and the
  predictors used are printed.
- `catmech_pool_lm()`, `catmech_pool_glm()` and `catmech_pool_coxph()` warn when the
  audited value is an analysis predictor but the outcome (for a Cox model, the event
  indicator) was not an imputation predictor. Leaving it out attenuated a logistic log
  odds ratio from 0.80 to 0.61, with 65% coverage of nominal 95% intervals, in 200
  simulated data sets. `catmech_pool_glm()` and `catmech_pool_coxph()` now also give the
  clustered-imputation warnings that `catmech_pool_lm()` already gave. The examples and
  the clinical audit vignette now put the outcome into the imputation model, and for Cox
  models the event indicator and the Nelson-Aalen cumulative hazard.
- The cluster-robust variance in `catmech_pool_glm()` used the dispersion-scaled
  covariance as its bread, so for gaussian, Gamma and quasi families its standard errors
  were multiplied by the estimated dispersion: too wide when it was above 1 (three-fold in
  an overdispersed count example) and too narrow when it was below 1, as in most Gamma
  fits (2.5-fold too narrow in a Gamma example with dispersion 0.41). Rerun such
  clustered fits. Binomial and Poisson results are unchanged.
- `catmech_audit_power()` is now labeled a rough, unvalidated planning aid. Its
  documentation states the assumption behind its default factor, and without `plan` it
  also reports the audit size that linear imputation needs under the double-sampling
  regression estimator at the same correlation, which is often several times larger.
  "Power" there is the share of plausible factors meeting the target, not the power of a
  test.
- `catmech_audit_trail()` no longer calls its record tamper-evident or says the pooled
  intervals are hashed. Its documentation now says which hashes are of serialized R data
  frames and which of canonical JSON text, and that compliance with FDA 21 CFR Part 11 or
  ICH E9 depends on the system around it. An unknown package version is recorded as
  "unknown" rather than "1.0.0", and the HTML title is escaped.
- The vignette, the README's method guide and the help pages no longer describe the
  package as production-ready, a compliance tool or a source of guarantees, and they
  follow the evidence on clustered audits: pool with `cluster =`, and use
  `catmech_reference(cluster = )` only as an experimental check when the audited value is
  the outcome. The validation records and the assessment that called 1.0.0
  production-ready now open with a dated correction.

# catmech 1.0.0

Some descriptions in this entry were corrected on 2026-11-07; see 1.0.1 above.

- First 1.0 release (not submitted to CRAN). The public API is stable and backward-compatible across
  model fitting, diagnostics, and downstream inference. The clustered imputation
  option (`catmech_reference(cluster = )`) and the joint-marginal engine remain
  experimental, and the audit-size planner keeps its NO-GO validation status.
- Add Generalized Linear Model pooling (`catmech_pool_glm()`) supporting logistic
  and Poisson regressions with Rubin's combining rules, Barnard-Rubin small-sample
  degrees of freedom, exponentiated effect reporting (odds and rate ratios), and
  Huber-White-Liang-Zeger cluster-robust sandwich variance.
- Add Cox proportional hazards survival pooling (`catmech_pool_coxph()`) with
  event-based Barnard-Rubin degrees of freedom, exponentiated hazard ratios, and
  Lin-Wei cluster-robust sandwich covariance.
- Add base R diagnostic `plot()` methods across the package:
  `plot.catmech_marginal_sensitivity_grid()`, `plot.catmech_audit_plan()`,
  `plot.catmech_holdout_check()`, `plot.catmech_recording_audit()`,
  `plot.catmech_audit_power()`, `plot.catmech_pooled_glm()`,
  `plot.catmech_pooled_coxph()`, and `plot.catmech_pooled_lm()`.
- Add multi-strategy audit comparison framework (`catmech_compare_strategies()`)
  with companion comparative forest plot (`plot.catmech_strategy_comparison()`).
- Add tidyverse and `broom` compatibility generics: `tidy()` and `glance()` methods
  for all pooled regression models (`catmech_pooled_glm`, `catmech_pooled_coxph`,
  `catmech_pooled_lm`).
- Add audit trail and provenance report generator (`catmech_audit_trail()`) producing
  tamper-evident SHA-256 data and completion checksums across Markdown, JSON,
  and standalone HTML formats.
- Add prospective audit power and sample size calculator (`catmech_audit_power()`)
  for sizing two-phase data audits before collecting gold-standard measurements.
- Add comprehensive clinical audit vignette (`vignettes/catmech-clinical-audit.Rmd`)
  demonstrating simulation, holdout calibration, audit planning, strategy comparisons,
  and LM/GLM/Cox pooling.
- Formalize S3 console output with dedicated `print()` and `summary()` methods.
- Pass all 1,754 testthat unit tests. `R CMD check --as-cran` reports 0 errors and
  0 warnings; the only notes are the new-submission note and local-environment notes.

# catmech 0.10.1

- Extend only the marginal engine's default maximum to 64,000 retained iterations
  per chain. Retain four chains, the same burn-in, doubling/restarts, every
  diagnostic monitor, R-hat cutoff and bulk/tail ESS requirements. Existing
  conditional-joint defaults and explicit user maxima retain their meaning.
- Both original v0.10.0 convergence refusals pass the selected numerical challenge
  at 64,000; earlier attempts reproduce the recorded history and diagnostics.
  The likelihood, prior and interval calculation are unchanged.
- The full same-cohort numerical replay returned all 12,000 MCMC fits and
  27,600 intervals; all four screens pass. All 27,596 previously successful
  analyses are unchanged. Original and interrupted evidence are preserved;
  this is not an independent scientific confirmation.
- Run validation workers in independent PSOCK R processes and reject missing
  worker results. Pin macOS CI to macos-15; all 1,526 assertions pass on all
  three platforms.

# catmech 0.10.0

- Add the explicit `cluster_model = "joint_marginal"` for a joint Gaussian model
  of the true measurement, recorded value and complete continuous predictors.
  Correlated site effects permit contextual associations. The new
  `catmech_marginal_summary()` reports population marginal coefficients from all
  parameter draws, including uncertainty in the predictor distribution.
- Factor covariance priors into predictor and conditional-measurement
  Huang--Wand blocks with an explicit proper matrix-normal regression prior.
  Add an original compiled blocked Gibbs sampler. An exact scale-interweaving move addresses slow mixing
  near zero site variance without changing the posterior or diagnostic screens.
- Preserve the existing conditional joint model and its priors. New marginal
  inference requires Gaussian continuous predictor associations, noninformative
  site size, exact references and ignorable selection; categorical, incomplete,
  site-constant or aliased predictors are unsupported.
- Add independent Gaussian conditioning and parameter-posterior oracles, prior
  distribution checks, matched-posterior comparison with `pan`, exact MH density
  and Jacobian checks, and tests of transformed marginal coefficients.
- Building from source now requires a C++17 toolchain, Rcpp and RcppArmadillo.
  Frozen v0.9.0 sources and the first numerical pilot remain reproducible.

- Validate on 3,600 fresh one-predictor cohorts. Contextual coverage was
  94.00%–95.75% versus 58.50%–67.00% for the conditional direct summary;
  measured slope sensitivity across both declared prior axes was 70.1% lower.
  Contextual, ordinary calibration and sensitivity screens passed. The strict
  availability screen failed: two of 3,600 candidate default fits refused at
  unchanged ESS thresholds, versus none for the baseline. Ordinary intervals
  were wider; the new model is an explicit experimental option.
- Qualify `stats::ave` for namespace-only input preparation. All 3,600 frozen
  confirmation inputs and prepared model objects are byte-identical after this
  engineering correction; no confirmation fit was replaced.

# catmech 0.9.0

- Add `catmech_joint_summary()` for population target coefficients in the exact
  fitted joint-model design. It uses all retained posterior draws, transforms
  the intercept and slopes to original units, and returns credible intervals,
  convergence diagnostics and Monte Carlo errors. It does not change
  `catmech_pool_lm()` or support arbitrary downstream formulas.
- Vectorize conditional-mean trace calculations and increase the joint sampler's
  initial retained run to 4,000 and maximum to 32,000 iterations per chain.
  Diagnostic thresholds remain unchanged; the original-scale target intercept
  is now monitored as well. Explicit historical controls retain their meaning.
- Add an independent known-covariance Gaussian coefficient oracle, transformation
  and refusal tests, numerical compatibility checks and a replay of the 22
  previously failed fits. All selected failures passed with the longer schedule;
  that selected challenge is not an estimate of future failure rates.
- Preserve v0.8.0 sources and evidence. A fresh 2,600-cohort comparison passed
  separate calibration, relative interval-width sensitivity and availability
  screens. Direct primary coverage ranged from 93.0% to 96.75%; all 3,800
  candidate fits returned. Covariance priors are unchanged and coefficient
  sensitivity is not uniformly reduced. Contextual-covariate stress produced
  63.5% coverage, limiting the direct summary to its declared model assumptions.
  Full results, refusals and Monte Carlo uncertainty are retained in
  `validation/JOINT_CALIBRATION_SYNTHESIS.md`.

# catmech 0.8.0

- Add explicit `cluster_model = "joint"` to normal reference imputation with
  recorded values and clusters. A bivariate Gaussian random-intercept model
  fitted by `pan` shares site uncertainty across recorded-present and
  recorded-missing rows. It also predicts entirely unobserved sites, while
  preserving observed gold references and raw input exactly.
- Expose and record covariance prior scales and multiple-chain controls;
  reject unsupported designs and failed convergence. Joint fitting does not
  use the moment model's deterministic all-agree shortcut.
- Run each `pan` chain in a fresh local R process to avoid hidden cached RNG
  state making same-seed fits depend on earlier unrelated `pan` calls.
- Add independent conditional-normal oracle tests, mixed-pattern integration
  tests, a runnable installed example and a separately planned simulation with
  prior sensitivity and misspecification stress. Historical real-audit results
  remain evidence for earlier models, not validation of this option.
- Preserve `cluster_model = "moment"` as the default. The joint model is an
  explicit change of likelihood and priors, not a general correction for
  rare-count undercoverage or predictor-slope attenuation.

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
