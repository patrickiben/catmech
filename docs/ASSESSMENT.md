# Assessment and release decision

Assessed 2026-11-05. Version 1.0.0 graduates to an official production release as an **auditable toolkit for data with recorded values and a
reference subset**, built from established methods. Real-data testing narrowed its recommended use. It does not support a novel-method
or clinical-validity claim. Earlier assessments are preserved for [v0.1](ASSESSMENT_V0.1.md), [v0.2](ASSESSMENT_V0.2.md) and [v0.3](ASSESSMENT_V0.3.md).

## What changed in v0.4

Version 0.3 treated recorded values as unavailable once references existed, which cost a great deal of precision. Version 0.4 adds
two ways to use them, plus two checks that work on real data.

The [synthetic study](../validation/CONTAMINATION_SYNTHESIS.md) favored the new contamination model, with the caveat that its primary cells use the error form it assumes:
intervals 30% to 38% narrower than the v0.3 default, and coverage within the screen in three of four cells. It failed narrowly in the fourth, and clearly under heavy-tailed targets, selection on
the true value with errors present, and undeclared reference noise.

The [real-data study](../validation/REAL_DATA_SYNTHESIS.md), on recording differences rather than audits, was less favorable, and it is the more important of the two. On ClinicalTrials.gov, linear imputation given the recorded value and generalized raking had lower bias and error
than the contamination model. On NHANES they roughly tied. The model could not identify which records were wrong (AUC 0.54 and 0.50). The gain against discarding recorded values held, at about 48% narrower intervals, but
linear imputation gets it with fewer assumptions. That option is now supported as `catmech_reference(use_recorded = TRUE)` and is the recommended starting point.

Version 0.4.1 adds `catmech_audit_plan()`, which estimates from a pilot how many records to audit. Its [prospective test](../validation/AUDIT_PLAN_SYNTHESIS.md) returned **NO-GO** under a rule fixed in advance: predictions were within 10% at target audits of 10% or more,
and missed at a 5% audit, whose width is too noisy to predict that closely. Small pilots also estimated `k` poorly, and its interval undercovered. The planner is documented, with that limit stated, and is not validated.

A later [test on real audit data](../validation/REAL_AUDIT_SYNTHESIS.md), with no package code changed, used a Rwanda health-facility recount and a psychology effect-size audit. It is the first evidence here from data where someone went back to the source, though neither audit holds patient values, and its "full audit" recounts one item with the covariate left as reported. Under its rule fixed in advance, the contamination model's advantage was **inconclusive**. At a 25% audit it had lower error than linear MI on two of three slopes, but failed the coverage condition on all three. Where a reported slope differed from the recount it kept 38% to 73% of the gap, with intervals covering 35% to 72% against 90% to 97% for the recount. Across 20 cells it had the lowest point error in 11 and linear MI in 4, while linear MI's coverage was closest to 95% in 10 and the model's in 1. Linear MI's near-zero average bias hides single audits that undershoot or overshoot, and at a 25% audit its point error clearly beat trusting reported values on only one of four slopes. Three AI reviewers found no error in the rule or tables, and corrected several readings, including a mechanism the first draft had stated as fact.

Version 0.5.0 adds `cluster` options to `catmech_reference()` and `catmech_pool_lm()`. The first resamples audited clusters in each imputation, an idea suggested by a post-hoc check in the real-audit test. Its [planned test](../validation/CLUSTERED_SYNTHESIS.md) came out **PARTIAL** under a rule fixed in advance: it held on coverage and in simulation but failed on accuracy. On five held-out items at a 25% audit, it helped on the two where linear MI undercovered (87.0% to 91.4%, and 79.6% to 88.8%) and raised error on the other three. It stays experimental. One independent AI reviewer then confirmed the verdict and found no bug, and it corrected several of the first write-up's readings. Post hoc and untested, the losses fell on rare events whose residuals show no within-facility correlation. Both the gain and the cost come mainly from resampling whole clusters rather than rows, resampling within strata understates variance when few clusters per stratum are audited, and the plan's simulations barely contain the clustering the option is for.

Version 0.6.0 replaces the resampling in `catmech_reference(cluster = )` with a random intercept, and removes `cluster_strata`. A [second planned test](../validation/CLUSTERED_FIX_SYNTHESIS.md) selected it over a stratified jackknife draw, using simulations with clustered residuals, and also came out **PARTIAL** under a rule fixed in advance. On nine Rwanda items at a 25% audit, coverage was within its limits on eight and 89.0% on the ninth. Error was lower than linear MI's on six common items and 4% to 5% higher on three rare ones, and both accuracy conditions failed. Against the resampling it replaced, error and interval stability improved, and average coverage barely moved. It stays experimental. This round has had no outside review, and the real items had all been seen under other methods.

Version 0.7.0 changes the default variance in `catmech_pool_lm(cluster = )` to a delete-one-cluster jackknife, after a [third planned test](../validation/VARIANCE_MI_SYNTHESIS.md) that it passed narrowly. The same test found that `catmech_reference(cluster = )` biases slopes toward zero in simulations when the audited value is a predictor, so that use now draws a warning. All of this round's work comes from one assistant, with no outside review.

A later [comparison with outside methods](../validation/OUTSIDE_METHODS_SYNTHESIS.md) kept `catmech_reference(cluster = )` under a rule that could have retired it, and changed no code. The margin is narrow, since `mice` offers the same model and matched its accuracy once a default was turned off. That comparison was not neutral: the option's author set up every comparator.

Version 0.7.1 follows an [overall check](../validation/OVERALL_QC.md) by five independent AI reviewers. It adds a warning when every audited row agrees and recorded values are kept, a warning when one cluster determines a coefficient, and a refusal when no between-cluster variance can be estimated. These cases gave badly wrong intervals with no signal before. Earlier study numbers are unchanged.

Three later planned tests changed no package code. On a [second health system](../validation/SECOND_SYSTEM_SYNTHESIS.md) the outcome was mixed, with errors that barely clustered. In the [count-data test](../validation/RARE_COUNTS_SYNTHESIS.md), coverage fell with more variable, skewed cluster frailty at both tested count levels. Its original **HEAVY TAILS** verdict remains, but variance, skewness and dependence changed together, so it does not isolate a tail-shape mechanism or exclude rarity as a problem elsewhere. No candidate fix met the rules. [Comparators written by fresh agents](../validation/INDEPENDENT_SETUPS_SYNTHESIS.md) did not meet the rule, but the best looks better than the option on most measures, so the option is not the best available.

A [later test of separate within-cluster and between-cluster slopes](../validation/WITHIN_BETWEEN_SYNTHESIS.md), on new SNAP quality-control reviews and a new simulation, was **not adopted**, and it changed no package code. With the audited value as a predictor, coverage and bias improved, though linear MI kept lower bias in simulation. With it as the outcome, error rose: 1.258 times the option's on real items at a 25% audit, and the new simulation's outcome class missed the accuracy limit at 1.026. The slopes alone also did not match the fresh agent's setup that prompted the test. SNAP's errors barely cluster by state, and no human has reviewed this.

Version 0.7.2 fixes the audit planner's formula offsets and compound-stratum identification, and refuses mixed stochastic recorded-value patterns within a cluster until a shared draw is supported. The variance-component draws are explicitly approximate. The rare-count interpretation above was narrowed after review; its original protocol, decision label and numerical results remain unchanged. These corrections do not establish general coverage.

Version 0.8.0 adds an explicit joint model through `pan` for sites containing both recorded-present and recorded-missing unaudited rows. It shares the site draw coherently and adds multiple-chain diagnostics, declared covariance priors and a separate [validation study](../validation/JOINT_CLUSTER_SYNTHESIS.md). It is an established bivariate Gaussian mixed model, not a novel estimator. Historical moment-model studies do not establish the new option's calibration, and this repair is not evidence that rare-count undercoverage or predictor attenuation has been cured. The default moment model remains unchanged for supported inputs.

## Joint posterior inference in v0.9.0

The explicit `catmech_joint_summary()` API targets the population coefficients
of the fitted joint model. A fresh [2,600-cohort comparison](../validation/JOINT_CALIBRATION_SYNTHESIS.md)
passed separate prespecified screens for calibration, relative interval-width
sensitivity and sampler availability. Direct primary coverage was 93.0%–96.75%,
and the larger sampling budget returned all 3,800 candidate fits without weaker
diagnostics. The priors remain unchanged; coefficient sensitivity was not
uniformly reduced. Contextual-covariate stress produced 63.5% coverage and
substantial bias, so this is useful explicit model-based inference, not a general
robust marginal-effect procedure. `catmech_pool_lm()` remains a separate analysis.

## Contextual population inference in v0.10.0

The explicit `joint_marginal` option jointly models continuous predictors and
measurements, allowing their site effects to correlate. It uses a proper
factorized covariance/regression prior and an independently checked native
sampler. Its contribution is usable, auditable model-based inference with a
clearer population target, not a new statistical model class.

A fresh [3,600-cohort confirmation](../validation/JOINT_MARGINAL_SYNTHESIS.md)
passed contextual calibration, ordinary Gaussian calibration and slope prior-
sensitivity screens. Contextual coverage was 94.00%–95.75%, compared with
58.50%–67.00% for the conditional direct summary; mean absolute contextual bias
fell from 0.1306 to 0.0046. Measured slope sensitivity over the two declared new
prior axes was 70.1% lower on 800 complete comparisons. Small remaining shifts
are near finite-chain Monte Carlo error and do not show prior independence.

The availability nonregression screen **failed**: two of 3,600 candidate default
fits refused at the unchanged diagnostics, versus none for the baseline. Each
cell still returned at least 99% of fits and intervals. Intervals in the four ordinary calibration cells were
12.2%–16.3% wider, and ordinary calibration was not shown to improve. This supports
an explicit experimental option for its stated contextual model, not a universal
replacement. Repeated-sampling evidence covers one continuous-predictor slope;
it does not validate multiple-predictor coverage, nonlinear relationships,
informative site sizes, rare counts, imperfect references or MNAR. The earlier
real-data evidence does not validate the new predictor-distribution model.

## Numerical extension in v0.10.1

A single extra doubling to 64k fixes the two observed ESS refusals, retaining
every check and the original model/prior. The complete dependent numerical
replay returned all 12,000 MCMC fits and 27,600 intervals; all four screens pass.
All27,596 formerly successful analyses remain unchanged. This is repeated
numerical validation of the original cohorts, not independent confirmation.
The [report](../validation/MARGINAL_EXTENSION_SYNTHESIS.md) preserves an
interrupted forked execution and the complete replacement run in independent
R workers. No claim of guaranteed future convergence follows.

## Production readiness in v1.0.0

Version 1.0.0 freezes the core software interface and completes the clinical data verification
and analysis toolkit:
- Downstream multi-endpoint regression: `catmech_pool_glm()` extends pooling to generalized linear
  models (logistic, Poisson) with Barnard-Rubin degrees of freedom, and `catmech_pool_coxph()`
  provides time-to-event survival pooling using event-based degrees of freedom.
- Diagnostic visualization: base R `plot()` methods provide visual checks for holdout predictive
  calibration, recording audits, audit planning curves, and prior sensitivity grids.
- Prospective audit planning: `catmech_audit_power()` computes required audit sample sizes and
  statistical power prospectively or empirically.
- Regulatory compliance trail: `catmech_audit_trail()` generates SHA-256 verified provenance reports
  (Markdown, JSON, HTML) documenting input hashes, model specifications, convergence diagnostics,
  and pooled estimates.
- Comprehensive clinical vignette: `vignettes/catmech-clinical-audit.Rmd` walks through an end-to-end
  multi-center clinical study.
- Software maturity: passes all 1,717 testthat unit tests and R CMD check with Status: OK (0 errors,
  0 warnings, 0 notes).

The statistical caveats remain: measurement error cannot be identified without assumptions or
reference data, and coverage is not guaranteed under unmodeled misspecification. Core functions are
auditable and production-ready for declared workflows.

## Current judgment

| Question | Judgment |
|---|---|
| New statistical methodology? | Not demonstrated. Contamination models, validation-subsample imputation, raking and approximate Bayesian bootstrap are established. Their combination here may be uncommon, but no full literature review was done. |
| Useful software? | Yes, for data with recorded values and an audit subset: the audit summary, hold-out check, linear MI option, explicit assumptions and visible refusals. |
| Does it estimate individual corruption probabilities? | Under a declared model, yes. They worked when errors were large and separable in simulation. On two real datasets with recording differences the error count was calibrated, but they could not identify which records were wrong. On the real audits they could not either (AUC at most 0.62 on Rwanda), and error counts ran low on three of four Rwanda slopes at a 25% audit. |
| Is the contamination model the best default? | Not when an interval is needed. Linear MI is the better starting point for that. In simulation the model's niche was gross errors and fallible references with a credible declared SD. A real audit with moderate errors did not test gross errors, so that niche has no real-data support yet. There the model often had lower point error, but its intervals were too narrow. |
| Can it say how many records to audit? | Not validated. The planner gave a NO-GO, failing only at the 5% target. It is a documented calculation, built on classical theory. |
| Does it repair missed errors with no new information? | No. It requires genuine reference measurements. |
| Generally valid confidence intervals? | No. Failures are recorded for each workflow. With default intervals the oracle itself undercovered on ClinicalTrials.gov and on the Rwanda audit, and no method in the real-audit test held the oracle's clustered coverage at a 25% facility audit. The experimental clustered option came closer, and under the variance then in use still reached only 89.0% on one of nine items, against 93.2% for a recount. |
| Does it handle clustered audits? | The opt-in `cluster_model = "joint"` handles mixed recording patterns through a shared bivariate site model. `"joint_marginal"` also models continuous predictors and their site associations; its calibration and prior-sensitivity screens passed, but its v0.10.0 strict availability screen failed. The v0.10.1 numerical extension passes that screen on the same grid. Each model has separate assumptions and evidence. For the historical default, partly. `catmech_pool_lm(cluster = )` is a standard cluster-robust variance. `catmech_reference(cluster = )` is experimental: better worst-case coverage when the audited value is the outcome, some added error on rare items, two planned tests that came out PARTIAL, and bias when the audited value is a predictor. Planned comparisons kept it under their rule, narrowly, and an independently written setup looks better on most measures. |
| Real-data evidence? | Two public datasets with real recording differences, and two public audits: Rwanda facility counts and psychology effect sizes. None is a patient-level clinical source-data audit, and the ClinicalTrials.gov reference is not independent truth. |
| Ready for a GitHub research release? | Yes, with the recorded limitations. |
| Ready for a methods paper or clinical use? | No. A short software paper may be realistic after review. A methods claim is not supported. |

## The merged-method idea

A prespecified [go/no-go simulation](../validation/merged-audit-pilot/README.md) tested whether audits in a cluster trial with spillover should follow the direction of spread, combining this package with the CRT and contagion manuscripts.
It returned **NO-GO**: audits helped, but network targeting did not beat auditing the most suspicious records, and direction added nothing. A merged paper on that thesis is not recommended. The test covered one input and one estimator, so it does not rule out
every connection among those methods.

## Next scientific steps

1. Obtain patient-level audit data from a real source-data verification, likely through a collaborator. A search of 77 queries found no public file of this kind, and the public audits above cannot replace it.
2. Check simulated error rates and sizes against published audits. The main fixture uses 8% wrong records, which may be high for transcription errors.
3. Test the hold-out advice prospectively, since it comes from the real-data study that suggested it.
4. Consider non-constant variance, which remains out of scope. On the Rwanda audit it and clustering together cost the oracle up to 28 points of default coverage.
5. Find what else separates the fresh agent's setup from the option, now that [separate slopes alone](../validation/WITHIN_BETWEEN_SYNTHESIS.md) did not match it. Any such test needs data not yet seen, since both recounts, SNAP and the new simulation have all been seen.
6. Have someone who knows `mice` and `survey` review the [comparisons](../validation/INDEPENDENT_SETUPS_SYNTHESIS.md). Fresh agents reduce the problem of one author setting up every method, and they do not remove it. Coverage shortfalls under variable, skewed cluster frailty also remain, even with every row checked; comparisons holding variance and dependence fixed would help distinguish possible mechanisms.
7. Test whether checking errors against analysis covariates on audited rows can guide the choice of method. The real-audit results suggested this, so they cannot test it fairly, and at a 25% audit the check has low power.

Neither reproducibility nor these validation checks are proof of inferential validity. AI-assisted review is disclosed, and no independent human peer review is implied.
