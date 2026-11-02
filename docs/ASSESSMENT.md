# Assessment and release decision

Assessed 2026-11-02. Version 0.7.1 is worth releasing as an **experimental, auditable toolkit for data with recorded values and a
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
| Does it handle clustered audits? | Partly. `catmech_pool_lm(cluster = )` is a standard cluster-robust variance. `catmech_reference(cluster = )` is experimental: better worst-case coverage when the audited value is the outcome, some added error on rare items, two planned tests that came out PARTIAL, and bias when the audited value is a predictor. A planned comparison kept it, narrowly, over `mice` and design-based methods. |
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
5. Test the clustered option on a second health system. Every item in the one recount available has now been seen, and two planned tests have run on it.
6. Have someone who knows `mice` and `survey` review the [comparison with outside methods](../validation/OUTSIDE_METHODS_SYNTHESIS.md), since one author set up every method. Rare counts also remain unsolved, even with every row checked.
7. Test whether checking errors against analysis covariates on audited rows can guide the choice of method. The real-audit results suggested this, so they cannot test it fairly, and at a 25% audit the check has low power.

Neither reproducibility nor these validation checks are proof of inferential validity. AI-assisted review is disclosed, and no independent human peer review is implied.
