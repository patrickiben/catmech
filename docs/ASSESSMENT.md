# Assessment and release decision

Assessed 2026-11-01. Version 0.4.0 is worth releasing as an **experimental, auditable toolkit for data with recorded values and a
reference subset**, built from established methods. Real-data testing narrowed its recommended use. It does not support a novel-method
or clinical-validity claim. Earlier assessments are preserved for [v0.1](ASSESSMENT_V0.1.md), [v0.2](ASSESSMENT_V0.2.md) and [v0.3](ASSESSMENT_V0.3.md).

## What changed in v0.4

Version 0.3 treated recorded values as unavailable once references existed, which cost a great deal of precision. Version 0.4 adds
two ways to use them, plus two checks that work on real data.

The [synthetic study](../validation/CONTAMINATION_SYNTHESIS.md) favored the new contamination model, with the caveat that its primary cells use the error form it assumes:
intervals 30% to 38% narrower than the v0.3 default, and coverage within the screen in three of four cells. It failed narrowly in the fourth, and clearly under heavy-tailed targets, selection on
the true value with errors present, and undeclared reference noise.

The [real-data audit](../validation/REAL_DATA_SYNTHESIS.md) was less favorable, and it is the more important result. On ClinicalTrials.gov, linear imputation given the recorded value and generalized raking had lower bias and error
than the contamination model. On NHANES they roughly tied. The model could not identify which records were wrong (AUC 0.54 and 0.50). The gain against discarding recorded values held, at about 48% narrower intervals, but
linear imputation gets it with fewer assumptions. That option is now supported as `catmech_reference(use_recorded = TRUE)` and is the recommended starting point.

## Current judgment

| Question | Judgment |
|---|---|
| New statistical methodology? | Not demonstrated. Contamination models, validation-subsample imputation, raking and approximate Bayesian bootstrap are established. Their combination here may be uncommon, but no full literature review was done. |
| Useful software? | Yes, for data with recorded values and an audit subset: the audit summary, hold-out check, linear MI option, explicit assumptions and visible refusals. |
| Does it estimate individual corruption probabilities? | Under a declared model, yes. They worked when errors were large and separable. On both real datasets the error count was calibrated, but they could not identify which records were wrong. |
| Is the contamination model the best default? | No. Linear MI is the better starting point. The model's niche is large or gross errors, and fallible references with a credible declared SD. |
| Does it repair missed errors with no new information? | No. It requires genuine reference measurements. |
| Generally valid confidence intervals? | No. Failures are recorded for each workflow, and real-data oracle coverage was itself below nominal on one dataset. |
| Real-data evidence? | Two public datasets with real recording differences. Neither is a clinical source-data audit, and the ClinicalTrials.gov reference is not independent truth. |
| Ready for a GitHub research release? | Yes, subject to Patrick's review and the recorded limitations. |
| Ready for a methods paper or clinical use? | No. A short software paper may be realistic after review. A methods claim is not supported. |

## The merged-method idea

A prespecified [go/no-go simulation](../validation/merged-audit-pilot/README.md) tested whether audits in a cluster trial with spillover should follow the direction of spread, combining this package with the CRT and contagion manuscripts.
It returned **NO-GO**: audits helped, but network targeting did not beat auditing the most suspicious records, and direction added nothing. A merged paper on that thesis is not recommended. The test covered one input and one estimator, so it does not rule out
every connection among those methods.

## Next scientific steps

1. Obtain record-level audit data from a real source-data verification, likely through a collaborator. Public data cannot replace this.
2. Check simulated error rates and sizes against published audits. The main fixture uses 8% wrong records, which may be high for transcription errors.
3. Test the hold-out advice prospectively, since it comes from the real-data study that suggested it.
4. Consider non-constant variance, which affected every method on ClinicalTrials.gov, and clustered records, which remain out of scope.

Neither reproducibility nor these validation checks are proof of inferential validity. AI-assisted review is disclosed, and no independent human peer review is implied.
