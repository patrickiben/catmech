# Assessment and release decision

Assessed 2026-10-31. Version 0.2.0 is worth releasing as an **experimental,
auditable implementation of established imputation models**, with explicit
assumptions and visible failures. It does not yet support a novel-methods or
broad clinical-validity claim.

## What changed after the original assessment

The [v0.1 assessment](ASSESSMENT_V0.1.md) preserves the original novelty review,
confirmed prototype defects and prior-art sources. Its starting judgment remains:
`mice`, interval imputation including `hmi`, and multiple overimputation already
cover much of this statistical territory. A targeted search cannot prove priority.
The new Gaussian likelihood and Gibbs augmentation do not change that conclusion.

The useful software contribution is a unified trust policy, explicit category
likelihood, preservation of inputs, machine-readable assumptions, sampler refusal,
analysis helpers and failure-aware sensitivity reporting. It is small enough to
inspect, test and extend. This is a defensible package niche even without a novel
estimator or a performance advantage over established normal imputation.

## The undercoverage diagnosis

Two failures in the historical benchmark motivated a change in the model:

- **Retained contamination:** q=.5 randomly retained some badly corrupted flagged
  values. A marginal fraction of corrupt flags is not a conditional posterior
  corruption probability. Coverage was .43 in the 100-replicate exploratory run,
  with large attenuation and already wide intervals.
- **Incorrect hard categories:** forced bin consistency treated miscoded categories
  as truth. Coverage was .37 in the historical category-error setting. Visibly
  consistent completions did not imply valid regression inference.

Version 0.2 masks every flagged value and adds a coherent Gaussian model. Ignoring
categories avoids the second premise; exact mode rejects trusted contradictions;
noisy mode uses a declared category likelihood. A fixed, wrong confusion matrix
can still fail. This is a change in assumptions and likelihood, not a post hoc
adjustment to confidence-interval width.

The [held-out protocol](../validation/COVERAGE_PLAN.md) uses disjoint seeds,
strong trusted-PMM and independent normal-MI comparators, a null control, oracle
bias decomposition and explicit failure denominators. See the
[synthetic results](../validation/coverage-confirmatory/RESULTS.md) and
[clinical-format stresses](../validation/coverage-clinical/RESULTS.md). The latter
expose severe residual bias from missed flags and site errors. No settings were
retuned to conceal those failures.

An [independent quadrature oracle](../validation/model-oracle-final/README.md)
checks posterior moments for one category example. It supports implementation
correctness there; coverage requires the separate repeated-simulation evidence.

## Novelty, usefulness and validity

| Question | Current judgment |
|---|---|
| New statistical methodology? | Not demonstrated. Normal MI, data augmentation, interval likelihoods and overimputation are established. |
| Useful software? | Yes: explicit assumptions, conflict handling, diagnostics, provenance and finite sensitivity comparisons. |
| Does it estimate individual corruption probabilities? | No. Flags define a masking policy; H is fixed externally or varied as an assumption. |
| Does it solve the original retain/replace bias mechanism? | It removes that policy from the principal model; held-out comparisons test its inferential consequences. |
| Generally valid confidence intervals? | No. Coverage is conditional on trust, selection, outcome-model compatibility and distributional assumptions. |
| Stronger than ordinary trusted normal MI? | The default is the same established model class. No novel-estimator advantage should be claimed. |
| Ready for a GitHub research release? | Yes, subject to the recorded package checks and prominent limitations. |
| Ready for a methods paper or clinical use? | No. External validation, broader model comparisons and independent scientific review remain needed. |

## Relevant research transfers

The carryover work contributes a decomposition into selection and retained error;
False Floor contributes strong alternatives and a null that cannot conceal bias
under nonzero effects; calibrated clinical ML contributes held-out evaluation and
known-answer controls; Anchor Audit contributes explicit unresolved assumption
paths. A local CDISC simulation workflow contributes frozen truth and separately
injected defects. [Transfer details](RESEARCH_TRANSFER_PLAN.md) document sources
and boundaries. No topology or interference estimator is imported.

## Next scientific steps

Before any claim beyond the tested settings, obtain validation measurements or
repeated measures that identify missed contamination and category-error rates.
Model uncertainty in H when estimated, and consider an explicit observation-error
likelihood for W rather than treating all unflagged W as exact. Test nonlinear,
heteroskedastic, informative-selection, censored and dependent-record designs.
An external interval/measurement-error implementation remains a valuable comparator;
`hmi` was reviewed but not executed in this assessment.

The [validation record](../validation/VALIDATION.md) separates software checks,
posterior verification and simulation evidence. Neither reproducibility nor the
research-validation process is itself proof of inferential validity. AI-assisted
review is disclosed; no independent human peer review is implied.
