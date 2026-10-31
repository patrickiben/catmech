# Assessment and release decision

Assessed 2026-10-31. Version 0.3.1 is worth releasing as an **experimental,
auditable implementation of established imputation models**, with explicit
assumptions and visible failures. It does not yet support a novel-methods or
broad clinical-validity claim.

## What changed after the original assessment

The [v0.1 assessment](ASSESSMENT_V0.1.md) preserves the original novelty review,
confirmed prototype defects and prior-art sources. Its starting judgment remains:
`mice`, interval imputation including `hmi`, and multiple overimputation already
cover much of this statistical territory. A targeted search cannot prove priority.
The Gaussian likelihood, Gibbs augmentation and reference-subset MI do not
change that conclusion. The [v0.2 assessment](ASSESSMENT_V0.2.md) is preserved
with its original scope.

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

## Reference information and fair comparisons in v0.3

`catmech_reference()` changes what is observed: independently measured reference
values are supplied on a subset, and every other true target is missing even
when recorded W is present. A random accurate reference subset can identify a
missing-data model under adequate regression assumptions. No reference data means
no such information; the API refuses that case. Selection and reference accuracy
remain substantive assumptions. This is established internal-validation MI, not
a new method for detecting corruption from flags.

`catmech_pmm()` adds an audited ordinary donor method that masks every flag,
without q retention or forced bins. The earlier heavy-tail comparator had also
received category covariates, so its difference from default Gaussian imputation
could not be attributed to PMM alone. The [new study plan](../validation/EXTENSION_PLAN.md)
compares the complete normal/PMM by category-information factorial, including
strong MAR and incorrect-category stresses. It also compares reference methods
with reference-only analysis and independent implementations given the same data.

The [reference method](REFERENCE_METHOD.md) documents the additional data contract,
primary sources and the distinction between known references and a second noisy
measurement. The new [predictive oracle](../validation/reference-oracle/README.md)
checks means, variances and cross-record covariance independently. These checks
verify implementation properties; the held-out study supplies separate operating
characteristics.

The [completed study](../validation/EXTENSION_SYNTHESIS.md) gives a qualified
positive result: with 25% random, accurate reference measurements, normal MI
covered 95.2% of nominal 95% intervals, versus 58.1% for flag masking with missed
errors and 12.7% under compound site/error stress. These paired comparisons use
the same 1,000 cohorts, and the reference method has additional measured
information. Reference-only analysis covered 95.3%; normal MI reduced mean
interval width from .3304 to .2782. Independent standard normal MI performed
similarly, supporting the established model rather than a proprietary advantage.

Reference PMM covered only 85.1%, exactly matching the independent mice PMM
comparator. Version 0.3.1 therefore makes partial-reference PMM pooled analysis
explicitly sensitivity-only. The normal reference method also failed with
nonignorable selection (5.5% coverage) and noisy references (67.8%). Ordinary
PMM helped the tested heavy-tail example but hurt strong-MAR coverage at equal
predictor information. These failures rule out a general robustness claim and
remain in the released evidence; they were not tuned away.

## Novelty, usefulness and validity

| Question | Current judgment |
|---|---|
| New statistical methodology? | Not demonstrated. Normal MI, data augmentation, interval likelihoods and overimputation are established. |
| Useful software? | Yes: explicit assumptions, conflict handling, diagnostics, provenance and finite sensitivity comparisons. |
| Does it estimate individual corruption probabilities? | No. Flags define a masking policy; H is fixed externally or varied as an assumption. |
| Does it solve the original retain/replace bias mechanism? | It removes that policy from the principal model; held-out comparisons test its inferential consequences. |
| Can it repair missed errors with no new information? | No. The reference path requires genuine reference measurements and an appropriate selection/model assumption. |
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

The new reference path makes independently obtained validation measurements usable;
it does not supply or scientifically validate them. Before claims beyond the tested
settings, obtain suitable real reference measurements or repeated measures that
identify missed contamination and category-error rates.
Model uncertainty in H when estimated, and consider an explicit observation-error
likelihood for W rather than treating all unflagged W as exact. Test nonlinear,
heteroskedastic, informative-selection, censored and dependent-record designs.
An external interval/measurement-error implementation remains a valuable comparator;
`hmi` was reviewed but not executed in this assessment.

The [validation record](../validation/VALIDATION.md) separates software checks,
posterior verification and simulation evidence. Neither reproducibility nor the
research-validation process is itself proof of inferential validity. AI-assisted
review is disclosed; no independent human peer review is implied.
