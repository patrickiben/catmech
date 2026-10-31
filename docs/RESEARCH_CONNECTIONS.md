# Connections to related work by the author

Reviewed against the author's GitHub repositories on 2026-10-31. These connections
suggest useful evaluation and application paths. They do not transfer validity
from one project to another or establish novelty for catmech.

| Existing project | Concrete contribution here | Boundary |
|---|---|---|
| [A False Floor](https://github.com/patrickiben/false-floor) | Test the strongest simple alternative before claiming a gain. The catmech benchmark includes flag-aware PMM and found no clear general advantage over it. | This is an evaluation principle, not evidence that either imputer is valid. The repository's main manuscript snapshot is historical; follow its current-paper links. |
| [Residual carryover study](https://github.com/patrickiben/biometrics-ai-portfolio/tree/main/sim_lifecycle/carryover_oc) | Separate bias caused by selecting retained observations from bias caused by residual error. Its oracle arms, paired comparisons, known-answer checks and MCSE reporting are a useful template for a larger catmech study. | Pharmacokinetic carryover is not the same mechanism as recording corruption. Its estimators and regulatory rules do not transfer directly. |
| [Calibrated clinical ML](https://github.com/patrickiben/calibrated-clinical-ml) | If q becomes a prediction model, evaluate corruption probabilities against independently adjudicated truth on held-out records. An anomaly ranking does not establish probability calibration. | The existing synthetic calibration experiments do not validate catmech's q. A corruption model needs new labels, assumptions and evidence. |
| [Contact-map Anchor Audit](https://github.com/patrickiben/contact-map-anchor-audit) | Report which conclusions change along a declared assumption path and retain unresolved outcomes. The supplied q-grid example is a first descriptive sensitivity view. | Persistence-diagram stability bounds cannot be reused as imputation guarantees. Observed stability over a finite grid is not proof over every admissible assumption. |
| [eCRF monitoring apps](https://github.com/patrickiben/ecrf-pipeline-apps) | A possible upstream source of documented quality flags, paired with synthetic laboratory data and per-run audit records. | Its anomaly scores are leads for review, not corruption probabilities. A category derived from the same corrupted value supplies no independent truth constraint. |
| [CRT design work](https://github.com/patrickiben/crt-design-demo) | Site-correlated error and contaminated donor groups would be useful future stress tests. Public R study code is linked from the demo's OSF deposit. | Spatial treatment interference differs from measurement error. Current catmech assumes an independent-record regression and does not implement network-aware donor validity. |

## Best next research study

Use synthetic laboratory data with a frozen original truth and independently
recorded interval category. Inject known value errors, false flags, missed errors,
rounding differences and laboratory-specific shifts. Compare flag-aware PMM,
catmech policies and an explicit measurement-error model. Decompose the effects of
donor selection and residual retained error. Predeclare the downstream estimand,
test analysis coverage, and report failure rates as well as successful fits.

Use only one analyte and visit for an initial independent-record experiment.
Clinical low/normal/high conventions may treat the lower limit differently from
catmech's right-closed intervals. Reconcile endpoint and rounding conventions before
using those categories. Below-quantification measurements need a separate censored
model: no observed low-region donors exist, and fitting only quantifiable values
introduces selection. An optional truncated-normal fallback does not solve that
identification or likelihood problem.

The immediate package therefore retains a small supported scope. The useful
connection across these projects is explicit assumptions, strong alternatives,
calibration, and failure reporting. Topology, interdiction and causal interpretation
would each require additional model development and validation.
