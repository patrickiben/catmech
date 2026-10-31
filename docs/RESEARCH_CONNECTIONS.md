# Research connections used in version 0.2

Reviewed 2026-10-31. These connections informed concrete changes; they do not
transfer a validity theorem from another project.

| Source | Implemented transfer | Boundary |
|---|---|---|
| [Residual carryover study](https://github.com/patrickiben/biometrics-ai-portfolio/tree/main/sim_lifecycle/carryover_oc) | Paired oracle decomposition of selection and retained measurement error; MCSEs and fixed streams | Carryover estimators and regulatory rules are not imported. |
| [A False Floor](https://github.com/patrickiben/false-floor) | Strong trusted PMM and independent normal-MI baselines; null and nonzero-effect controls | A calibrated null does not establish power or remove attenuation. |
| [Calibrated clinical ML](https://github.com/patrickiben/calibrated-clinical-ml) | Separate engineering and held-out streams; independent known-answer checks | No anomaly score is converted into a corruption probability. H is not learned here. |
| [Contact-map Anchor Audit](https://github.com/patrickiben/contact-map-anchor-audit) | `catmech_sensitivity()` retains errors, missing planned cases and unresolved conclusions | Finite-grid stability is descriptive; topological stability bounds do not apply. |
| Local CDISC evaluation work | Standalone synthetic laboratory fixture, frozen truth, independent defects, missed flags, rounding and site shifts | No real patient data or private source files are bundled. Clinical categories are ignored until their endpoint/rounding model is supported. |
| [CRT design work](https://github.com/patrickiben/crt-design-demo) | Motivation for site-error stress | The current fixed laboratory shift is not a random-cluster outcome model or an interference analysis. |

The [transfer plan](RESEARCH_TRANSFER_PLAN.md) records inspected source versions
and the [prospective coverage plan](../validation/COVERAGE_PLAN.md) specifies the
executed study. The [clinical results](../validation/coverage-clinical/RESULTS.md)
show that missed corruption remains a major limitation, despite this better
validation workflow.

An upstream app such as [eCRF monitoring](https://github.com/patrickiben/ecrf-pipeline-apps)
could provide documented review flags, but its alerts would still need validation.
Below-quantification records need a censored likelihood, not an ad hoc donor fallback.
Repeated visits, random cluster errors and validation-sample uncertainty remain
future work. No unimplemented topology, causal or calibration feature is advertised.
