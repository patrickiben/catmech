# Applying the author's related research to catmech

Prepared 2026-10-31. This plan turns related research controls into a concrete
validation fixture and specifies how to test changes intended to reduce
undercoverage. It does not transfer validity between projects. The fixture is
synthetic; it contains no actual trial or laboratory observations.

## What has been adopted

| Source examined | Adopted control | Implementation or evaluation use |
|---|---|---|
| `LLM_CDISC_Eval/R/01_simulate_sdtm.R` | Deterministic synthetic laboratory records, participant identifiers, analyte units, explicit reference limits. | `tools/lab-fixture.R` generates one CREAT result at baseline per participant, with truth separated from analysis data. |
| `LLM_CDISC_Eval/R/04_freeze_truth.R` | Separate artifact existence from value correctness; compute truth independently of the evaluated pipeline. | The fixture writer records row counts, hashes and independently checked values. Imputation code does not create the reference category or oracle coefficient. |
| `LLM_CDISC_Eval/R/05_seed_defects.R` | A claimed defect must visibly alter the data or output it is meant to test. | Self-checks require every injected recording-error scenario to change observed values and verify missed/false flags against separate truth. |
| Residual carryover operating-characteristics study | Separate selection from retained-error effects; use oracle controls and paired repetitions. | Retain matched clean and contaminated versions of each cohort. Report oracle, selection-only and observed/error-retaining comparisons rather than attributing all error to imputation. |
| Residual carryover operating-characteristics study | Random-number streams by scenario group, shared draws within comparisons, MCSE-sized experiments. | Separate calibration, development and held-out streams; paired scenarios share latent truth, outcome noise and missingness. |
| False-Floor evaluation principle | Compare against the strongest simple available alternative. | Keep flag-aware PMM, complete-case exclusion and oracle analysis alongside any new catmech policy. |

The clinical source files were read locally; the new fixture is an independent
implementation of their controls. It takes the synthetic CREAT mean scale,
units and reference limits from that generator. Its larger sample, covariates,
sites, errors, outcome model and residual SD are newly declared stress-test
choices. They are not estimates from real participants. Source SHA-256 values
are recorded in the fixture manifest so the provenance can be checked without
bundling private source directories or data.

The public carryover source inspected was
[`R/analyze.R`](https://github.com/patrickiben/biometrics-ai-portfolio/blob/04a6d4fcd097ab0caf75b585e0daa073ad86081a/sim_lifecycle/carryover_oc/R/analyze.R),
with its [runner](https://github.com/patrickiben/biometrics-ai-portfolio/blob/04a6d4fcd097ab0caf75b585e0daa073ad86081a/sim_lifecycle/carryover_oc/R/run_oc.R).
Its identity `total = retained residual + selection`, paired MCSE reporting,
and known-answer tests motivate the controls here. Its pharmacokinetic
estimators and regulatory decision rules are not applied to laboratory data.

## Clinical endpoints and rounding must be separate tests

The local clinical laboratory generator defines LOW below the lower normal
limit, HIGH above the upper normal limit, and NORMAL at either limit. For this
fixture, LLN is 0.60 mg/dL and ULN is 1.30 mg/dL. These are synthetic fixture
settings, not clinical recommendations.

| Interval | Required convention |
|---|---|
| LOW | `x < LLN` |
| NORMAL | `LLN <= x <= ULN` |
| HIGH | `x > ULN` |

The initial catmech convention puts a value exactly at LLN in LOW. It therefore
cannot represent this clinical convention merely by changing category names.
An explicit clinical-closure option or equivalent interval-specific inclusion
flags would be needed before applying category constraints to this fixture.
The new model can use it with `category_mode = "ignore"`, which imposes no
category boundary convention. Do not move LLN by an arbitrary epsilon.
The support file supplies a separate clinical oracle and endpoint negative
control; it does not silently change package policy.

There is an additional issue in the source generator: `LBNRIND` is calculated
from unrounded values while `LBSTRESN` is rounded to two decimal places. An
underlying value 0.596 is LOW but appears numerically as 0.60, and 1.304 is HIGH
but appears as 1.30. Inclusive clinical endpoints alone do not reconcile these
records. The observation contains rounding information, not an exact equality
between recorded and underlying values. The fixture includes both cases as
fixed controls outside the random cohorts, so deliberately inserted endpoints
cannot distort the coverage experiment.

In an exact-target release, reject or explicitly identify this incompatible
rounding scenario. A later rounding-aware model must treat the reported value
as an interval and intersect it with the independent category information,
with an explicit rounding/tie convention. Neither clipping values to limits
nor re-deriving the reference category from corrupted measurements is a valid
repair of the original evidence.

## Frozen one-analyte fixture

Source `tools/lab-fixture.R` and call `make_lab_fixture()`. The return value has
four parts: `data`, separate `truth`, fixed `boundaries`, and `manifest`.
Only `data` may enter an ordinary imputation procedure. The additional truth
table is available exclusively to the scorer and clearly labeled oracle arms.
Use `write_lab_fixture()` with a new or empty directory to freeze all four
parts and their artifact hashes. Reusing an output directory is refused.

The target is the coefficient of the uncorrupted continuous laboratory value
in `outcome ~ LBSTRESN + age + sex + inflammation + SITEID`. The generative
coefficient is 0.8 synthetic outcome units per mg/dL. The generator adds
independent outcome noise after the true laboratory value and covariates are
fixed; its deterministic zero-noise oracle recovers the coefficient exactly.
This is a known-answer wiring check, not evidence of confidence-interval coverage.

Each participant contributes one baseline CREAT measurement. True values have
a mean near 0.95 mg/dL and conditional residual SD 0.20, broadened from the
local generator's 0.13 to provide useful tail support at moderate sample sizes.
Nonpositive draws are rejected. Observed age, sex, inflammation and fixed site
effects enter the generator. The category is obtained from an independent
synthetic reference before injected errors. Such independent categorical
information must not be presumed to exist in real electronic records.

| Scenario | Recording mechanism | Flag mechanism |
|---|---|---|
| `clean` | MCAR missingness of 20%; no recording error. | No flags. |
| `flagged_errors` | Additive errors on 8% of observed records; shift mean 0.35, SD 0.07 mg/dL. | Sensitivity 1; false-positive probability 0.05 among uncorrupted observed records. |
| `missed_flags` | Same latent errors and missingness as the paired flagged-error cohort. | Sensitivity 0.60; false-positive probability 0.05. |
| `site_shift` | Same individual errors, plus 0.20 mg/dL shift on every observed record at site 003. | Sensitivity 0.20 for any corruption; false-positive probability 0.05 otherwise. |
| `rounded_records` | No additive error; all observed values rounded to two decimals. | No flags; reference categories remain based on unrounded truth. |

Site labels are observed. The shift is a shared recording mechanism; this
fixture does not model random site effects or establish cluster-robust inference.
The site scenario tests whether unrecognized, in-bound errors can enter the
trusted donor pool even when site is included as a predictor. Diagnostic
comparisons can hold out the shifted site or use independent site calibration;
they must not quietly replace observed flags with oracle corruption labels.

## Development, calibration and held-out execution

`phase` selects a full L'Ecuyer-CMRG stream; `replicate` selects a substream.
The three named phases are `calibration`, `development`, and `heldout`.
The same phase/repetition across scenarios uses common latent values, outcome
noise, missingness and error/flag uniforms. Different phases use distinct
streams. The generator restores the caller's RNG state.

1. Use development repetitions for implementation checks and to diagnose the
   current undercoverage. Record all changes and do not call these data held-out.
2. If fitting a corruption model or variance adjustment, estimate it only on
   calibration repetitions. Specify its predictors and tuning rule in advance.
   Real calibration would require independent adjudication or repeated/reference
   measurements; the synthetic labels only test whether the procedure is wired
   correctly under a known model.
3. Freeze policy choices, source hashes, analysis formula, interval construction,
   nominal confidence level, failure definitions, and MCSE target before running
   held-out repetitions. Capture the plan hash before reading results.
4. Run all selected methods on the same held-out cohorts. Report failures and
   reasons alongside coverage conditional on successful fits. Report all declared
   scenarios, including the ones that weaken the proposed method.
5. A change after viewing held-out results requires a new versioned experiment
   with new held-out streams. Repeating a selected favorable seed is development,
   not confirmation.

This document reserves a mechanism and workflow, not an irrevocable statistical
preregistration. It does not claim that the held-out performance study has been
executed. Its final model policies remain the main implementation task's decision.

## Diagnose undercoverage before changing the interval

Compare empirical across-repetition SD with the mean reported SE and inspect
bias relative to interval width. Distinguish point-estimate bias, underestimated
variance, finite-imputation Monte Carlo noise, and donor/model failures.
Widening intervals cannot by itself repair biased retained values or an incorrect
category premise.

For exclusion-only analyses, use the same observed selection mask in two arms:
one analyzes recorded values, the other replaces them by frozen true values.
Their contrast isolates retained recording error conditional on that selection;
the selected-truth arm versus the full-truth arm isolates selection. The three
estimates satisfy an exact within-repetition additive identity. This is a
diagnostic decomposition, not a claim that a selected-truth estimator is usable
on real data.

For imputation, do not attach the same simple decomposition without defining
counterfactual donor pools and completing rules. Instead compare paired oracle
donor and observed donor variants, and label changes to both donor selection
and recipient replacement. These diagnostic variants cannot count as practical
competitors.

The held-out report should include bias, empirical SD, mean reported SE, RMSE,
coverage, interval width, failure proportion, and paired method differences
with MCSEs. At 95% coverage, 1,900 independent successful repetitions give an
MCSE of approximately 0.005; a 100-repetition pilot is not precise enough to
certify nominal coverage. Preserve nominal and observed levels, including
Wilson intervals when plug-in binomial MCSEs are zero.

Success means an improvement survives an appropriate strong comparator and
new repetitions under its stated assumptions. It does not mean that all
misclassification, missed-flag, rounding, or site-error scenarios must pass.
Explicit refusal or a correctly exposed failure can be the valid outcome for
an unsupported measurement mechanism.
