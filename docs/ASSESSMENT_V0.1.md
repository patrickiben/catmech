# Assessment of the catmech prototype

Assessed: 2026-10-31. Scope: the supplied `catmech-local` v0.9 source,
README, example specification, saved example outputs, and QC report, followed
by a targeted review of prior literature and package documentation. The original
draft was preserved. This document assesses the starting idea and its defensible
release scope; it does not certify the final package build or report its check results.

## Recommendation

Continue as a small, experimental R package for auditable sensitivity completion.
The useful contribution is an explicit workflow for combining trusted interval
categories with uncertain measurement-quality flags. The evidence does not
establish a new imputation methodology, recovery of individual latent truths,
or generally calibrated downstream confidence intervals.

| Question | Judgment |
|---|---|
| Is the statistical idea new? | Its main components are established. An original statistical contribution has not been demonstrated. |
| Is a package useful? | Yes, if it validates assumptions, handles contradictions explicitly, preserves input values, and records what was changed and why. |
| Does the original QC establish validity? | It establishes limited example behavior. It misses counterexamples and overstates some conclusions. |
| Is a GitHub research release worthwhile? | Yes, with a narrow experimental scope and reproducible checks. |
| Is a methods-paper claim ready? | No. It needs a coherent target, stronger comparisons, and repeated-simulation evidence. |

A suitable description is: *An R package for assumption-based sensitivity
completion of continuous measurements with trusted interval categories and
quality flags, using established imputation components with explicit checks and
audit records.*

## Prior art and the actual contribution

The following is a targeted scan, not a systematic review or proof of priority.
Searches covered interval/coarsened data imputation, measurement error,
multiple overimputation, predictive mean matching, and sensitivity analysis.

| Existing work | Established capability | Consequence for catmech |
|---|---|---|
| `mice` | Custom imputation methods, selected-cell overimputation through `where`, passive consistency, analysis, and pooling. | Integration with these facilities is useful engineering, not a new MI framework. |
| `mice` predictive mean matching | Regression-based donor selection and random choice among nearby donors. | Restricted donor eligibility is an adaptation of existing PMM. |
| `hmi` | Imputation of interval responses, rounding-related measurement error, and a workflow compatible with `mice`. | Broad claims to introduce interval-aware MI or MI for measurement error are untenable. |
| Multiple overimputation and `Amelia::moPrep()` | Joint treatment of missingness and measurement error; priors for selected mismeasured cells, with validation or error-size information. | Replacing suspect observed measurements is established. |
| `inlamemi` | Hierarchical models for missingness and classical/Berkson measurement error, using error-variance information or priors. | A fuller measurement model is an existing alternative when inferential correction is the goal. |
| `mice` NARFCS | User-specified MNAR sensitivity offsets through its imputation interface. | A future delta sensitivity feature would extend convenience, not establish novelty by itself. |

Sources: [mice](https://amices.org/mice/reference/mice.html),
[PMM](https://amices.org/mice/reference/mice.impute.pmm.html),
[hmi paper](https://doi.org/10.18637/jss.v095.i09),
[multiple overimputation](https://doi.org/10.1177/0049124115585360),
[Amelia documentation](https://search.r-project.org/CRAN/refmans/Amelia/html/moPrep.html),
[inlamemi](https://emmaskarstein.github.io/inlamemi/), and
[NARFCS](https://amices.org/mice/reference/mice.impute.mnar.html).

`hmi` deserves particular attention as a strong comparator. Its published work
covers interval responses with covariates, not merely imputation of missing
ordinal labels. Inspection of the official archived version 1.0.0 confirms that
`imp_interval()` accepts lower/upper bounds, including infinite endpoints, and
uses parameter draws and truncated-normal completion. This is close to the
core data problem here. Its rounding-error model is not identical to catmech's
arbitrary flag-based retain/replace rule. Neither paper nor source inspection
establishes that it offers the same assumption ledger and conflict policies.
See the [paper](https://doi.org/10.18637/jss.v095.i09) and
[archived source](https://cran.r-project.org/src/contrib/Archive/hmi/hmi_1.0.0.tar.gz).

CRAN records that `hmi` was archived on 2022-05-10 because it required the
archived `linLIR` package. That creates a possible maintenance and usability
niche, but does not erase the prior art or establish superior statistical
performance. No comparative execution of `hmi` was performed in this review.
[CRAN status](https://cran.r-project.org/web/packages/hmi/index.html).

The plausible contribution is the combination of a small validated specification,
an explicit rule for which source takes precedence when observations disagree,
category-constrained completion, sensitivity over flag treatment, and machine-readable
provenance. A hash is a useful identity check for an artifact; it is not evidence
that the assumptions or imputations are correct. No claim that this combination
is the first of its kind is warranted by this scan.

## Target and validity boundary

For a participant, let `T` be a latent continuous measurement, `W` its recorded
value, `C` an observed category, `F` a quality flag, and `Z` the conditioning
variables. A hard-category declaration assumes `T` lies in the interval named
by `C`. It does not identify the particular value of `T` within that interval.
The imputation procedure supplies additional distributional assumptions.

The implemented quantity is therefore a distribution of completed datasets
under an analyst-selected rule. A downstream regression coefficient is a
separate estimand and must be specified explicitly. Reliable estimation of that
coefficient requires adequate imputation models and relevant conditioning
information, including the outcome and needed interactions where appropriate.

Several assumptions do most of the work:

1. Categories are accurate measurements of the same underlying quantity and
   refer to the same time and units. A category derived from a corrupted value
   is not independent corroborating information.
2. The trusted observations provide adequate information about the conditional
   distribution in the rows being completed. Sparse categories and selection
   into the trusted subset can defeat this assumption.
3. The flag policy describes how uncertainty should be explored. An unflagged
   or missing flag does not empirically prove that its measurement is correct.
4. The missingness and regression assumptions support the intended analysis.
   Repeated plausible draws alone do not establish proper MI or valid pooling.

In the draft, `corruption_probability` is an imposed reconstruction probability
for eligible flagged values. There is no specified and fitted observation-error
likelihood that turns it into a posterior probability of corruption. The
retain/replace rule can be useful for sensitivity analysis, but should be
described that way.

If an observed value conflicts with a hard category, the pair is inconsistent.
Replacing the value resolves the inconsistency only under the assumption that
the category takes precedence. The conflict alone does not determine whether
the category, the value, or their linkage was wrong.

For a flagged observed subset, let `f` be the fraction outside its declared
interval and let `q` be the reconstruction probability for the remainder.
If conflicts force reconstruction, the expected reconstruction fraction is
`f + (1 - f) * q`. A realized fraction need not exceed `q`. Also, reconstructing
a cell can draw its original value again. Reconstruction decisions, numerical
changes, and threshold violations are different diagnostics.

## Confirmed problems in the supplied v0.9 draft

These are findings about the original source, not a statement that they remain
in the refined package. Small direct probes sourced the original utility file
without modifying it or the saved example outputs.

| Finding | Evidence and implication |
|---|---|
| Singleton donor selection can escape the eligible pool | `pick_pmm_donor()` uses `sample(nearest, 1)`. When the sole eligible row is 8, R interprets the scalar as `1:8`. With seed 1, 20 direct draws returned only rows 1 through 7. Use indexed `sample.int()` selection. |
| Threshold endpoints are inconsistent | `cut(right = TRUE)` assigns 3.5 to Low, while the original bounds reject 3.5 for both Low and Medium. Every finite boundary needs one explicit inclusion convention used throughout. |
| Hash canonicalization loses distinctions | Direct probes show identical hashes for `list(x = Inf)`, `list(x = -Inf)`, and `list(x = NA_real_)`. The normalizer converts both infinities to NA. A hash over that representation cannot establish assumption identity. Numeric rounding also discards sufficiently small changes. |
| Trusted observations can violate the claimed hard constraint | The runner marks missing and flagged measurements for processing, but does not check all observed unflagged measurements against categories. The original zero-violation example does not prove a general invariant. |
| A diagnostic has the wrong meaning | Saved `threshold_violations.csv` reports zero violations with `violation_rate = 1` in every imputation. The rate is a consistency rate, despite its name. |
| The specification overstates implemented behavior | `both_missing: joint_impute` is declared but the runner rejects both-missing rows. `conflict_policy` is not enforced for all rows; `allow_true_extreme` and `min_donors` do not implement their apparent promises. |
| Fallback changes the trust model | With too few trusted rows, the imputer falls back to observed rows, which can include suspect measurements. That weakens the stated donor protection and should not occur silently within an assumption-based workflow. |
| Missing predictors remain problematic | Dropping incomplete predictor rows from fitting does not make predictions for incomplete recipient rows finite. Input validation or an explicit predictor-imputation workflow is needed. |
| The simulation does not test uncertain flags | Every flagged observed measurement is actually corrupted, and categories are generated from uncorrupted truth. The uncertain-flag model and hard-category assumption are not challenged. |
| The runner is not a general specification compiler | It hard-codes column names and generates the example rather than reading arbitrary data from the declared data-file field. The public interface must state its supported scope. |

The original QC report also contains unsupported API-history explanations.
Returning a full-length vector was a contract error, not an established change
introduced in `mice` 3.19. The official `mice` 3.0.0 source already has a `wy`
argument and returns only `sum(wy)` draws; version 2.25 already returns
missing-cell draws rather than a full vector. These archives were inspected
directly: [mice 3.0.0](https://cran.r-project.org/src/contrib/Archive/mice/mice_3.0.0.tar.gz),
[mice 2.25](https://cran.r-project.org/src/contrib/Archive/mice/mice_2.25.tar.gz).
The accompanying stories about R 4.0 and `quickpred()` API drift were not
historically verified and should not be repeated as facts. State the observed
contract errors and their fixes instead.

Passing the supplied assertions is not enough to validate the method: the
singleton-donor and hash counterexamples pass outside their coverage. Likewise,
the one-seed regression example is a demonstration, not evidence of general
bias reduction or nominal confidence-interval coverage. The reported
coefficient comparison should remain labeled as an original smoke test.

## Useful release scope and next research steps

For the first release, support one continuous measurement with known category
intervals and optional quality flags. Reject unsupported settings, preserve the
original data, expose conflicts, make the replacement policy explicit, and
record the specification, software versions, random seed, and diagnostics.
Keep speculative topology, contamination-ring, and interdiction extensions
outside the advertised implementation.

Before claiming an inferential improvement, pre-specify a simulation study
using the ADEMP structure: aims, data-generating mechanisms, estimands, methods,
and performance measures. Report Monte Carlo standard errors and unsuccessful
runs, rather than selecting favorable examples.
[Morris, White, and Crowther (2019)](https://doi.org/10.1002/sim.8086).

The study should compare ordinary PMM, flag-aware PMM after masking suspicious
values, interval-constrained completion without uncertain retention, and an
oracle analysis using the simulated truth. Where feasible, include a fitted
interval or measurement-error model such as the relevant `hmi` approach.
Evaluate bias, RMSE, coverage, interval width, failure rates, and constraint
violations. Use paired repetitions for comparative Monte Carlo uncertainty.

Stress conditions should include no corruption, false-positive and false-negative
flags, within-interval corruption, incorrect categories, rare categories,
nonlinear relationships, heavy-tailed errors, and selection into missingness.
Choose the number of repetitions from the desired MCSE. A small development
pilot can reveal failures but cannot settle small performance differences.

A further methods contribution would require deriving a coherent observation
model for the category and corrupted measurement, explaining what identifies
its parameters, and assessing how its uncertainty propagates to analysis.
Validation measurements, repeated measurements, or justified external error
information would strengthen that work. Until then, the appropriate claim is
transparent sensitivity analysis under declared assumptions.

The assessment and refinement were AI-assisted. The author remains responsible
for confirming the intended scientific use, reviewing the assumptions, and
deciding whether future evidence supports a publication claim. Package release
checks and executed simulation results should be read in their own versioned
validation records.
