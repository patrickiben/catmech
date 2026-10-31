# Prospective v0.3 extension evaluation

Written 2026-10-31 before running the new engines. Earlier v0.1 and v0.2 results
remain versioned evidence. This study adds information and model comparisons;
it does not relabel the earlier failures as repaired without new measurements.

## Questions and identification

1. In the v0.2 heavy-tail comparison, PMM also received category covariates that
   the default normal model did not. Separate the effect of method from that
   additional information in a complete 2x2 comparison.
2. Can a genuinely measured, randomly selected internal reference sample recover
   inference when unflagged recorded measurements are corrupted? What precision
   is lost, and how does performance compare with analyzing only reference records?

The reference workflow treats true measurement as observed ONLY where the reference
column is measured. Every other target is imputed, including unflagged recorded
values. Recorded W is retained as evidence, not silently promoted to truth. The
primary reference models use complete covariates, outcome and site, and ignore W.
Random internal validation makes this a conventional missing-true-value problem
under the imputation/analysis model. It does not identify truth when no reference
measurements exist. Targeted validation requires an appropriate selection model;
an inaccurate reference assay is not a gold standard.

This is established measurement-error MI, not a new identification theorem. See
[STRATOS Part 2](https://doi.org/10.1002/sim.8531). PMM is an alternative to test,
not a universal robustness guarantee; a [2025 primary preprint](https://arxiv.org/abs/2506.22981)
provides MAR counterexamples motivating the strong-MAR negative control below.

## Distribution and information ablation

Use n=300, 20 imputations, five donors for PMM, and 1,000 replicates per setting.
Use the existing simulator and its coefficient 0.8 (zero in the null setting),
thresholds (-Inf,5.5,7,Inf), and the following settings:

| Setting | Latent error | Recorded missingness | Flags/corruption | Category error |
|---|---|---|---|---:|
| missing_only | normal | 20% MCAR | none | 0 |
| perfect_flags | normal | 20% MCAR | 10% flags, all corrupted | 0 |
| false_positive_flags | normal | 20% MCAR | 10% flags, half corrupted | 0 |
| category_error | normal | 20% MCAR | 10% flags, all corrupted | 10% |
| null_effect | normal | 20% MCAR | 10% flags, half corrupted | 0 |
| heavy_tails | standardized t(3) | 20% MCAR | 10% flags, all corrupted | 0 |
| heavy_category_error | standardized t(3) | 20% MCAR | 10% flags, all corrupted | 10% |
| strong_mar | normal | logistic selection on inflammation | none | 0 |

For strong MAR, P(missing)=logistic(3*(inflammation-2)); missingness depends only
on a complete predictor and averages about 50%. Categories and outcome remain
observed. The normal conditional model remains correctly specified; PMM must
handle limited donor overlap. This stress is not retuned to make either method win.
Adding an exact category as an ordinary factor makes the latent conditional law
truncated and its mean nonlinear even in a Gaussian base DGM. Consequently the
normal-plus-category cell is an information/model ablation, not a correctly
specified Gaussian conditional model.

Cross normal versus PMM with category-as-covariate absent versus present. All four
methods mask all flags, and all receive age, sex, inflammation and outcome. The
category is an ordinary factor predictor, not a hard constraint. Include the same
four combinations implemented independently with mice.impute.norm/pmm, plus the
full-truth oracle and trusted observed-record regression. Compare methods within
fixed information sets and category inclusion within fixed methods.

## Internal-reference study

Use tools/lab-fixture.R unchanged, n=600, its heldout stream under NEW master seed
390173. Primary settings are clean, flagged_errors, missed_flags, site_shift and
rounded_records, each with 1,000 replicates and independent Bernoulli(.25) reference
sampling. Reference values are the exact simulated truth only on selected rows;
all nonselected values are genuinely NA in the supplied reference column.

The same cohort and reference-selection uniforms are reused across recording-error
scenarios. Therefore comparisons isolate recording errors, and these scenarios are
not five independent replications of the same reference model. Ignoring W makes
reference-model outputs invariant to changes in W when cohort, reference subset,
predictors and imputation seed are unchanged. This is intended and checked.

Comparators: full-truth oracle; the v0.2 flag-masking normal workflow without gold;
reference-only regression; new reference-normal and reference-PMM workflows; and
independent mice normal/PMM imputers given exactly the same reference rows and
predictors. Include outcome and SITEID. Do not give hidden truth to any method
except the oracle or the selected reference cells. No category covariates or raw-W
auxiliaries enter the primary reference imputation models.

The site_shift fixture is a compound setting: a fixed +.20 shift at one site,
8% random errors and only 20% flag sensitivity. SITEID is in the analysis; this
study does not isolate the shift's incremental effect.

## Support and assumption-violation stresses

Use 250 replicates each for 10% and 50% random reference fractions in missed_flags
and site_shift. These are precision/support sensitivity analyses, not chosen after
seeing the results.

Use 1,000 replicates for a nonignorable reference-selection control in clean data.
Fit the simulation-only full-truth regression on the selected complete predictors
and outcome to construct a standardized latent residual. Reference probability is
logistic(a + 2*residual), with a fixed by root-finding to give mean probability .25.
The hidden residual drives selection beyond supplied predictors; it never enters
the imputer. This preserves the expected reference fraction while violating the
reference missingness assumption. No method is required to pass this control.

Use 500 replicates of missed_flags with a 25% reference sample whose supplied
reference measurement has independent N(0,.08^2) error. This explicitly violates
reference accuracy. The estimand remains the coefficient of error-free truth.

All-reference and no-reference controls use 20 replicates each. All-reference
completed targets, coefficient, standard error, degrees of freedom and confidence
interval must agree with the direct full-truth linear analysis. No-reference calls
must refuse to produce imputed truth; the presence of recorded W is not a substitute.

## Seeds, engineering and freezing

Distribution data seeds: 60,000,000 + 100,000*setting_index + replicate. Pilot seeds
use 50,000,000 instead. Null outcomes use a separate +125,000 stream; strong-MAR
masking uses +42,000. Reference fixtures use master 390173, independent of the
v0.2 master. Reference selection uses seed 80,100,000+replicate (pilot 70,100,000),
with shared uniforms across recording scenarios and reference fractions. Reference
assay noise uses 86,100,000+replicate. Reference imputation uses 90,000,000+replicate
(pilot 75,000,000), shared across scenarios to test W-invariance.

Five-replicate engineering pilots check API, runtime, source stability and worker
schedule independence. They do not select scientific settings using coverage or
bias. All held-out settings and seeds are frozen before execution. Failures are
recorded without retrying until favorable. The default four independent workers
have explicit per-replicate seeds; a serial/parallel pilot comparison is required.

## Measures and interpretation

Report all requests and statuses, returned intervals, conditional coverage, interval
yield, Wilson intervals and failure rates. For successful estimates report bias,
RMSE, empirical SD, mean model SE, mean interval width, and MCSEs for headline
performance measures. Report paired differences with complete-pair counts.
Separately assess the 2x2 method/information contrasts.

For reference studies, paired oracle decomposition separates selection of reference
records from error in supplied reference measurements. A retrospective bias-centered
coverage diagnostic may be reported with bootstrap uncertainty; it uses simulation
truth and is not a real-data correction. All-reference numerical equivalence and
no-reference refusal are functional controls, not coverage claims.

Primary Gaussian coverage criterion: the 95% Wilson interval is contained in
[.925,.975] at 1,000 replicates. This criterion applies to normal-without-category
in the six Gaussian distribution settings and to normal-reference MI with random
accurate gold at 25%; the clinical Gaussian approximation includes rare positivity
rejection. Added-category, PMM and deliberate stress cells are descriptive. Report failure to meet it without tuning. Smaller
support studies and deliberate assumption violations receive descriptive results,
not forced pass labels. PMM gains must be compared at equal information and cannot
be generalized to arbitrary MAR selection or sparse reference samples. Reference
improvement over v0.2 receives additional measured information; do not present it
as an algorithm solving missed corruption without that information.
