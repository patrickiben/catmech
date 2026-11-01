# Prospective v0.4 contamination evaluation

Written 2026-11-01 before any run on the held-out streams below. Earlier v0.1 to v0.3
results remain versioned evidence and are not relabeled.

## Questions

The v0.3 default treats every unreferenced target as missing and ignores the recorded
value W. In the v0.3 fixture most recorded values are exactly right, so this discards
a great deal of information: reference-normal interval width was 0.278 against 0.161
for the full-truth oracle.

1. Given the same reference rows, does an explicit recording-error model that uses W
   keep nominal coverage while shortening intervals, relative to the v0.3 default?
2. Does it beat the strongest simple alternatives that also use W: a linear
   validation-sample imputation (Cole, Chu and Greenland, 2006) and generalized raking
   (Lumley, Shaw and Dai, 2011; Oh et al., 2021)?
3. How much does it depend on its own assumptions? It assumes additive normal errors
   mixed with exact recordings, a linear normal target model, nondifferential errors
   and ignorable selection given W. The primary cells use additive normal errors, so
   they favor this model by construction and cannot by themselves show robustness.
   The wrong-model and credibility cells are the deliberate counterweight.
4. Does bootstrapping PMM donors repair the 85.1% partial-reference coverage in v0.3?
5. Does `catmech_holdout_check()` catch wrong models on data alone, and how often does
   it alarm when the model is right?

This is an established contamination-model idea (Di Zio and Guarnera, 2013) applied with
validation data. It is not presented as a new identification result.

## Development disclosure

Before this plan was frozen, 60-replicate development runs used the `development` phase of
fixture master 410273 and streams 170,100,000/175,000,000/176,100,000/178,100,000, which are
disjoint from the held-out streams. Earlier single-replicate smoke runs used master 510001.
Those runs changed the following, and these changes are part of the frozen procedure:

- The first sampler mixed slowly when errors were rare. It was replaced by slice updates of
  the recording parameters with status and target integrated out, plus a scale move across
  the prior's funnel when fewer than five errors are observed.
- The error prior was tightened from `kappa0 = 0.01, nu0 = 1` to `kappa0 = 0.1, nu0 = 2`,
  because a nearly improper prior caused convergence refusals on data with no errors.
- A restart rule doubles chain length after a failed screen, up to `max_iterations`. Default
  length became `burnin = 200, iterations = 300, max_iterations = 1200` to bound run time.
- The Cole comparator now drops aliased columns and needs ten flagged calibration rows before
  adding flag terms. Without this it was unstable at sparse references, which would have
  flattered the new method.
- PMM variants were screened at 10%, 25% and 50% references. Type-1 matching on bootstrapped
  donors was chosen over type-0 and over `mice.impute.midastouch`.

These development results are not reported as evidence. Cell definitions and comparators were
written before the first development run and were not removed because of its results.

## Design

n = 600, 20 imputations, analysis `outcome ~ lab_value + age + sex + inflammation + SITEID`,
estimand the coefficient of error-free truth (0.8). Cohorts come from `tools/lab-fixture.R`,
unchanged, with master seed 410273, `heldout` phase. Mechanisms absent from that fixture are in
`tools/contamination-fixtures.R`. The default corruption stratum is the observed QC flag.

| Family | Cell | Mechanism and declared model | Replicates |
|---|---|---|---:|
| primary | clean | No errors | 1000 |
| primary | flagged_errors | Flag sensitivity 1, 5% false positives | 1000 |
| primary | missed_flags | Flag sensitivity 0.6, 5% false positives | 1000 |
| primary | site_shift_declared | +0.20 at site 003 plus random errors; corruption and error strata by site | 1000 |
| support | reference_05pct, 10pct, 50pct | missed_flags at other random fractions | 150 each |
| support | rounded_records | W rounded to 0.01; `agreement_tol = 0.005` | 150 |
| support | site_shift_undeclared | Same data, flag stratum only | 150 |
| credibility | targeted_audit | Selection probability 0.80 if flagged, 0.65 if unflagged and \|z(W)\| > 1.5, otherwise 0.15; known to the analyst | 150 |
| credibility | nonignorable_clean, nonignorable_missed | Selection on latent target residual, slope 2, mean 25%, as in v0.3 | 150 each |
| credibility | noisy_declared | Reference error SD 0.08 declared correctly | 150 |
| credibility | noisy_wrong_sd | Same data, declared 0.04 | 150 |
| credibility | noisy_undeclared | Same data, treated as exact | 150 |
| wrong_model | heavy_tail_target | Target residuals standardized t(3) | 150 |
| wrong_model | nonlinear_target | Target mean quadratic in inflammation | 150 |
| wrong_model | heteroskedastic_target | Target SD depends on inflammation | 150 |
| wrong_model | error_bimodal | Errors ±(0.35, SD 0.07) | 150 |
| wrong_model | error_gross_and_small | 25% of errors are W = 10T; the rest small, mean 0.10 | 150 |
| wrong_model | error_multiplicative | W = 1.4T | 150 |
| wrong_model | differential_outcome | Error probability depends on outcome | 150 |
| wrong_model | value_dependent_error | Error probability depends on truth | 150 |
| wrong_model | flag_from_recorded | Flag is a range check on W itself | 150 |
| wrong_model | assay_noise_tolerance | Clean records have SD 0.02 noise; `agreement_tol = 0.06` | 150 |
| wrong_model | assay_noise_no_tolerance | Same data, no tolerance declared | 150 |
| wrong_model | recorded_missing_mnar | W missing depends on truth | 150 |

Stress cells are descriptive at 150 replicates (MCSE about 1.8 points at 95% coverage). This
size keeps the whole run near 30 minutes, and a larger rerun would use new streams.

Methods in each cell, all given the same data: full-truth `oracle`; `gold_only` regression;
`v02_mask_normal`; `reference_normal` (v0.3 default); `cole_mi`; `raking`; `contamination`.
Variant cells that change only the model declaration (noisy_wrong_sd, noisy_undeclared,
assay_noise_no_tolerance) run `contamination` alone. `reference_pmm` and `reference_pmm_abb`
also run in missed_flags, reference_05pct, 10pct, 50pct, heavy_tail_target and nonlinear_target.
Raking receives the true selection probabilities where the analyst would know them; under
nonignorable selection it receives only the average rate.

Seeds: reference selection 180,100,000 + replicate, noise 186,100,000 + replicate, mechanisms
188,100,000 + replicate, imputation 190,000,000 + replicate. Uniforms are shared across cells
for pairing, so cells sharing a fixture are not independent replications.

`catmech_holdout_check()` runs with five folds on the first 40 replicates of clean,
missed_flags, site_shift_undeclared, heavy_tail_target, nonlinear_target, error_bimodal,
differential_outcome and flag_from_recorded.

## Measures and criteria

Report all requests and statuses, returned intervals, conditional coverage, Wilson intervals,
bias, RMSE, empirical SD, mean SE and width, each with MCSE, and paired differences against each
comparator. Report corruption-probability calibration against simulated truth, which is available
only because these data are synthetic.

Primary screen, carried over from v0.3: for `contamination` in each primary cell, the 95% Wilson
interval for coverage lies within [0.925, 0.975]. A failed screen is reported, not retuned.

Claims are limited in advance as follows.

- A precision gain is claimed only where coverage holds, against comparators given the same data.
- A wrong-model cell is called a failure when coverage falls below 0.90 or absolute bias exceeds
  0.03, about half the oracle standard error. Failures stay in every summary.
- If `reference_pmm_abb` falls below 0.90 in missed_flags or reference_50pct, it receives the same
  sensitivity gate as `reference_pmm`, added as a scoped post-study guard with its numerics unchanged.
- Hold-out detection rates are descriptive at 40 replicates; alarms under the correct model are
  reported as the false-alarm rate.

Failures are recorded without retrying until favorable. A change after viewing held-out results
requires a new versioned experiment with new streams. All data are synthetic; nothing here is
clinical validation or a general coverage guarantee.
