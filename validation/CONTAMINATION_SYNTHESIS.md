# Frozen v0.4 contamination evidence and interpretation

The substantive improvement is that recorded values are now used without being trusted. With a
random 25% reference sample, `catmech_contamination()` kept coverage near nominal in three of four
primary cells and shortened intervals by 30% to 38% relative to the v0.3 default, which discards
recorded values. It failed its screen narrowly in the fourth, and it failed under heavy-tailed targets,
selection on the true value with errors present, and undeclared reference noise. These primary cells use
the error form the model assumes, so they favor it by construction. The stress cells are the counterweight.

This report follows the [plan written first](CONTAMINATION_PLAN.md), with the development runs disclosed there.
All data are synthetic. There were **52,950 requests: 52,900 returned intervals, 16 expected refusals, 33
convergence refusals and 1 other error** (a singular raking fit). The run took 2,396 seconds at 14 workers, longer
than my estimate of about 30 minutes. Every source and plan hash matched before and after the run, and the
[independent verifier](../tools/verify-artifacts.py) recomputes every summary from replicate records. A pilot was byte-identical
between 1 and 14 workers ([record](contamination-pilot-check.txt)).

Each entry below is **bias, coverage, mean interval width**, conditional on a returned interval. The estimand
is the coefficient of error-free truth (0.8), with n = 600 and 20 imputations. Primary cells have 1,000 cohorts
and MCSE about 0.7 percentage points. Stress cells have 150 cohorts and MCSE about 1.8 points. Cells sharing a fixture
reuse cohorts and are not independent replications. "Returned" counts the contamination model.

See the [evidence figure](contamination-figures/contamination-evidence.pdf).

## Where the model holds

| Cell | Returned | Contamination | Reference normal | Linear MI | Raking | Flag masking | Oracle |
|---|---:|---:|---:|---:|---:|---:|---:|
| Clean records | 1000/1000 | -0.000, 96.5%, 0.170 | -0.010, 95.4%, 0.277 | -0.008, 95.7%, 0.173 | -0.006, 94.5%, 0.196 | -0.003, 95.3%, 0.173 | -0.002, 95.7%, 0.161 |
| Flagged errors | 1000/1000 | +0.001, 96.3%, 0.175 | -0.010, 95.4%, 0.277 | -0.014, 94.9%, 0.193 | -0.010, 94.6%, 0.209 | -0.004, 95.6%, 0.181 | -0.002, 95.7%, 0.161 |
| Missed flags | 1000/1000 | -0.005, 94.9%, 0.184 | -0.010, 95.4%, 0.277 | -0.019, 92.4%, 0.210 | -0.015, 92.9%, 0.233 | -0.082, 55.9%, 0.175 | -0.002, 95.7%, 0.161 |
| Site shift, declared | 1000/1000 | -0.008, 93.7%, 0.193 | -0.010, 95.4%, 0.277 | -0.018, 91.4%, 0.219 | -0.012, 92.8%, 0.243 | -0.144, 11.8%, 0.173 | -0.002, 95.7%, 0.161 |

| Primary cell | Coverage | 95% Wilson interval | Screen [0.925, 0.975] | Bias (MCSE) |
|---|---:|---:|---|---:|
| Clean records | 0.965 | [0.9517, 0.9747] | Pass | -0.0004 (0.0013) |
| Flagged errors | 0.963 | [0.9494, 0.9730] | Pass | +0.0012 (0.0013) |
| Missed flags | 0.949 | [0.9336, 0.9610] | Pass | -0.0048 (0.0014) |
| Site shift, declared | 0.937 | [0.9202, 0.9505] | **Fail** | -0.0085 (0.0016) |

The site cell's failure stands and was not retuned. Its errors combine a fixed shift with occasional larger errors, which one
normal error distribution per site does not fully describe.

In missed flags, paired against the v0.3 default, mean width fell by 0.0927 (MCSE 0.0013) and the coverage difference was -0.005
(MCSE 0.008). Against linear MI using the recorded value, width fell by 0.0259 (MCSE 0.0007) and coverage rose by 0.025 (MCSE 0.009).
Against generalized raking, width fell by 0.0488 (MCSE 0.0010) and coverage rose by 0.020 (MCSE 0.009). The full-truth oracle
remains narrower, by 0.0229 (MCSE 0.0003). Paired contrasts are in the [study folder](contamination-study/paired-comparisons.csv).

These gains are expected where the assumed error form is true. They do not show the model is better in general.

## Fewer references, rounding and undeclared structure

| Cell | Returned | Contamination | Reference normal | Linear MI | Raking |
|---|---:|---:|---:|---:|---:|
| References 5% | 133/150 | -0.010, 91.7%, 0.207 | -0.049, 96.0%, 0.769 | -0.095, 80.0%, 0.481 | -0.058, 75.3%, 0.525 |
| References 10% | 148/150 | -0.010, 93.9%, 0.200 | -0.027, 95.3%, 0.461 | -0.041, 88.7%, 0.315 | -0.033, 87.3%, 0.332 |
| References 50% | 150/150 | -0.006, 94.0%, 0.174 | -0.006, 96.7%, 0.205 | -0.010, 96.0%, 0.180 | -0.011, 95.3%, 0.190 |
| Rounded records | 150/150 | -0.003, 93.3%, 0.171 | -0.011, 95.3%, 0.282 | -0.010, 94.0%, 0.174 | -0.012, 95.3%, 0.199 |
| Site shift, undeclared | 150/150 | +0.007, 93.3%, 0.201 | -0.011, 95.3%, 0.282 | -0.019, 89.3%, 0.221 | -0.019, 91.3%, 0.246 |

Precision barely depends on reference count in this setting, because correct recordings carry most of the information. That has a cost:
at 5% references, 17 of 150 requests returned no interval (14 unsupported strata, 3 convergence refusals), and coverage among the
rest was 91.7%. Linear MI and raking were unstable at low fractions, likely because a small calibration sample leaves their proxies poorly
fitted. Reference normal held coverage there by being very wide.

## Selective or inaccurate references

| Cell | Returned | Contamination | Reference normal | Linear MI | Raking |
|---|---:|---:|---:|---:|---:|
| Targeted audit | 150/150 | -0.007, 94.0%, 0.181 | -0.075, 80.0%, 0.231 | -0.024, 92.0%, 0.204 | -0.003, 96.7%, 0.239 |
| Selective, clean records | 133/150 | -0.014, 92.5%, 0.176 | +0.263, 7.3%, 0.279 | -0.038, 85.3%, 0.170 | -0.063, 78.7%, 0.198 |
| Selective, missed flags | 150/150 | -0.051, 79.3%, 0.174 | +0.263, 7.3%, 0.279 | +0.076, 74.0%, 0.250 | +0.007, 87.3%, 0.256 |
| Noisy, SD declared | 147/150 | -0.007, 91.8%, 0.190 | -0.124, 52.7%, 0.277 | -0.134, 38.0%, 0.227 | -0.132, 43.3%, 0.248 |
| Noisy, wrong SD | 150/150 | +0.010, 92.0%, 0.204 | not run | not run | not run |
| Noisy, not declared | 150/150 | -0.026, 88.7%, 0.200 | not run | not run | not run |

A targeted audit is now supported, because selection depends on recorded values the model includes. Reference normal lost coverage (80.0%).
Raking also held coverage there, with wider intervals, and it is the design-consistent comparator. Linear MI's fit for records with no
recorded value uses all references, which a targeted design selects, so its 92.0% partly reflects that simple form.

Selection on the true value is not ignorable given recorded values when some recordings are wrong. With no errors present the model held at 92.5%,
but 17 of 150 fits were refused for slow mixing. With errors present it **failed: bias -0.051, coverage 79.3%**. That is still far better than
the v0.3 default (7.3%), but raking did better (+0.007, 87.3%) even with only the average selection rate.

A declared reference SD repaired the v0.3 noisy-reference failure: bias -0.007 and 91.8% against -0.124 and 52.7%. Declaring half the true SD still gave 92.0%,
but expected error counts rose to 186 against 39 actual, so individual error probabilities should not be trusted when that SD is wrong. Leaving the noise
undeclared gave 88.7% and raised the near-universal-disagreement warning in all 150 replicates.

## Assumptions deliberately broken

| Cell | Returned | Contamination | Reference normal | Linear MI | Raking |
|---|---:|---:|---:|---:|---:|
| Heavy-tailed target | 150/150 | +0.033, 87.3%, 0.212 | +0.013, 93.3%, 0.315 | -0.000, 92.0%, 0.239 | -0.006, 92.0%, 0.273 |
| Nonlinear target | 150/150 | -0.001, 94.0%, 0.172 | -0.008, 94.7%, 0.255 | -0.019, 91.3%, 0.194 | -0.022, 84.7%, 0.215 |
| Unequal target variance | 150/150 | +0.013, 93.3%, 0.193 | -0.001, 93.3%, 0.289 | -0.016, 90.0%, 0.218 | -0.024, 85.3%, 0.246 |
| Two-sided errors | 150/150 | -0.029, 90.0%, 0.186 | -0.011, 95.3%, 0.282 | -0.023, 90.7%, 0.217 | -0.025, 92.7%, 0.241 |
| Gross and small errors | 150/150 | -0.012, 94.0%, 0.171 | -0.011, 95.3%, 0.282 | -0.235, 64.7%, 0.256 | -0.284, 53.7%, 0.274 |
| Multiplicative errors | 150/150 | -0.008, 95.3%, 0.183 | -0.011, 95.3%, 0.282 | -0.020, 88.0%, 0.215 | -0.021, 92.7%, 0.248 |
| Errors depend on outcome | 150/150 | +0.002, 92.7%, 0.172 | -0.011, 95.3%, 0.282 | -0.015, 92.7%, 0.205 | -0.015, 92.7%, 0.234 |
| Errors depend on value | 150/150 | +0.011, 94.0%, 0.182 | -0.011, 95.3%, 0.282 | -0.016, 88.7%, 0.203 | -0.018, 87.3%, 0.240 |
| Flag computed from value | 140/150 | -0.011, 94.3%, 0.191 | -0.011, 95.3%, 0.282 | -0.020, 90.0%, 0.219 | -0.021, 92.7%, 0.243 |
| Assay noise, tolerance | 150/150 | -0.014, 94.0%, 0.184 | -0.011, 95.3%, 0.282 | -0.018, 91.3%, 0.215 | -0.018, 93.3%, 0.239 |
| Assay noise, no tolerance | 150/150 | -0.007, 90.7%, 0.202 | not run | not run | not run |
| Recorded missing on value | 150/150 | +0.010, 92.0%, 0.196 | -0.011, 95.3%, 0.282 | +0.011, 92.7%, 0.225 | -0.017, 92.7%, 0.241 |

The plan's failure rule is coverage below 0.90 or absolute bias above 0.03.

- **Heavy-tailed target fails:** bias +0.033 and coverage 87.3%, the worst of the four methods. Genuine extremes get treated as errors and shrunk.
- **Two-sided errors are borderline:** bias -0.029 and coverage exactly 90.0%.
- **No failure by the rule:** nonlinear target, unequal variance, gross errors, multiplicative errors, outcome- or value-dependent errors, range-check flags,
  assay noise and informatively missing recorded values. Each result is one set of 150 cohorts at one mechanism and strength, so this is not a robustness claim.
- **Gross errors favor the model strongly:** linear MI (64.7%) and raking (53.7%) break down because their proxies absorb values ten times too large, while the model treats
  them as errors. Robust proxies might narrow that gap, and they were not tested.
- **Refusals:** range-check flags had 10 of 150 convergence refusals.

The reference normal default remains the better choice when the target's residual distribution is doubtful and recorded values are not needed for precision.

## Per-record error probabilities

These checks use simulated truth, which real data do not have. In missed flags, the model expected 30.1 errors per cohort against 29.0 actual among unreferenced recorded rows.
AUC was 0.962 and the Brier score 0.0305, against 0.0738 for the prevalence alone. Mean probability was 0.629 for wrong records and 0.036 for correct ones, so many errors stay uncertain.
Posterior means had RMSE 0.077 against 0.101 for recorded values. Calibration broke in the gross-and-small cell (8.6 expected against 29.7 actual), where small errors resemble correct values.

## Hold-out check on data alone

| Cell | Is the model wrong? | Alarms | Unresolved (fold failed) | Cohorts |
|---|---|---:|---:|---:|
| Clean records | No | 0 | 0 | 40 |
| Missed flags | No | 0 | 1 | 40 |
| Site shift, undeclared | Yes, strongly | 18 | 0 | 40 |
| Heavy-tailed target | Yes, and coverage failed | 1 | 6 | 40 |
| Nonlinear target | Yes | 0 | 1 | 40 |
| Two-sided errors | Yes | 1 | 1 | 40 |
| Errors depend on outcome | Yes | 1 | 0 | 40 |
| Flag computed from value | Yes | 1 | 6 | 36 |

False alarms were absent under the correct model in this small sample. Power was low: it caught strong structural misfit about half the time and barely registered the heavy-tail failure,
because hold-out rows test prediction rather than the regression coefficient. **A clean hold-out result is weak evidence.** It is a screen, not validation.

## Repairing reference matching

| References | Standard PMM | Bootstrap-donor PMM | Reference normal |
|---|---:|---:|---:|
| 5% | -0.030, 78.0%, 0.530 | +0.047, 97.3%, 1.216 | -0.049, 96.0%, 0.769 |
| 10% | -0.000, 80.7%, 0.295 | +0.029, 93.3%, 0.572 | -0.027, 95.3%, 0.461 |
| 25% | +0.001, 86.7%, 0.213 | +0.008, 95.9%, 0.306 | -0.010, 95.4%, 0.277 |
| 50% | -0.004, 94.0%, 0.188 | -0.001, 96.0%, 0.213 | -0.006, 96.7%, 0.205 |
| 25%, heavy tails | +0.021, 83.3%, 0.268 | +0.030, 96.0%, 0.384 | +0.013, 93.3%, 0.315 |
| 25%, nonlinear | -0.000, 85.3%, 0.194 | +0.012, 97.3%, 0.279 | -0.008, 94.7%, 0.255 |

Standard PMM conditions on one fixed donor set. Resampling donors for each imputation raised coverage from 86.7% to 95.9% in the 25% cell and held at 96.0% in the 50% cell, so the plan's
gate rule was not triggered. The repair is confined to coverage: intervals are wider than the normal method's, and bias reached +0.047 at 5% references and +0.029 at 10%. It is useful
when the normal model is doubtful and references are not sparse. It adds no precision.

## Implementation and scope

An [independent Metropolis oracle](contamination-oracle/README.md), sharing no code with the package sampler, matched all 38 posterior moments and 36 row-level quantities, and it detected all three
deliberately wrong likelihoods. [All 18 older workflow runs](v0.3-unchanged-check.txt) are numerically identical to frozen v0.3.1.

**A bookkeeping bug is disclosed.** The frozen runner overwrote corruption scores and hold-out records in the six cells that also ran reference PMM. Main estimates were unaffected. A
[supplement](contamination-supplement/README.md) recovered those records from the same seeds and refused to write unless every estimate matched the main study, which it did to within 1e-9.

This work supports a reference-calibrated contamination model in the tested synthetic settings. It does not establish a new identification result, general robustness, clinical validity
or reliable individual error labels. [Real-data checks](REAL_DATA_SYNTHESIS.md) are reported separately and are less favorable: there, simple linear MI and raking did as well or better.
