# Development review: substantive gains, unresolved positive-context bias

The first scored development stage is complete and independently verified. It
supports a substantial reduction in the declared coefficient scale-sensitivity
measure and improved contextual coverage, but it does **not** pass the contextual
or ordinary calibration screens. Investigate the residual positive-context bias
before deciding whether to freeze this unchanged candidate for confirmation.
This document is a development assessment, not a confirmation or release claim.

The [prospective plan](JOINT_MARGINAL_PLAN.md) fixed nine cells, targets, priors,
controls and fresh seeds before scoring. The [completion record](joint-marginal-development/COMPLETE.json)
contains 270 cohorts, 1,020 MCMC fits plus 270 complete-truth oracles and 2,310
analysis requests. All requests returned. The [independent verifier](joint-marginal-development/verification.json)
checked 35 source files, 19,380 diagnostic rows, 2,926 summary values, 408 paired
contrasts, every one of 2,000 cohort bootstrap replicates and all four scientific
screens; nine deliberately corrupted copies were rejected. No diagnostic
threshold, width multiplier, nominal level or prior scale was tuned to these
accuracy results.

## Contextual dependence: improvement with remaining bias

All methods retain the same marginal targets: 0.9, 0.9 and zero. The structural
within-site coefficients are 0.7, 0.7 and -0.2 and are never substituted as the
scoring target.

| Contextual cell | Old direct coverage | New marginal coverage | New bias (MCSE) | New mean width |
| --- | ---: | ---: | ---: | ---: |
| Balanced | 19/30 (63.3%) | 29/30 (96.7%) | -0.0491 (0.0180) | 0.4207 |
| Unbalanced | 17/30 (56.7%) | 28/30 (93.3%) | -0.0593 (0.0203) | 0.4184 |
| Marginal null | 19/30 (63.3%) | 30/30 (100%) | -0.0129 (0.0160) | 0.3920 |

The equal-cell coverage-per-request improvement is 35.56 percentage points
(MCSE 5.13; 95% Monte Carlo interval 25.50 to 45.61 points). Equal-cell mean
absolute bias decreases from 0.13164 to 0.04040; its paired bootstrap difference
is -0.09124 (95% interval -0.10412 to -0.07545). These are substantial development
improvements. [Raw and paired results](joint-marginal-development/paired.csv).

Nevertheless, both positive-context cells exceed the prespecified absolute-bias
limit of 0.02; the contextual-null cell also exceeds the upper coverage screen.
As an explicitly post-hoc diagnostic, the equal mean signed bias in the two
positive-context cells is about -0.0542 (MCSE 0.0136; normal 95% Monte Carlo
interval approximately -0.0808 to -0.0276). This is a reason to investigate,
not evidence that the remaining bias can safely be dismissed as 30-cohort
noise. The post-hoc calculation does not replace the three-cell planned screen.
The model's corrected marginal estimand, wider uncertainty and reduced
scale sensitivity do not themselves establish absence of covariance-prior
shape effects or small-sample attenuation.

## Ordinary cells and distributional stress

| Cell | Old direct coverage | New marginal coverage | New bias (MCSE) |
| --- | ---: | ---: | ---: |
| Gaussian balanced | 90.0% | 100% | -0.0188 (0.0159) |
| Gaussian unbalanced | 100% | 100% | -0.0105 (0.0137) |
| Near-zero ICC | 90.0% | 93.3% | -0.0019 (0.0126) |
| Observed-data MAR | 96.7% | 96.7% | -0.0150 (0.0143) |
| Null slope | 90.0% | 90.0% | 0.0075 (0.0155) |
| Student-t5 | 93.3% | 96.7% | -0.0273 (0.0163) |

All four primary ordinary bias bounds pass, but the primary coverage screen
fails: the two 100% cells and the 90% null cell fall outside 92.5--97.5%. With
only 30 cohorts these proportions are imprecise. The equal-cell absolute
distance from 95% is 4.17 percentage points for both procedures; the paired
bootstrap difference interval is -3.33 to 0.83 points. The larger planned
confirmation is needed to estimate calibration more precisely, but should not
be used to choose a new prior or an interval correction. Full Wilson intervals
and all secondary CV3/oracle comparisons remain in the [summary table](joint-marginal-development/summary.csv).

## Coefficient scale-sensitivity and simulation noise

The sensitivity metric is the average absolute shift in the posterior mean at
variance-scale settings 0.25 and four relative to one. All 120 paired prior
triples returned.

| Cell | Old mean coefficient sensitivity | New mean coefficient sensitivity | Paired difference (MCSE) |
| --- | ---: | ---: | ---: |
| Gaussian balanced | 0.01134 | 0.00187 | -0.00947 (0.00203) |
| Near-zero ICC | 0.00865 | 0.00175 | -0.00690 (0.00133) |
| Observed-data MAR | 0.01267 | 0.00247 | -0.01021 (0.00144) |
| Contextual X | 0.01201 | 0.00637 | -0.00564 (0.00212) |

The equal-cell means are 0.0111693 and 0.0031154, a 72.1% reduction in the
measured perturbation response. The paired difference is -0.0080539 (MCSE
0.0008828; 95% Monte Carlo interval -0.0097841 to -0.0063238), with all cell
point differences negative. The planned coefficient sensitivity screen passes.
The secondary relative-width sensitivity also decreases. These compare two
models and prior families over declared multiplicative scale ranges; they do
not isolate a single causal contribution of the prior family.

Finite-chain Monte Carlo error matters for interpreting the very small new
shifts. Averaged over the three prior fits, the new posterior-mean MCSE is about
0.00219 in balanced Gaussian, 0.00141 near zero ICC, 0.00274 in MAR and 0.00286
in contextual X. Some measured shifts are comparable to this simulation error.
The old MCSEs are about 0.00171, 0.00156, 0.00234 and 0.00194. Shared seeds can
correlate simulation errors, so these single-summary MCSEs cannot simply be
subtracted from absolute shifts. The evidence supports a reduced measured
sensitivity for the declared fitting procedure, not exact prior independence
or a noise-free estimate of limiting-posterior sensitivity.

## Numerical availability and the next decision

Both old and new models returned all 510 requested fits per model, including
prior alternatives. Every direct and CV3 analysis returned. The unchanged
R-hat/ESS screens passed and the availability nonregression screen passes.
The first masked pilot's four refusals remain saved separately; the exact
interweaving update removed those selected numerical failures without weakening
criteria. Such a selected rerun is not an estimate of future refusal frequency.

Before a confirmation freeze, examine whether the remaining contextual bias is
an implementation issue, a finite-sample feature of the joint covariance prior,
or a limitation of the chosen model. An unchanged-candidate confirmation can
measure that limitation more precisely; it cannot be described as already
repairing it. Any scientific candidate change requires a dated amendment,
preserved current evidence, another distinct development seed set and a new
source freeze before its accuracy run. The current confirmation seeds remain
untouched.
