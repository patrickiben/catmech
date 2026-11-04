# D2 decision: freeze the factorized statistical candidate for fresh confirmation

Proceed to the prespecified fresh confirmation with the conditional-regression
covariance prior, marginal target, default scales and numerical controls fixed
in Amendment 2. This is a decision to test the candidate more precisely, not a
claim that its calibration screens have already passed. The separately
identified one-paired-site support correction and its input-equivalence checks
are complete; freeze the corrected source before confirmation. That correction
does not change the statistical specification or the planned scoring.

D2 used source commit `288febf` on 270 new cohorts, with 1,260 MCMC fits plus
270 complete-truth oracle workflows and 2,790 analyses. All fits and analyses
returned. The [independent verification](joint-marginal-development-factorized/verification.json)
checked all 35 source files, 32,880 diagnostic rows, 3,534 summary values,
504 paired contrasts, all 2,000 cohort bootstrap replicates and all four
screens; nine corrupted copies were rejected. The [prospective plan](JOINT_MARGINAL_PLAN.md)
and [earlier failed candidate](JOINT_MARGINAL_DEVELOPMENT_REVIEW.md) remain
part of the evidence. No interval multiplier, nominal level or prior-scale
grid was selected to make these development proportions equal 95%.

## Contextual evidence and limits

| Cell and marginal target | Old direct coverage | New marginal coverage | New bias (MCSE) |
| --- | ---: | ---: | ---: |
| Contextual balanced, 0.9 | 18/30 (60.0%) | 29/30 (96.7%) | 0.0247 (0.0175) |
| Contextual unbalanced, 0.9 | 14/30 (46.7%) | 30/30 (100%) | -0.0336 (0.0182) |
| Contextual null, 0 | 17/30 (56.7%) | 28/30 (93.3%) | 0.0101 (0.0183) |

The equal-cell coverage-per-request gain is 42.22 percentage points (MCSE 5.70;
95% Monte Carlo interval 31.06 to 53.39 points). Equal-cell mean absolute bias
decreases from 0.14231 to 0.02280. Its paired bootstrap difference is -0.11951
(95% interval -0.13928 to -0.08975). The marginal targets are identical for both
procedures; the within-site coefficients are not used as substitute targets.

The contextual screen nevertheless **fails in D2**: one cell covers 100%,
and the two positive-context point biases exceed the absolute 0.02 bound.
These should remain explicit warnings for confirmation. Unlike the first
candidate's consistently negative positive-context biases, the D2 positive-
context biases have opposite signs and each is imprecise. As a post-hoc
orientation check only, their equal mean signed bias is about -0.00444 (MCSE
0.01264). This does not prove that bias is absent or permit changing the
prespecified per-cell screen. D1 and D2 use different cohorts, so differences
between those stages do not isolate a causal prior-factorization effect.

## Ordinary calibration also requires confirmation

| Cell | Old direct coverage | New marginal coverage | New bias (MCSE) |
| --- | ---: | ---: | ---: |
| Gaussian balanced | 96.7% | 96.7% | 0.0057 (0.0151) |
| Gaussian unbalanced | 93.3% | 93.3% | -0.0220 (0.0158) |
| Near-zero ICC | 96.7% | 96.7% | -0.0018 (0.0124) |
| Observed-data MAR | 100% | 100% | -0.0124 (0.0138) |
| Null slope | 93.3% | 100% | 0.0043 (0.0154) |
| Student-t5 | 100% | 100% | 0.0013 (0.0116) |

The four-cell ordinary screen **fails in D2**: MAR and null coverage exceed
97.5%, and the unbalanced bias slightly exceeds 0.02. The baseline's unbalanced
bias is also negative, -0.0164 (MCSE 0.0154). The equal-cell absolute coverage
distance from 95% is 3.33 percentage points for the candidate versus 2.50 for
the baseline; the paired bootstrap difference interval is -1.67 to 0.83 points.
With 30 cohorts, a 100% proportion has substantial uncertainty. The planned
400-cohort cells must resolve the remaining calibration uncertainty; their
outcomes cannot be used to choose narrower intervals or a new prior.
[Full summaries, including CV3 and the extra-information oracle](joint-marginal-development-factorized/summary.csv).

## The strengthened coefficient-sensitivity screen passes

For every sensitivity cohort, the old model contributes three inverse-Wishart
scale fits. The candidate contributes three HW-scale fits and two additional
regression-scale fits. Its score is the **larger** of its two mean absolute
coefficient shifts, calculated within each cohort before averaging. All 120
common old-three/new-five sets are available.

| Cell | Old coefficient sensitivity | New strongest-axis sensitivity | Difference (MCSE) |
| --- | ---: | ---: | ---: |
| Gaussian balanced | 0.00801 | 0.00362 | -0.00439 (0.00116) |
| Near-zero ICC | 0.00809 | 0.00092 | -0.00716 (0.00123) |
| Observed-data MAR | 0.01036 | 0.00432 | -0.00603 (0.00145) |
| Contextual X | 0.01089 | 0.00299 | -0.00791 (0.00166) |

Equal-cell means are 0.0093372 and 0.0029631: a 68.3% reduction in the measured
perturbation response. The paired difference is -0.0063741 (MCSE 0.0006943;
95% Monte Carlo interval -0.0077350 to -0.0050133), with every cell point
contrast negative. Each axis and its signed shifts remain in the
[prior results](joint-marginal-development-factorized/prior-summary.csv).
The secondary maximum-axis width sensitivity also decreases; it does not
substitute for the coefficient result.

Finite-chain uncertainty remains important. Mean candidate posterior-mean MCSEs
are approximately 0.00241, 0.00128, 0.00268 and 0.00255 in these four cells,
respectively. Several new measured shifts are of comparable size. Shared
random-number streams can correlate their Monte Carlo errors; no simulation
error is subtracted from absolute shifts. This supports reduced sensitivity
for the declared numerical procedure and scale ranges, not zero prior influence
or robustness to every prior family.

## Availability and confirmation freeze

All 510 baseline fits and all 750 candidate fits returned, including every
prior alternative. The unchanged R-hat/ESS screen was enforced, and the
availability nonregression screen passes. Candidate retained iterations were
4k for 419 fits, 8k for 302, 16k for 27 and 32k for two. Zero observed refusals
are not a guarantee of zero future refusal probability.

Retain the Gaussian joint likelihood, factorized HW/matrix-normal priors,
`prior_scale = 1`, `regression_scale = 10`, both declared sensitivity axes,
posterior-mean/equal-tailed summary, numerical schedule and all separate gates.
The exact 35-file D2 source tree is preserved under
`validation/joint-marginal-factorized-development-source` at archive commit
`18a0542`. Amendment 3 documents the one-paired-site input correction. The
[equivalence audit](v0.10.0-development/single-paired-site-equivalence.json)
matched all 270 public input hashes and verified identical prepared model
objects and serialized bytes. Every D2 cohort had at least 12 paired sites and
never entered the corrected branch. No fitting or accuracy reevaluation was
performed by that audit; no D2 cohort or fitted result is replaced.

Then run the untouched confirmation grid: 400 cohorts in each of nine cells,
with sensitivity fits on the first 200 in each of four fixed cells. This is
3,600 cohorts, 12,000 MCMC fits plus 3,600 oracle workflows and 27,600 analyses.
The data/imputation bases remain 1710000000/1720000000; the cohort-bootstrap
seed remains 1740000000. Any failed confirmation component must narrow release
claims; it cannot be hidden by success in another component or repaired by
retuning against these confirmation outcomes.
