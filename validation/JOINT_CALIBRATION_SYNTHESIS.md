# Joint posterior inference: fresh v0.9.0 confirmation

The three prospective statistical screens passed: average calibration across the four primary cells, relative interval-width sensitivity to the covariance prior, and fit availability. This supports the explicit `catmech_joint_summary()` API under its joint-model assumptions. It does **not** establish uniformly nominal coverage: the null-slope cell covered 93.0%, and the contextual-X stress covered only 63.5% for the marginal slope. The covariance priors remain unchanged and materially affect interval width.

The calibration change requires calling the new posterior-summary API. `catmech_pool_lm()` still performs the existing CV3/Rubin analysis; its intervals did not acquire the direct posterior's calibration merely because the numerical sampler improved.

## Design and integrity

The [prospective plan](JOINT_CALIBRATION_PLAN.md) preceded accuracy runs. The [development decision](JOINT_CALIBRATION_DECISION.md) froze the unchanged-prior candidate before fresh confirmation. Development's null result, 25/30 covering intervals, remains in that record; no interval multiplier, nominal level, prior, diagnostic threshold or confirmation seed was tuned to repair it. The prior-development branch was not activated because the prespecified development width-sensitivity screen passed.

Confirmation used 400 independent cohorts in each of six original cells and 200 in the contextual-X stress: 2,600 cohorts, 6,400 MCMC fits, 2,600 full-truth oracle analyses and 12,800 analysis requests. Alternative priors used the first 200 cohorts in each of three fixed cells. These repeated fits are paired analyses, not additional independent cohorts. All 12,755 returned intervals and 45 refusals are retained. The candidate's CV3 and direct analyses use the same fitted posterior in each cohort. The old comparator loads exact frozen v0.8.0 sources, rather than an unrecorded installed version.

The [completed run](joint-calibration-confirmation/COMPLETE.json) and [settings](joint-calibration-confirmation/settings.json) record the new seed streams, 16 outer workers, `pan` 2.0 and 1,721.4 seconds elapsed. All 27 hashed source files agreed before and after the run, fingerprint `a260c005f1ddda9fb27f63225d7c312bd9e98a384a30383f75e089922b372b6a`; the candidate study source was frozen at `b0d4211`, and the confirmation decision at `a9034d8`. The [independent verifier](joint-calibration-confirmation/verification.json) checked the exact grid, all refusals, 80,600 diagnostic rows, raw arithmetic, 1,520 summary fields, 198 paired contrasts, all 2,000 paired bootstrap draws, prior triples and all gates. Nine temporary-copy tamper controls were rejected. It uses base R only to reproduce the predeclared bootstrap indices, then computes the bootstrap statistics independently in Python.

## Separate prospective screens

The accompanying four-panel evidence figure is available as [PNG](joint-calibration-figures/joint-calibration-evidence.png) and [PDF](joint-calibration-figures/joint-calibration-evidence.pdf).

Values below are Monte Carlo comparisons across simulated cohorts. Percentage-point (pp) differences are not differences in individual-study confidence levels. The complete machine-readable decisions are in [adoption.json](joint-calibration-confirmation/adoption.json).

| Screen | Confirmation result | Interpretation |
|---|---|---|
| Primary calibration | Direct mean absolute distance from 95%: 1.625 pp; old CV3: 3.168 pp; candidate CV3: 3.125 pp. Direct-minus-old difference −1.543 pp, paired bootstrap 95% interval [−2.602, −0.356]; direct-minus-candidate −1.500 pp [−2.563, −0.313]. | Passed the fixed average-improvement, per-cell 92.5–97.5% coverage and absolute-bias ≤0.02 screens. Improvement was not uniform across cells. |
| Relative width sensitivity | Equal-cell mean `S_direct − S_CV3 = −0.022691`, paired MCSE 0.001449; 95% Monte Carlo interval [−0.025531, −0.019851]. Every cell contrast was negative. | Passed for the specified joint perturbation of both covariance-prior scales. This is a width-specific claim. |
| Fit availability | Original six cells: 0/2,400 candidate refusals versus 33/2,400 old-profile refusals. Paired equal-cell failure-rate difference −1.375 pp, MCSE 0.237 pp; 95% Monte Carlo interval [−1.840, −0.910] pp. | Passed with unchanged diagnostic thresholds. Zero observed failures is not a guarantee of future convergence. |

## Coverage, bias and precision

All candidate default-prior analyses returned. Therefore their conditional and per-request coverage coincide. Old CV3 coverage below is conditional on return, with the returned/requested denominator shown; [summary.csv](joint-calibration-confirmation/summary.csv) separately records old per-request coverage, counting refusals as no covering interval. Wilson intervals quantify simulation uncertainty. The oracle uses otherwise unavailable full truth and is a complete-data CV3 comparator, not an oracle for the Bayesian posterior.

| Cell | Direct covered/requested | Direct coverage: MCSE; Wilson 95% interval | Candidate CV3 coverage | Old CV3 coverage (returned/requested) | Full-truth CV3 coverage |
|---|---:|---:|---:|---:|---:|
| Balanced Gaussian | 387/400 | 96.75%: 0.89 pp; 94.52–98.09% | 99.00% | 98.99% (396/400) | 95.75% |
| Unbalanced Gaussian | 385/400 | 96.25%: 0.95 pp; 93.91–97.71% | 99.00% | 98.98% (392/400) | 96.00% |
| Observed-data MAR | 386/400 | 96.50%: 0.92 pp; 94.21–97.90% | 97.75% | 97.96% (392/400) | 96.25% |
| Null slope | 372/400 | 93.00%: 1.28 pp; 90.07–95.11% | 96.75% | 96.74% (399/400) | 94.50% |
| Near-zero ICC | 376/400 | 94.00%: 1.19 pp; 91.23–95.94% | 96.50% | 96.75% (400/400) | 93.50% |
| Variance-matched t5 | 380/400 | 95.00%: 1.09 pp; 92.40–96.74% | 97.75% | 97.68% (388/400) | 95.00% |
| Contextual-X stress | 127/200 | 63.50%: 3.40 pp; 56.63–69.86% | 96.50% | 96.28% (188/200) | 93.00% |

The null cell is a remaining calibration concern, not a result to dismiss because the screen passed. Direct inference excluded zero in 28/400 null cohorts. Its posterior RMS spread was 0.07750 against empirical estimate SD 0.08425. Candidate CV3 covered 15 more of these same cohorts: direct-minus-CV3 coverage −3.75 pp, paired MCSE 1.13 pp. Direct coverage was within the prospectively allowed band, but its point distance from 95% was slightly larger than candidate CV3's in this cell. Fresh confirmation is more informative than the small development sample; it does not turn the observed null undercoverage into evidence of exact calibration.

For precision, entries in parentheses below are MCSEs. Width differences use paired cohorts, not independent-method standard errors. Sources: [raw records](joint-calibration-confirmation/records.csv), [summaries](joint-calibration-confirmation/summary.csv) and [paired contrasts](joint-calibration-confirmation/paired.csv).

| Cell | Direct bias (MCSE) | Direct mean width (MCSE) | Candidate CV3 mean width | Paired direct-minus-CV3 width (MCSE) |
|---|---:|---:|---:|---:|
| Balanced Gaussian | 0.00731 (0.00371) | 0.31378 (0.00171) | 0.41212 | −0.09834 (0.00229) |
| Unbalanced Gaussian | −0.00070 (0.00379) | 0.31200 (0.00185) | 0.40927 | −0.09727 (0.00218) |
| Observed-data MAR | 0.00004 (0.00414) | 0.34038 (0.00219) | 0.44770 | −0.10732 (0.00257) |
| Null slope | −0.00493 (0.00421) | 0.30267 (0.00170) | 0.39239 | −0.08972 (0.00212) |
| Near-zero ICC | −0.00580 (0.00376) | 0.28455 (0.00154) | 0.33970 | −0.05516 (0.00183) |
| Variance-matched t5 | −0.00572 (0.00391) | 0.31219 (0.00220) | 0.41359 | −0.10140 (0.00265) |
| Contextual-X stress | −0.12789 (0.00508) | 0.30655 (0.00236) | 0.46270 | −0.15615 (0.00528) |

The interval score improved in the four primary cells; paired direct-minus-CV3 differences ranged from −0.04095 to −0.09347. The near-zero-ICC difference was −0.00988 (MCSE 0.01207), which does not demonstrate an improvement there. In contextual-X it worsened by 0.54985 (MCSE 0.08767). Narrower intervals alone are therefore not the success criterion.

The contextual failure limits the API's interpretation. Its marginal population slope is 0.9, whereas the within-site structural coefficient is 0.7; the data violate the fitted independent-site-effect assumption. All procedures were scored against the same marginal 0.9 target. Direct bias was −0.12789 versus candidate CV3 bias −0.03719 (MCSE 0.00657). Posterior convergence did not diagnose this estimand/model mismatch. The direct summary must not replace an arbitrary marginal or cluster-robust analysis, and this study does not establish performance for a downstream regression having the imputed target as a predictor.

## Prior sensitivity remains visible

For each cohort, `S` is the mean absolute log width ratio at scales 0.25 and 4 relative to scale 1. Each comparison below uses all six returned analyses from the same 200 cohorts: three scales times two analysis methods. No alternative-prior subset is compared with the full 400-cohort default result. Both the row and site covariance-prior scales change together; these contrasts do not isolate either prior component.

| Cell | Common triples | Mean direct S | Mean CV3 S | Paired difference (MCSE) |
|---|---:|---:|---:|---:|
| Balanced Gaussian | 200 | 0.09623 | 0.10688 | −0.01065 (0.00267) |
| Near-zero ICC | 200 | 0.09801 | 0.13045 | −0.03244 (0.00191) |
| Observed-data MAR | 200 | 0.09848 | 0.12346 | −0.02498 (0.00285) |

Mean absolute point-estimate movement across the same alternative scales does not improve in every cell. The tiny increases in the first two cells are uncertain, but preclude describing the result as a demonstrated reduction in every aspect of prior sensitivity. Sources: [prior triples](joint-calibration-confirmation/prior-triples.csv) and [prior summary](joint-calibration-confirmation/prior-summary.csv).

| Cell | Direct mean absolute estimate shift | CV3 mean absolute estimate shift | Paired difference (MCSE) |
|---|---:|---:|---:|
| Balanced Gaussian | 0.009380 | 0.009319 | 0.000061 (0.000271) |
| Near-zero ICC | 0.008512 | 0.008339 | 0.000173 (0.000244) |
| Observed-data MAR | 0.011916 | 0.013271 | −0.001355 (0.000317) |

All 200 fits returned at every scale in each cell. Thus the following covered counts also give per-request coverage. The default rows are recalculated on the **first 200** cohorts in [records.csv](joint-calibration-confirmation/records.csv). Width MCSEs are in parentheses.

| Cell / prior scale | Direct covered/200 | CV3 covered/200 | Direct width (MCSE) | CV3 width (MCSE) |
|---|---:|---:|---:|---:|
| Balanced / 0.25 | 190 | 194 | 0.29842 (0.00225) | 0.38847 (0.00425) |
| Balanced / 1 | 193 | 198 | 0.31327 (0.00233) | 0.40988 (0.00445) |
| Balanced / 4 | 197 | 200 | 0.36183 (0.00280) | 0.47865 (0.00524) |
| Near-zero / 0.25 | 186 | 191 | 0.27055 (0.00187) | 0.31578 (0.00327) |
| Near-zero / 1 | 188 | 194 | 0.28602 (0.00200) | 0.34042 (0.00359) |
| Near-zero / 4 | 192 | 197 | 0.32919 (0.00230) | 0.41056 (0.00453) |
| MAR / 0.25 | 190 | 194 | 0.32701 (0.00310) | 0.42797 (0.00531) |
| MAR / 1 | 194 | 196 | 0.34405 (0.00329) | 0.45465 (0.00566) |
| MAR / 4 | 198 | 199 | 0.39872 (0.00409) | 0.54903 (0.00868) |

Scale four still substantially widens direct intervals and makes these particular cells more conservative. The result supports reduced *relative width sensitivity against this CV3 comparator*, not insensitivity to the retained inverse-Wishart priors or a generally calibrated choice among the three scales.

## Numerical availability and scope

All 3,800 candidate fits returned, including the 1,200 alternative-prior fits; all 7,600 associated analyses returned. The old default profile refused 4 balanced, 8 unbalanced, 0 near-zero, 8 MAR, 1 null and 12 t5 requests, plus 12 contextual-X requests. All were convergence refusals, preserved in [errors.csv](joint-calibration-confirmation/errors.csv). There were no manual successful-fit replacements. The six-cell primary availability comparison is in [availability.csv](joint-calibration-confirmation/availability.csv); the contextual stress is reported separately rather than improving the gate by adding another favorable cell.

Four chains, burn-in 1,000, initial 4,000 retained iterations and the fixed 4k→8k→16k→32k schedule were used. R-hat ≤1.05 and bulk/tail ESS ≥400 were unchanged, with an additional original-scale intercept monitor. Among the 2,600 default-prior candidate fits, 1,598 stopped at 4k, 957 at 8k and 45 at 16k. None required 32k in this run. Zero refusals in a 400-cohort cell has a Wilson 95% upper failure-rate bound of 0.951%; it is not evidence of zero population failure probability. The separate selected 22-case numerical challenge is supportive engineering evidence, not an independent estimate of refusal frequency.

Recorded mean per-fit elapsed time was 1.481 seconds for the candidate default and 5.080 seconds for the old profile in this 16-worker study. These are descriptive concurrent-run timings; this paired accuracy study is not a controlled hardware benchmark. The source-level vectorization/equivalence checks provide separate evidence that faster trace construction preserved the old-control sampling results.

The direct summary uses all retained coefficient draws; the CV3 comparisons here use 40 imputations. An [end-to-end test](../tests/testthat/test-joint-cluster.R) confirms identical full posterior traces and direct summaries at `m = 20` and `m = 40` for the same tested fixture and seed. This does not extend the pooled CV3 calibration evidence to other imputation counts. Calibration and availability are distinct from the historical [v0.8.0 200-cohort study](JOINT_CLUSTER_SYNTHESIS.md). This new evidence concerns 24-site synthetic cohorts, exact references, complete predictors, the specified observed-data selection mechanisms and a target-as-outcome regression. The unbalanced cell also has 35% rather than 30% missing recorded values because of within-site rounding, so it does not isolate imbalance. The variance-matched t5 stress is one particular departure, not general robustness evidence. Contextual dependence, inaccurate references, nonignorable reference selection and different downstream estimands require separate justification and validation.
