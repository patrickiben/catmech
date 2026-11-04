# Stage D decision: freeze the unchanged candidate for confirmation

The development study is complete and independently verified. **Proceed to the prespecified fresh confirmation with the current direct-posterior candidate, unchanged inverse-Wishart priors and numerical controls. This is permission to test the frozen candidate, not an adoption or calibration claim.** No interval multiplier, nominal level, degrees-of-freedom adjustment or prior was selected to make development coverage approach 95%.

The [prospective plan](JOINT_CALIBRATION_PLAN.md) governs the decision. Development used frozen source `b0d4211`: 210 cohorts, 600 MCMC fits and 210 complete-truth oracle analyses, yielding 1,200 analysis requests. There were 1,198 returned intervals and two retained refusals, both from the old v0.8.0 sampler profile. All 390 candidate fits returned, including 210 default-prior fits and 180 alternative-prior fits; both candidate analyses returned from each. The independent verifier checked the complete grid, raw calculations, 198 paired comparisons, the exact 2,000-replicate cohort bootstrap, all three screens and source hashes. See the development [completion record](joint-calibration-development/COMPLETE.json), [summaries](joint-calibration-development/summary.csv), [paired results](joint-calibration-development/paired.csv), [screens](joint-calibration-development/adoption.json) and [raw requests](joint-calibration-development/records.csv).

## Calibration remains unresolved

| Primary cell | Direct posterior covers / returns | Candidate CV3 | Old CV3 |
| --- | ---: | ---: | ---: |
| Gaussian balanced | 30 / 30 | 30 / 30 | 30 / 30 |
| Gaussian unbalanced | 29 / 30 | 30 / 30 | 30 / 30 |
| Observed-data MAR | 29 / 30 | 30 / 30 | 29 / 29 |
| Null slope | 25 / 30 | 27 / 30 | 27 / 30 |

The calibration screen **did not pass**. The mean absolute coverage distance from 95% was 0.05 for each analysis, and the paired bootstrap interval for the direct-minus-CV3 difference included zero. Development cannot support a claim of improved calibration.

The null result is a substantive warning, not an observation to discard: direct coverage was 83.3% (MCSE 6.8 percentage points; Wilson interval 66.4%–92.7%). Its bias was small, 0.0044 (MCSE 0.0164), but empirical estimate SD was 0.0901 versus RMS posterior SD 0.0787. Candidate CV3 had similar empirical SD, 0.0904, and RMS reported SE 0.0884. The complete-truth oracle covered 29/30 with empirical SD 0.0590 and RMS SE 0.0635; its extra truth makes it a different-information comparator.

Direct intervals missed in null replicates 4, 8, 9, 15 and 17. The CV3 analyses missed in 7, 8 and 9, and the oracle missed in 9. The direct-minus-candidate-CV3 coverage difference was −6.7 percentage points with paired MCSE 6.7 points. The missed direct endpoints were more than five reported endpoint Monte Carlo SEs from zero; there is no evidence that a small extra sampling run would repair these particular misses. Independent transformation, conditional-Gaussian and source checks found no material implementation defect. This motivates the already planned 400 fresh null cohorts, with the unfavorable development result retained, rather than a numerical width adjustment or a claim that the warning has been explained away.

## Relative prior sensitivity passed its development screen

All 90 common prior triples were available. The prespecified sensitivity statistic averages absolute log-width changes at scales 0.25 and 4 relative to scale one.

| Cell | Direct sensitivity | Candidate CV3 sensitivity | Paired difference (MCSE) |
| --- | ---: | ---: | ---: |
| Gaussian balanced | 0.0892 | 0.1086 | −0.0194 (0.0047) |
| Near-zero ICC | 0.0966 | 0.1351 | −0.0385 (0.0042) |
| Observed-data MAR | 0.1038 | 0.1193 | −0.0155 (0.0079) |

The equal-cell mean difference was −0.02448 (MCSE 0.00338; 95% Monte Carlo interval −0.03110 to −0.01786), and no cell's point contrast was positive. Thus the prespecified Huang–Wand development branch is **not activated**. The current priors remain substantive assumptions; this result concerns relative interval-width sensitivity within the tested scale grid. Fresh confirmation is required, and estimate shifts and failures remain part of the assessment. [Prior triples](joint-calibration-development/prior-triples.csv), [prior summaries](joint-calibration-development/prior-summary.csv).

## Numerical readiness and scientific limits

All candidate fits returned without weaker R-hat or ESS thresholds. The two old-profile refusals were MAR replicate 5 and t5 replicate 15. The equal-cell candidate-minus-old refusal difference was −0.01111, but its 95% Monte Carlo interval, −0.02651 to 0.00429, included zero. The demonstrated-availability screen therefore **did not pass in development**; its population improvement remains to be tested. The separately selected 22 historical failures all passed the enlarged budget, but that challenge set cannot estimate a future refusal rate.

The planned contextual-X stress is an important limitation: direct posterior inference covered the marginal slope in only 18/30 cohorts, with bias −0.1233; candidate CV3 covered 29/30, with bias −0.0243. Here the marginal population target is 0.9, while the structural within-site coefficient is 0.7, and the fitted independent-site-effect assumption fails. This result rules out presenting the direct summary as a drop-in robust marginal-effect analysis. It remains an explicit summary of the fitted joint model's coefficients. [Contextual and other results](joint-calibration-development/summary.csv).

## Freeze for Stage C

Retain the same likelihood, inverse-Wishart priors and scale perturbations, direct equal-tail summary, CV3 comparator, diagnostic thresholds and numerical controls: four chains, burn-in 1,000, initial 4,000 retained iterations and maximum 32,000. Do not add the hierarchical prior candidate or change the calibration screens after this development result.

Proceed with the untouched confirmation seeds and grid already declared in the plan: 400 cohorts in each original cell, 200 contextual-X cohorts, and prior alternatives on the first 200 cohorts in the three fixed sensitivity cells. The confirmation contains 2,600 cohorts and 12,800 analysis requests. The original implementation, calibration, sensitivity and availability screens remain separate and unchanged. A failed confirmation screen must be reported and narrow the claim; it does not authorize tuning against confirmation data.
